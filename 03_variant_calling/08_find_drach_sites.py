#!/usr/bin/env python3
# [[03_variant_calling/08_find_drach_sites.py]]
# Goal: Scan TRESS consensus peak sequences for all DRACH motif occurrences
#       and convert each hit to genomic coordinates (BED format).

# Steps:
#   1. Parse consensus_peaks.bed to build a coordinate lookup per peak.
#   2. Extract peak sequences from consensus_peaks.fa (bedtools getfasta output).
#   3. Scan each sequence for DRACH regex: [AGT][AG]AC[ACT] (+ and - strand).
#   4. Convert local motif offsets to genomic chromosome:start-end coordinates.
#   5. Write all DRACH sites to _data/peaks/drach_sites.bed.
#   6. Assert output is non-empty.

import os
import re
import subprocess
import sys

# --- Dynamic Pathing ---
script_dir   = os.path.dirname(os.path.abspath(__file__))
project_root = os.environ.get("PROJECT_ROOT") or os.path.normpath(
    os.path.join(script_dir, "..")
)

PEAKS_BED = os.path.join(project_root, "_data", "peaks", "consensus_peaks.bed")
PEAKS_FA  = os.path.join(project_root, "_data", "peaks", "consensus_peaks.fa")
OUT_BED   = os.path.join(project_root, "_data", "peaks", "drach_sites.bed")
LOG_DIR   = os.path.join(project_root, "_data", "logs", "drach")
os.makedirs(LOG_DIR, exist_ok=True)

# DRACH regex: D=[AGT], R=[AG], A=A, C=C, H=[ACT]
DRACH_FWD = re.compile(r"(?=([AGT][AG]AC[ACT]))")
# Reverse complement: DRACH on minus strand → look for H'CA'C'R'D' in the fwd string
# i.e., complement of DRACH = [ACT][CT]GT[ACG], reversed = [GCA][TC]GT[TCА]
# Simplified: detect DRACH on both strands by running on the sequence and its RC
COMPLEMENT = str.maketrans("ACGTagct", "TGCAtcga")

def rev_comp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]

# --- Step 1: Parse BED to build coordinate lookup ---
# BED format: chr  start  end  peak_name  (0-based start, end exclusive)
print("[INFO] Loading consensus peaks BED...")
peak_coords = {}   # name → (chrom, start, end, strand)

with open(PEAKS_BED) as fh:
    for line in fh:
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
        chrom  = parts[0]
        start  = int(parts[1])
        end    = int(parts[2])
        strand = parts[5] if len(parts) >= 6 else "+"
        # bedtools getfasta -name uses field 4 as sequence ID; fall back to coord key
        name   = parts[3] if len(parts) >= 4 else f"{chrom}:{start}-{end}"
        peak_coords[name] = (chrom, start, end, strand)

print(f"[INFO] Loaded {len(peak_coords)} peaks.")

# --- Step 2: Parse FASTA ---
print("[INFO] Scanning peak FASTA for DRACH motifs...")
if not os.path.isfile(PEAKS_FA):
    print(f"[INFO] consensus_peaks.fa not found. Generating with bedtools getfasta...")
    REF_FA = os.path.join(project_root, "_data", "genome", "Homo_sapiens.GRCh38.dna.primary_assembly.fa")
    if not os.path.isfile(REF_FA):
        print(f"Assertion Error: Reference genome missing: {REF_FA}", file=sys.stderr)
        sys.exit(1)
    result = subprocess.run(
        ["bedtools", "getfasta", "-fi", REF_FA, "-bed", PEAKS_BED, "-name", "-fo", PEAKS_FA],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"Assertion Error: bedtools getfasta failed:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)
    print(f"[INFO] Generated {PEAKS_FA}")

def iter_fasta(path):
    """Yield (header, sequence) pairs from a FASTA file."""
    header, seqparts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seqparts)
                header   = line[1:].split()[0]
                seqparts = []
            else:
                seqparts.append(line.upper())
    if header is not None:
        yield header, "".join(seqparts)

# --- Step 3 & 4: Scan sequences and convert offsets ---
records = []
unresolved = 0

for seq_id, seq in iter_fasta(PEAKS_FA):
    # bedtools getfasta -name writes header as:  peakname::chr:start-end (0-based coords in header)
    # We strip the ::coord suffix to recover the peak name for lookup
    peak_name = seq_id.split("::")[0]
    coord_key  = seq_id  # fallback: try the full ID

    entry = peak_coords.get(peak_name) or peak_coords.get(coord_key)
    if entry is None:
        unresolved += 1
        continue

    chrom, p_start, p_end, p_strand = entry

    # Scan forward strand
    for m in DRACH_FWD.finditer(seq):
        motif_start = m.start()
        motif_end   = motif_start + 5
        genomic_start = p_start + motif_start
        genomic_end   = p_start + motif_end
        # Methylation site = position 3 (0-indexed) of the motif = A
        methyl_pos = p_start + motif_start + 2
        records.append((chrom, genomic_start, genomic_end,
                        f"{peak_name}:DRACH_fwd:{m.group(1)}", 0, "+", methyl_pos))

    # Scan reverse complement (minus-strand DRACH sites)
    rc_seq = rev_comp(seq)
    for m in DRACH_FWD.finditer(rc_seq):
        motif_local_start = m.start()
        motif_local_end   = motif_local_start + 5
        # Convert RC offset back to genomic coords (BED 0-based)
        genomic_end   = p_end - motif_local_start
        genomic_start = genomic_end - 5
        methyl_pos    = genomic_end - 3
        records.append((chrom, genomic_start, genomic_end,
                        f"{peak_name}:DRACH_rev:{m.group(1)}", 0, "-", methyl_pos))

print(f"[INFO] Found {len(records)} DRACH sites across all peaks.")
if unresolved > 0:
    print(f"[WARNING] {unresolved} sequences could not be matched to a peak coordinate.")

# --- Assertion: at least some sites found ---
if len(records) == 0:
    print("Assertion Error: Zero DRACH sites found. Check FASTA and regex logic.")
    sys.exit(1)

# --- Step 5: Write BED output ---
print(f"[INFO] Writing DRACH sites to: {OUT_BED}")
records.sort(key=lambda r: (r[0], r[1]))

with open(OUT_BED, "w") as out:
    out.write("# chrom\tstart\tend\tname\tscore\tstrand\tmethyl_genomic_pos\n")
    for r in records:
        out.write("\t".join(str(x) for x in r) + "\n")

# --- Assertion: output non-empty ---
if os.path.getsize(OUT_BED) == 0:
    print("Assertion Error: Output DRACH BED file is empty.")
    sys.exit(1)

print(f"[DONE] {len(records)} DRACH sites written to {OUT_BED}")
print("       Next: run 09_intersect_drach_snps.sh")
