# [[04_asm/02_verify_drach_reference.py]]
# Goal: Independently verify DRACH motif geometry from peak sequence before any disruption annotation is trusted.
#
# Steps:
#   1. Load consensus_peaks.bed coordinates and consensus_peaks.fa sequence.
#   2. Re-derive each motif 5-mer from the peak FASTA (reverse complement for '-' strands).
#   3. Assert IUPAC validity and methyl_pos == genomic_start + 2.
#   4. When the SNP x DRACH intersection exists, assert REF alleles match donor VCF REF.
#   5. Write the QC summary and exit non-zero on any hard failure.
#
# Notes: deliberately avoids Biostrings so this gate stays independent of the
# producer. Verified to fail on injected corruption. See docs/04_asm_methods.md.

import csv
import gzip
import os
import sys

# --- Dynamic Pathing ---
script_dir   = os.path.dirname(os.path.abspath(__file__))
project_root = os.environ.get("PROJECT_ROOT") or os.path.normpath(
    os.path.join(script_dir, "..")
)

PEAKS_BED = os.path.join(project_root, "_data", "peaks", "consensus_peaks.bed")
PEAKS_FA  = os.path.join(project_root, "_data", "peaks", "consensus_peaks.fa")
DRACH_BED = os.path.join(project_root, "_data", "peaks", "drach_sites.bed")
VCF_DIR   = os.path.join(project_root, "_data", "variant", "filtered_vcf")
INTERSECT = os.path.join(
    project_root, "_data", "logs", "drach", "tmp", "snp_drach_intersect.bed"
)
OUT_DIR   = os.path.join(project_root, "results", "04_asm", "tables")
OUT_CSV   = os.path.join(OUT_DIR, "02_drach_reference_qc.csv")
OUT_FAIL_CSV = os.path.join(OUT_DIR, "02_drach_reference_failures.csv")

# --- IUPAC (transcript orientation, 5'->3') ---
# D = A/G/T (not C) | R = A/G | A = A (core m6A) | C = C | H = A/C/T (not G)
IUPAC_ALLOWED = {
    0: set("AGT"),
    1: set("AG"),
    2: set("A"),
    3: set("C"),
    4: set("ACT"),
}

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rev_comp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def is_valid_drach(motif):
    """Assert motif matches [AGT][AG]AC[ACT] and central base is adenine."""
    if len(motif) != 5:
        return False, "length!=5"
    for offset, base in enumerate(motif):
        if base not in IUPAC_ALLOWED[offset]:
            return False, f"offset{offset}={base} violates IUPAC"
    if motif[2] != "A":
        return False, "central base is not adenine"
    return True, ""


def iter_fasta(path):
    header, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.upper())
    if header is not None:
        yield header, "".join(parts)


# --- Input assertions ---
for label, path in (("consensus_peaks.bed", PEAKS_BED),
                    ("consensus_peaks.fa", PEAKS_FA),
                    ("drach_sites.bed", DRACH_BED)):
    if not os.path.isfile(path):
        print(f"Assertion Error: {label} not found: {path}", file=sys.stderr)
        print("Run 01_find_drach_sites.py first.", file=sys.stderr)
        sys.exit(1)

os.makedirs(OUT_DIR, exist_ok=True)

# --- Step 1: Load peak coordinates and sequences ---
print("[INFO] Loading consensus peak coordinates and sequences...")
peak_coords = {}
with open(PEAKS_BED) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < 4:
            continue
        peak_coords[p[3]] = (p[0], int(p[1]), int(p[2]))

peak_seq = {}
for seq_id, seq in iter_fasta(PEAKS_FA):
    peak_seq[seq_id.split("::")[0]] = seq

print(f"[INFO] Loaded {len(peak_coords):,} peaks "
      f"({len(peak_seq):,} with sequence).")

# --- Step 2: Load DRACH motif BED ---
motifs = []
with open(DRACH_BED) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 7:
            continue
        motifs.append({
            "chrom":   p[0],
            "gs":      int(p[1]),
            "ge":      int(p[2]),
            "name":    p[3],
            "strand":  p[5],
            "methyl":  int(p[6]),
        })

if not motifs:
    print("Assertion Error: drach_sites.bed contains no motifs.", file=sys.stderr)
    sys.exit(1)
print(f"[INFO] Loaded {len(motifs):,} DRACH motifs.")

# --- Step 3: Stage A - motif geometry ---
print("[STEP A] Verifying motif geometry against peak FASTA...")
qc_rows = []
n_checked = n_bad_geom = n_bad_iupac = n_bad_methyl = n_missing_peak = 0

for m in motifs:
    parts = m["name"].split(":")
    peak_name = parts[0]
    recorded = parts[-1].strip().upper()
    strand = m["strand"]

    if peak_name not in peak_seq or peak_name not in peak_coords:
        n_missing_peak += 1
        qc_rows.append({**m, "motif_seq": recorded, "derived_seq": "",
                        "check": "peak_missing", "detail": peak_name})
        n_bad_geom += 1
        continue

    _, p_start, _ = peak_coords[peak_name]
    off = m["gs"] - p_start
    span = peak_seq[peak_name][off:off + 5]

    if len(span) != 5:
        n_bad_geom += 1
        qc_rows.append({**m, "motif_seq": recorded, "derived_seq": span,
                        "check": "span_truncated",
                        "detail": f"off={off} len={len(span)}"})
        continue

    derived = span if strand == "+" else rev_comp(span)

    ok_iupac, why_iupac = is_valid_drach(derived)
    if not ok_iupac:
        n_bad_iupac += 1
    if derived != recorded:
        n_bad_geom += 1
    if m["methyl"] != m["gs"] + 2:
        n_bad_methyl += 1

    n_checked += 1
    if derived != recorded or not ok_iupac or m["methyl"] != m["gs"] + 2:
        qc_rows.append({
            **m, "motif_seq": recorded, "derived_seq": derived,
            "check": "geometry_mismatch" if derived != recorded else why_iupac or "methyl_pos_mismatch",
            "detail": f"off={off} strand={strand} methyl={m['methyl']} expected={m['gs'] + 2}",
        })

print(f"[STEP A] checked={n_checked:,}  geom_mismatch={n_bad_geom:,}  "
      f"iupac_invalid={n_bad_iupac:,}  methyl_pos_wrong={n_bad_methyl:,}  "
      f"missing_peak={n_missing_peak:,}")

# --- Step 4: Stage B - allele agreement (optional) ---
n_b_rows = n_b_bad = n_b_checked = 0
b_detail = ""
if os.path.isfile(INTERSECT) and os.path.getsize(INTERSECT) > 0:
    print("[STEP B] Verifying REF alleles against donor VCFs...")
    vcf_ref = {}
    vcf_files = sorted(
        os.path.join(VCF_DIR, f) for f in os.listdir(VCF_DIR)
        if f.endswith("_input.filtered.vcf.gz")
    ) if os.path.isdir(VCF_DIR) else []

    for vcf_path in vcf_files:
        donor = os.path.basename(vcf_path).replace("_input.filtered.vcf.gz", "")
        with gzip.open(vcf_path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 5:
                    continue
                vcf_ref[(donor, p[0].replace("chr", ""), int(p[1]))] = p[3]

    print(f"       indexed {len(vcf_ref):,} VCF records from {len(vcf_files)} donors")

    with open(INTERSECT) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            n_b_rows += 1
            donor = os.path.basename(p[3]).replace("_input", "").replace("_m6A", "")
            key = (donor, p[0].replace("chr", ""), int(p[2]))
            ref_bed = p[4].upper()
            if key not in vcf_ref:
                continue
            n_b_checked += 1
            if vcf_ref[key] != ref_bed:
                n_b_bad += 1
                if len(qc_rows) < 5000:
                    qc_rows.append({
                        "chrom": p[0], "gs": int(p[2]), "ge": "", "name": donor,
                        "strand": "", "methyl": "", "motif_seq": ref_bed,
                        "derived_seq": vcf_ref[key],
                        "check": "ref_allele_mismatch",
                        "detail": f"bed_REF={ref_bed} vcf_REF={vcf_ref[key]}",
                    })
    b_detail = f" rows={n_b_rows:,} matched={n_b_checked:,} mismatch={n_b_bad:,}"
    print(f"[STEP B] {b_detail}")
else:
    print("[STEP B] SKIPPED - SNP x DRACH intersect absent or empty "
          "(run 03_intersect_drach_genotypes.sh to enable).")

# --- Step 5: Report ---
n_failures = n_bad_geom + n_bad_iupac + n_bad_methyl + n_missing_peak + n_b_bad

summary_items = [
    ("motifs_checked",          n_checked),
    ("geometry_mismatches",     n_bad_geom),
    ("iupac_invalid",           n_bad_iupac),
    ("methyl_pos_incorrect",    n_bad_methyl),
    ("missing_peak_sequence",   n_missing_peak),
    ("ref_allele_mismatches",   n_b_bad),
    ("total_hard_failures",     n_failures),
]

# Summary is always written, so a clean run still leaves an auditable record.
with open(OUT_CSV, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["metric", "value"])
    for k, v in summary_items:
        w.writerow([k, v])
    w.writerow([])
    w.writerow(["VERDICT", "PASS" if n_failures == 0
                else f"FAIL ({n_failures} hard failure(s))"])

# Failure detail is written only when something actually failed.
if n_failures:
    with open(OUT_FAIL_CSV, "w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["chrom", "gs", "ge", "name", "strand", "methyl",
                        "motif_seq", "derived_seq", "check", "detail"],
        )
        w.writeheader()
        w.writerows(qc_rows)

print("=" * 72)
print("   DRACH REFERENCE VERIFICATION")
print(f"   Motifs checked        : {n_checked:,}")
print(f"   Geometry mismatches   : {n_bad_geom:,}")
print(f"   IUPAC invalid         : {n_bad_iupac:,}")
print(f"   methyl_pos incorrect  : {n_bad_methyl:,}")
print(f"   Missing peak sequence : {n_missing_peak:,}")
print(f"   REF allele mismatches : {n_b_bad:,}")
print(f"   Verdict               : {'PASS' if n_failures == 0 else 'FAIL'}")
print(f"   QC summary            : {OUT_CSV}")
if n_failures:
    print(f"   Failure detail        : {OUT_FAIL_CSV}")
print("=" * 72)

if n_failures:
    print(f"VERIFICATION FAILED: {n_failures:,} hard failure(s). "
          f"Do NOT trust downstream disruption annotations.", file=sys.stderr)
    sys.exit(1)

print("VERIFICATION PASSED: all motifs independently re-derived from peak FASTA.")