#!/usr/bin/env python3
# [[03_variant_calling/11_asm_build_matrix.py]]
# Goal: Build the tidy per-(locus, subject) allele-count matrix for the ASM group-comparison.
#       DATA MUNGING ONLY (no statistics) — all hypothesis testing happens in
#       03_variant_calling/12_asm_compare.R using established R packages.
#
# Counts come from GATK ASEReadCounter (identical pipeline for IP and Input fractions),
#   *_ip_ase_counts.csv / *_input_ase_counts.csv -> (subject, chr, pos) ref/alt counts.
# M6Allele *_asmPeak.txt provide the m6A-peak -> gene mapping, the per-subject peak q-value,
#   and the SNP positions (majorAlleleBase) that define which sites belong to a peak.
#
# Steps:
#   1. Load metadata (sample_sheet.csv) -> subject -> (condition, sex).
#   2. Parse *_asmPeak.txt -> per (subject, peak): chr, start, end, gene, q, SNP positions.
#      Dedup identical (subject, chr, start, end) rows (keep min q, union SNP positions).
#   3. If _data/peaks/consensus_peaks.bed exists, use it as authoritative loci; else
#      cluster per-subject peak intervals by overlap (union-merge) into consensus loci.
#   4. Load IP (and, if present, Input) ASEReadCounter counts.
#   5. Aggregate per (locus, subject): ip_ref/ip_total, in_ref/in_total (Input NA if absent),
#      min_q, n_snps, gene. Allele A = REF allele (fixed genome base, consistent across subjects).
#   6. Write tidy TSV to _data/RData/03_variant_calling/asm_matrix.tsv for the R analysis.

import csv
import os
import pickle
import re

# --- Dynamic Pathing ---
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.environ.get("PROJECT_ROOT") or os.path.normpath(os.path.join(SCRIPT_DIR, ".."))

ASM_DIR    = os.path.join(PROJECT_ROOT, "_data", "asm")
ASE_DIR    = os.path.join(PROJECT_ROOT, "_data", "ase")
PEAKS_BED  = os.path.join(PROJECT_ROOT, "_data", "peaks", "consensus_peaks.bed")
META_FILE  = os.path.join(PROJECT_ROOT, "_data", "_metadata", "sample_sheet.csv")
OUT_DIR    = os.path.join(PROJECT_ROOT, "_data", "RData", "03_variant_calling")
LOG_DIR    = os.path.join(PROJECT_ROOT, "_data", "logs", "asm")
for d in (OUT_DIR, LOG_DIR):
    os.makedirs(d, exist_ok=True)
OUT_TSV = os.path.join(OUT_DIR, "asm_matrix.tsv")
PICKLE  = os.path.join(OUT_DIR, "asm_matrix.pkl")
LOG_FILE = os.path.join(LOG_DIR, "11_asm_build_matrix.log")

CHR_ORDER = {str(i): i for i in range(1, 23)} | {"X": 23, "Y": 24, "MT": 25}


def log(msg):
    import datetime
    line = "[%s] %s" % (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), msg)
    print(line, flush=True)
    with open(LOG_FILE, "a") as fh:
        fh.write(line + "\n")


def chr_key(c):
    return CHR_ORDER.get(c, 99)


def load_metadata():
    subjects = {}
    with open(META_FILE) as fh:
        for row in csv.DictReader(fh):
            sample = re.sub(r"_[12]$", "", row["File1"])
            subj = re.sub(r"_(input|m6A)$", "", sample)
            subjects[subj] = {"condition": row["SampleType"], "sex": row["Gender"]}
    return subjects


def parse_bases(s):
    """'16600576:G;16600733:C' -> [16600576, 16600733]"""
    out = []
    for item in s.split(";"):
        item = item.strip()
        if ":" in item:
            pos = item.rsplit(":", 1)[0]
            try:
                out.append(int(pos))
            except ValueError:
                pass
    return out


def load_asm_peaks(subjects):
    """Parse *_asmPeak.txt -> dict[(subject, chr, start, end)] = peak record.
       Dedup identical coords (keep min q), union SNP positions."""
    peaks = {}
    nfiles = 0
    nrows = 0
    ndup = 0
    for fn in sorted(os.listdir(ASM_DIR)):
        if not fn.endswith("_asmPeak.txt"):
            continue
        subj = fn.replace("_asmPeak.txt", "")
        if subj not in subjects:
            continue
        nfiles += 1
        with open(os.path.join(ASM_DIR, fn)) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                nrows += 1
                try:
                    q = float(r["qValue"])
                    start = int(r["peakStart"])
                    end = int(r["peakEnd"])
                except (KeyError, ValueError):
                    continue
                pos = parse_bases(r.get("majorAlleleBase", ""))
                key = (subj, r["chr"], start, end)
                rec = peaks.get(key)
                if rec is None:
                    peaks[key] = {"subject": subj, "chr": r["chr"], "start": start, "end": end,
                                  "gene_id": r.get("geneId", ""), "gene_name": r.get("geneName", ""),
                                  "q": q, "pos": set(pos)}
                else:
                    ndup += 1
                    if q < rec["q"]:
                        rec["q"] = q
                        rec["gene_id"] = r.get("geneId", "")
                        rec["gene_name"] = r.get("geneName", "")
                    rec["pos"] |= set(pos)
    log("asmPeak: %d files, %d rows, %d duplicate coords merged" % (nfiles, nrows, ndup))
    return peaks


def load_consensus_loci():
    """Return list of {chr,start,end} from TRESS consensus_peaks.bed (0-based start)."""
    loci = []
    with open(PEAKS_BED) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            try:
                loci.append({"chr": p[0], "start": int(p[1]), "end": int(p[2])})
            except ValueError:
                continue
    return loci


def cluster_loci(peaks):
    """Union-merge peak intervals per chromosome into consensus loci."""
    by_chr = {}
    for rec in peaks.values():
        by_chr.setdefault(rec["chr"], []).append(rec)
    loci = []
    for chrm in sorted(by_chr, key=chr_key):
        for rec in sorted(by_chr[chrm], key=lambda x: (x["start"], x["end"])):
            if loci and loci[-1]["chr"] == chrm and rec["start"] <= loci[-1]["end"]:
                loci[-1]["end"] = max(loci[-1]["end"], rec["end"])
            else:
                loci.append({"chr": chrm, "start": rec["start"], "end": rec["end"]})
    return loci


def build_locus_index(loci):
    """Map (chr, start, end) -> locus record for exact matches."""
    return {(lo["chr"], lo["start"], lo["end"]): lo for lo in loci}


def find_locus_by_pos(loci, chrm, pos):
    """First locus on chr containing pos (loci are disjoint after union-merge)."""
    for lo in loci:
        if lo["chr"] == chrm and lo["start"] <= pos <= lo["end"]:
            return lo
    return None


def load_ase_counts(ase_dir, subjects, fraction):
    """Parse '*_<fraction>_ase_counts.csv' -> {(subject, chr, pos): (ref, alt, total)}."""
    counts = {}
    nfiles = 0
    for fn in sorted(os.listdir(ase_dir)):
        if not fn.endswith("_%s_ase_counts.csv" % fraction):
            continue
        subj = fn.replace("_%s_ase_counts.csv" % fraction, "")
        if subj not in subjects:
            continue
        nfiles += 1
        with open(os.path.join(ase_dir, fn)) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                try:
                    pos = int(row["position"])
                    ref = int(row["refCount"])
                    alt = int(row["altCount"])
                except (KeyError, ValueError):
                    continue
                counts[(subj, row["contig"], pos)] = (ref, alt, ref + alt)
    log("ASE counts [%s]: %d files, %d sites" % (fraction, nfiles, len(counts)))
    return counts


def main():
    log("=== 11_asm_build_matrix.py ===")
    subjects = load_metadata()
    log("subjects: %d" % len(subjects))

    peaks = load_asm_peaks(subjects)

    # --- Loci: authoritative consensus BED if available, else overlap clustering ---
    if os.path.isfile(PEAKS_BED):
        loci = load_consensus_loci()
        log("loci: %d from %s" % (len(loci), PEAKS_BED))
        by_chr_loci = {}
        for lo in loci:
            by_chr_loci.setdefault(lo["chr"], []).append(lo)
        def find_locus(chrm, pos):
            for lo in by_chr_loci.get(chrm, []):
                if lo["start"] <= pos <= lo["end"]:
                    return lo
            return None
    else:
        loci = cluster_loci(peaks)
        log("loci: %d by overlap clustering of %d unique peaks (consensus_peaks.bed not present)"
            % (len(loci), len(peaks)))
        find_locus = lambda chrm, pos: find_locus_by_pos(loci, chrm, pos)

    for lo in loci:
        lo["id"] = "%s:%d:%d" % (lo["chr"], lo["start"], lo["end"])

    # --- Assign each peak to a locus (exact containment) ---
    locus_idx = build_locus_index(loci)
    peak_locus = {}
    for key, rec in peaks.items():
        lo = locus_idx.get((rec["chr"], rec["start"], rec["end"]))
        if lo is None:
            lo = find_locus(rec["chr"], rec["start"] + (rec["end"] - rec["start"]) // 2)
        if lo is None:
            continue
        peak_locus[key] = lo

    log("peaks assigned to loci: %d / %d" % (len(peak_locus), len(peaks)))

    # --- Load counts ---
    ip_counts = load_ase_counts(ASE_DIR, subjects, "ip")
    input_counts = load_ase_counts(ASE_DIR, subjects, "input")

    # --- Aggregate per (locus, subject) ---
    # row: locus_id chr start end gene subject condition sex ip_ref ip_total in_ref in_total in_avail min_q n_snps
    agg = {}
    for key, lo in peak_locus.items():
        subj, chrm, start, end = key
        rec = peaks[key]
        gk = (lo["id"], subj)
        row = agg.setdefault(gk, {
            "locus_id": lo["id"], "chr": lo["chr"], "start": lo["start"], "end": lo["end"],
            "subject": subj, "condition": subjects[subj]["condition"], "sex": subjects[subj]["sex"],
            "ip_ref": 0, "ip_total": 0, "in_ref": 0, "in_total": 0, "in_avail": 0,
            "min_q": rec["q"], "n_snps": 0, "n_peaks": 0, "genes": {},
        })
        row["min_q"] = min(row["min_q"], rec["q"])
        row["n_peaks"] += 1
        row["genes"][rec["gene_name"]] = row["genes"].get(rec["gene_name"], 0) + 1
        for pos in rec["pos"]:
            ip = ip_counts.get((subj, chrm, pos))
            if ip is None:
                continue
            row["ip_ref"] += ip[0]
            row["ip_total"] += ip[2]
            row["n_snps"] += 1
            inn = input_counts.get((subj, chrm, pos)) if input_counts else None
            if inn is not None:
                row["in_ref"] += inn[0]
                row["in_total"] += inn[2]
                row["in_avail"] = 1

    log("aggregated (locus, subject) rows: %d" % len(agg))

    # --- Write TSV ---
    fields = ["locus_id", "chr", "start", "end", "gene", "subject", "condition", "sex",
              "ip_ref", "ip_total", "in_ref", "in_total", "in_avail", "min_q", "n_snps", "n_peaks"]
    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in agg.values():
            gene = max(row["genes"], key=lambda g: (row["genes"][g], -len(g))) if row["genes"] else ""
            row["gene"] = gene
            if row["in_avail"]:
                row["in_ref"], row["in_total"] = row["in_ref"], row["in_total"]
            else:
                row["in_ref"], row["in_total"] = "", ""
            w.writerow(row)

    with open(PICKLE, "wb") as fh:
        pickle.dump({"loci": loci, "n_loci": len(loci), "n_rows": len(agg),
                     "input_available": bool(input_counts), "n_input_files": len(set(k[0] for k in input_counts))}, fh)

    log("written: %s (%d rows)" % (OUT_TSV, len(agg)))
    log("input counts available: %s" % bool(input_counts))
    log("done.")


if __name__ == "__main__":
    main()
