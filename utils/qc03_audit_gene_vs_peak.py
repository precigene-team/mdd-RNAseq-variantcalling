#!/usr/bin/env python3
# [[utils/audit_gene_vs_peak.py]]
# Goal: Audit gene-level vs peak-level resolution by mapping narrowPeak coordinates to genes
#       and checking overlap/concordance with DESeq2 differential methylation results.

import os
import sys
import re
import argparse
import pandas as pd
from collections import defaultdict
import bisect

def parse_gtf(gtf_path):
    """
    Parses GTF file to extract exonic coordinates for genes.
    Optimized to deduplicate and merge overlapping exons per gene to minimize memory footprint.
    Uses pickle caching and grep streaming for high performance.
    Returns:
        gene_intervals: dict mapping chrom -> list of tuples (start, end, gene_id, gene_name) sorted by start
        gene_id_to_name: dict mapping gene_id -> gene_name
    """
    import pickle
    import subprocess

    cache_path = gtf_path + ".exons.pkl"
    if os.path.exists(cache_path):
        print(f"[INFO] Loading cached exon annotations from: {cache_path} ...")
        try:
            with open(cache_path, 'rb') as f:
                gene_intervals, gene_id_to_name = pickle.load(f)
            # After loading, each chrom value is {'intervals': [...], 'starts': [...]}.
            # Count exon tuples from the nested list, not the dict keys.
            n_exons = sum(len(v['intervals']) for v in gene_intervals.values())
            if n_exons < 100_000:
                print(f"  [WARN] Cache appears corrupt ({n_exons} exons < 100,000 expected). Re-parsing GTF...")
                raise ValueError("corrupt cache")
            print(f"  [SUCCESS] Loaded {n_exons} merged exons from cache.")
            return gene_intervals, gene_id_to_name
        except Exception as e:
            print(f"  [WARN] Cache load failed: {e}. Re-parsing GTF...")

    print(f"[INFO] Parsing GTF annotation file: {gtf_path} ...")
    gene_exons = defaultdict(list)
    gene_id_to_name = {}
    
    # Regular expressions for attributes
    id_re = re.compile(r'gene_id "([^"]+)"')
    name_re = re.compile(r'gene_name "([^"]+)"')
    
    # Attempt to use grep to pre-filter exon lines for maximum speed
    use_grep = True
    try:
        proc = subprocess.Popen(["grep", "-P", r"\texon\t", gtf_path], stdout=subprocess.PIPE, text=True, bufsize=1024*1024)
        infile = proc.stdout
    except Exception:
        use_grep = False
        infile = open(gtf_path, 'r')
        
    try:
        for line in infile:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            # If we didn't use grep, we must check the feature type here
            if not use_grep and fields[2] != 'exon':
                continue
                
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            attributes = fields[8]
            
            id_match = id_re.search(attributes)
            name_match = name_re.search(attributes)
            
            if id_match:
                gene_id = id_match.group(1)
                gene_name = name_match.group(1) if name_match else gene_id
                gene_id_to_name[gene_id] = gene_name
                gene_exons[(chrom, gene_id, gene_name)].append((start, end))
    finally:
        infile.close()
        if use_grep:
            proc.wait()
                
    # Merge overlapping intervals per gene to minimize memory usage
    gene_intervals = defaultdict(list)
    total_loaded_exons = 0
    
    for (chrom, gene_id, gene_name), intervals in gene_exons.items():
        if not intervals:
            continue
        # Sort and merge intervals for this gene
        intervals.sort(key=lambda x: x[0])
        merged = [intervals[0]]
        for current in intervals[1:]:
            prev_start, prev_end = merged[-1]
            curr_start, curr_end = current
            if curr_start <= prev_end:
                merged[-1] = (prev_start, max(prev_end, curr_end))
            else:
                merged.append(current)
                
        for start, end in merged:
            gene_intervals[chrom].append((start, end, gene_id, gene_name))
            total_loaded_exons += 1
            
    # Structure for high-performance lookup (precompute starts)
    final_intervals = {}
    for chrom in gene_intervals:
        intervals = gene_intervals[chrom]
        intervals.sort(key=lambda x: x[0])
        final_intervals[chrom] = {
            'intervals': intervals,
            'starts': [x[0] for x in intervals]
        }
    gene_intervals = final_intervals
        
    print(f"  [SUCCESS] Loaded {total_loaded_exons} merged exons for {len(gene_id_to_name)} genes.")
    
    # Save cache
    try:
        print(f"[INFO] Caching parsed exons to: {cache_path} ...")
        with open(cache_path, 'wb') as f:
            pickle.dump((dict(gene_intervals), gene_id_to_name), f, protocol=pickle.HIGHEST_PROTOCOL)
    except Exception as e:
        print(f"  [WARN] Failed to write cache: {e}")
        
    return gene_intervals, gene_id_to_name

def find_overlapping_genes(chrom, start, end, gene_intervals):
    """
    Finds all genes whose exons overlap with the given interval.
    Uses binary search for high performance.
    """
    if chrom not in gene_intervals:
        # Try prepending/removing 'chr' prefix
        if chrom.startswith('chr'):
            chrom = chrom[3:]
        else:
            chrom = 'chr' + chrom
        if chrom not in gene_intervals:
            return set()
            
    chrom_data = gene_intervals[chrom]
    starts = chrom_data['starts']
    intervals = chrom_data['intervals']
    
    # Binary search to find candidate exons
    idx = bisect.bisect_left(starts, start)
    
    overlapping = set()
    # Search backwards from the bisect index
    i = idx - 1
    while i >= 0:
        ex_start, ex_end, gid, gname = intervals[i]
        # Since intervals are sorted by start, if ex_start is far back,
        # we can't stop immediately because of large exons, but we can limit search distance
        if start - ex_start > 2000000: # 2MB limit
            break
        if ex_end >= start and ex_start <= end:
            overlapping.add((gid, gname))
        i -= 1
        
    # Search forwards from the bisect index
    i = idx
    while i < len(intervals):
        ex_start, ex_end, gid, gname = intervals[i]
        if ex_start > end:
            break
        if ex_end >= start and ex_start <= end:
            overlapping.add((gid, gname))
        i += 1
        
    return overlapping

def audit_peaks(peak_path, gene_intervals, dmg_path, output_report_path=None, quiet=False):
    """Runs the audit comparing narrowPeak calls against DMG csv results."""
    if not quiet:
        print(f"\n[AUDIT] Auditing peaks: {peak_path}")
        print(f"        Against DMG results: {dmg_path}")
    
    if not os.path.exists(peak_path):
        if not quiet:
            print(f"  [ERROR] Peak file not found: {peak_path}")
        return None
    if not os.path.exists(dmg_path):
        if not quiet:
            print(f"  [ERROR] DMG file not found: {dmg_path}")
        return None
        
    # 1. Load DMGs
    df_dmg = pd.read_csv(dmg_path)
    # Check column names
    gene_col = 'gene_name' if 'gene_name' in df_dmg.columns else 'gene_id'
    dmgs = set(df_dmg[gene_col].dropna().unique())
    if not quiet:
        print(f"  Loaded {len(df_dmg)} significant DMGs from gene-level analysis.")
    
    # 2. Map peaks to genes
    peak_to_genes = defaultdict(list)
    gene_to_peaks = defaultdict(list)
    
    total_peaks = 0
    mapped_peaks = 0
    
    with open(peak_path, 'r') as f:
        # TRESS peaks.xls: tab-separated, 10 columns:
        # 0=chr, 1=start, 2=end, 3=strand, 4=summit, 5=input_counts,
        # 6=ip_counts, 7=pvals, 8=p.adj, 9=lg.fc
        # The first line is always a quoted header — skip it unconditionally.
        first_line = f.readline()
        is_header = (
            first_line.startswith('"chr"')
            or first_line.startswith('chr')
            or first_line.startswith('"seqnames"')
        )
        if not is_header:
            # Not a header line — rewind so this data row is not lost
            f.seek(0)

        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            total_peaks += 1
            chrom = fields[0].replace('"', '')
            start = int(fields[1])
            end = int(fields[2])
            peak_name = f"peak_{total_peaks}"
            # lg.fc is at index 9 in the TRESS peaks.xls format
            lg_fc = float(fields[9]) if len(fields) > 9 else 0.0
            fold_enrich = 2 ** abs(lg_fc)
            
            overlapping = find_overlapping_genes(chrom, start, end, gene_intervals)
            if overlapping:
                mapped_peaks += 1
                for gid, gname in overlapping:
                    # Match by gene_name or gene_id depending on how DMGs are listed
                    key = gname if gene_col == 'gene_name' else gid
                    peak_to_genes[peak_name].append(key)
                    gene_to_peaks[key].append((peak_name, fold_enrich))
                    
    if not quiet:
        print(f"  Processed {total_peaks} peaks. Mapped {mapped_peaks} ({mapped_peaks/total_peaks*100:.1f}%) to exonic gene intervals.")
    
    # 3. Calculate Sensitivity Check
    dmgs_with_peaks = dmgs.intersection(gene_to_peaks.keys())
    sensitivity = len(dmgs_with_peaks) / len(dmgs) * 100 if dmgs else 0.0
    
    # 4. Calculate Conflict/Cancellation Check
    genes_with_multiple_peaks = 0
    
    for g, peaks in gene_to_peaks.items():
        if len(peaks) > 1:
            genes_with_multiple_peaks += 1
            # Check if fold enrichment directions conflict
            # In MACS2, fold enrichment > 1 is typical. Direction conflicts occur in differential peaks.
            # But for baseline validation, we look at the spread of enrichment.
            enrichments = [p[1] for p in peaks]
            max_enrich = max(enrichments)
            min_enrich = min(enrichments)
            
    # Print Report
    report = []
    report.append("="*80)
    report.append(f"🧠 METHODOLOGICAL AUDIT (QC-03): GENE-LEVEL VS PEAK-LEVEL")
    report.append("="*80)
    report.append(f"  Peak file:            {os.path.basename(peak_path)}")
    report.append(f"  DMG file:             {os.path.basename(dmg_path)}")
    report.append(f"  Total Peaks:          {total_peaks}")
    report.append(f"  Exon-Mapped Peaks:    {mapped_peaks} ({mapped_peaks/total_peaks*100:.1f}%)")
    report.append(f"  Unique Genes w/ Peaks: {len(gene_to_peaks)}")
    report.append("-"*80)
    report.append(f"📌 SENSITIVITY CHECK:")
    report.append(f"  Total DMGs:           {len(dmgs)}")
    report.append(f"  DMGs with Peaks:      {len(dmgs_with_peaks)} ({sensitivity:.1f}%)")
    if sensitivity >= 85:
        report.append("  Result: PASS. Most gene-level calls have corresponding narrowPeak coordinates.")
    else:
        report.append("  Result: WARNING. Significant proportion of DMGs lack defined peak coordinates.")
    report.append("-"*80)
    report.append(f"📌 RESOLUTION DENSITY:")
    report.append(f"  Genes w/ Multi-peaks: {genes_with_multiple_peaks}")
    report.append(f"  Avg Peaks per Gene:   {mapped_peaks/len(gene_to_peaks) if gene_to_peaks else 0:.2f}")
    report.append("="*80)
    
    report_text = "\n".join(report)
    if not quiet:
        print(report_text)
    
    if output_report_path:
        with open(output_report_path, 'w') as f:
            f.write(report_text)
        if not quiet:
            print(f"[SUCCESS] Saved audit report to: {output_report_path}")
            
    return {
        'total_peaks': total_peaks,
        'mapped_peaks': mapped_peaks,
        'mapping_rate': mapped_peaks / total_peaks * 100 if total_peaks else 0.0,
        'unique_genes_w_peaks': len(gene_to_peaks),
        'total_dmgs': len(dmgs),
        'dmgs_with_peaks': len(dmgs_with_peaks),
        'sensitivity': sensitivity,
        'genes_w_multi_peaks': genes_with_multiple_peaks,
        'avg_peaks_per_gene': mapped_peaks / len(gene_to_peaks) if gene_to_peaks else 0.0
    }

def run_cohort_audit(peaks_dir, dmg_male, dmg_female, gene_intervals, output_base):
    print(f"\n[COHORT AUDIT] Starting cohort audit in: {peaks_dir}")
    print(f"               Male DMG: {dmg_male}")
    print(f"               Female DMG: {dmg_female}")
    
    # 1. Locate all peaks.xls files under peaks_dir
    peak_files = []
    for root, dirs, files in os.walk(peaks_dir):
        for f in files:
            if f.endswith('_peaks.xls'):
                peak_files.append(os.path.join(root, f))
                
    if not peak_files:
        print(f"[ERROR] No _peaks.xls files found in: {peaks_dir}")
        return
        
    print(f"Found {len(peak_files)} samples to audit.")
    
    results = []
    
    for idx, peak_path in enumerate(sorted(peak_files)):
        # Determine metadata from directory or file name
        sample_name = os.path.basename(os.path.dirname(peak_path))
        if not sample_name or sample_name == 'peaks':
            sample_name = os.path.basename(peak_path).replace('_peaks.xls', '')
            
        name_lower = sample_name.lower()
        # Check 'female' with startswith to prevent 'female' matching the
        # 'male' substring check when order is reversed.
        if name_lower.startswith('female') or '_female_' in name_lower:
            sex = 'female'
            dmg_path = dmg_female
        elif name_lower.startswith('male') or '_male_' in name_lower:
            sex = 'male'
            dmg_path = dmg_male
        else:
            print(f"[WARN] Could not determine sex for sample {sample_name}, defaulting to male DMG.")
            sex = 'unknown'
            dmg_path = dmg_male
            
        if 'mdd' in name_lower:
            condition = 'mdd'
        elif 'ctrl' in name_lower or 'control' in name_lower:
            condition = 'ctrl'
        else:
            condition = 'unknown'
            
        group = f"{sex.upper()} {condition.upper()}"
        
        print(f"  [{idx+1}/{len(peak_files)}] Auditing {sample_name} ({group}) ...")
        stats = audit_peaks(peak_path, gene_intervals, dmg_path, quiet=True)
        if stats:
            stats['sample_name'] = sample_name
            stats['sex'] = sex
            stats['condition'] = condition
            stats['group'] = group
            results.append(stats)
            
    if not results:
        print("[ERROR] No successful audit results generated.")
        return
        
    df = pd.DataFrame(results)
    
    # Reorder columns to put metadata first
    cols = ['sample_name', 'sex', 'condition', 'group'] + [c for c in df.columns if c not in ['sample_name', 'sex', 'condition', 'group']]
    df = df[cols]
    
    # Calculate group summaries
    summary_lines = []
    summary_lines.append("="*80)
    summary_lines.append("🧠 COHORT-WIDE METHODOLOGICAL AUDIT (QC-03) SUMMARY")
    summary_lines.append("="*80)
    summary_lines.append(f"Total Samples Processed: {len(df)}")
    summary_lines.append("-"*80)
    
    group_summary = df.groupby('group').agg({
        'total_peaks': 'mean',
        'mapping_rate': 'mean',
        'unique_genes_w_peaks': 'mean',
        'total_dmgs': 'first',
        'dmgs_with_peaks': 'mean',
        'sensitivity': 'mean',
        'avg_peaks_per_gene': 'mean'
    })
    
    summary_lines.append("Group Averages:")
    summary_lines.append(group_summary.to_string())
    summary_lines.append("="*80)
    
    summary_text = "\n".join(summary_lines)
    print(summary_text)
    
    # Write outputs
    if output_base:
        # Save summary report text
        report_path = output_base if output_base.endswith('.txt') else f"{output_base}_summary.txt"
        with open(report_path, 'w') as f:
            f.write(summary_text)
        print(f"[SUCCESS] Saved cohort summary report to: {report_path}")
        
        # Save detailed CSV
        csv_path = report_path.replace('.txt', '.csv')
        df.to_csv(csv_path, index=False)
        print(f"[SUCCESS] Saved detailed sample-level table to: {csv_path}")

def main():
    parser = argparse.ArgumentParser(description="Audit gene-level vs peak-level resolution.")
    parser.add_argument("--peaks", help="Path to narrowPeak file.")
    parser.add_argument("--peaks-dir", help="Path to directory containing all sample peak folders.")
    parser.add_argument("--dmg", help="Path to DESeq2 significant results CSV (required for single sample mode).")
    parser.add_argument("--dmg-male", help="Path to Male significant results CSV (used in cohort mode).")
    parser.add_argument("--dmg-female", help="Path to Female significant results CSV (used in cohort mode).")
    parser.add_argument("--gtf", required=True, help="Path to GTF genome annotation file.")
    parser.add_argument("--out", help="Path to write the report txt/csv.")
    
    args = parser.parse_args()
    
    # Run parsing
    gene_intervals, _ = parse_gtf(args.gtf)
    
    if args.peaks_dir:
        # Check required cohort mode args
        if not args.dmg_male or not args.dmg_female:
            print("[ERROR] Cohort mode (--peaks-dir) requires both --dmg-male and --dmg-female.")
            sys.exit(1)
        run_cohort_audit(args.peaks_dir, args.dmg_male, args.dmg_female, gene_intervals, args.out)
    elif args.peaks:
        if not args.dmg:
            print("[ERROR] Single-sample mode (--peaks) requires --dmg.")
            sys.exit(1)
        audit_peaks(args.peaks, gene_intervals, args.dmg, args.out)
    else:
        print("[ERROR] Either --peaks or --peaks-dir must be provided.")
        sys.exit(1)

if __name__ == "__main__":
    main()
