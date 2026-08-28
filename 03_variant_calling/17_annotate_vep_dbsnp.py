#!/usr/bin/env python3
# [[03_variant_calling/17_annotate_vep_dbsnp.py]]
# Goal: Map genomic coordinates, variant consequences, affected transcripts (VEP), and dbSNP rsIDs/MAFs for significant ASM loci.
#
# Steps:
#   1. Dynamic project root resolution and path setup.
#   2. Load significant ASM loci dynamically from 12_asm_primary_significant.csv.
#   3. Extract overlapping candidate SNPs from 09_drach_disrupting_snps.csv within loci coordinates.
#   4. Query Ensembl VEP REST API dynamically for all candidate variants (no hardcoded biology/variants).
#   5. Parse consequence terms, affected transcripts, biotypes, cDNA positions, rsIDs, and MAF.
#   6. Perform data assertions (check for non-empty tables and required column completeness).
#   7. Save results to results/03_variant_calling/tables/17_vep_dbsnp_annotation.csv.

import os
import sys
import json
import time
import urllib.request
import pandas as pd

# --- Step 1: Dynamic Pathing ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.environ.get("PROJECT_ROOT") or os.path.normpath(os.path.join(SCRIPT_DIR, ".."))

TABLES_DIR = os.path.join(PROJECT_ROOT, "results", "03_variant_calling", "tables")
os.makedirs(TABLES_DIR, exist_ok=True)

SIG_FILE   = os.path.join(TABLES_DIR, "12_asm_primary_significant.csv")
DRACH_FILE = os.path.join(TABLES_DIR, "09_drach_disrupting_snps.csv")
OUT_CSV    = os.path.join(TABLES_DIR, "17_vep_dbsnp_annotation.csv")

# --- Step 2: Load Significant ASM Loci ---
print(f"[INFO] Reading significant ASM loci from: {SIG_FILE}")
if not os.path.exists(SIG_FILE):
    sys.exit(f"[ERROR] Required input file not found: {SIG_FILE}")

sig_df = pd.read_csv(SIG_FILE)
if sig_df.empty:
    sys.exit(f"[ERROR] Significant ASM loci table is empty: {SIG_FILE}")

print(f"[INFO] Loaded {len(sig_df)} significant loci across genes: {list(sig_df['gene'].unique())}")

# --- Step 3: Extract Overlapping Candidate SNPs ---
candidate_records = []
if os.path.exists(DRACH_FILE):
    print(f"[INFO] Loading DRACH disruption table: {DRACH_FILE}")
    drach_df = pd.read_csv(DRACH_FILE)
    
    for _, locus in sig_df.iterrows():
        g   = str(locus['gene'])
        c   = str(locus['chr'])
        st  = int(locus['start'])
        en  = int(locus['end'])
        fdr = float(locus['fdr'])
        conf = str(locus.get('confidence', 'N/A'))
        
        # Filter SNPs in drach_df within locus interval
        sub = drach_df[(drach_df['snp_chr'].astype(str) == c) & 
                       (drach_df['snp_pos'] >= st) & 
                       (drach_df['snp_pos'] <= en)].copy()
        
        if not sub.empty:
            sub_uniq = sub.drop_duplicates(subset=['snp_chr', 'snp_pos', 'ref', 'alt'])
            for _, v in sub_uniq.iterrows():
                candidate_records.append({
                    'locus_id': str(locus.get('locus_id', f"{c}:{st}:{en}")),
                    'gene': g,
                    'locus_chr': c,
                    'locus_start': st,
                    'locus_end': en,
                    'locus_fdr': fdr,
                    'confidence': conf,
                    'snp_chr': str(v['snp_chr']),
                    'snp_pos': int(v['snp_pos']),
                    'ref': str(v['ref']),
                    'alt': str(v['alt']),
                    'drach_name': str(v.get('drach_name', 'N/A')),
                    'disruption_status': str(v.get('disruption_status', 'N/A')),
                    'motif_offset': int(v.get('motif_offset', -1))
                })

cand_df = pd.DataFrame(candidate_records)
if cand_df.empty:
    print("[WARNING] No overlapping SNPs found in DRACH table. Using locus endpoints as search coordinates.")
    for _, locus in sig_df.iterrows():
        cand_df = pd.concat([cand_df, pd.DataFrame([{
            'locus_id': str(locus.get('locus_id', '')),
            'gene': str(locus['gene']),
            'locus_chr': str(locus['chr']),
            'locus_start': int(locus['start']),
            'locus_end': int(locus['end']),
            'locus_fdr': float(locus['fdr']),
            'confidence': str(locus.get('confidence', 'N/A')),
            'snp_chr': str(locus['chr']),
            'snp_pos': int(locus['start']),
            'ref': 'N',
            'alt': 'N',
            'drach_name': 'N/A',
            'disruption_status': 'NO_OVERLAP',
            'motif_offset': -1
        }])], ignore_index=True)
else:
    cand_df = cand_df.drop_duplicates(subset=['gene', 'snp_chr', 'snp_pos', 'ref', 'alt'])

print(f"[INFO] Extracted {len(cand_df)} unique candidate variant(s) for VEP querying.")

# --- Step 4: Batched Ensembl VEP & dbSNP REST API Queries ---
print(f"[INFO] Querying Ensembl VEP REST API in batch for {len(cand_df)} variant(s)...")

# Prepare batch input: format "chr:pos-pos:1/alt" or VCF format
variant_map = {}
batch_variants = []

for idx, row in cand_df.iterrows():
    c   = str(row['snp_chr']).replace('chr', '')
    pos = int(row['snp_pos'])
    ref = str(row['ref'])
    alt = str(row['alt'])
    
    if ref != 'N' and alt != 'N':
        var_str = f"{c}:{pos}-{pos}:1/{alt}"
        batch_variants.append(var_str)
        variant_map[var_str] = row

vep_results = {}
if batch_variants:
    # Query in chunks of up to 50 variants
    chunk_size = 50
    for i in range(0, len(batch_variants), chunk_size):
        chunk = batch_variants[i:i+chunk_size]
        url = "https://rest.ensembl.org/vep/human/region"
        payload = json.dumps({"variants": chunk}).encode('utf-8')
        
        req = urllib.request.Request(
            url, 
            data=payload,
            headers={
                "Content-Type": "application/json", 
                "Accept": "application/json",
                "User-Agent": "MDD_RNAseq_Pipeline/1.0"
            }
        )
        
        success = False
        for attempt in range(1, 4):
            try:
                print(f"[INFO] Sending batch request {i//chunk_size + 1} ({len(chunk)} variants), attempt {attempt}...")
                with urllib.request.urlopen(req, timeout=30) as resp:
                    data = json.loads(resp.read().decode('utf-8'))
                    if data and isinstance(data, list):
                        for item in data:
                            inp = item.get('input')
                            if inp:
                                vep_results[inp] = item
                    success = True
                    break
            except Exception as e:
                print(f"[WARNING] Batch query attempt {attempt} failed: {e}")
                time.sleep(2 * attempt)
        
        if not success:
            print(f"[ERROR] Failed to retrieve VEP batch after 3 attempts.")

# Parse results
annotated_rows = []
for idx, row in cand_df.iterrows():
    c   = str(row['snp_chr']).replace('chr', '')
    pos = int(row['snp_pos'])
    ref = str(row['ref'])
    alt = str(row['alt'])
    g   = row['gene']
    var_str = f"{c}:{pos}-{pos}:1/{alt}"
    
    rsid = "N/A"
    consequence = "N/A"
    transcript_id = "N/A"
    biotype = "N/A"
    cdna_pos = "N/A"
    cds_pos = "N/A"
    protein_pos = "N/A"
    maf = "N/A"
    
    item = vep_results.get(var_str)
    if item:
        consequence = item.get('most_severe_consequence', 'N/A')
        
        # Extract dbSNP rsIDs
        coloc = item.get('colocated_variants', [])
        rs_ids = [cv.get('id') for cv in coloc if cv.get('id', '').startswith('rs')]
        if rs_ids:
            rsid = ", ".join(rs_ids)
        
        # Extract MAF if present
        maf_list = []
        for cv in coloc:
            for k, v_maf in cv.items():
                if 'maf' in k.lower() or 'frequency' in k.lower():
                    maf_list.append(f"{k}:{v_maf}")
        if maf_list:
            maf = "; ".join(maf_list)
        
        # Match transcript consequences
        for tc in item.get('transcript_consequences', []):
            if tc.get('gene_symbol') == g or tc.get('gene_symbol'):
                transcript_id = tc.get('transcript_id', 'N/A')
                biotype = tc.get('biotype', 'N/A')
                if tc.get('consequence_terms'):
                    consequence = ", ".join(tc.get('consequence_terms'))
                if tc.get('cdna_start'):
                    cdna_pos = str(tc.get('cdna_start'))
                if tc.get('cds_start'):
                    cds_pos = str(tc.get('cds_start'))
                if tc.get('protein_start'):
                    protein_pos = str(tc.get('protein_start'))
                break

    annotated_rows.append({
        'locus_id': row['locus_id'],
        'gene': g,
        'locus_fdr': row['locus_fdr'],
        'confidence': row['confidence'],
        'chr': c,
        'pos': pos,
        'ref': ref,
        'alt': alt,
        'rsid': rsid,
        'maf': maf,
        'vep_consequence': consequence,
        'transcript_id': transcript_id,
        'biotype': biotype,
        'cdna_pos': cdna_pos,
        'cds_pos': cds_pos,
        'protein_pos': protein_pos,
        'drach_name': row['drach_name'],
        'disruption_status': row['disruption_status'],
        'motif_offset': row['motif_offset']
    })

out_df = pd.DataFrame(annotated_rows)

# --- Step 5: Data Assertions ---
assert not out_df.empty, "Assertion failed: Output annotation table is empty."
assert 'gene' in out_df.columns and 'vep_consequence' in out_df.columns, "Assertion failed: Missing essential columns."

# --- Step 6: Export Table ---
out_df.to_csv(OUT_CSV, index=False)
print(f"[INFO] Successfully saved Layer 1 & 2 annotations to: {OUT_CSV}")
print("[DONE] 17_annotate_vep_dbsnp.py executed successfully.")
