#!/usr/bin/env python3
# [[utils/ame_plot_enrichment.py]]
# Goal: Generate an enrichment figure from AME --text TSV output.
# Usage: python3 utils/ame_plot_enrichment.py <input_tsv> <output_png>

import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} <input_tsv> <output_png>", file=sys.stderr)
    sys.exit(1)

in_tsv = sys.argv[1]
out_png = sys.argv[2]

if not pd.io.common.file_exists(in_tsv):
    print(f"Assertion Error: Input TSV not found: {in_tsv}", file=sys.stderr)
    sys.exit(1)

df = pd.read_csv(in_tsv, sep='\t', comment='#')
if df.empty:
    print("Assertion Error: Input TSV is empty.", file=sys.stderr)
    sys.exit(1)

motif = df['motif_ID'].values[0]
tp    = int(df['TP'].values[0])
fp    = int(df['FP'].values[0])
n     = int(df['pos'].values[0])
pval  = df['p-value'].values[0]
odds_ratio = (tp / (n - tp)) / (fp / (n - fp))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))

x = np.arange(2)
width = 0.35
bars1 = ax1.bar(x - width/2, [tp, n - tp], width, label='Peak sequences',
                 color='#2ecc71', edgecolor='white')
bars2 = ax1.bar(x + width/2, [fp, n - fp], width, label='Shuffled (control)',
                 color='#e74c3c', edgecolor='white')
for bars in [bars1, bars2]:
    for bar in bars:
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 500,
                f'{int(bar.get_height()):,}', ha='center', va='bottom', fontsize=9)
ax1.set_xticks(x)
ax1.set_xticklabels(['Has GGACH motif', 'No GGACH motif'])
ax1.set_ylabel('Sequence count')
ax1.legend(fontsize=8)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

ax2.barh([''], [odds_ratio], color='#3498db', height=0.4, edgecolor='white')
ax2.axvline(1.0, color='gray', linestyle='--', linewidth=1.2)
ax2.text(odds_ratio + 0.02, 0, f'OR = {odds_ratio:.2f}  (p = {pval:.2e})',
         va='center', fontsize=11, fontweight='bold')
ax2.set_xlim(0, max(odds_ratio * 1.7, 1.7))
ax2.set_ylim(-0.5, 0.5)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.tick_params(left=False, labelleft=False)
ax2.set_xlabel('Odds ratio')
ax2.text(1.0, 0.38, 'OR = 1 (no enrichment)', ha='center', va='bottom',
         fontsize=7.5, color='gray', style='italic')

fig.suptitle(f'{motif} Motif Enrichment in TRESS Consensus Peaks', fontsize=13, y=1.02)
plt.tight_layout()
plt.savefig(out_png, dpi=150, bbox_inches='tight')
print(f'    Figure saved: {out_png}')
