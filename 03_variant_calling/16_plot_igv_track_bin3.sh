#!/usr/bin/env bash
# [[03_variant_calling/16_plot_igv_track_bin3.sh]]
# Goal: Render pyGenomeTracks multi-track figure for BIN3 locus.

export PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

mkdir -p "$PROJECT_ROOT/results/03_variant_calling/figures"

cd "$PROJECT_ROOT"

echo ">>> Rendering pyGenomeTracks figure for BIN3 (chr8:22,619,800-22,621,500)..."

pyGenomeTracks \
  --tracks 03_variant_calling/16_bin3_tracks.ini \
  --region 8:22619800-22621500 \
  --outFileName results/03_variant_calling/figures/16_igv_track_bin3.png \
  --dpi 300

pyGenomeTracks \
  --tracks 03_variant_calling/16_bin3_tracks.ini \
  --region 8:22619800-22621500 \
  --outFileName results/03_variant_calling/figures/16_igv_track_bin3.pdf

echo ">>> pyGenomeTracks rendering complete for BIN3."
