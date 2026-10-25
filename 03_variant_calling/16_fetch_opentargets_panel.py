#!/usr/bin/env python3
# [[03_variant_calling/16_fetch_opentargets_panel.py]]
# Goal: Fetch the primary gene panel -- Open Targets MDD, GENETIC evidence only
#       -- and cache it locally. NETWORK STEP ONLY; no statistics.
#
# Why genetic evidence only:
#   Filtering to `genetic_association` excludes literature and clinical
#   evidence, which would otherwise echo the same publications this project
#   cites. What remains is derived from large-scale GWAS and is independent of
#   these 39 donors, so using it is not circular and not hand-picked.
#   Source: Ochoa et al. 2023, NAR 51, doi:10.1093/nar/gkac1046
#
# Why cache:
#   m6AVar, RMVar and m6A-TSHub all refused connection from this host during
#   this project. A cached, versioned panel file means 17 never breaks because
#   a remote host is down. Re-running this script is opt-in via FORCE=1.
#
# Output:
#   16_opentargets_mdd_panel.tsv   gene, genetic_score, ensembl_id, retrieved
#
# Usage:
#   python3 03_variant_calling/16_fetch_opentargets_panel.py
#   FORCE=1 python3 03_variant_calling/16_fetch_opentargets_panel.py

from __future__ import annotations

import json
import os
import sys
import urllib.request
from datetime import date

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.environ.get("PROJECT_ROOT") or os.path.normpath(os.path.join(SCRIPT_DIR, ".."))
TABLE_DIR = os.path.join(PROJECT_ROOT, "results", "03_variant_calling", "tables")
OUT_PANEL = os.path.join(TABLE_DIR, "16_opentargets_mdd_panel.tsv")

API = "https://api.platform.opentargets.org/api/v4/graphql"
DISEASE = "MONDO_0002009"  # major depressive disorder
THRESHOLDS = (0.50, 0.35, 0.20, 0.10)

QUERY = """
{ disease(efoId:"%s") {
    id
    name
    associatedTargets(page:{index:0,size:2000}) {
      count
      rows { score datatypeScores { id score } target { id approvedSymbol } }
    }
  } }
""" % DISEASE


def log(msg: str) -> None:
    print(f"[16] {msg}", flush=True)


def fetch() -> list[dict]:
    req = urllib.request.Request(
        API,
        data=json.dumps({"query": QUERY}).encode(),
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=120) as resp:
        payload = json.loads(resp.read().decode())
    if payload.get("errors"):
        raise RuntimeError(f"Open Targets API error: {payload['errors']}")
    return payload["data"]["disease"]["associatedTargets"]["rows"]


def main() -> int:
    log("=" * 78)
    log("FETCH OPENTARGETS MDD PANEL (genetic evidence only)")
    log("=" * 78)

    if os.path.isfile(OUT_PANEL) and os.environ.get("FORCE", "") not in ("1", "true"):
        with open(OUT_PANEL) as fh:
            n = sum(1 for _ in fh) - 1
        log(f"cached panel already present: {OUT_PANEL} ({n:,} genes)")
        log("re-run with FORCE=1 to refresh")
        return 0

    os.makedirs(TABLE_DIR, exist_ok=True)
    log(f"querying {API} for {DISEASE} ...")
    rows = fetch()
    log(f"associated targets returned: {len(rows):,}")

    out = []
    for r in rows:
        scores = {d["id"]: d["score"] for d in r.get("datatypeScores", [])}
        genetic = float(scores.get("genetic_association", 0.0))
        if genetic <= 0:
            continue
        t = r["target"]
        out.append({
            "gene": t["approvedSymbol"],
            "ensembl_id": t["id"],
            "genetic_score": round(genetic, 6),
            "overall_score": round(float(r["score"]), 6),
        })
    out.sort(key=lambda d: -d["genetic_score"])

    retrieved = date.today().isoformat()
    with open(OUT_PANEL, "w") as fh:
        fh.write("gene\tensembl_id\tgenetic_score\toverall_score\tretrieved\n")
        for d in out:
            fh.write(f"{d['gene']}\t{d['ensembl_id']}\t{d['genetic_score']}\t"
                     f"{d['overall_score']}\t{retrieved}\n")

    log(f"genes with genetic evidence > 0: {len(out):,}")
    for thr in THRESHOLDS:
        log(f"  score >= {thr:.2f}: {sum(1 for d in out if d['genetic_score'] >= thr):>4,} genes")
    log(f"top 10: " + ", ".join(f"{d['gene']}({d['genetic_score']:.3f})" for d in out[:10]))
    log(f"panel -> {OUT_PANEL}")
    return 0


if __name__ == "__main__":
    sys.exit(main())