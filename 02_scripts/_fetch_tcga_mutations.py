# -*- coding: utf-8 -*-
"""_fetch_tcga_mutations.py — cBioPortal brca_tcga 突变数据抓取。

产出 01_data/raw/brca_tcga/brca_tcga_mutations.csv:
  sampleId, mutation_count (全外显子组非同义计数, 来自 clinical MUTATION_COUNT),
  TP53, BRCA1, BRCA2, PIK3CA (0/1, 任一非 Silent 突变)

API:
  GET  /studies/brca_tcga/clinical-data?attributeId=MUTATION_COUNT  (sample)
  POST /molecular-profiles/brca_tcga_mutations/mutations/fetch      (分页)
"""
from __future__ import annotations

import json
import ssl
import time
import urllib.request
from pathlib import Path

import pandas as pd

OUT = Path(r"D:/嗜神经浸润/01_data/raw/brca_tcga/brca_tcga_mutations.csv")
API = "https://www.cbioportal.org/api"
CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE
FOCUS = {"TP53", "BRCA1", "BRCA2", "PIK3CA"}


def req(url: str, data: bytes | None = None, retries: int = 4) -> bytes:
    for i in range(retries):
        try:
            r = urllib.request.Request(url, data=data, method="POST" if data else "GET")
            r.add_header("accept", "application/json")
            if data:
                r.add_header("content-type", "application/json")
            with urllib.request.urlopen(r, timeout=180, context=CTX) as resp:
                return resp.read()
        except Exception as ex:
            print(f"  retry {i+1}: {ex}")
            time.sleep(5)
    raise RuntimeError(f"API 失败: {url}")


def main() -> None:
    # ---- 1. MUTATION_COUNT (sample 级)
    raw = req(f"{API}/studies/brca_tcga/clinical-data?attributeId=MUTATION_COUNT"
              "&clinicalDataType=SAMPLE")
    d = json.loads(raw)
    mcount = {r["sampleId"]: int(r["value"]) for r in d if r.get("value", "").isdigit()}
    print(f"[MUTATION_COUNT] n={len(mcount)}, median={pd.Series(list(mcount.values())).median():.0f}")

    # ---- 2. mutations/fetch 分页
    body = json.dumps({"sampleListId": "brca_tcga_all"}).encode()
    records: list[dict] = []
    page = 0
    while True:
        url = (f"{API}/molecular-profiles/brca_tcga_mutations/mutations/fetch"
               f"?direction=ASC&pageNumber={page}&pageSize=10000&projection=DETAILED")
        raw = req(url, data=body)
        chunk = json.loads(raw)
        if not chunk:
            break
        records.extend(chunk)
        print(f"  page {page}: 累计 {len(records)} 条")
        page += 1
        if page > 40:
            break

    rows = {}
    nonsilent_total: dict[str, int] = {}
    for m in records:
        gene = (m.get("gene") or {}).get("hugoGeneSymbol", "")
        smp = m.get("sampleId", "")
        mtype = m.get("mutationType", "")
        if mtype in ("Silent", "Intron", "3'UTR", "5'UTR", "3'Flank", "5'Flank",
                     "RNA", "IGR", "LINCRNA"):
            continue
        nonsilent_total[smp] = nonsilent_total.get(smp, 0) + 1
        if gene in FOCUS:
            rows.setdefault(smp, {})[gene] = 1

    out = []
    for s, c in mcount.items():
        r = rows.get(s, {})
        out.append({"sampleId": s, "mutation_count": c,
                    "n_nonsilent_fetch": nonsilent_total.get(s, 0),
                    "TP53": r.get("TP53", 0), "BRCA1": r.get("BRCA1", 0),
                    "BRCA2": r.get("BRCA2", 0), "PIK3CA": r.get("PIK3CA", 0)})
    df = pd.DataFrame(out)
    df.to_csv(OUT, index=False)
    print(f"\n[输出] {OUT}  n={len(df)}")
    print(df[["TP53", "BRCA1", "BRCA2", "PIK3CA"]].sum())
    j = df[df.n_nonsilent_fetch > 0]
    corr = j[["mutation_count", "n_nonsilent_fetch"]].corr().iloc[0, 1]
    print(f"[交叉验证] MUTATION_COUNT vs fetch 非同义计数 ρ={corr:.3f} (n={len(j)})")


if __name__ == "__main__":
    main()
