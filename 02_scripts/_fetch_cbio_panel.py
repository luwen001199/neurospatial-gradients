# -*- coding: utf-8 -*-
"""通用 cBioPortal 面板抓取：任选 study，抓「模块/程序/反卷积」基因 + 生存临床。

用法:
    python 02_scripts/_fetch_cbio_panel.py brca_tcga
    python 02_scripts/_fetch_cbio_panel.py paad_tcga hnsc_tcga

输出:
    01_data/raw/<study>/<study>_panel.csv       sample × gene (linear, RSEM/TPM)
    01_data/raw/<study>/<study>_clin.csv        sample × 临床字段
"""
from __future__ import annotations

import json
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, r"D:/嗜神经浸润/02_scripts")
sys.path.insert(0, r"D:/嗜神经浸润/00_pipeline")
import bulk_lib as BL  # noqa: E402
from nsgrad.bulk import symbols_to_entrez  # noqa: E402

CBIO = "https://www.cbioportal.org/api"
UA = {"User-Agent": "Mozilla/5.0"}

# 各 study 的表达 profile 与样本列表（TCGA PanCan 与 legacy 命名不同）
STUDY_CFG = {
    "brca_tcga": ("brca_tcga_rna_seq_v2_mrna", "brca_tcga_rna_seq_v2_mrna"),
    "paad_tcga": ("paad_tcga_rna_seq_v2_mrna", "paad_tcga_rna_seq_v2_mrna"),
    "hnsc_tcga": ("hnsc_tcga_rna_seq_v2_mrna", "hnsc_tcga_rna_seq_v2_mrna"),
}

# 反卷积面板额外需要的区室标记（不在模块/PROGRAMS 里的）
DECONV_EXTRA = sorted({g for v in BL.DECONV_PANEL.values() for g in v})

CLIN_ATTRS = ["OS_MONTHS", "OS_STATUS", "AGE", "AJCC_PATHOLOGIC_TUMOR_STAGE",
              "AJCC_NODES_PATHOLOGIC_PN", "GRADE"]


def get(url: str, retries=3, timeout=180):
    last = None
    for i in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers=UA)
            return json.loads(urllib.request.urlopen(req, timeout=timeout).read())
        except Exception as e:
            last = e
            time.sleep(1.5 * i)
    raise last


def fetch_gene(profile, sample_list, sym, eid):
    url = (f"{CBIO}/molecular-profiles/{profile}/molecular-data"
           f"?sampleListId={sample_list}&entrezGeneId={eid}&pageSize=100000")
    for attempt in range(1, 4):
        try:
            recs = get(url, retries=1)
            return sym, {r["sampleId"]: r.get("value") for r in recs}
        except Exception:
            if attempt == 3:
                return sym, {}
            time.sleep(1.5 * attempt)
    return sym, {}


def fetch_clin(study, attr):
    url = (f"{CBIO}/studies/{study}/clinical-data"
           f"?clinicalDataType=PATIENT&attributeId={attr}&pageSize=100000")
    try:
        recs = get(url, retries=2)
    except Exception:
        return attr, {}
    return attr, {r["patientId"]: r.get("value") for r in recs}


def run(study: str):
    if study not in STUDY_CFG:
        raise SystemExit(f"未配置 {study}")
    profile, sample_list = STUDY_CFG[study]
    outdir = BL.RAW / study
    outdir.mkdir(parents=True, exist_ok=True)

    mods = BL.load_modules()
    genes = set(DECONV_EXTRA) | set(BL.PURITY_PROXY) | set(CLIN_ATTRS[:0])
    for v in mods.values():
        genes |= set(v)
    for v in BL.PROGRAMS.values():
        genes |= set(v)
    # 纤维化"质"、EMT、神经、增殖（与 BRCA 那一轮保持一致）
    genes |= {"COL11A1", "COL12A1", "COL5A2", "LOXL1", "LOXL2", "PLOD1", "PLOD2",
              "P4HA1", "P4HA2", "SERPINH1", "ADAMTS2", "BMP1", "MMP2", "MMP11",
              "MMP14", "TIMP1", "YAP1", "WWTR1", "VIM", "SNAI1", "SNAI2", "TWIST1",
              "ZEB1", "ZEB2", "SMAD3", "SKIL", "TGFB3", "TGFBR3", "MKI67",
              "NGF", "NTF3", "NGFR", "S100B", "PLP1", "MPZ", "SOX10", "SCN7A",
              "KRT5", "KRT14", "MUC1", "CD68", "HLA-A", "B2M", "PDCD1", "CD274", "LAG3"}
    genes = sorted(genes)
    print(f"[{study}] 目标基因 {len(genes)}")

    path = outdir / f"{study}_panel.csv"
    df = None
    if path.exists():
        df = pd.read_csv(path, index_col=0)
        miss = [g for g in genes if g not in df.columns]
        print(f"cache {df.shape}, 缺 {len(miss)}")
    else:
        miss = genes

    if miss:
        id_map, unmatched = symbols_to_entrez(miss)
        print(f"Entrez 映射 {len(id_map)}/{len(miss)}; 未匹配 {unmatched[:12]}")
        new, t0 = {}, time.time()
        with ThreadPoolExecutor(max_workers=10) as ex:
            futs = [ex.submit(fetch_gene, profile, sample_list, s, e)
                    for s, e in sorted(id_map.items())]
            for k, f in enumerate(futs, 1):
                sym, d = f.result()
                if d:
                    new[sym] = d
                if k % 40 == 0:
                    print(f"  {k}/{len(futs)}  {time.time()-t0:.0f}s", flush=True)
        add = pd.DataFrame(new)
        if df is not None:
            add = add.reindex(df.index)
            df = pd.concat([df, add[[c for c in add.columns if add[c].notna().any()]]], axis=1)
        else:
            df = add
        df.index.name = "sample"
        df.to_csv(path)
        print(f"矩阵 {df.shape} -> {path}")

    # ---- 临床
    clin_path = outdir / f"{study}_clin.csv"
    cols = {}
    with ThreadPoolExecutor(max_workers=6) as ex:
        for attr, d in ex.map(lambda a: fetch_clin(study, a), CLIN_ATTRS):
            if d:
                cols[attr] = d
    clin = pd.DataFrame(cols)
    # patient -> sample：TCGA barcode 前 12 位；同一病人多次取样时只保留第一个
    clin.index = clin.index.astype(str)
    samples = pd.Series(df.index.astype(str))
    pid = samples.str[:12]
    first = samples[~pid.duplicated()].to_numpy()
    pid2sample = pd.Series(first, index=pd.Series(first).str[:12])
    mapped = pid2sample.reindex(clin.index)
    clin = clin.loc[mapped.notna().to_numpy()]
    clin.index = mapped.dropna().to_numpy()
    clin = clin.loc[~clin.index.duplicated()]
    df = df.loc[clin.index]
    df.to_csv(path)
    clin.to_csv(clin_path)
    os_t = pd.to_numeric(clin.get("OS_MONTHS"), errors="coerce")
    ev = clin.get("OS_STATUS", pd.Series(index=clin.index, dtype=object)).astype(str)
    print(f"[{study}] 临床 {clin.shape}; OS 可算 {int(os_t.notna().sum())}; "
          f"事件 {int(ev.str.startswith('1').sum())}")
    return df, clin


if __name__ == "__main__":
    targets = sys.argv[1:] or ["brca_tcga"]
    for s in targets:
        run(s)
