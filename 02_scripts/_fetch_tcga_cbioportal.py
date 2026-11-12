# -*- coding: utf-8 -*-
"""抓取 TCGA-BRCA 目标基因表达 + 临床（cBioPortal per-gene API + Xena 临床表）。

为什么不用整矩阵: Xena HiSeqV2 64 MB 单连接限速，且我们只需要 ~190 个基因。
cBioPortal 的 per-gene 接口按 Entrez ID 取，直接返回 symbol 级矩阵。

输出:
    01_data/raw/tcga_brca/tcga_brca_target_genes.csv   sample × gene (RSEM)
    01_data/raw/tcga_brca/tcga_brca_clin.csv           sample × 临床字段
"""
from __future__ import annotations

import json
import sys
import time
import urllib.parse
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
STUDY = "brca_tcga"
PROFILE = "brca_tcga_rna_seq_v2_mrna"
SAMPLE_LIST = "brca_tcga_rna_seq_v2_mrna"
UA = {"User-Agent": "Mozilla/5.0"}
RAW = BL.RAW / "tcga_brca"


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


def fetch_gene(sym: str, eid: int):
    url = (f"{CBIO}/molecular-profiles/{PROFILE}/molecular-data"
           f"?sampleListId={SAMPLE_LIST}&entrezGeneId={eid}&pageSize=100000")
    for attempt in range(1, 4):
        try:
            recs = get(url, retries=1)
            return sym, {r["sampleId"]: r.get("value") for r in recs}
        except Exception:
            if attempt == 3:
                return sym, {}
            time.sleep(1.5 * attempt)
    return sym, {}


# 下游链（纤维化程序 / 纯度代理 / 神经-免疫轴）额外需要的基因。
# 放在这里而不是 bulk_lib，是因为它们不参与任何模块评分，只做协变量与程序定量。
EXTRA_GENES = [
    # 上皮/肿瘤纯度代理
    "EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1", "MUC1", "KRT5", "KRT14",
    # 胶原成熟与交联（纤维化"质"而非"量"）
    "COL11A1", "COL12A1", "COL5A2", "LOXL1", "LOXL2", "PLOD1", "PLOD2",
    "P4HA1", "P4HA2", "SERPINH1", "ADAMTS2", "BMP1",
    # 基质刚度 / 机械转导
    "THBS1", "TNC", "CTHRC1", "MMP2", "MMP11", "MMP14", "TIMP1", "YAP1", "WWTR1",
    # EMT / TGF-β
    "VIM", "SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2", "SMAD3", "SKIL", "TGFB3", "TGFBR3",
    # 神经（受体/施万/嗜铬粒）
    "NGF", "NTF3", "NGFR", "S100B", "PLP1", "MPZ", "SOX10", "SCN7A",
    # 免疫区室（用于纯度校正后的独立读数）
    "CD3D", "CD3E", "CD4", "CD68", "KLRD1", "HLA-A", "B2M", "PDCD1", "CD274", "LAG3",
    # 增殖
    "MKI67", "TOP2A", "AURKA",
]


def main():
    mods = BL.load_modules()
    genes = set(EXTRA_GENES) | set(BL.PURITY_PROXY)
    for v in mods.values():
        genes |= set(v)
    for v in BL.PROGRAMS.values():
        genes |= set(v)
    genes = sorted(genes)
    print(f"目标基因 {len(genes)} 个")

    path = RAW / "tcga_brca_target_genes.csv"
    df = None
    if path.exists():
        df = pd.read_csv(path, index_col=0)
        miss = [g for g in genes if g not in df.columns]
        if not miss:
            print(f"cache hit {df.shape}")
        else:
            # 增量补抓：只请求缺失基因，与已有矩阵按 sample 对齐后拼接
            print(f"cache {df.shape} 缺 {len(miss)} 基因 → 增量抓取")
            id_map, unmatched = symbols_to_entrez(miss)
            print(f"Entrez 映射 {len(id_map)}/{len(miss)}; 未匹配 {unmatched[:10]}")
            new = {}
            t0 = time.time()
            with ThreadPoolExecutor(max_workers=10) as ex:
                futs = [ex.submit(fetch_gene, s, e) for s, e in sorted(id_map.items())]
                for k, f in enumerate(futs, 1):
                    sym, d = f.result()
                    if d:
                        new[sym] = d
                    if k % 30 == 0:
                        print(f"  {k}/{len(futs)}  {time.time()-t0:.0f}s", flush=True)
            add = pd.DataFrame(new).reindex(df.index)
            got = [c for c in add.columns if add[c].notna().any()]
            df = pd.concat([df, add[got]], axis=1)
            df.index.name = "sample"
            df.to_csv(path)
            print(f"合并后 {df.shape}（新增 {len(got)} 基因）-> {path}")

    if df is None:
        id_map, unmatched = symbols_to_entrez(genes)
        print(f"Entrez 映射 {len(id_map)}/{len(genes)}; 未匹配 {unmatched[:10]}")
        out = {}
        t0 = time.time()
        with ThreadPoolExecutor(max_workers=10) as ex:
            futs = [ex.submit(fetch_gene, s, e) for s, e in sorted(id_map.items())]
            for k, f in enumerate(futs, 1):
                sym, d = f.result()
                if d:
                    out[sym] = d
                if k % 30 == 0:
                    print(f"  {k}/{len(futs)}  {time.time()-t0:.0f}s", flush=True)
        df = pd.DataFrame(out)
        df.index.name = "sample"
        df.to_csv(path)
        print(f"表达矩阵 {df.shape} -> {path}")

    # ---- 临床（Xena）
    clin = pd.read_csv(RAW / "BRCA_clinicalMatrix", sep="\t", index_col=0, low_memory=False)
    keep = ["ER_Status_nature2012", "PR_Status_nature2012", "HER2_Final_Status_nature2012",
            "PAM50Call_RNAseq", "Age_at_Initial_Pathologic_Diagnosis_nature2012",
            "Node_nature2012", "AJCC_Stage_nature2012", "OS_Time_nature2012", "OS_event_nature2012"]
    keep = [c for c in keep if c in clin.columns]
    clin = clin[keep].copy()
    clin.columns = ["ER", "PR", "HER2", "PAM50", "AGE", "NODE", "STAGE", "OS_TIME", "OS_EVENT"]
    clin = clin.reindex(df.index)
    clin.to_csv(RAW / "tcga_brca_clin.csv")
    print(f"临床 {clin.shape}; ER 分布 {clin['ER'].value_counts(dropna=False).to_dict()}")
    neg = lambda s: s.astype(str).str.lower().str.startswith(("neg", "n"))
    tn = neg(clin["ER"]) & neg(clin["PR"]) & neg(clin["HER2"])
    print(f"TNBC n={int(tn.sum())} / {clin['ER'].notna().sum()} 有 ER 记录")
    print(f"PAM50: {clin['PAM50'].value_counts(dropna=False).to_dict()}")


if __name__ == "__main__":
    main()
