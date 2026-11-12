# -*- coding: utf-8 -*-
"""bulk 队列工具: (1) 从 cBioPortal 取 METABRIC 指定基因表达; (2) 空间签名打分。

复用 fetch_metabric_mrna.py 踩过的坑:
  * symbol->Entrez 必须走 NCBI E-utilities (cBioPortal /api/genes 返回的 ID 不可靠)
  * /molecular-data 的 entrezGeneId 是单值, 一次只能查一个基因 -> 并发拉
"""
from __future__ import annotations
import json
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

CBIO = "https://www.cbioportal.org/api"
PROFILE = "brca_metabric_mrna"
SAMPLE_LIST = "brca_metabric_mrna"
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)",
      "Accept": "application/json"}
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _get(url, timeout=120, retries=3):
    last = None
    for i in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers=UA)
            return json.loads(urllib.request.urlopen(req, timeout=timeout).read())
        except Exception as e:
            last = e
            time.sleep(1.5 * i)
    raise last


def symbols_to_entrez(symbols, batch=40):
    """NCBI E-utilities 批量 symbol -> Entrez GeneID(严格校验回传 name)。"""
    symbols = [s for s in dict.fromkeys(symbols)]
    mapping, unmatched = {}, []
    for s in range(0, len(symbols), batch):
        chunk = symbols[s:s + batch]
        term = "(" + " OR ".join(f"{x}[sym]" for x in chunk) + ") AND 9606[taxid]"
        url = (f"{EUTILS}/esearch.fcgi?db=gene&retmax=200&retmode=json"
               f"&term={urllib.parse.quote(term)}")
        try:
            ids = _get(url).get("esearchresult", {}).get("idlist", [])
        except Exception:
            unmatched.extend(chunk)
            continue
        if not ids:
            unmatched.extend(chunk)
            continue
        try:
            res = _get(f"{EUTILS}/esummary.fcgi?db=gene&retmode=json&id={','.join(ids)}").get("result", {})
        except Exception:
            unmatched.extend(chunk)
            continue
        got = {}
        for gid in ids:
            rec = res.get(gid) or {}
            name = (rec.get("name") or "").upper()
            if name in {c.upper() for c in chunk}:
                got[name] = int(gid)
        for c in chunk:
            if c.upper() in got:
                mapping[c] = got[c.upper()]
            else:
                unmatched.append(c)
        time.sleep(0.35)
    return mapping, sorted(set(unmatched))


def fetch_metabric_genes(symbols, workers=8, cache: Path | None = None, verbose=True):
    """返回 DataFrame: sampleId x symbol (mRNA, Illumina HT-12)。"""
    if cache is not None and Path(cache).exists():
        df = pd.read_csv(cache, index_col=0)
        missing = [g for g in symbols if g not in df.columns]
        if not missing:
            if verbose:
                print(f"    cache hit: {cache.name} {df.shape}")
            return df[list(dict.fromkeys(symbols))]
    id_map, unmatched = symbols_to_entrez(symbols)
    if verbose:
        print(f"    mapped {len(id_map)}/{len(set(symbols))} genes; "
              f"unmatched={unmatched[:8]}{'...' if len(unmatched) > 8 else ''}")

    def one(sym, eid):
        url = (f"{CBIO}/molecular-profiles/{PROFILE}/molecular-data"
               f"?sampleListId={SAMPLE_LIST}&entrezGeneId={eid}&pageSize=100000")
        for attempt in range(1, 4):
            try:
                req = urllib.request.Request(url, headers=UA)
                return sym, json.loads(urllib.request.urlopen(req, timeout=180).read())
            except Exception:
                if attempt == 3:
                    return sym, []
                time.sleep(1.5 * attempt)
        return sym, []

    order = sorted(id_map.items())
    out = {}
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(one, s, e) for s, e in order]
        for k, f in enumerate(futs, 1):
            sym, recs = f.result()
            if recs:
                out[sym] = {r["sampleId"]: r.get("value") for r in recs}
            if verbose and k % 25 == 0:
                print(f"      {k}/{len(futs)}", flush=True)
    df = pd.DataFrame(out)
    df.index.name = "sampleId"
    if cache is not None:
        Path(cache).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(cache)
        if verbose:
            print(f"    cached -> {cache} {df.shape}")
    return df


def zscore_df(df: pd.DataFrame) -> pd.DataFrame:
    """按列(基因)在样本间 z-score。"""
    return (df - df.mean()) / df.std(ddof=0).replace(0, np.nan)


def signature_score(expr: pd.DataFrame, up: list, down: list | None = None,
                    min_genes: int = 5) -> pd.Series:
    """NPPscore = mean z(up 基因) - mean z(down 基因)。

    expr: 样本 x 基因 (raw 或 TPM 均可, 内部做 z-score)。
    """
    Z = zscore_df(expr)
    up = [g for g in up if g in Z.columns]
    down = [g for g in (down or []) if g in Z.columns]
    if len(up) < min_genes:
        raise ValueError(f"up 基因在表达矩阵中仅剩 {len(up)} 个 (<{min_genes})")
    s = Z[up].mean(axis=1)
    if len(down) >= min_genes:
        s = s - Z[down].mean(axis=1)
    elif down:
        s = s - Z[down].mean(axis=1)
    return s.rename("NPPscore")
