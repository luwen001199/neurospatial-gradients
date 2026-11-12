# -*- coding: utf-8 -*-
"""下载 METABRIC mRNA 表达 —— 限定目标基因 (修正版).

踩坑记录 (写进注释以免重复浪费时间):
1. cBioPortal `POST /api/genes/fetch` 在本题环境下返回 400, 不可用;
   `GET /api/genes?hugoGeneSymbol=` 返回的 entrezGeneId 也是错的 (NGF 返回 1)。
   => 改用 **NCBI E-utilities db=gene** 做 symbol -> Entrez 映射, 已验证准确:
      NGF->4803, CD8A->925, CALCA->796, ACTA2->59, UCHL1->7345
2. `/molecular-data` 必须同时提供两个参数, 缺任一即 HTTP 400:
      sampleListId  = brca_metabric_mrna   (Samples with mRNA data, Illumina HT-12 v3)
      entrezGeneIds = <逗号分隔>
   `hugoGeneSymbols` 参数不被该接口接受。
3. 只取方案的神经/CAF/免疫/TLS/PAM50 签名基因, 而非全转录组, 把数据量压到几 MB。
"""
import json
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

ROOT = Path(r"D:/嗜神经浸润")
OUT = ROOT / "01_data" / "raw" / "metabric"
OUT.mkdir(parents=True, exist_ok=True)

CBIO = "https://www.cbioportal.org/api"
PROFILE = "brca_metabric_mrna"
SAMPLE_LIST = "brca_metabric_mrna"
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)", "Accept": "application/json"}
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

TARGET_GENES = sorted({g.strip() for g in """
NGF UCHL1 TRPV1 CALCA CALCB S100B TH DBH SLC18A3 CHAT TUBB3 MAP2 NEFL NEFM SYP PLP1 MPZ
SCN9A SCN10A NTRK1 NTRK2 NTRK3 RAMP1 CALCRL CALCR ADCYAP1 VIP TAC1 GAL NPY SST PENK
SEMA3A SEMA3B SEMA3C SEMA3D SEMA3E SEMA3F SEMA3G NRP1 NRP2 PLXNA1 PLXNA2 PLXNA3 PLXNA4 PLXNB1 PLXND1
SLIT1 SLIT2 ROBO1 ROBO2 GAD1 GAD2 ABAT SLC6A1 GABBR1 GABBR2 GABRA1 GABRG2 BDNF NTF3 NTF4 SORT1 NGFR
ACTA2 TAGLN MYL9 CNN1 FAP PDGFRB PDGFRA PDPN S100A4 COL1A1 COL1A2 COL3A1 COL5A1 COL11A1 COL10A1
FN1 POSTN THBS2 INHBA CTHRC1 TGFB1 MMP2 MMP9 MMP11 MMP14 LOX LOXL2 SPARC VCAN TNC SERPINE1
LUM DCN GSN ITGA1 CXCL12 IL6 CXCL1 CXCL2 CXCL8 CCL2 LIF
CD8A CD8B CD3D CD3E CD2 GZMA GZMB GZMK PRF1 NKG7 IFNG TNF IL2 CCL5 CXCL9 CXCL10 PDCD1 CTLA4
LAG3 TIGIT HAVCR2 TOX TCF7 LEF1 SELL CD4 FOXP3 IL2RA CCR7 MKI67
CCL21 CCL19 CXCL13 LAMP3 MS4A1 CD79A CD79B CXCR5 CR2 AICDA POU2AF1
B2M HLA-A HLA-B HLA-C TAP1 TAP2 NLRC5 PSMB8 PSMB9 STAT1 IRF1 GBP1 IDO1 CD274 PDCD1LG2
HIF1A VEGFA SLC2A1 ADM EGLN3 CA9 LDHA PGK1 TWIST1 SNAI1 SNAI2 ZEB1 ZEB2 CDH1 CDH2 VIM
ESR1 PGR ERBB2 GRB7 AURKA FOXA1 MYBL2 CCNB1 MYC CCNE1 CDC20 CEP55 ANLN EXO1 KIF2C MELK
NDC80 NUF2 PTTG1 RRM2 TYMS UBE2C UBE2T CENPF KIF11 KPNA2 EGFR KRT5 KRT14 KRT17 MDM2 MLPH
MMP11 NAT1 PHGDH SCUBE2 SFRP1 SLC39A6 STARD3 GPR160 FOXC1 PRAME ESR1 PGR FOXA1 MKI67
""".split()})


def _get(url, timeout=90, retries=4):
    last = None
    for i in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers=UA)
            return json.loads(urllib.request.urlopen(req, timeout=timeout).read())
        except Exception as e:
            last = e
            print(f"    retry{i} {type(e).__name__} {str(e)[:70]}")
            time.sleep(2 * i)
    raise last


def symbols_to_entrez(symbols, batch=40):
    """NCBI E-utilities 批量 symbol -> Entrez GeneID. 严格校验回传 name 等于输入 symbol."""
    mapping, unmatched = {}, []
    for s in range(0, len(symbols), batch):
        chunk = symbols[s:s + batch]
        term = "(" + " OR ".join(f"{x}[sym]" for x in chunk) + ") AND 9606[taxid]"
        url = (f"{EUTILS}/esearch.fcgi?db=gene&retmax=200&retmode=json"
               f"&term={urllib.parse.quote(term)}")
        try:
            ids = _get(url).get("esearchresult", {}).get("idlist", [])
        except Exception as e:
            print(f"    esearch batch failed: {str(e)[:70]}")
            unmatched.extend(chunk)
            continue
        if not ids:
            unmatched.extend(chunk)
            continue
        url2 = f"{EUTILS}/esummary.fcgi?db=gene&retmode=json&id={','.join(ids)}"
        try:
            res = _get(url2).get("result", {})
        except Exception as e:
            print(f"    esummary batch failed: {str(e)[:70]}")
            unmatched.extend(chunk)
            continue
        got = {}
        for gid in ids:
            rec = res.get(gid)
            if not rec:
                continue
            name = (rec.get("name") or "").upper()
            if name in {c.upper() for c in chunk}:
                got[name] = int(gid)
        for c in chunk:
            if c.upper() in got:
                mapping[c] = got[c.upper()]
            else:
                unmatched.append(c)
        print(f"    mapped {len(mapping)}/{len(symbols)}", flush=True)
        time.sleep(0.4)
    return mapping, sorted(set(unmatched))


def fetch_matrix(id_map, workers=8):
    """逐基因拉取 molecular-data。

    踩坑: 该接口的 entrezGeneId 是 **单值** (java Long), 传入 "4803,925,796" 会抛
    NumberFormatException。因此不能一次批多个基因, 必须一基因一请求 (~2.3s/次),
    换来 8 路并发把总耗时从 ~10 分钟压到 ~1.5 分钟。
    """
    import urllib.request
    from concurrent.futures import ThreadPoolExecutor, as_completed

    def one(sym, eid):
        url = (f"{CBIO}/molecular-profiles/{PROFILE}/molecular-data"
               f"?sampleListId={SAMPLE_LIST}&entrezGeneId={eid}&pageSize=100000")
        for attempt in range(1, 4):
            try:
                req = urllib.request.Request(url, headers=UA)
                return sym, json.loads(urllib.request.urlopen(req, timeout=180).read())
            except Exception as e:
                if attempt == 3:
                    print(f"    !! {sym}({eid}) failed: {str(e)[:70]}")
                    return sym, []
                time.sleep(1.5 * attempt)
        return sym, []

    records, done = [], 0
    items = sorted(id_map.items())
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(one, s, e): s for s, e in items}
        for fu in as_completed(futs):
            sym, recs = fu.result()
            # 记录里带上 gene symbol, 便于后续透视
            for r in recs:
                r.setdefault("gene", {})
                if not r["gene"].get("hugoGeneSymbol"):
                    r["gene"]["hugoGeneSymbol"] = sym
            records.extend(recs)
            done += 1
            if done % 25 == 0 or done == len(items):
                print(f"    {done}/{len(items)} genes, {len(records)} records", flush=True)
    return records


def main():
    print(f"[1/4] target genes: {len(TARGET_GENES)}")
    cache = OUT / "metabric_genes_target.tsv"
    if cache.exists():
        id_map = {}
        for line in cache.read_text(encoding="utf8").strip().split("\n")[1:]:
            a, b = line.split("\t")
            id_map[a] = int(b)
        print(f"      reuse cached map ({len(id_map)} genes)")
    else:
        id_map, missing = symbols_to_entrez(TARGET_GENES)
        print(f"      mapped={len(id_map)} unmatched={len(missing)}")
        if missing:
            print(f"      unmatched symbols: {missing[:40]}")
        cache.write_text("hugoGeneSymbol\tentrezGeneId\n" +
                         "\n".join(f"{k}\t{v}" for k, v in sorted(id_map.items())),
                         encoding="utf8")

    print(f"[2/4] fetch expression ({len(id_map)} genes) ...")
    records = fetch_matrix(id_map)
    with open(OUT / f"{PROFILE}_targetgenes.json", "w", encoding="utf8") as f:
        json.dump(records, f)
    print(f"      raw records: {len(records)}")

    print("[3/4] pivot ...")
    import pandas as pd
    mat = defaultdict(dict)
    for r in records:
        sid = r.get("sampleId")
        gsym = (r.get("gene") or {}).get("hugoGeneSymbol")
        if sid and gsym:
            try:
                mat[sid][gsym] = float(r.get("value"))
            except (TypeError, ValueError):
                pass
    df = pd.DataFrame(mat).T
    df.index.name = "sampleId"
    df.to_csv(OUT / "metabric_mrna_targetgenes.csv")
    print(f"      matrix: {df.shape[0]} samples x {df.shape[1]} genes")

    if df.shape[0] > 0:
        na = df.isna().mean().sort_values(ascending=False)
        print("[4/4] top missing genes:", dict((na.head(8) * 100).round(1)))
    print(f"      -> {OUT/'metabric_mrna_targetgenes.csv'}")


if __name__ == "__main__":
    main()
