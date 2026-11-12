# -*- coding: utf-8 -*-
"""GEO / cBioPortal 数据集可达性与元数据探查脚本.

项目: TNBC感觉神经驱动的空间免疫架构分型与PNI进展梯度研究
用途: 在执行正式分析前, 确认每个公开数据集的原始文件 URL、样本量、平台与补充文件格式。
"""
import gzip
import json
import os
import re
import urllib.request
from pathlib import Path

ROOT = Path(r"D:/嗜神经浸润")
OUT = ROOT / "01_data" / "availability"
OUT.mkdir(parents=True, exist_ok=True)

UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}

# acc -> GEO FTP 目录规则
TARGETS = {
    "GSE12276": "GSE12nnn",
    "GSE256279": "GSE256nnn",
    "GSE325196": "GSE325nnn",
    "GSE240212": "GSE240nnn",
}


def fetch(url, timeout=90, decode=True):
    req = urllib.request.Request(url, headers=UA)
    raw = urllib.request.urlopen(req, timeout=timeout).read()
    return raw.decode("utf8", "ignore") if decode else raw


def fetch_gzip_text(url, timeout=120):
    return gzip.decompress(fetch(url, timeout, decode=False)).decode("utf8", "ignore")


def probe_geo(acc, dirn):
    """返回该 GSE 的 soft 元数据摘要与可下载文件清单."""
    info = {"acc": acc, "soft_url": None, "title": None, "n_samples": None,
            "platforms": [], "supp": [], "matrix": None, "error": None}
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{dirn}/{acc}"
    # 1) family.soft.gz
    try:
        u = f"{base}/soft/{acc}_family.soft.gz"
        txt = fetch_gzip_text(u)
        info["soft_url"] = u
        m = re.search(r"!Series_title\t(.*)", txt)
        info["title"] = m.group(1).strip() if m else None
        info["platforms"] = list(dict.fromkeys(re.findall(r"!Sample_platform_id\t(\S+)", txt)))
        info["n_samples"] = len([l for l in txt.split("\n") if l.startswith("^SAMPLE")])
        summ = re.search(r"!Series_summary\t(.*)", txt)
        info["summary"] = (summ.group(1).strip()[:400] if summ else None)
        pub = re.findall(r"!Series_pubmed_id\t(\S+)", txt)
        info["pubmed"] = pub
    except Exception as e:
        info["error"] = f"soft:{type(e).__name__}:{e}"[:150]

    # 2) suppl file list via GEO acc.cgi API (lighter)
    for ftype in ["file", "brief"]:
        pass
    return info


def probe_matrix(acc, dirn):
    """尝试 series matrix. 有些 GSE 只有 GSE***_series_matrix.txt.gz."""
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{dirn}/{acc}/matrix"
    res = {}
    for name in [f"{acc}_series_matrix.txt.gz"]:
        try:
            req = urllib.request.Request(f"{base}/{name}", headers=UA, method="HEAD")
            r = urllib.request.urlopen(req, timeout=60)
            res[name] = {"status": r.status, "bytes": r.headers.get("Content-Length")}
        except Exception as e:
            res[name] = {"status": "FAIL", "error": str(e)[:90]}
    return res


def main():
    report = {"geo": {}, "matrix": {}, "cbioportal": {}}
    for acc, dirn in TARGETS.items():
        try:
            report["geo"][acc] = probe_geo(acc, dirn)
            print(f"[OK] {acc}: n_samples={report['geo'][acc]['n_samples']} "
                  f"plat={report['geo'][acc]['platforms']} title={str(report['geo'][acc]['title'])[:70]}")
        except Exception as e:
            report["geo"][acc] = {"error": str(e)[:120]}
            print(f"[FAIL] {acc}: {e}")
        report["matrix"][acc] = probe_matrix(acc, dirn)
        print(f"     matrix -> {report['matrix'][acc]}")

    # cBioPortal METABRIC
    try:
        u = "https://www.cbioportal.org/api/studies/brca_metabric"
        d = json.loads(fetch(u, 60))
        report["cbioportal"]["brca_metabric"] = {
            "name": d.get("name"),
            "description": (d.get("description") or "")[:200],
            "allSampleCount": d.get("allSampleCount"),
            "sequencedSampleCount": d.get("sequencedSampleCount"),
            "cnaSampleCount": d.get("cnaSampleCount"),
            "pmid": d.get("pmid"),
        }
        print("[OK] METABRIC:", report["cbioportal"]["brca_metabric"])
    except Exception as e:
        report["cbioportal"]["brca_metabric"] = {"error": str(e)[:150]}
        print("[FAIL] METABRIC:", e)

    p = OUT / "availability_report.json"
    p.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf8")
    print(f"\nSaved -> {p}")


if __name__ == "__main__":
    main()
