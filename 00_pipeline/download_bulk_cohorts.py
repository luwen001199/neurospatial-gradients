# -*- coding: utf-8 -*-
"""下载本方案中不依赖争议的 bulk 队列:

1) METABRIC   —— cBioPortal REST API (2509 例, Illumina HT-12)
2) GSE12276   —— NCBI GEO series matrix (204 例, GPL570 Affymetrix HG-U133 Plus 2.0)

这两个队列均已被 tool-use 阶段的 GEO/cBioPortal 元数据核实为真实存在且内容匹配,
因此先行落地, 为 Step4 的临床验证环节提供数据底座。
网络较慢, 全部请求带重试与超时控制。
"""
import gzip
import json
import os
import shutil
import time
import urllib.request
from pathlib import Path

ROOT = Path(r"D:/嗜神经浸润")
RAW = ROOT / "01_data" / "raw"
RAW.mkdir(parents=True, exist_ok=True)

UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}
CBIO = "https://www.cbioportal.org/api"


def http_get(url, timeout=90, retries=3, decode=True, stream_to=None):
    """带重试的 GET。stream_to 给定 Path 时流式落盘。"""
    last = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                if stream_to:
                    tmp = str(stream_to) + ".part"
                    with open(tmp, "wb") as f:
                        shutil.copyfileobj(r, f, length=1 << 20)
                    os.replace(tmp, stream_to)
                    return stream_to.stat().st_size
                raw = r.read()
                return raw.decode("utf8", "ignore") if decode else raw
        except Exception as e:
            last = e
            print(f"    retry {attempt}/{retries} after {type(e).__name__}: {str(e)[:80]}")
            time.sleep(2 * attempt)
    raise last


# ---------------------------------------------------------------- METABRIC
def download_metabric():
    print("[METABRIC] query molecular profiles ...")
    prof_url = f"{CBIO}/studies/brca_metabric/molecular-profiles?projection=SUMMARY&pageSize=200"
    try:
        profs = json.loads(http_get(prof_url, timeout=120))
    except Exception as e:
        print("[METABRIC] FAIL profiles:", e)
        return {"status": "fail_profiles", "error": str(e)[:150]}

    ids = [p.get("molecularProfileId") for p in profs]
    print("[METABRIC] profiles:", ids)

    out = {"profiles": ids, "downloaded": {}}
    (RAW / "metabric").mkdir(exist_ok=True)

    # Save profile list for later selection
    (RAW / "metabric" / "molecular_profiles.json").write_text(
        json.dumps(profs, ensure_ascii=False, indent=2), encoding="utf8")

    # 1) Clinical patient-level
    for ctype, fname in [("PATIENT", "clinical_patient.txt"),
                         ("SAMPLE", "clinical_sample.txt")]:
        u = (f"{CBIO}/studies/brca_metabric/clinical-data"
             f"?clinicalDataType={ctype}&projection=DETAILED&pageSize=100000")
        try:
            txt = http_get(u, timeout=180)
            p = RAW / "metabric" / fname
            p.write_text(txt, encoding="utf8")
            out["downloaded"][fname] = len(txt)
            print(f"[METABRIC] {fname}: {len(txt)} chars")
        except Exception as e:
            out["downloaded"][fname] = f"FAIL {e}"
            print(f"[METABRIC] {fname} FAIL: {e}")

    # 2) Pick mRNA expression profile (Illumina HT-12 => usually *_mrna)
    cand = [i for i in ids if "rna_seq" not in i and ("mrna" in i or "mrna_median" in i)]
    for pid in cand:
        u = f"{CBIO}/molecular-profiles/{pid}/molecular-data?projection=SUMMARY&pageSize=100000"
        try:
            txt = http_get(u, timeout=600)
            p = RAW / "metabric" / f"{pid}.json"
            p.write_text(txt, encoding="utf8")
            out["downloaded"][pid] = len(txt)
            print(f"[METABRIC] {pid}: {len(txt)} chars")
        except Exception as e:
            out["downloaded"][pid] = f"FAIL {str(e)[:120]}"
            print(f"[METABRIC] {pid} FAIL: {str(e)[:120]}")
    return out


# ---------------------------------------------------------------- GSE12276
def download_gse12276():
    print("\n[GSE12276] try multiple mirrors for series matrix ...")
    urls = [
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE12nnn/GSE12276/matrix/GSE12276_series_matrix.txt.gz",
        "http://ftp.ncbi.nlm.nih.gov/geo/series/GSE12nnn/GSE12276/matrix/GSE12276_series_matrix.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE12nnn/GSE12276/soft/GSE12276_family.soft.gz",
    ]
    dest_dir = RAW / "gse12276"
    dest_dir.mkdir(exist_ok=True)
    for u in urls:
        name = u.rsplit("/", 1)[-1]
        dest = dest_dir / name
        try:
            n = http_get(u, timeout=420, retries=2, stream_to=dest)
            print(f"[GSE12276] OK {name} -> {n} bytes")
            return {"status": "ok", "file": str(dest), "bytes": n, "url": u}
        except Exception as e:
            print(f"[GSE12276] FAIL {name}: {str(e)[:100]}")
    return {"status": "fail"}


def main():
    rep = {"metabric": {}, "gse12276": {}}
    try:
        rep["metabric"] = download_metabric()
    except Exception as e:
        rep["metabric"] = {"status": "fail", "error": str(e)[:200]}
        print("[METABRIC] fatal:", e)
    rep["gse12276"] = download_gse12276()

    p = ROOT / "01_data" / "availability" / "download_report.json"
    p.write_text(json.dumps(rep, ensure_ascii=False, indent=2), encoding="utf8")
    print(f"\nSaved -> {p}")


if __name__ == "__main__":
    main()
