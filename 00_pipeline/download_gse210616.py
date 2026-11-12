# -*- coding: utf-8 -*-
"""下载 GSE210616 (43 例人源 TNBC Visium) —— 逐样本直取, 绕开 46 GB 的 RAW.tar.

关键发现: GEO 为每个 GSM 提供独立 suppl 目录
    https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6433nnn/{GSM}/suppl/{file}
因此只抓取每个样本真正需要的三类文件, 跳过体积占比 99% 的 .cloupe 文件:
    - *_filtered_feature_bc_matrix.h5   (~12 MB/样本, 表达矩阵, SpaceRanger 标准输出)
    - *_tissue_positions_list.csv.gz    (~60 KB/样本, spot 坐标)
    - *_scalefactors_json.json.gz       (~170 B/样本, 图中 spot 直径换算)
    - *_tissue_hires_image.png.gz       (~5 MB/样本, 可选, 病理图)

预计总下载量 ~700 MB, 而非 RAW.tar 的 46.65 GB.

同时落地 GSE210616 的 family.soft.gz 以获得样本级临床注释 (种族分层/年龄/分期等)。
"""
import gzip
import json
import os
import re
import shutil
import time
import urllib.request
from pathlib import Path

ROOT = Path(r"D:/嗜神经浸润")
RAW = ROOT / "01_data" / "raw" / "gse210616"
RAW.mkdir(parents=True, exist_ok=True)

UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}
ACC = "GSE210616"
FILELIST = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE210nnn/{ACC}/suppl/filelist.txt"

WANT_KINDS = {"H5", "CSV", "JSON", "PNG"}
SKIP_KINDS = {"CLOUPE"}


def http_stream(url, dest, expected_size=None, timeout=300, retries=3):
    """流式下载到 dest, 带重试。

    dest 存在且**字节数与 filelist 声明一致**才跳过 —— 只查 size>0 会让
    truncated 文件永久混入(GSM6433614/618 的教训)。
    """
    if dest.exists() and dest.stat().st_size > 0 and \
            (expected_size is None or dest.stat().st_size == expected_size):
        return dest.stat().st_size, "cached"
    last = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                tmp = str(dest) + ".part"
                with open(tmp, "wb") as f:
                    shutil.copyfileobj(r, f, length=1 << 20)
                os.replace(tmp, dest)
                return dest.stat().st_size, "ok"
        except Exception as e:
            last = e
            print(f"      retry {attempt}/{retries}: {type(e).__name__} {str(e)[:70]}")
            time.sleep(3 * attempt)
    return 0, f"FAIL {last}"


def gsm_dir(gsm):
    """GSM6433585 -> GSM6433nnn"""
    return f"{gsm[:7]}nnn"


def parse_filelist():
    txt = urllib.request.urlopen(
        urllib.request.Request(FILELIST, headers=UA), timeout=120).read().decode("utf8", "ignore")
    rows = []
    for line in txt.split("\n"):
        parts = line.rstrip("\r").split("\t")
        if len(parts) < 5:
            continue
        kind, name, size = parts[0], parts[1], parts[3]
        if kind == "File":
            m = re.match(r"(GSM\d+)_", name)
            if m:
                rows.append({"gsm": m.group(1), "name": name,
                             "size": int(size), "kind": parts[4]})
    return rows


def main():
    print("[1/3] parse filelist ...")
    rows = parse_filelist()
    gsms = sorted({r["gsm"] for r in rows})
    print(f"      {len(rows)} files / {len(gsms)} samples")

    want = [r for r in rows if r["kind"] in WANT_KINDS and r["kind"] not in SKIP_KINDS]
    skip = [r for r in rows if r["kind"] in SKIP_KINDS]
    print(f"      download {len(want)} files ({sum(r['size'] for r in want)/1e6:.1f} MB), "
          f"skip {len(skip)} cloupe ({sum(r['size'] for r in skip)/1e9:.1f} GB)")

    log = {}
    n = 0
    total = len(want)
    for r in want:
        n += 1
        gsm = r["gsm"]
        outdir = RAW / gsm
        outdir.mkdir(exist_ok=True)
        dest = outdir / (r["name"][:-3] if r["name"].endswith(".gz") else r["name"])
        url = (f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm_dir(gsm)}/{gsm}/suppl/{r['name']}")
        size, status = http_stream(url, dest, expected_size=int(r["size"]))
        log[r["name"]] = {"size": size, "status": status, "gsm": gsm}
        print(f"      [{n}/{total}] {r['name'][:58]:<58} {size/1e6:8.2f} MB {status}")
        if status.startswith("FAIL"):
            print(f"      !! failed: {url}")

    print("[2/3] download family.soft.gz for clinical annotation ...")
    soft_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE210nnn/{ACC}/soft/{ACC}_family.soft.gz"
    dest = RAW / f"{ACC}_family.soft.gz"
    sz, st = http_stream(soft_url, dest, timeout=600)
    print(f"      soft: {sz/1e6:.2f} MB {st}")

    report = {
        "acc": ACC,
        "n_samples": len(gsms),
        "downloaded_files": sum(1 for v in log.values() if v["status"] in ("ok", "cached")),
        "failed_files": [k for k, v in log.items() if str(v["status"]).startswith("FAIL")],
        "total_mb": round(sum(v["size"] for v in log.values()) / 1e6, 2),
        "detail": log,
    }
    p = ROOT / "01_data" / "availability" / "download_gse210616.json"
    p.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf8")
    print(f"[3/3] saved -> {p}")
    print(f"      OK={report['downloaded_files']} FAIL={len(report['failed_files'])} "
          f"total={report['total_mb']} MB")


if __name__ == "__main__":
    main()
