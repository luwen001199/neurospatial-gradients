# -*- coding: utf-8 -*-
"""GSE210616 并行下载器 (多线程 + Range 断点续传).

为什么必须并行: 实测 NCBI GEO FTP 单条连接吞吐约 0.16 MB/s, 是**连接级**限速而非账号级。
43 个样本串行需要 ~1.1 小时; 开 N 条并行连接, 总吞吐近似 N 倍叠加, 可压到 5-10 分钟。

同时实现断点续传: 已下载的 .part 文件用 HTTP Range 头从断点继续, 网络中断不会浪费进度。

下载内容优先级:
  P1 (必须) *_filtered_feature_bc_matrix.h5  表达矩阵
     (必须) *_tissue_positions_list.csv.gz   spot 坐标
     (必须) *_scalefactors_json.json.gz      图像换算参数
  P2 (后续) *_tissue_hires_image.png.gz      病理全切片图, 用于空间图背景
"""
import json
import os
import shutil
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(r"D:/嗜神经浸润")
RAW = ROOT / "01_data" / "raw" / "gse210616"
RAW.mkdir(parents=True, exist_ok=True)

UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}
ACC = "GSE210616"
FILELIST = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE210nnn/{ACC}/suppl/filelist.txt"
N_WORKERS = 12          # 并行连接数
CHUNK = 1 << 18         # 256 KB

P2_KINDS = {"PNG"}      # 病理图放第二轮


def gsm_dir(gsm):
    return f"{gsm[:7]}nnn"


def parse_filelist():
    txt = urllib.request.urlopen(
        urllib.request.Request(FILELIST, headers=UA), timeout=120).read().decode("utf8", "ignore")
    rows = []
    for line in txt.split("\n"):
        p = line.rstrip("\r").split("\t")
        if len(p) >= 5 and p[0] == "File" and p[1].startswith("GSM"):
            gsm = p[1].split("_")[0]
            rows.append({"gsm": gsm, "name": p[1], "size": int(p[3]), "kind": p[4]})
    return rows


def download_one(row, include_p2=False, retries=4):
    """单文件下载, 支持 Range 续传。返回 (name, bytes_done, status)."""
    if row["kind"] in P2_KINDS and not include_p2:
        return row["name"], 0, "skip-p2"
    gsm = row["gsm"]
    outdir = RAW / gsm
    outdir.mkdir(exist_ok=True)
    final = outdir / (row["name"][:-3] if row["name"].endswith(".gz") else row["name"])
    # 关键修复: 早前的单线程下载器没有做 size 校验, 中断时会留下 *非零但被截断* 的文件。
    # 若仅判断 size>0 就当作 cached, 这些残帧会永久混进数据里 (h5py 打开时报 truncated file)。
    # 因此必须要求字节数与 GEO filelist 声明的 size **完全相等**才算完成, 否则删除重下。
    if final.exists():
        if final.stat().st_size == row["size"]:
            return row["name"], final.stat().st_size, "cached"
        final.unlink()   # 大小不符 = 残帧, 重来

    part = Path(str(final) + ".part")
    have = part.stat().st_size if part.exists() else 0
    if have == row["size"]:           # 已下完但没改名
        os.replace(part, final)
        return row["name"], have, "completed"

    url = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm_dir(gsm)}/{gsm}/suppl/{row['name']}"
    last = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers=UA)
            if have > 0:
                req.add_header("Range", f"bytes={have}-")
            ifhave = have
            with urllib.request.urlopen(req, timeout=120) as r:
                mode = "ab" if (ifhave > 0 and r.status == 206) else "wb"
                with open(part, mode) as f:
                    while True:
                        b = r.read(CHUNK)
                        if not b:
                            break
                        f.write(b)
            sz = part.stat().st_size
            if sz >= row["size"]:
                os.replace(part, final)   # 移动后 part 即消失, 必须先用变量接住 size
                return row["name"], sz, "ok"
            have = sz
        except Exception as e:
            last = e
            time.sleep(1.5 * attempt)
            have = part.stat().st_size if part.exists() else 0
    return row["name"], (part.stat().st_size if part.exists() else 0), f"FAIL:{type(last).__name__}"


def main(include_p2=False):
    rows = parse_filelist()
    want = [r for r in rows if r["size"] < 200 * 1024 * 1024
            and (include_p2 or r["kind"] not in P2_KINDS)]
    total_bytes = sum(r["size"] for r in want)
    gsms = sorted({r["gsm"] for r in want})
    print(f"[plan] {len(want)} files / {len(gsms)} samples / "
          f"{total_bytes/1e6:.1f} MB  workers={N_WORKERS}")

    results = {}
    t0 = time.time()
    done = 0
    with ThreadPoolExecutor(max_workers=N_WORKERS) as ex:
        futs = {ex.submit(download_one, r, include_p2): r for r in want}
        for fu in as_completed(futs):
            r = futs[fu]
            name, size, status = fu.result()
            results[name] = {"size": size, "status": status, "gsm": r["gsm"]}
            done += 1
            el = time.time() - t0
            got = sum(v["size"] for v in results.values())
            spd = got / 1e6 / max(el, 1)
            if status == "skip-p2":
                continue
            eta = (total_bytes - got) / 1e6 / spd if spd > 0 else -1
            print(f"  [{done}/{len(want)}] {r['gsm']} {r['kind']:<5} "
                  f"{size/1e6:7.2f}MB {status:<9} "
                  f"| total {got/1e6:7.1f}MB {spd:5.2f}MB/s eta {eta:6.0f}s", flush=True)

    el = time.time() - t0
    got = sum(v["size"] for v in results.values())
    fails = [k for k, v in results.items() if str(v["status"]).startswith("FAIL")]
    samples_done = len({v["gsm"] for v in results.values()
                        if v["status"] in ("ok", "cached", "completed")})
    rep = {
        "acc": ACC, "seconds": round(el, 1),
        "files_ok": len(results) - len(fails), "failed_files": fails,
        "downloaded_mb": round(got / 1e6, 2), "target_mb": round(total_bytes / 1e6, 2),
        "samples_with_files": samples_done, "n_samples": len(gsms),
    }
    p = ROOT / "01_data" / "availability" / "download_gse210616_parallel.json"
    p.write_text(json.dumps(rep, ensure_ascii=False, indent=2), encoding="utf8")
    print(f"\n[done] {rep['files_ok']} files OK, {len(fails)} failed, "
          f"{rep['downloaded_mb']} MB in {el/60:.1f} min")
    print(f"       report -> {p}")
    if fails:
        print("[failed sample names]:", fails[:10])


if __name__ == "__main__":
    import sys
    main(include_p2=("--with-images" in sys.argv))
