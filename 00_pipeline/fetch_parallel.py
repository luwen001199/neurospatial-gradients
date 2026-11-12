# -*- coding: utf-8 -*-
"""Parallel range-request downloader —— 单连接被 NCBI 限速到 ~10 KB/s 时,
开 N 个并发连接各取一段, 再拼回完整文件。每段独立断点续传。

用法: python fetch_parallel.py <url> <out> [n_chunks=8]
"""
import os
import sys
import threading
import time
import urllib.request

UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}


def head_size(url):
    req = urllib.request.Request(url, headers=UA, method="HEAD")
    with urllib.request.urlopen(req, timeout=30) as r:
        return int(r.headers.get("Content-Length", -1))


def fetch_range(url, start, end, path, idx, errors):
    while True:
        pos = os.path.getsize(path) if os.path.exists(path) else 0
        want = end - start + 1
        if pos >= want:
            return
        hdr = dict(UA)
        hdr["Range"] = f"bytes={start + pos}-{end}"
        try:
            req = urllib.request.Request(url, headers=hdr)
            with urllib.request.urlopen(req, timeout=30) as r:
                while True:
                    b = r.read(1 << 16)
                    if not b:
                        break
                    with open(path, "ab") as f:
                        f.write(b)
        except Exception as e:
            errors.append(e)
            time.sleep(2)


def main(url, out, n_chunks=8):
    total = head_size(url)
    print(f"total {total/1e6:.1f} MB, {n_chunks} parallel chunks", flush=True)
    edges = [int(total * i / n_chunks) for i in range(n_chunks + 1)]
    threads, errors = [], []
    for i in range(n_chunks):
        start, end = edges[i], edges[i + 1] - 1
        p = f"{out}.p{i:02d}"
        t = threading.Thread(target=fetch_range,
                             args=(url, start, end, p, i, errors), daemon=True)
        t.start()
        threads.append(t)
    t0 = time.time()
    last = -1
    while any(t.is_alive() for t in threads):
        done = sum(os.path.getsize(f"{out}.p{i:02d}")
                   if os.path.exists(f"{out}.p{i:02d}") else 0
                   for i in range(n_chunks))
        if done // (4 << 20) != last:
            last = done // (4 << 20)
            rate = done / 1e6 / max(time.time() - t0, 1e-9)
            print(f"  {done/1e6:6.1f}/{total/1e6:.1f} MB  {rate:5.2f} MB/s",
                  flush=True)
        time.sleep(3)
    for t in threads:
        t.join()
    sizes = [os.path.getsize(f"{out}.p{i:02d}") if os.path.exists(f"{out}.p{i:02d}")
             else 0 for i in range(n_chunks)]
    if sum(sizes) != total or min(sizes) != max(sizes) and False:
        pass
    expect = [edges[i + 1] - edges[i] for i in range(n_chunks)]
    if sizes != expect:
        raise SystemExit(f"chunk mismatch: {sizes} vs {expect}; rerun to resume")
    with open(out, "wb") as fo:
        for i in range(n_chunks):
            with open(f"{out}.p{i:02d}", "rb") as fi:
                fo.write(fi.read())
            os.remove(f"{out}.p{i:02d}")
    print(f"done {total/1e6:.1f} MB in {time.time()-t0:.0f}s -> {out}", flush=True)


if __name__ == "__main__":
    url, out = sys.argv[1], sys.argv[2]
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    main(url, out, n)
