"""Resumable HTTP(S) downloader with per-chunk retries (GEO ftp throttling tolerant)."""
import os
import sys
import time
import urllib.request
import urllib.error

UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}
CHUNK = 1 << 18  # 256 KB read blocks
MAX_TRIES = 40
READ_TIMEOUT = 45


def head_size(url):
    req = urllib.request.Request(url, headers=UA, method="HEAD")
    with urllib.request.urlopen(req, timeout=30) as r:
        return int(r.headers.get("Content-Length", -1))


def download(url, out, chunk=CHUNK):
    os.makedirs(os.path.dirname(out), exist_ok=True)
    total = head_size(url)
    part = out + ".part"
    pos = os.path.getsize(part) if os.path.exists(part) else 0
    if os.path.exists(out) and total > 0 and os.path.getsize(out) == total:
        print("already complete", out, flush=True)
        return out
    if total > 0 and pos > total:
        pos = 0
    print("total=%.1f MB, resume at %.1f MB" % (total / 1e6, pos / 1e6), flush=True)
    t0 = time.time()
    tries = 0
    while total < 0 or pos < total:
        end = "" if total < 0 else str(min(pos + chunk - 1, total - 1))
        hdr = dict(UA)
        hdr["Range"] = "bytes=%d-%s" % (pos, end)
        try:
            req = urllib.request.Request(url, headers=hdr)
            with urllib.request.urlopen(req, timeout=READ_TIMEOUT) as r:
                buf = r.read()
            if not buf:
                raise urllib.error.URLError("empty chunk")
            with open(part, "ab") as f:
                f.write(buf)
            pos += len(buf)
            tries = 0
            if total > 0 and pos % (4 << 20) < chunk:
                el = time.time() - t0
                rate = pos / 1e6 / max(el, 1e-9)
                eta = (total - pos) / 1e6 / max(rate, 1e-9)
                print("  %.1f/%.1f MB  %.2f MB/s  eta %.0fs"
                      % (pos / 1e6, total / 1e6, rate, eta), flush=True)
        except Exception as e:
            tries += 1
            if tries > MAX_TRIES:
                print("!! give up after %d tries: %s" % (tries, e), flush=True)
                raise
            time.sleep(min(30, 1.5 * tries))
    if total > 0 and os.path.getsize(part) != total:
        raise IOError("size mismatch: %d vs %d" % (os.path.getsize(part), total))
    os.replace(part, out)
    print("done %.1f MB in %.0fs -> %s" % (total / 1e6, time.time() - t0, out), flush=True)
    return out


if __name__ == "__main__":
    download(sys.argv[1], sys.argv[2])
