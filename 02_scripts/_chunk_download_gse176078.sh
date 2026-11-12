#!/bin/bash
# GSE176078 剩余 tar.gz 分块并行下载 (与 h5 修复同一策略: 8 块并行 + 校验 + 拼接)
cd "D:/嗜神经浸润/01_data/external/gse176078"
LOG="D:/嗜神经浸润/03_results/gse176078_chunk.log"
: > "$LOG"
NCHUNK=8

chunk_worker() {
  url=$1; out=$2; s=$3; e=$4
  want=$((e - s + 1))
  for round in $(seq 1 60); do
    curl -sL --ssl-no-revoke --connect-timeout 30 --max-time 900 -r "${s}-${e}" -o "$out" "$url"
    got=$(stat -c%s "$out" 2>/dev/null || echo 0)
    [ "$got" -eq "$want" ] && return 0
    sleep 3
  done
  return 1
}

mapfile -t files < _dl.txt
for f in "${files[@]}"; do
  [ -z "$f" ] && continue
  exp=$(awk -F'\t' -v n="$f" '$2==n {print $4}' filelist.txt)
  sz=$(stat -c%s "$f" 2>/dev/null || echo 0)
  if [ "$sz" -ge "$exp" ]; then echo "SKIP $f (done)" >> "$LOG"; continue; fi
  gsm=${f%%_*}
  URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5354nnn/${gsm}/suppl/${f}"
  rm -f _chunk.* 2>/dev/null
  CSZ=$(( (exp + NCHUNK - 1) / NCHUNK ))
  pids=()
  for i in $(seq 0 $((NCHUNK-1))); do
    s=$((i * CSZ)); e=$(( (i+1)*CSZ - 1 )); [ $e -ge $exp ] && e=$((exp-1))
    [ $s -gt $e ] && continue
    chunk_worker "$URL" "_chunk.$i" "$s" "$e" &
    pids+=($!)
  done
  fail=0
  for p in "${pids[@]}"; do wait "$p" || fail=1; done
  if [ $fail -eq 1 ]; then echo "$f CHUNK_FAIL" >> "$LOG"; continue; fi
  cat $(for i in $(seq 0 $((NCHUNK-1))); do [ -f "_chunk.$i" ] && echo "_chunk.$i"; done) > "$f"
  got=$(stat -c%s "$f" 2>/dev/null || echo 0)
  if [ "$got" -eq "$exp" ]; then echo "$f OK $got $(date +%H:%M:%S)" >> "$LOG"; else echo "$f SIZE_MISMATCH $got/$exp" >> "$LOG"; fi
  rm -f _chunk.* 2>/dev/null
done
echo "GSE176078_ALL_DONE $(date)" >> "$LOG"
