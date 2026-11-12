#!/bin/bash
# 分块并行下载器: 每个文件切成 8 个字节区间并行拉取, 每块独立重试,
# 全部到位后按序拼接 -> 整体大小校验 -> h5py 结构校验。
# 适用于 NCBI 单连接限速 (~13KB/s) 但多连接可扩速的网络环境。
cd "D:/嗜神经浸润"
BASE="01_data/raw/gse210616"
LOG="03_results/fix13_chunk.log"
: > "$LOG"

declare -A EXP=(
  [GSM6433585]=11971283  [GSM6433588]=11154472 [GSM6433590]=10250753
  [GSM6433591]=5947059   [GSM6433592]=9585572  [GSM6433594]=13058530
  [GSM6433596]=10991372  [GSM6433610]=6358567  [GSM6433614]=9042318
  [GSM6433618]=11364941  [GSM6433622]=6928440
)

NCHUNK=8

chunk_worker() {
  # $1 url $2 dest.partN $3 start $4 end ; 循环直到块完整
  url=$1; out=$2; s=$3; e=$4
  want=$((e - s + 1))
  for round in $(seq 1 40); do
    curl -s --ssl-no-revoke --connect-timeout 30 --max-time 600 -r "${s}-${e}" -o "$out" "$url"
    got=$(stat -c%s "$out" 2>/dev/null || echo 0)
    [ "$got" -eq "$want" ] && return 0
    sleep 3
  done
  return 1
}

dl_one() {
  GSM=$1
  WANT=${EXP[$GSM]}
  DIR="$BASE/$GSM"
  H5=$(ls "$DIR"/*_filtered_feature_bc_matrix.h5 2>/dev/null | head -1)
  [ -z "$H5" ] && { echo "$GSM NO_H5" >> "$LOG"; return 1; }
  URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/${GSM:0:7}nnn/$GSM/suppl/$(basename "$H5")"
  rm -f "$DIR"/_chunk.*
  CSZ=$(( (WANT + NCHUNK - 1) / NCHUNK ))
  pids=()
  for i in $(seq 0 $((NCHUNK-1))); do
    s=$((i * CSZ)); e=$(( (i+1)*CSZ - 1 )); [ $e -ge $WANT ] && e=$((WANT-1))
    [ $s -gt $e ] && continue
    chunk_worker "$URL" "$DIR/_chunk.$i" "$s" "$e" &
    pids+=($!)
  done
  fail=0
  for p in "${pids[@]}"; do wait "$p" || fail=1; done
  if [ $fail -eq 1 ]; then echo "$GSM CHUNK_FAIL" >> "$LOG"; rm -f "$DIR"/_chunk.*; return 1; fi
  cat $(for i in $(seq 0 $((NCHUNK-1))); do [ -f "$DIR/_chunk.$i" ] && echo "$DIR/_chunk.$i"; done) > "$H5"
  got=$(stat -c%s "$H5" 2>/dev/null || echo 0)
  if [ "$got" -ne "$WANT" ]; then echo "$GSM SIZE_MISMATCH $got/$WANT" >> "$LOG"; rm -f "$DIR"/_chunk.*; return 1; fi
  rm -f "$DIR"/_chunk.*
  if python -c "
import h5py, sys
try:
    with h5py.File(r'''$H5''', 'r') as f:
        assert list(f['matrix/shape'][:])[0] == 36601
    sys.exit(0)
except Exception as e:
    print(e); sys.exit(1)
" >>"$LOG" 2>&1; then
    echo "$GSM OK $got bytes $(date +%H:%M:%S)" >> "$LOG"; return 0
  else
    echo "$GSM H5PY_FAIL" >> "$LOG"; return 1
  fi
}

echo "start $(date)" >> "$LOG"
for GSM in GSM6433585 GSM6433588 GSM6433590 GSM6433591 GSM6433592 GSM6433594 GSM6433596 GSM6433610 GSM6433614 GSM6433618 GSM6433622; do
  dl_one "$GSM" >> "$LOG" 2>&1
done
rm -f "$BASE/GSM6433624"/*.png.part 2>/dev/null
echo "ALL_DONE $(date)" >> "$LOG"
