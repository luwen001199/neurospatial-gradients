#!/bin/bash
# MSR scan 并行驱动: 4 个 worker lane, 每样本独立 tag 输出, 避免互相覆盖。
# 用法: bash _msr_parallel.sh <bin_target> <tag_prefix>   例如:
#   bash _msr_parallel.sh 100 ""      -> msr_per_sample.csv.gz (主)
#   bash _msr_parallel.sh 40  fine    -> msr_per_sample.fine.csv.gz
cd "D:/嗜神经浸润"
BIN=${1:-100}
TAG=${2:-}
NW=${MSR_WORKERS:-12}
LOG="03_results/msr_par_${TAG:-main}.log"
: > "$LOG"

mapfile -t GSMS < <(ls 01_data/processed/spot_level/*.h5ad | sed 's|.*/||; s|\.h5ad||' | sort)

worker() {
  w=$1; n=${#GSMS[@]}
  for ((i=w; i<n; i+=NW)); do
    g=${GSMS[$i]}
    echo "[w$w] $g start $(date +%H:%M:%S)" >> "$LOG"
    # 每样本独立 tag, 防止同 worker 内后一样本覆盖前一样本
    OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 \
    python 02_scripts/run_msr_spatial.py scan --only "$g" \
      --bin "$BIN" --tag "s_${g}_${TAG}" >> "03_results/msr_${TAG:-main}_$g.log" 2>&1
    echo "[w$w] $g done $(date +%H:%M:%S)" >> "$LOG"
  done
}
for ((w=0; w<NW; w++)); do worker $w & done
wait

# 合并所有 lane 输出 -> 正式文件名
python - "$TAG" <<'EOF'
import glob, sys, pandas as pd
tag = sys.argv[1]
parts = sorted(glob.glob(rf"D:/嗜神经浸润/01_data/processed/msr_per_sample.s_GSM*_{tag}.csv.gz"))
print("merging", len(parts), "parts")
assert len(parts) >= 40, f"only {len(parts)} parts, expected 43"
df = pd.concat([pd.read_csv(p) for p in parts], ignore_index=True)
out = rf"D:/嗜神经浸润/01_data/processed/msr_per_sample{'.' + tag if tag else ''}.csv.gz"
df.to_csv(out, index=False, compression="gzip")
print("->", out, len(df), "rows")
EOF
echo "PARALLEL_SCAN_DONE $(date)" >> "$LOG"
