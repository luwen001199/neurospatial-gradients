#!/bin/bash
# 重下 11 个截断的 GSE210616 h5: 全新下载 + 大小校验 + h5py 可读性验证 + 指数退避
cd "D:/嗜神经浸润"
BASE="01_data/raw/gse210616"
LOG="03_results/fix13_download.log"
: > "$LOG"

declare -A EXP=(
  [GSM6433585]=11971283  [GSM6433588]=11154472 [GSM6433590]=10250753
  [GSM6433591]=5947059   [GSM6433592]=9585572  [GSM6433594]=13058530
  [GSM6433596]=10991372  [GSM6433610]=6358567  [GSM6433614]=9042318
  [GSM6433618]=11364941  [GSM6433622]=6928440
)

for GSM in GSM6433585 GSM6433588 GSM6433590 GSM6433591 GSM6433592 GSM6433594 GSM6433596 GSM6433610 GSM6433614 GSM6433618 GSM6433622; do
  DIR="$BASE/$GSM"
  H5=$(ls "$DIR"/*_filtered_feature_bc_matrix.h5 2>/dev/null | head -1)
  [ -z "$H5" ] && { echo "$GSM NO_H5_FILE" >> "$LOG"; continue; }
  WANT=${EXP[$GSM]}
  URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/${GSM:0:7}nnn/$GSM/suppl/$(basename "$H5")"

  ok=0
  for attempt in 1 2 3 4 5 6; do
    rm -f "$H5"
    curl -s --ssl-no-revoke --connect-timeout 60 --max-time 900 -o "$H5" "$URL"
    sz=$(stat -c%s "$H5" 2>/dev/null || echo 0)
    if [ "$sz" -eq "$WANT" ]; then
      # h5py 结构校验
      if python -c "
import h5py, sys
try:
    with h5py.File(r'''$H5''', 'r') as f:
        assert list(f['matrix/shape'][:])[0] == 36601
    sys.exit(0)
except Exception as e:
    print(e); sys.exit(1)
" 2>>"$LOG"; then ok=1; break; fi
    fi
    echo "$GSM attempt $attempt: size $sz / want $WANT, retrying" >> "$LOG"
    sleep $((attempt * 20))
  done
  if [ $ok -eq 1 ]; then echo "$GSM OK ($(stat -c%s "$H5") bytes)" >> "$LOG"; else echo "$GSM FAILED after 6 attempts" >> "$LOG"; fi
done

# 清掉 GSM6433624 的 .part 残留(它会让 discover_samples 跳过该样本)
rm -f "$BASE/GSM6433624"/*.png.part
echo "DONE" >> "$LOG"
