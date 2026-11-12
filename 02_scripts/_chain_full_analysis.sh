#!/bin/bash
# 全量 GSE176078 (26/26) 到位后自动执行: 来源归属 -> 解卷积(主+no-module) -> 共聚集 -> S9图 -> S18/S19表
cd "D:/嗜神经浸润"
SRC=01_data/external/gse176078

echo "[chain] waiting for 26/26 GSE176078 downloads + extraction..."
for i in $(seq 1 720); do   # 最多等 24 小时 (每 2 分钟一查)
  ok=0
  while read f; do
    [ -z "$f" ] && continue
    exp=$(awk -F'\t' -v n="$f" '$2==n {print $4}' "$SRC/filelist.txt")
    sz=$(stat -c%s "$SRC/$f" 2>/dev/null || echo 0)
    cid=${f#*_}; cid=${cid%.tar.gz}
    if [ "$sz" -ge "$exp" ] && [ -f "$SRC/$cid/count_matrix_sparse.mtx" ]; then ok=$((ok+1)); fi
  done < "$SRC/_dl.txt"
  if [ "$ok" -ge 26 ]; then echo "[chain] all 26 ready at $(date +%H:%M)"; break; fi
  [ $((i % 15)) -eq 0 ] && echo "[chain] progress $ok/26 at $(date +%H:%M)"
  sleep 120
done
if [ "$ok" -lt 26 ]; then echo "[chain] TIMEOUT: only $ok/26"; exit 1; fi

echo "[chain] STEP1 scrna_cell_origin (26-sample reference)"
python 02_scripts/run_scrna_cell_origin.py > 03_results/scrna_origin_full.log 2>&1 || { echo "[chain] FAIL origin"; exit 1; }

echo "[chain] STEP2 deconvolution primary"
python 02_scripts/run_deconvolution.py > 03_results/deconv_full.log 2>&1 || { echo "[chain] FAIL deconv"; exit 1; }

echo "[chain] STEP3 deconvolution no-module sensitivity"
python 02_scripts/run_deconvolution.py --no-module > 03_results/deconv_nomodule.log 2>&1 || { echo "[chain] FAIL deconv nomodule"; exit 1; }

echo "[chain] STEP4 colocalization primary"
python 02_scripts/run_colocalization.py > 03_results/coaggr_full.log 2>&1 || { echo "[chain] FAIL coaggr"; exit 1; }

echo "[chain] STEP5 colocalization no-module"
python 02_scripts/run_colocalization.py --no-module > 03_results/coaggr_nomodule.log 2>&1 || { echo "[chain] FAIL coaggr nomodule"; exit 1; }

echo "[chain] STEP6 figure S9"
python 02_scripts/make_figure_scrna_anchor.py > 03_results/figS9_full.log 2>&1 || { echo "[chain] FAIL figS9"; exit 1; }

echo "[chain] STEP7 tables S18/S19 + xlsx rebuild"
python 02_scripts/export_scrna_tables.py > 03_results/tabS18_full.log 2>&1 && \
python 02_scripts/build_supp_tables_xlsx.py >> 03_results/tabS18_full.log 2>&1 || { echo "[chain] FAIL tables"; exit 1; }

echo "[chain] ALL DONE at $(date +%H:%M)"
