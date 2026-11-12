#!/bin/bash
# GSE350571 (Xenium Prime 5K, 6 例乳腺 FFPE) 必需文件下载编排
# 阶段 1: 仅拉 6 个 gene_panel.json.gz (每个 ~100-200KB) -> 验证神经标记 go/no-go
# 阶段 2 (由 run_gse350571_nerve.sh 或人工触发): cells/h5/transcripts
set -u
ROOT="D:/嗜神经浸润/01_data/raw/gse350571"
mkdir -p "$ROOT"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM10116nnn"

declare -A GSMS=(
  [GSM10116352]=P17  [GSM10116353]=P22  [GSM10116354]=P24
  [GSM10116355]=P26  [GSM10116356]=P29  [GSM10116357]=P32
)

# ---- 阶段 1: panel json (小文件, 直接 curl 重试) ----
for gsm in GSM10116352 GSM10116353 GSM10116354 GSM10116355 GSM10116356 GSM10116357; do
  p="${GSMS[$gsm]}"
  out="$ROOT/${p}_panel.json.gz"
  [ -f "$out" ] && gunzip -t "$out" 2>/dev/null && { echo "SKIP $p"; continue; }
  url="$BASE/$gsm/suppl/${gsm}_${p}_Xe5k_gene_panel.json.gz"
  for i in $(seq 1 60); do
    curl -sS --ssl-no-revoke -f --connect-timeout 30 --max-time 120 -o "$out" "$url" 2>/dev/null
    if [ -f "$out" ] && gunzip -t "$out" 2>/dev/null; then echo "OK $p ($(stat -c%s "$out") bytes)"; break; fi
    rm -f "$out"; sleep 8
  done
  [ -f "$out" ] || echo "FAIL $p"
done
echo "=== 阶段 1 完成 ==="
ls -la "$ROOT"/*_panel.json.gz 2>/dev/null

# ---- 阶段 2 (仅当 STAGE2=1): 三个核心数据文件, 用分块下载器 ----
if [ "${STAGE2:-0}" = "1" ]; then
  DL="D:/嗜神经浸润/02_scripts/_chunkdl_url.sh"
  for gsm in GSM10116352 GSM10116353 GSM10116354 GSM10116355 GSM10116356 GSM10116357; do
    p="${GSMS[$gsm]}"
    for kind in cells.parquet.gz cell_feature_matrix.h5 transcripts.parquet.gz; do
      out="$ROOT/${p}_Xe5k_${kind}"
      # 跳过前先做轻量完整性校验: gzip 流校验 / HDF5 magic 校验
      ok=0
      if [ -f "$out" ] && [ "$(stat -c%s "$out")" -gt 1000000 ]; then
        case "$kind" in
          *.gz) gunzip -t "$out" 2>/dev/null && ok=1 ;;
          *)    [ "$(head -c 8 "$out" | od -An -tx1 | tr -d ' \n')" = "894844460d0a1a0a" ] && ok=1 ;;
        esac
      fi
      [ "$ok" = "1" ] && { echo "SKIP $p $kind"; continue; }
      [ -f "$out" ] && rm -f "$out"
      url="$BASE/$gsm/suppl/${gsm}_${p}_Xe5k_${kind}"
      echo ">>> downloading $p $kind"
      bash "$DL" "$url" "$out" "${CHUNKS:-4}" || echo "FAIL $p $kind"
    done
  done
fi
echo "=== ALL DONE ==="
