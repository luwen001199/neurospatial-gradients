#!/usr/bin/env bash
# 续传下载 bulk 队列（支持 -C - 断点续传；Windows 需 --ssl-no-revoke）
set -u
cd "D:/嗜神经浸润/01_data/raw" || exit 1
CURL="curl -sS -g --ssl-no-revoke --retry 8 --retry-delay 3 --retry-all-errors -C -"

dl () {  # url out want_bytes
  local url="$1" out="$2" want="$3"
  for i in $(seq 1 30); do
    "$CURL" -m 600 -o "$out" "$url" && break
    sleep 2
  done
  local have
  have=$(stat -c %s "$out" 2>/dev/null || echo 0)
  echo "[$(basename "$out")] $have / $want bytes"
}

echo ">>> TCGA HiSeqV2"
dl "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FHiSeqV2.gz" \
   "tcga_brca/HiSeqV2.gz" 64286337

echo ">>> TCGA clinical"
dl "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix" \
   "tcga_brca/BRCA_clinicalMatrix" 1698413

echo ">>> GSE163882 TPM"
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE163nnn/GSE163882/suppl/GSE163882_all.data.tpms_222Samples.csv.gz" \
   "gse163882/GSE163882_all.data.tpms_222Samples.csv.gz" 8752814

echo ">>> GSE25065 series matrix"
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE25nnn/GSE25065/matrix/GSE25065_series_matrix.txt.gz" \
   "gse25065/GSE25065_series_matrix.txt.gz" 0

echo "=== DONE ==="
ls -la tcga_brca gse163882 gse25065
