#!/bin/bash
# GSE75688 带退避重试: TPM 矩阵 + 样本信息 (NCBI 403 限流后自动重试)
cd "D:/嗜神经浸润/01_data/external/gse75688"
D=01_data/external/gse75688

want() {  # url, out, expected_min_bytes
  local url=$1 out=$2 min=$3
  for i in $(seq 1 30); do
    sz=$(stat -c%s "$out" 2>/dev/null || echo 0)
    if [ "$sz" -ge "$min" ] && gzip -t "$out" 2>/dev/null; then echo "OK $out ($sz)"; return 0; fi
    # 文件存在但不是 gzip 魔数 -> 删除重来, 避免续传拼接损坏
    if [ -f "$out" ]; then
      magic=$(head -c 2 "$out" | od -An -tx1 | tr -d ' \n')
      [ "$magic" != "1f8b" ] && rm -f "$out"
    fi
    curl -sL --ssl-no-revoke --http1.1 --max-time 900 -C - "$url" -o "$out"
    sleep $((i * 120))   # 递增退避: 2,4,6,...分钟
  done
  echo "FAIL $out"
  return 1
}

# 以 gzip 完整性为准; TPM 至少 1MB, sample_info 至少 2KB
want "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75688/suppl/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt.gz" TPM.txt.gz 1000000
want "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75688/suppl/GSE75688_final_sample_information.txt.gz" sample_info.txt.gz 2000
zcat TPM.txt.gz | wc -l
zcat sample_info.txt.gz | head -3
echo GSE75688_DONE
