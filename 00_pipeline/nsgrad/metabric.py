# -*- coding: utf-8 -*-
"""METABRIC 临床数据加载 (cBioPortal 导出的是 JSON 数组, 不是 TSV)。"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np
import pandas as pd


def _load_json_records(path: Path) -> pd.DataFrame:
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    rows = []
    for rec in raw:
        attr = rec.get("clinicalAttribute") or {}
        rows.append({
            "patientId": rec.get("patientId"),
            "sampleId": rec.get("sampleId"),
            "attrId": rec.get("clinicalAttributeId") or attr.get("clinicalAttributeId"),
            "value": rec.get("value"),
        })
    return pd.DataFrame(rows)


def wide(df: pd.DataFrame, key: str) -> pd.DataFrame:
    """长表 -> 宽表 (key = patientId 或 sampleId)。"""
    d = df.dropna(subset=[key]).copy()
    d = d.drop_duplicates(subset=[key, "attrId"])
    out = d.pivot(index=key, columns="attrId", values="value")
    out.columns.name = None
    return out


def load_clinical(dir_raw: Path):
    dir_raw = Path(dir_raw)
    pat = wide(_load_json_records(dir_raw / "clinical_patient.txt"), "patientId")
    smp = wide(_load_json_records(dir_raw / "clinical_sample.txt"), "sampleId")
    return pat, smp


def build_outcome_table(pat: pd.DataFrame, smp: pd.DataFrame) -> pd.DataFrame:
    """构造生存分析表: 每个 PATIENT_ID 一行, 含 RFS / OS / 亚型 / 分期等。"""
    p = pat.copy()
    # sample -> patient 映射
    if "PATIENT_ID" in smp.columns:
        s2p = smp["PATIENT_ID"]
    else:
        s2p = pd.Series(smp.index, index=smp.index)

    keep = [c for c in [
        "OS_MONTHS", "OS_STATUS", "RFS_MONTHS", "RFS_STATUS",
        "AGE_AT_DIAGNOSIS", "TUMOR_SIZE", "TUMOR_STAGE", "GRADE",
        "LYMPH_NODES_EXAMINED_POSITIVE",
        "ER_IHC", "ER_STATUS", "PR_STATUS", "HER2_SNP6", "HER2_STATUS",
        "CLAUDIN_SUBTYPE", "PAM50_SUBTYPE", "INTCLUST", "THREEGENE",
        "CHEMOTHERAPY", "HORMONE_THERAPY", "RADIO_THERAPY", "RADIOTHERAPY",
        "BREAST_SURGERY", "INFERRED_MENOPAUSAL_STATE", "MENOPAUSAL_STATUS",
        "NPI", "CELLULARITY", "COHORT", "SEX", "VITAL_STATUS",
        "HISTOLOGICAL_SUBTYPE", "LATERALITY", "SAMPLE_COUNT",
        "FIRST_RELAPSE_SITE", "RELAPSE_SITE",
    ] if c in p.columns]
    t = p[keep].copy()
    t.index.name = "PATIENT_ID"
    if "RADIO_THERAPY" in t and "RADIOTHERAPY" not in t:
        t["RADIOTHERAPY"] = t["RADIO_THERAPY"]

    def num(x):
        return pd.to_numeric(x, errors="coerce")

    t["OS_MONTHS"] = num(t.get("OS_MONTHS"))
    t["RFS_MONTHS"] = num(t.get("RFS_MONTHS"))
    # 1:0 编码: cBioPortal 的 STATUS 形如 "1:DECEASED" / "0:LIVING"
    def ev(x):
        s = str(x)
        return 1 if s.startswith("1") else (0 if s.startswith("0") else np.nan)
    t["OS_EVENT"] = t.get("OS_STATUS").map(ev) if "OS_STATUS" in t else np.nan
    t["RFS_EVENT"] = t.get("RFS_STATUS").map(ev) if "RFS_STATUS" in t else np.nan
    t["AGE"] = num(t.get("AGE_AT_DIAGNOSIS"))
    t["SIZE"] = num(t.get("TUMOR_SIZE"))
    t["NODES_POS"] = num(t.get("LYMPH_NODES_EXAMINED_POSITIVE"))
    t["GRADE"] = num(t.get("GRADE"))
    t["STAGE"] = num(t.get("TUMOR_STAGE"))
    t["NPI"] = num(t.get("NPI"))
    for c in ["ER_IHC", "ER_STATUS", "PR_STATUS", "HER2_SNP6", "HER2_STATUS",
              "CHEMOTHERAPY", "HORMONE_THERAPY", "RADIOTHERAPY",
              "INFERRED_MENOPAUSAL_STATE", "MENOPAUSAL_STATUS"]:
        if c in t:
            t[c] = (t[c].astype(str).str.upper()
                    .replace({"NAN": np.nan, "NONE": np.nan, "": np.nan}))
    return t
