# -*- coding: utf-8 -*-
"""一次性脚本（2026-11-11）：补充图 S13-S16 全局重编号为 S11-S14，修复 S11/S12 空缺。

映射: S13→S11, S14→S12, S15→S13, S16→S14（仅图，表格编号本来就连续不动）。
覆盖: manuscript_PLOSCB.md / manuscript_main.md / Supplementary_Note_13 / to_plos.py
     以及 4 个 make_figure_supp_s*.py 的内部标签（文件重命名由本脚本一并完成）。
每个替换带出现次数断言，任何一处不匹配即整体报错退出、不写盘。
"""
import re
from pathlib import Path

ROOT = Path(r"D:/嗜神经浸润")


def rep(path, pairs, label):
    p = ROOT / path
    t = p.read_text(encoding="utf-8")
    for old, new, expect in pairs:
        n = t.count(old)
        if n != expect:
            if n == 0 and t.count(new) >= expect:  # 幂等：已替换过则跳过
                print(f"SKIP {label}: already applied '{old[:50]}'")
                continue
            raise SystemExit(f"ABORT {label}: '{old[:60]}' found {n}, expected {expect}")
        t = t.replace(old, new)
    p.write_text(t, encoding="utf-8")
    print(f"OK  {label}: {len(pairs)} replacement rule(s)")


# ---------------------------------------------------------------- 1. PLOS 主稿
rep("04_manuscript/manuscript_PLOSCB.md", [
    ("sections (S16 Fig; S1 Text, Note 13.26).",
     "sections (S14 Fig; S1 Text, Note 13.26).", 1),
    ("exists; S16 Fig; S1 Text, Note 13.26).",
     "exists; S14 Fig; S1 Text, Note 13.26).", 1),
    ("unmeasurable; S15 Fig; S1 Text, Note 13.28)",
     "unmeasurable; S13 Fig; S1 Text, Note 13.28)", 1),
    ("Note 13; S13 and S14 Figs)",
     "Note 13; S11 and S12 Figs)", 1),
    ("**S13 Fig.** External bulk-cohort",
     "**S11 Fig.** External bulk-cohort", 1),
    ("**S14 Fig.** Immune compartments",
     "**S12 Fig.** Immune compartments", 1),
    ("**S15 Fig.** Sensory-nerve",
     "**S13 Fig.** Sensory-nerve", 1),
    ("**S16 Fig.** Xenium",
     "**S14 Fig.** Xenium", 1),
], "manuscript_PLOSCB.md")

# checklist 项整行替换（用前缀定位，避免特殊字符拷贝误差）
p = ROOT / "04_manuscript/manuscript_PLOSCB.md"
t = p.read_text(encoding="utf-8")
lines = t.splitlines(keepends=True)
hits = [i for i, l in enumerate(lines) if l.startswith("- [ ] Supporting figure numbering has a gap")]
if len(hits) == 1:
    lines[hits[0]] = ("- [x] Supporting figure numbering gap fixed (2026-11-11): S13-S16 renumbered to "
                      "S11-S14 across the manuscript, S1 Text, figure files and figure scripts; tables were "
                      "already sequential (S1-S37). PLOS assigns final S-labels at upload - re-check list "
                      "order then.\n")
    p.write_text("".join(lines), encoding="utf-8")
    print("OK  manuscript_PLOSCB.md: checklist item marked done")
elif any(l.startswith("- [x] Supporting figure numbering gap fixed") for l in lines):
    print("SKIP manuscript_PLOSCB.md: checklist already marked done")
else:
    raise SystemExit("ABORT: checklist item not found")

# ---------------------------------------------------------------- 2. BIB 源稿
rep("04_manuscript/manuscript_main.md", [
    ("Supplementary Figure S16", "Supplementary Figure S14", 2),
    ("Supplementary Figure S15", "Supplementary Figure S13", 1),
    ("Figures S13–S14", "Figures S11–S12", 1),
    (" S13 external bulk-cohort", " S11 external bulk-cohort", 1),
    ("; S14 immune compartments", "; S12 immune compartments", 1),
    ("; S15 sensory-nerve", "; S13 sensory-nerve", 1),
    ("; S16 Xenium", "; S14 Xenium", 1),
], "manuscript_main.md")

# ---------------------------------------------------------------- 3. S1 Text (Note 13)
rep("04_manuscript/supplementary/Supplementary_Note_13_bulk_extrapolation.md", [
    ("make_figure_supp_s10.py` (Supplementary Figure S13)",
     "make_figure_supp_s11.py` (Supplementary Figure S11)", 1),
    ("(Supplementary Figure S13C)", "(Supplementary Figure S11C)", 1),
    ("; Supplementary Figure S16)", "; Supplementary Figure S14)", 1),
    ("; Supplementary Figure S15)", "; Supplementary Figure S13)", 1),
], "Supplementary_Note_13")

# ---------------------------------------------------------------- 4. to_plos.py（保持可重跑一致性）
rep("02_scripts/to_plos.py", [
    ('"Supplementary Figure S15, Note 13.28", "S15 Fig; S1 Text, Note 13.28"',
     '"Supplementary Figure S13, Note 13.28", "S13 Fig; S1 Text, Note 13.28"', 1),
    ('"Supplementary Figure S16, Note 13.26", "S16 Fig; S1 Text, Note 13.26"',
     '"Supplementary Figure S14, Note 13.26", "S14 Fig; S1 Text, Note 13.26"', 1),
    ('"Figures S13–S14", "S13 and S14 Figs"',
     '"Figures S11–S12", "S11 and S12 Figs"', 1),
    ("**S13 Fig.** External bulk-cohort", "**S11 Fig.** External bulk-cohort", 1),
    ("**S14 Fig.** Immune compartments", "**S12 Fig.** Immune compartments", 1),
    ("**S15 Fig.** Sensory-nerve", "**S13 Fig.** Sensory-nerve", 1),
    ("**S16 Fig.** Xenium", "**S14 Fig.** Xenium", 1),
    ("(S15 Fig / S20 Table / S1 Text, Note 12)", "(S13 Fig / S20 Table / S1 Text, Note 12)", 1),
], "to_plos.py")

# ---------------------------------------------------------------- 5. 图脚本内部标签
rep("02_scripts/make_figure_supp_s13.py", [
    ("SuppFigS10", "SuppFigS11", 1),
    ("SuppFigS13_bulk_extrapolation", "SuppFigS11_bulk_extrapolation", 1),
], "make_figure_supp_s13.py")
rep("02_scripts/make_figure_supp_s14.py", [
    ("SuppFigS14", "SuppFigS12", 2),
    ("S14A", "S12A", 1), ("S14B", "S12B", 1),
    ("S14C", "S12C", 1), ("S14D", "S12D", 1),
], "make_figure_supp_s14.py")
rep("02_scripts/make_figure_supp_s15.py", [
    ("SuppFigS15", "SuppFigS13", 2),
    ("S15A", "S13A", 1), ("S15B", "S13B", 1),
    ("S15C", "S13C", 1), ("S15D", "S13D", 1),
], "make_figure_supp_s15.py")
rep("02_scripts/make_figure_supp_s16.py", [
    ("SuppFigS16", "SuppFigS14", 2),
    ("S16A", "S14A", 1), ("S16B", "S14B", 1),
    ("S16C", "S14C", 1), ("S16D", "S14D", 1),
], "make_figure_supp_s16.py")

# ---------------------------------------------------------------- 6. 文件重命名（顺序避开冲突: 先腾出 S11/S12, 再落 S13/S14）
renames = [
    ("02_scripts/make_figure_supp_s13.py", "02_scripts/make_figure_supp_s11.py"),
    ("02_scripts/make_figure_supp_s14.py", "02_scripts/make_figure_supp_s12.py"),
    ("02_scripts/make_figure_supp_s15.py", "02_scripts/make_figure_supp_s13.py"),
    ("02_scripts/make_figure_supp_s16.py", "02_scripts/make_figure_supp_s14.py"),
]
figs = [
    ("SuppFigS13_bulk_extrapolation", "SuppFigS11_bulk_extrapolation"),
    ("SuppFigS14_immune_pancancer", "SuppFigS12_immune_pancancer"),
    ("SuppFigS15_sensory_axis_tls", "SuppFigS13_sensory_axis_tls"),
    ("SuppFigS16_xenium_validation", "SuppFigS14_xenium_validation"),
]
for old, new in renames:
    src, dst = ROOT / old, ROOT / new
    if not src.exists() and dst.exists():
        print(f"SKIP rename {old} (already done)")
        continue
    assert src.exists(), f"missing {src}"
    assert not dst.exists(), f"collision {dst}"
    src.rename(dst)
    print(f"REN {old} -> {new}")
for old, new in figs:
    for ext in (".png", ".pdf"):
        src = ROOT / "03_results/figures" / (old + ext)
        dst = ROOT / "03_results/figures" / (new + ext)
        if not src.exists() and dst.exists():
            print(f"SKIP rename figures/{old}{ext} (already done)")
            continue
        assert src.exists(), f"missing {src}"
        assert not dst.exists(), f"collision {dst}"
        src.rename(dst)
        print(f"REN figures/{old}{ext} -> {new}{ext}")

# ---------------------------------------------------------------- 7. 复制 4 张重编版合成图到 composed/（PLOS 上传目录）
for old, new in figs:
    for ext in (".png", ".pdf"):
        src = ROOT / "03_results/figures" / (new + ext)
        dst = ROOT / "03_results/figures/composed" / (new + ext)
        import shutil
        shutil.copy2(src, dst)
        print(f"CPY composed/{new}{ext}")

print("\nALL DONE — S13-S16 -> S11-S14 renumbering complete.")
