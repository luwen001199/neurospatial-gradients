# -*- coding: utf-8 -*-
"""把 04_manuscript/manuscript_main.md 转成投稿用 Word (.docx), BIB (format-free) 版式。

支持: #/##/### 标题、**粗体**、*斜体*、无序/有序/复选列表、分隔线、行内代码。
 Times New Roman 10pt, 1.5 倍行距, 双倍行距可切换 (BIB 要求 double-spaced 提交)。
输出: 04_manuscript/manuscript_main_BIB.docx
"""
import re
from pathlib import Path

from docx import Document
from docx.enum.text import WD_LINE_SPACING
from docx.oxml.ns import qn
from docx.shared import Pt

ROOT = Path(r"D:/嗜神经浸润")
SRC = ROOT / "04_manuscript/manuscript_main.md"
OUT = ROOT / "04_manuscript/manuscript_main_BIB.docx"

LINE_SPACING = 1.5   # 改为 2.0 即为期刊要求的 double-spaced


def set_font(run, size=10, bold=None, italic=None, mono=False):
    run.font.name = "Consolas" if mono else "Times New Roman"
    run._element.rPr.rFonts.set(qn("w:eastAsia"), "宋体")
    run.font.size = Pt(size)
    if bold is not None:
        run.bold = bold
    if italic is not None:
        run.italic = italic


INLINE = re.compile(r"(\*\*\*.+?\*\*\*|\*\*.+?\*\*|\*(?!\s)[^*]+?\*|`[^`]+`)")


def add_inline(par, text, size=10, base_bold=False):
    for tok in INLINE.split(text):
        if not tok:
            continue
        if tok.startswith("***") and tok.endswith("***") and len(tok) > 6:
            set_font(par.add_run(tok[3:-3]), size, True, True)
        elif tok.startswith("**") and tok.endswith("**") and len(tok) > 4:
            set_font(par.add_run(tok[2:-2]), size, True, base_bold)
        elif tok.startswith("*") and tok.endswith("*") and len(tok) > 2:
            set_font(par.add_run(tok[1:-1]), size, base_bold, True)
        elif tok.startswith("`") and tok.endswith("`") and len(tok) > 2:
            set_font(par.add_run(tok[1:-1]), size - 0.5, mono=True)
        else:
            set_font(par.add_run(tok), size, base_bold)


def para(doc, text, style=None, size=10, bold=False, space_before=0, space_after=6,
         indent=0.0, align_left=True):
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.line_spacing_rule = WD_LINE_SPACING.MULTIPLE
    pf.line_spacing = LINE_SPACING
    pf.space_before = Pt(space_before)
    pf.space_after = Pt(space_after)
    if indent:
        pf.left_indent = Pt(indent)
    if not align_left:
        p.alignment = 1
    add_inline(p, text, size=size, base_bold=bold)
    return p


def main():
    lines = SRC.read_text(encoding="utf-8").splitlines()
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(10)

    for raw in lines:
        line = raw.rstrip()
        if not line.strip():
            continue
        if line.startswith("### "):
            para(doc, line[4:], size=11, bold=True, space_before=10, space_after=4)
        elif line.startswith("## "):
            para(doc, line[3:], size=12, bold=True, space_before=14, space_after=6)
        elif line.startswith("# "):
            para(doc, line[2:], size=14, bold=True, space_after=10)
        elif re.match(r"^-{3,}$", line.strip()):
            continue
        elif re.match(r"^\s*-\s*\[[ x]\]", line):
            m = re.match(r"^\s*-\s*\[([ x])\]\s*(.*)$", line)
            para(doc, ("[x] " if m.group(1) == "x" else "[ ] ") + m.group(2),
                 indent=18, space_after=2)
        elif re.match(r"^\s*-\s+", line):
            para(doc, "• " + re.sub(r"^\s*-\s+", "", line), indent=18, space_after=2)
        elif re.match(r"^\s*\d+\.\s+", line):
            para(doc, line.strip(), indent=18, space_after=2)
        else:
            para(doc, line.strip())

    doc.save(OUT)
    print("saved ->", OUT)


if __name__ == "__main__":
    main()
