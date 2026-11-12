# -*- coding: utf-8 -*-
"""把 04_manuscript/manuscript_PLOSCB.md 转成 PLOS Computational Biology 投稿版 docx。

版式: Times New Roman 12pt, 双倍行距, 全文连续行号 (PLOS 审稿要求),
内部提交清单段(## Submission checklist 之后)不进入 docx。
输出: 04_manuscript/manuscript_PLOSCB.docx
"""
import re
from pathlib import Path

from docx import Document
from docx.enum.text import WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt

ROOT = Path(r"D:/嗜神经浸润")
SRC = ROOT / "04_manuscript/manuscript_PLOSCB.md"
OUT = ROOT / "04_manuscript/manuscript_PLOSCB.docx"

LINE_SPACING = 2.0
BASE_SIZE = 12


def set_font(run, size=BASE_SIZE, bold=None, italic=None, mono=False):
    run.font.name = "Consolas" if mono else "Times New Roman"
    run._element.rPr.rFonts.set(qn("w:eastAsia"), "宋体")
    run.font.size = Pt(size)
    if bold is not None:
        run.bold = bold
    if italic is not None:
        run.italic = italic


INLINE = re.compile(r"(\*\*\*.+?\*\*\*|\*\*.+?\*\*|\*(?!\s)[^*]+?\*|`[^`]+`)")


def add_inline(par, text, size=BASE_SIZE, base_bold=False):
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
            set_font(par.add_run(tok[1:-1]), size - 1, mono=True)
        else:
            set_font(par.add_run(tok), size, base_bold)


def para(doc, text, size=BASE_SIZE, bold=False, space_before=0, space_after=6,
         indent=0.0, italic=False):
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.line_spacing_rule = WD_LINE_SPACING.MULTIPLE
    pf.line_spacing = LINE_SPACING
    pf.space_before = Pt(space_before)
    pf.space_after = Pt(space_after)
    if indent:
        pf.left_indent = Pt(indent)
    add_inline(p, text, size=size, base_bold=bold)
    if italic:
        for r in p.runs:
            r.italic = True
    return p


def enable_line_numbers(doc):
    sectPr = doc.sections[0]._sectPr
    ln = OxmlElement("w:lnNumType")
    ln.set(qn("w:countBy"), "1")
    ln.set(qn("w:restart"), "continuous")
    sectPr.append(ln)


def main():
    lines = SRC.read_text(encoding="utf-8").splitlines()
    # 内部清单不进 docx
    try:
        cut = next(i for i, ln in enumerate(lines) if ln.startswith("## Submission checklist"))
        lines = lines[:cut]
    except StopIteration:
        pass

    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(BASE_SIZE)

    for raw in lines:
        line = raw.rstrip()
        if not line.strip():
            continue
        if line.startswith("### "):
            para(doc, line[4:], size=12, bold=True, space_before=12, space_after=4)
        elif line.startswith("## "):
            para(doc, line[3:], size=13, bold=True, space_before=16, space_after=6)
        elif line.startswith("# "):
            para(doc, line[2:], size=14, bold=True, space_after=10)
        elif re.match(r"^-{3,}$", line.strip()):
            continue
        elif re.match(r"^\s*-\s*\[[ x]\]", line):
            continue  # 检查框行不进投稿版
        elif re.match(r"^\s*-\s+", line):
            para(doc, "• " + re.sub(r"^\s*-\s+", "", line), indent=18, space_after=2)
        elif re.match(r"^\s*\d+\.\s+", line):
            para(doc, line.strip(), indent=18, space_after=2)
        else:
            para(doc, line.strip())

    enable_line_numbers(doc)
    doc.save(OUT)
    print("saved ->", OUT)


if __name__ == "__main__":
    main()
