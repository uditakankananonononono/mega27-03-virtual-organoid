"""Build Times New Roman research papers (docx) with blue page borders,
numbered equations, tables and figures from real result files."""
from __future__ import annotations

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


def _page_border(section, color="1F4E79", size=12):
    sectPr = section._sectPr
    pb = OxmlElement("w:pgBorders"); pb.set(qn("w:offsetFrom"), "page")
    for side in ("top", "left", "bottom", "right"):
        el = OxmlElement(f"w:{side}")
        el.set(qn("w:val"), "single"); el.set(qn("w:sz"), str(size))
        el.set(qn("w:space"), "24"); el.set(qn("w:color"), color)
        pb.append(el)
    sectPr.append(pb)


class Paper:
    def __init__(self, title: str, authors: str):
        self.doc = Document(); self.eq = 0; self.fig = 0; self.tab = 0
        st = self.doc.styles["Normal"]; st.font.name = "Times New Roman"; st.font.size = Pt(12)
        st.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
        for s in ("Heading 1", "Heading 2", "Heading 3", "Title"):
            self.doc.styles[s].font.name = "Times New Roman"
            self.doc.styles[s].font.color.rgb = RGBColor(0x1F, 0x4E, 0x79)
        _page_border(self.doc.sections[0])
        t = self.doc.add_paragraph(title, style="Title"); t.alignment = WD_ALIGN_PARAGRAPH.CENTER
        a = self.doc.add_paragraph(authors); a.alignment = WD_ALIGN_PARAGRAPH.CENTER

    def h(self, text, level=1):
        self.doc.add_heading(text, level=level)

    def p(self, text):
        para = self.doc.add_paragraph(text); para.paragraph_format.space_after = Pt(6)
        para.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        return para

    def equation(self, text: str) -> int:
        self.eq += 1
        para = self.doc.add_paragraph(); para.alignment = WD_ALIGN_PARAGRAPH.CENTER
        r = para.add_run(text); r.italic = True
        para.add_run(f"\t\t({self.eq})")
        return self.eq

    def table(self, header: list[str], rows: list[list], caption: str):
        self.tab += 1
        self.doc.add_paragraph(f"Table {self.tab}. {caption}").runs[0].bold = True
        t = self.doc.add_table(rows=1, cols=len(header)); t.style = "Light Grid Accent 1"
        for i, hdr in enumerate(header):
            t.rows[0].cells[i].text = str(hdr)
        for r in rows:
            cells = t.add_row().cells
            for i, v in enumerate(r):
                cells[i].text = f"{v:.3f}" if isinstance(v, float) else str(v)

    def figure(self, path: str, caption: str, width: float = 6.0):
        self.fig += 1
        self.doc.add_picture(path, width=Inches(width))
        self.doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
        c = self.doc.add_paragraph(f"Figure {self.fig}. {caption}"); c.runs[0].italic = True

    def page_break(self):
        self.doc.add_page_break()

    def save(self, path: str):
        self.doc.save(path)
