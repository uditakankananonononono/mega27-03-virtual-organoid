"""Transparent rendered-PDF body-page inventory; diagnostics, not a page-rule substitute."""
from __future__ import annotations
import json
import re
from pathlib import Path
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
PDF = ROOT / 'paper' / 'mega27-03-virtual-organoid-paper.pdf'
reader = PdfReader(PDF)
rows = []
for index, page in enumerate(reader.pages, 1):
    text = page.extract_text() or ''
    words = re.findall(r"\b[\w-]+\b", text)
    rows.append({'page': index, 'word_tokens': len(words),
                 'contains_references_heading': bool(re.search(r'(?m)^References\s*$', text)),
                 'contains_appendix_heading': bool(re.search(r'(?m)^Appendix\b', text)),
                 'images': len(page.images),
                 'excerpt': ' '.join(text.split())[:145]})
assert len(rows) == 91, 'Unexpected physical page count; re-establish body boundary'
assert rows[61]['contains_references_heading'], 'References did not begin on p62'
body = rows[1:61]  # pp2-61; title/abstract p1 and references p62 onward excluded
assert len(body) == 60
result = {
    'source_pdf': str(PDF.relative_to(ROOT)),
    'sha256': __import__('hashlib').sha256(PDF.read_bytes()).hexdigest(),
    'physical_pages': len(rows),
    'body_candidate_pages': '2-61',
    'candidate_count_before_qualitative_exclusions': len(body),
    'diagnostic_counts_not_a_rule': {str(n): sum(p['word_tokens'] >= n for p in body) for n in (200, 300, 350, 400)},
    'body_candidate_word_tokens': sum(p['word_tokens'] for p in body),
    'per_page': rows,
    'caveat': 'Word-token thresholds are diagnostics only. A heading/table/figure-heavy page is not automatically a substantive text-body page; manual inspection decides the user page gate.'
}
path = ROOT / 'results' / 'rendered_body_page_inventory.json'
path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n')
print(json.dumps({k:v for k,v in result.items() if k != 'per_page'}, indent=2))
