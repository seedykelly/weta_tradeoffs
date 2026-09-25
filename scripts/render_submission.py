#!/usr/bin/env python3
"""Render all maintained QMDs without Quarto replacing earlier Word outputs."""
from pathlib import Path
import shutil
import re
from html import escape
from zipfile import ZipFile
from xml.dom import minidom
import subprocess
import tempfile

def preserve_figure_alt_text(document, source):
    """Preserve QMD descriptions when a Quarto version omits Word alt text."""
    descriptions = re.findall(r'^#\| fig-alt: "(.*)"$', source.read_text(), re.M)
    with ZipFile(document) as archive:
        entries = [(item, archive.read(item.filename)) for item in archive.infolist()]
    count = 0
    def describe(match):
        nonlocal count
        if count >= len(descriptions):
            raise RuntimeError("More Word figures than QMD descriptions")
        description = escape(descriptions[count], quote=True)
        count += 1
        tag = match.group()
        if re.search(r'\bdescr="[^"]*"', tag):
            return re.sub(r'\bdescr="[^"]*"', lambda _: 'descr="' + description + '"', tag)
        end = '/>' if tag.endswith('/>') else '>'
        return tag[:-len(end)] + ' descr="' + description + '"' + end
    modified = []
    for info, data in entries:
        if info.filename == 'word/document.xml':
            text = re.sub(r'<wp:docPr\b[^>]*>', describe, data.decode('utf-8'))
            dom = minidom.parseString(text)
            for paragraph in dom.getElementsByTagName('w:p'):
                caption = ''.join(node.firstChild.data if node.firstChild else ''
                                  for node in paragraph.getElementsByTagName('w:t'))
                is_caption = re.match(r'^(?:Table|Figure)\s+(?:S)?\d+[a-z]?[.:]', caption) or caption.startswith('Supplementary Figure S')
                if not is_caption and not paragraph.getElementsByTagName('wp:docPr'):
                    continue
                properties = next((node for node in paragraph.childNodes
                                   if node.nodeType == node.ELEMENT_NODE and node.tagName == 'w:pPr'), None)
                if properties is None:
                    properties = dom.createElement('w:pPr')
                    paragraph.insertBefore(properties, paragraph.firstChild)
                below_figure = bool(re.match(r'^Figure\s+\d+:', caption))
                is_picture = bool(paragraph.getElementsByTagName('wp:docPr'))
                keep_next = (source.name == 'manuscript.qmd') if is_picture else not below_figure
                for setting in ('w:keepNext', 'w:keepLines'):
                    value = '1' if setting == 'w:keepLines' or keep_next else '0'
                    existing = properties.getElementsByTagName(setting)
                    if existing:
                        existing[0].setAttribute('w:val', value)
                    else:
                        element = dom.createElement(setting)
                        element.setAttribute('w:val', value)
                        properties.appendChild(element)
            data = dom.toxml(encoding='utf-8')
        modified.append((info, data))
    if count != len(descriptions):
        raise RuntimeError("Word figure count does not match QMD descriptions")
    with ZipFile(document, 'w') as archive:
        for info, data in modified:
            archive.writestr(info, data)

root = Path(__file__).resolve().parents[1]
quarto = shutil.which("quarto")
if quarto is None:
    raise SystemExit("Quarto is required to render the submission documents.")
outputs = root / "submission-review"
with tempfile.TemporaryDirectory(prefix=".submission-build-", dir=root) as scratch:
    build = Path(scratch)
    documents = ("manuscript", "supplement", "title-page")
    for name in documents:
        # Standalone Quarto rendering can clear an output directory. A distinct
        # directory for every source prevents deletion of earlier documents.
        relative = str((build / name).relative_to(root))
        subprocess.run([quarto, "render", name + ".qmd", "--to", "docx",
                        "--output-dir", relative], cwd=root, check=True)
        if not (build / name / (name + ".docx")).is_file():
            raise RuntimeError("Expected document not generated: " + name)
    for name in documents:
        preserve_figure_alt_text(build / name / (name + ".docx"), root / (name + ".qmd"))
    outputs.mkdir(exist_ok=True)
    for name in documents:
        shutil.copy2(build / name / (name + ".docx"), outputs / (name + ".docx"))
print("Rendered all three documents in", outputs)
