#!/usr/bin/env python3
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
from xml.etree import ElementTree as ET
import io
import re

SRC = Path('/Users/Aamir/Desktop/rUTIs/outputs/manual-20260527-current-review/presentations/ruti-longitudinal-vf-pipeline-review-v3/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-27_REPAIRED.pptx')
OUT = Path('/Users/Aamir/Desktop/rUTIs/outputs/manual-20260527-current-review/presentations/ruti-longitudinal-vf-pipeline-review-v3/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-27_REPAIRED_IDS.pptx')

ET.register_namespace('a', 'http://schemas.openxmlformats.org/drawingml/2006/main')
ET.register_namespace('r', 'http://schemas.openxmlformats.org/officeDocument/2006/relationships')
ET.register_namespace('p', 'http://schemas.openxmlformats.org/presentationml/2006/main')
ET.register_namespace('p14', 'http://schemas.microsoft.com/office/powerpoint/2010/main')
ET.register_namespace('p15', 'http://schemas.microsoft.com/office/powerpoint/2012/main')
ET.register_namespace('mc', 'http://schemas.openxmlformats.org/markup-compatibility/2006')
ET.register_namespace('a16', 'http://schemas.microsoft.com/office/drawing/2014/main')
ET.register_namespace('adec', 'http://schemas.microsoft.com/office/drawing/2017/decorative')
ET.register_namespace('asvg', 'http://schemas.microsoft.com/office/drawing/2016/SVG/main')
ET.register_namespace('ma14', 'http://schemas.microsoft.com/office/mac/powerpoint/2011/main')

slide_re = re.compile(r'^ppt/slides/slide\d+\.xml$')

def patch_slide(data: bytes) -> bytes:
    root = ET.fromstring(data)
    n = 1
    for el in root.iter():
        if el.tag.endswith('}cNvPr'):
            el.set('id', str(n))
            n += 1
    return ET.tostring(root, encoding='utf-8', xml_declaration=True)

with ZipFile(SRC, 'r') as zin, ZipFile(OUT, 'w', ZIP_DEFLATED) as zout:
    for info in zin.infolist():
        data = zin.read(info.filename)
        if slide_re.match(info.filename):
            data = patch_slide(data)
        zout.writestr(info, data)

print(OUT)
