#!/usr/bin/env python3
from __future__ import annotations

import re
import sys
import posixpath
from pathlib import Path
from zipfile import ZipFile
from xml.etree import ElementTree as ET


REL_NS = {"r": "http://schemas.openxmlformats.org/package/2006/relationships"}


def main() -> int:
    if len(sys.argv) != 2:
        print("Usage: qa_pptx_package.py deck.pptx", file=sys.stderr)
        return 2

    pptx_path = Path(sys.argv[1])
    if not pptx_path.exists():
        print(f"Missing PPTX: {pptx_path}", file=sys.stderr)
        return 2

    failures: list[str] = []
    with ZipFile(pptx_path, "r") as zf:
        names = zf.namelist()
        name_set = set(names)
        slide_xml = sorted(
            [n for n in names if re.match(r"^ppt/slides/slide\d+\.xml$", n)],
            key=lambda n: int(re.search(r"slide(\d+)\.xml", n).group(1)),
        )
        media = [n for n in names if n.startswith("ppt/media/") and not n.endswith("/")]

        if len(slide_xml) != 17:
            failures.append(f"expected 17 slides, found {len(slide_xml)}")

        bin_media = [n for n in media if n.lower().endswith(".bin")]
        if bin_media:
            failures.append(f"found .bin media: {bin_media}")

        non_png_media = [n for n in media if not n.lower().endswith(".png")]
        if non_png_media:
            failures.append(f"found non-PNG media: {non_png_media}")

        external = []
        missing = []
        for rel_name in [n for n in names if n.endswith(".rels")]:
            root = ET.fromstring(zf.read(rel_name))
            rel_base = Path(rel_name).parent
            # Convert from e.g. ppt/slides/_rels/slide1.xml.rels to ppt/slides/
            if rel_name == "_rels/.rels":
                source_dir = ""
            else:
                source_dir = str(rel_base).replace("/_rels", "")
            for rel in root.findall("r:Relationship", REL_NS):
                target = rel.attrib.get("Target", "")
                mode = rel.attrib.get("TargetMode", "")
                if mode == "External":
                    external.append((rel_name, target))
                    continue
                if not target or target.startswith("#"):
                    continue
                if target.startswith("/"):
                    normalized = posixpath.normpath(target.lstrip("/"))
                else:
                    normalized = posixpath.normpath(posixpath.join(source_dir, target))
                if normalized and normalized not in name_set:
                    # Some package-level relationship targets are outside ppt/ and still valid.
                    if not normalized.startswith(("ppt/", "docProps/", "_rels/", "customXml/")):
                        continue
                    missing.append((rel_name, target, normalized))

        if external:
            failures.append(f"found external relationships: {external}")
        if missing:
            failures.append(f"found missing relationship targets: {missing[:10]}")

        xml_text = "\n".join(
            zf.read(n).decode("utf-8", errors="ignore")
            for n in names
            if n.endswith(".xml")
        )
        forbidden = [
            "How to say it live",
            "Do not compare denominators",
            "ASB-versus-UTI",
            "ASB vs UTI",
            "significant lpf",
            "global AMR association",
            "stable VF proves",
        ]
        for phrase in forbidden:
            if phrase.lower() in xml_text.lower():
                failures.append(f"forbidden deck phrase present: {phrase}")

    print(f"PPTX: {pptx_path}")
    print(f"slides: {len(slide_xml)}")
    print(f"media: {len(media)}")
    print("media_types:", ", ".join(sorted({Path(n).suffix.lower() for n in media})) or "none")
    if failures:
        print("FAIL")
        for failure in failures:
            print(f"- {failure}")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
