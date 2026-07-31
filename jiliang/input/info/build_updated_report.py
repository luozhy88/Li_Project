#!/usr/bin/env python3
"""Append the info-folder reference material to the supplied Word report."""

from __future__ import annotations

import html
import os
import shutil
import tempfile
import zipfile
from pathlib import Path


SOURCE = Path("/Users/apple/Desktop/github.local/Li_Project/jiliang/output/report/False_Localizing_Sign_文献计量分析报告.docx")
INFO = Path("/Users/apple/Desktop/github.local/Li_Project/jiliang/input/info")
DESTINATION = Path("/Users/apple/Desktop/github.local/Li_Project/jiliang/input/info/False_Localizing_Sign_文献计量分析报告_含info材料补充版.docx")

IMAGES = [
    ("50ce9df587dde4e95152e0087a36f24a.jpg", "图 B1  原始文献中的参考文献分析图题截取。"),
    ("b498c0cefb30543535e8bf3a8bf8a321.jpg", "图 B2  拔毛癖研究的高被引文献（前 10 位）示例。"),
    ("3e420714a25417eb1ba3753577a76152.jpg", "图 B3  参考文献共被引网络的聚类可视化示例。"),
    ("5351faf8f800c6345e3e67d088150a94.png", "图 B4  拔毛癖研究的高产期刊（前 10 位）示例。"),
    ("e4734bba5eecc72e4f4c9f1fbc14f931.jpg", "图 B5  拔毛癖研究的高产作者（前 10 位）示例。"),
    ("446138e9a4f53277ff509aad6397de4d.jpg", "图 B6  参考文献突现（citation bursts）分析示例。"),
]

P = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"


def paragraph(text: str, *, bold: bool = False, size: int | None = None, center: bool = False) -> str:
    properties = ""
    if bold or size:
        properties = "<w:rPr>" + ("<w:b/>" if bold else "") + (f'<w:sz w:val="{size}"/>' if size else "") + "</w:rPr>"
    align = '<w:pPr><w:jc w:val="center"/></w:pPr>' if center else ""
    return f"<w:p>{align}<w:r>{properties}<w:t>{html.escape(text)}</w:t></w:r></w:p>"


def image_paragraph(rel_id: str, width: int, height: int, name: str) -> str:
    # EMUs: scale every image to a maximum width of 6.1 inches.
    max_width = 5577840
    if width > max_width:
        height = round(height * max_width / width)
        width = max_width
    return f'''<w:p><w:pPr><w:jc w:val="center"/></w:pPr><w:r><w:drawing><wp:inline distT="0" distB="0" distL="0" distR="0" xmlns:wp="http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"><wp:extent cx="{width}" cy="{height}"/><wp:docPr id="{rel_id[3:]}" name="{name}"/><a:graphic xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main"><a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture"><pic:pic xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture"><pic:nvPicPr><pic:cNvPr id="0" name="{name}"/><pic:cNvPicPr/></pic:nvPicPr><pic:blipFill><a:blip r:embed="{rel_id}"/><a:stretch><a:fillRect/></a:stretch></pic:blipFill><pic:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="{width}" cy="{height}"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr></pic:pic></a:graphicData></a:graphic></wp:inline></w:drawing></w:r></w:p>'''


def image_size(path: Path) -> tuple[int, int]:
    # JPEG / PNG dimensions are read without third-party dependencies.
    data = path.read_bytes()
    if data.startswith(b"\x89PNG"):
        return int.from_bytes(data[16:20], "big"), int.from_bytes(data[20:24], "big")
    i = 2
    while i < len(data):
        if data[i] != 0xFF:
            i += 1
            continue
        marker = data[i + 1]
        length = int.from_bytes(data[i + 2:i + 4], "big")
        if marker in range(0xC0, 0xC4):
            return int.from_bytes(data[i + 5:i + 7], "big"), int.from_bytes(data[i + 7:i + 9], "big")
        i += 2 + length
    raise ValueError(f"Unsupported image: {path}")


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        temp = Path(tmp)
        with zipfile.ZipFile(SOURCE) as archive:
            archive.extractall(temp)

        rels_path = temp / "word/_rels/document.xml.rels"
        rels = rels_path.read_text(encoding="utf-8")
        existing_ids = [int(value.split('"')[0]) for value in rels.split('Id="rId')[1:]]
        next_id = max(existing_ids) + 1
        drawing_parts: list[str] = []
        for index, (filename, caption) in enumerate(IMAGES, start=1):
            image_path = INFO / filename
            target_name = f"info_material_{index}{image_path.suffix.lower()}"
            shutil.copy2(image_path, temp / "word/media" / target_name)
            rel_id = f"rId{next_id}"
            next_id += 1
            rels = rels.replace("</Relationships>", f'<Relationship Id="{rel_id}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="media/{target_name}"/></Relationships>')
            px_w, px_h = image_size(image_path)
            # 96 dpi gives a clear but compact page layout.
            drawing_parts.append(image_paragraph(rel_id, px_w * 9525, px_h * 9525, target_name))
            drawing_parts.append(paragraph(caption, center=True))
        rels_path.write_text(rels, encoding="utf-8")

        types_path = temp / "[Content_Types].xml"
        types = types_path.read_text(encoding="utf-8")
        if 'Extension="jpg"' not in types:
            types = types.replace("</Types>", '<Default Extension="jpg" ContentType="image/jpeg"/></Types>')
        types_path.write_text(types, encoding="utf-8")

        appendix = [
            '<w:p><w:r><w:br w:type="page"/></w:r></w:p>',
            paragraph("附录 B  信息文件夹材料：可借鉴的文献计量分析框架", bold=True, size=32),
            paragraph("材料说明", bold=True, size=26),
            paragraph("本附录整合 info 文件夹中的 6 张图像材料。图像出自《Trichotillomania: A perspective from bibliometric analysis》所示的拔毛癖文献计量分析；其研究对象与本报告的假定位体征不同，故其中的发文量、作者、期刊及被引数据不作为本研究结果，也不与本研究样本合并。它们在本文中的作用是为假定位体征研究后续深化分析提供图表表达与方法学参考。"),
            paragraph("可迁移的方法学启示", bold=True, size=26),
            paragraph("（1）共被引网络与聚类：可用于识别假定位体征领域的知识基础、经典文献及主题群。"),
            paragraph("（2）高被引文献、期刊与作者清单：可补充当前报告的发文量分析，呈现学术影响力与核心传播渠道。"),
            paragraph("（3）引文突现分析：可识别某一时期迅速受到关注的文献或议题，从而辅助判断研究前沿。"),
            paragraph("建议后续从 Web of Science 导出包含 Cited References 的完整记录，并以 CiteSpace 或 VOSviewer 开展共被引、聚类、时间线与突现检测；分析时应保持检索式、时间范围、去重规则和阈值参数的可复现记录。"),
            paragraph("原始图像材料", bold=True, size=26),
            *drawing_parts,
            paragraph("注：图 B1—B6 保留原始英文图像，以便核对来源信息。图中涉及的疾病、期刊、作者与引文指标均属于拔毛癖研究示例，不应解读为假定位体征领域的统计结果。"),
        ]
        document_path = temp / "word/document.xml"
        document = document_path.read_text(encoding="utf-8")
        document = document.replace("</w:body>", "".join(appendix) + "</w:body>")
        document_path.write_text(document, encoding="utf-8")

        with zipfile.ZipFile(DESTINATION, "w", zipfile.ZIP_DEFLATED) as output:
            for file in temp.rglob("*"):
                if file.is_file():
                    output.write(file, file.relative_to(temp))
    print(DESTINATION)


if __name__ == "__main__":
    main()
