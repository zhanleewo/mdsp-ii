"""从 mdsp-ii.pdf 中提取 1/2/3/4 级标题，含精确位置(bbox)"""
import fitz
import sys
import io
import re

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

pdf_path = r"mdsp-ii.pdf"
doc = fitz.open(pdf_path)

SIZE_TO_LEVEL = {
    22.4: 1,
    14.0: 2,
    10.5: 3,
    8.2:  4,
}

all_headings = []

for page_num in range(doc.page_count):
    page = doc[page_num]
    page_rect = page.rect  # 页面尺寸
    blocks = page.get_text("dict")["blocks"]

    page_spans = []  # (level, size, text, bbox)

    for block in blocks:
        if "lines" not in block:
            continue
        for line in block["lines"]:
            for span in line["spans"]:
                size = round(span["size"], 1)
                text = span["text"].strip()
                if not text:
                    continue
                level = SIZE_TO_LEVEL.get(size)
                if level is not None:
                    page_spans.append((level, size, text, span["bbox"]))

    # 合并同一行内、同字号、同level的相邻span（只合并同行，不同行不合并）
    merged = []
    i = 0
    while i < len(page_spans):
        level, size, text, bbox = page_spans[i]
        x0, y0, x1, y1 = bbox
        j = i + 1
        while j < len(page_spans) and page_spans[j][0] == level and page_spans[j][1] == size:
            bx0, by0, bx1, by1 = page_spans[j][3]
            # 只在同一行内合并（y0 差距小于 6pt 视为同一行）
            if abs(by0 - y0) >= 6:
                break
            text += page_spans[j][2]
            x1 = max(x1, bx1)
            y1 = max(y1, by1)
            j += 1
        text = re.sub(r'\s+', ' ', text).strip()
        merged.append((page_num + 1, level, size, text, (x0, y0, x1, y1)))
        i = j

    # 过滤
    for page, level, size, text, bbox in merged:
        if len(text) < 2:
            continue
        if text.replace(" ", "").replace(".", "").replace("-", "").isdigit():
            continue
        if re.match(r'^[\s\.\-·•/\\]+$', text):
            continue
        if level == 4 and not re.match(r'^[\d]+\.[\d]+', text):
            continue
        all_headings.append((page, level, size, text, bbox))

# 输出
print(f"PDF 总页数: {doc.page_count} | 标题总数: {len(all_headings)}")
print(f"坐标单位: points (1pt = 1/72 inch), 原点: 页面左上角")
print("=" * 90)

for page, level, size, text, bbox in all_headings:
    x0, y0, x1, y1 = bbox
    indent = "  " * (level - 1)
    prefix = {1: "H1", 2: "H2", 3: "H3", 4: "H4"}.get(level, "??")
    # 输出: 页码, 位置(x0,y0)-(x1,y1), 标题
    print(f"[p.{page:>3}] ({x0:>5.0f},{y0:>5.0f})-({x1:>5.0f},{y1:>5.0f}) {indent}{prefix} {text}")

doc.close()
