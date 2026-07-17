"""为 LASSO.pdf 提取目录并添加书签（精确跳转至标题位置）"""
import fitz
import sys
import io
import re

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

pdf_path = r"04-03-bayesian-signal-process-uprng.pdf"
output_path = r"04-03-bayesian-signal-process-uprng-bookmarked.pdf"

doc = fitz.open(pdf_path)

# LASSO.pdf 字号 → 标题层级映射
SIZE_TO_LEVEL = {
    21.0: 1,   # H1 讲标题
    15.8: 2,   # H2 章标题
    12.3: 3,   # H3 节标题
    10.5: 4,   # H4 小节标题
}

all_headings = []

# ---- 第一步：提取所有标题及其精确位置（来自 extract_toc.py）----
for page_num in range(doc.page_count):
    page = doc[page_num]
    page_rect = page.rect
    blocks = page.get_text("dict")["blocks"]

    page_spans = []

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
        # H4 必须匹配三级编号 N.N.N，排除表格中的伪标题（如 "3.1 惩罚项角度"）
        if level == 4 and not re.match(r'^[\d]+\.[\d]+\.[\d]+', text):
            continue
        all_headings.append((page, level, text, bbox))

# ---- 第二步：构建 set_toc 书签列表 ----
toc_entries = []

for page, level, text, (x0, y0, x1, y1) in all_headings:
    dest = {
        "kind": fitz.LINK_GOTO,
        "page": page - 1,
        "to": fitz.Point(x0, y0),
    }
    toc_entries.append([level, text, page, dest])

doc.set_toc(toc_entries)
doc.save(output_path, garbage=4, deflate=True)
doc.close()

# ---- 第三步：输出结果 ----
print(f"原文件: {pdf_path}")
print(f"输出文件: {output_path}")
print(f"书签总数: {len(toc_entries)}")
print(f"  H1: {sum(1 for e in toc_entries if e[0]==1)}")
print(f"  H2: {sum(1 for e in toc_entries if e[0]==2)}")
print(f"  H3: {sum(1 for e in toc_entries if e[0]==3)}")
print(f"  H4: {sum(1 for e in toc_entries if e[0]==4)}")
print(f"\n每个书签均包含精确的页面坐标，点击后可跳转到标题所在位置。")
