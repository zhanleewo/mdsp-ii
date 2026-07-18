"""从 mdsp-ii.pdf 提取目录并添加书签（精确跳转至标题位置）"""
import fitz
import sys
import io
import re

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

pdf_path = r"mdsp-ii.pdf"
output_path = r"现代数字信号处理 II 笔记.pdf"

doc = fitz.open(pdf_path)

SIZE_TO_LEVEL = {
    29.6: 1,   # H1: 第a单元
    18.5: 2,   # H2: 第b讲
    13.9: 3,   # H3: c.
    10.8: 4,   # H4: c.d
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
        if level == 4 and not re.match(r'^[\d]+\.[\d]+', text):
            continue
        all_headings.append((page, level, text, bbox))

# ---- 第二步：构建 set_toc 书签列表 ----
toc_entries = []

# 插入虚拟 H1 作为引言等前置内容的父级
first_page = all_headings[0][0]
first_x0, first_y0 = all_headings[0][3][0], all_headings[0][3][1]
toc_entries.append([1, "课程导论 / 前言", first_page, {
    "kind": fitz.LINK_GOTO,
    "page": first_page - 1,
    "to": fitz.Point(first_x0, first_y0),
}])

parent_stack = [1]  # 虚拟 H1 已加入

for page, level, text, (x0, y0, x1, y1) in all_headings:
    dest = {
        "kind": fitz.LINK_GOTO,
        "page": page - 1,
        "to": fitz.Point(x0, y0),
    }

    # 弹出所有 >= 当前 level 的父级（同级或更深）
    while parent_stack and parent_stack[-1] >= level:
        parent_stack.pop()

    # 插入缺失的中间层级（虚拟父级）
    parent_level = parent_stack[-1] if parent_stack else 0
    for missing_level in range(parent_level + 1, level):
        toc_entries.append([missing_level, "(节点)", page, dest])
        parent_stack.append(missing_level)

    toc_entries.append([level, text, page, dest])
    parent_stack.append(level)

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
