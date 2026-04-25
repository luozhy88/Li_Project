import markdown
import subprocess
import os

md_path = "/Users/apple/Desktop/kimi/litao/test1/针灸穴位-基因关联数据库构建可行性分析报告.md"
html_path = "/Users/apple/Desktop/kimi/litao/test1/report.html"
pdf_path = "/Users/apple/Desktop/kimi/litao/test1/针灸穴位-基因关联数据库构建可行性分析报告.pdf"

with open(md_path, "r", encoding="utf-8") as f:
    md_text = f.read()

html_body = markdown.markdown(md_text, extensions=["tables", "fenced_code", "toc"])

html_full = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8">
<title>针灸穴位-基因关联数据库构建可行性分析报告</title>
<style>
  @page {{ size: A4; margin: 2cm; }}
  body {{
    font-family: "PingFang SC", "Hiragino Sans GB", "Microsoft YaHei", "WenQuanYi Micro Hei", sans-serif;
    font-size: 11pt;
    line-height: 1.7;
    color: #333;
    max-width: 100%;
  }}
  h1 {{
    font-size: 20pt;
    color: #1a5276;
    border-bottom: 3px solid #1a5276;
    padding-bottom: 8px;
    margin-top: 30px;
  }}
  h2 {{
    font-size: 15pt;
    color: #2874a6;
    border-bottom: 1px solid #aed6f1;
    padding-bottom: 5px;
    margin-top: 25px;
  }}
  h3 {{
    font-size: 13pt;
    color: #2e86c1;
    margin-top: 20px;
  }}
  h4 {{
    font-size: 12pt;
    color: #5dade2;
    margin-top: 15px;
  }}
  table {{
    border-collapse: collapse;
    width: 100%;
    margin: 15px 0;
    font-size: 10pt;
  }}
  th, td {{
    border: 1px solid #bbb;
    padding: 6px 8px;
    text-align: left;
    vertical-align: top;
  }}
  th {{
    background-color: #eaf2f8;
    font-weight: bold;
    color: #1a5276;
  }}
  tr:nth-child(even) {{
    background-color: #f8f9f9;
  }}
  code {{
    background-color: #f4f4f4;
    padding: 2px 5px;
    border-radius: 3px;
    font-family: "Courier New", monospace;
    font-size: 10pt;
  }}
  pre {{
    background-color: #f4f4f4;
    padding: 12px;
    border-radius: 5px;
    overflow-x: auto;
    font-size: 9.5pt;
    line-height: 1.4;
  }}
  blockquote {{
    border-left: 4px solid #5dade2;
    margin: 10px 0;
    padding: 8px 15px;
    background-color: #eaf2f8;
    color: #1a5276;
  }}
  ul, ol {{
    margin: 8px 0;
    padding-left: 25px;
  }}
  li {{
    margin: 4px 0;
  }}
  hr {{
    border: none;
    border-top: 1px solid #ddd;
    margin: 20px 0;
  }}
  strong {{
    color: #1a5276;
  }}
  p {{
    margin: 8px 0;
  }}
</style>
</head>
<body>
{html_body}
</body>
</html>"""

with open(html_path, "w", encoding="utf-8") as f:
    f.write(html_full)

chrome_path = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
cmd = [
    chrome_path,
    "--headless",
    "--disable-gpu",
    "--no-pdf-header-footer",
    "--print-to-pdf-no-header",
    "--run-all-compositor-stages-before-draw",
    f"--print-to-pdf={pdf_path}",
    f"file://{html_path}"
]

result = subprocess.run(cmd, capture_output=True, text=True)
if result.returncode != 0:
    print("STDERR:", result.stderr)
    print("STDOUT:", result.stdout)
else:
    print(f"PDF generated successfully: {pdf_path}")
    print(f"File size: {os.path.getsize(pdf_path)} bytes")
