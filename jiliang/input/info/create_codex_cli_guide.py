from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
from reportlab.platypus import (
    KeepTogether, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table,
    TableStyle,
)
from reportlab.pdfbase.pdfmetrics import stringWidth

OUT = "output/pdf/Codex_终端使用快速指南.pdf"
pdfmetrics.registerFont(UnicodeCIDFont("STSong-Light"))

NAVY = colors.HexColor("#17365D")
BLUE = colors.HexColor("#2F75B5")
LIGHT_BLUE = colors.HexColor("#EAF2F8")
TEXT = colors.HexColor("#1F2937")
MUTED = colors.HexColor("#5B6573")

styles = getSampleStyleSheet()
styles.add(ParagraphStyle(name="TitleCN", fontName="STSong-Light", fontSize=27, leading=38,
                          textColor=NAVY, alignment=TA_CENTER, spaceAfter=14))
styles.add(ParagraphStyle(name="SubtitleCN", fontName="STSong-Light", fontSize=13, leading=21,
                          textColor=MUTED, alignment=TA_CENTER))
styles.add(ParagraphStyle(name="H1CN", fontName="STSong-Light", fontSize=18, leading=28,
                          textColor=NAVY, spaceBefore=10, spaceAfter=10))
styles.add(ParagraphStyle(name="H2CN", fontName="STSong-Light", fontSize=13, leading=20,
                          textColor=BLUE, spaceBefore=10, spaceAfter=5))
styles.add(ParagraphStyle(name="BodyCN", fontName="STSong-Light", fontSize=10.5, leading=18,
                          textColor=TEXT, spaceAfter=6))
styles.add(ParagraphStyle(name="SmallCN", fontName="STSong-Light", fontSize=8.8, leading=14,
                          textColor=MUTED))
styles.add(ParagraphStyle(name="CodeCN", fontName="Courier", fontSize=8.5, leading=14,
                          textColor=colors.HexColor("#112A46"), backColor=colors.HexColor("#F2F6FA"),
                          borderColor=colors.HexColor("#D4E1EF"), borderWidth=.5, borderPadding=7,
                          spaceBefore=3, spaceAfter=9))

def p(text, style="BodyCN"):
    return Paragraph(text, styles[style])

def code(text):
    return Paragraph(text.replace("&", "&amp;").replace("\n", "<br/>"), styles["CodeCN"])

def footer(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor("#D9E2F0"))
    canvas.line(doc.leftMargin, 1.35*cm, A4[0]-doc.rightMargin, 1.35*cm)
    canvas.setFont("STSong-Light", 8)
    canvas.setFillColor(MUTED)
    canvas.drawString(doc.leftMargin, .88*cm, "Codex 终端使用快速指南")
    canvas.drawRightString(A4[0]-doc.rightMargin, .88*cm, f"第 {doc.page} 页")
    canvas.restoreState()

doc = SimpleDocTemplate(OUT, pagesize=A4, leftMargin=2.1*cm, rightMargin=2.1*cm,
                        topMargin=1.8*cm, bottomMargin=1.8*cm, title="Codex 终端使用快速指南")
story = []

story += [Spacer(1, 3.0*cm), p("Codex 终端使用快速指南", "TitleCN"),
          p("从安装到完成第一个代码任务", "SubtitleCN"), Spacer(1, .7*cm)]
intro = [[p("适用对象", "H2CN"), p("第一次在 macOS、Linux 或 Windows 终端中使用 Codex 的开发者。", "BodyCN")],
         [p("你将学会", "H2CN"), p("安装 Codex CLI、登录账号、在项目目录中发起任务，并理解常见的安全确认。", "BodyCN")],
         [p("核心理念", "H2CN"), p("把 Codex 当作会阅读项目、提出方案、修改文件并执行验证命令的协作式编程助手。", "BodyCN")]]
t = Table(intro, colWidths=[3.0*cm, 12.0*cm], hAlign="LEFT")
t.setStyle(TableStyle([("BACKGROUND", (0,0), (-1,-1), colors.white), ("GRID", (0,0), (-1,-1), .35, colors.HexColor("#D6E4F0")), ("BACKGROUND", (0,0), (0,-1), LIGHT_BLUE), ("VALIGN", (0,0), (-1,-1), "TOP"), ("LEFTPADDING", (0,0), (-1,-1), 10), ("RIGHTPADDING", (0,0), (-1,-1), 10), ("TOPPADDING", (0,0), (-1,-1), 8), ("BOTTOMPADDING", (0,0), (-1,-1), 8)]))
story += [t, Spacer(1, .8*cm), p("使用前，请确认你的电脑能够访问互联网，并已准备好可登录的 OpenAI / ChatGPT 账号。", "SmallCN"), PageBreak()]

story += [p("1. 安装 Codex CLI", "H1CN"), p("在终端中选择一种安装方式即可。macOS / Linux 推荐使用官方安装脚本，它不要求你预先安装 Node.js。", "BodyCN"), p("macOS / Linux", "H2CN"), code("curl -fsSL https://chatgpt.com/codex/install.sh | sh"), p("Windows PowerShell", "H2CN"), code('powershell -ExecutionPolicy ByPass -c "irm https://chatgpt.com/codex/install.ps1 | iex"'), p("也可选择包管理器", "H2CN"), code("# npm\nnpm install -g @openai/codex\n\n# macOS Homebrew\nbrew install --cask codex"), p("安装完成后，关闭并重新打开终端，检查是否成功：", "BodyCN"), code("codex --version"), p("提示：npm 安装时包名必须是 <font name=\"Courier\">@openai/codex</font>，不要使用 <font name=\"Courier\">npm install -g codex</font>。", "SmallCN")]

story += [p("2. 启动与登录", "H1CN"), p("先进入你要处理的项目目录，再启动 Codex。它会将当前目录作为工作范围，因此在正确的项目根目录启动很重要。", "BodyCN"), code("cd /你的/项目路径\ncodex"), p("首次启动会打开或提示你完成登录。按屏幕指引使用 OpenAI / ChatGPT 账号授权即可。登录完成后，直接在对话输入框写下任务。", "BodyCN"), p("如果只是希望执行一条明确任务，也可以在命令后直接附上文字：", "BodyCN"), code('codex "分析这个项目的结构，并告诉我如何运行它"')]

section3 = [p("3. 怎样提出高质量任务", "H1CN"), p("把目标、范围与验收标准说清楚。复杂任务可要求它先分析再修改；涉及生产环境、密钥或删除数据时，应额外明确边界。", "BodyCN")]
examples = [[p("场景", "H2CN"), p("可直接输入的示例", "H2CN")],
            [p("了解项目"), p("分析这个项目的目录结构、启动方式和测试命令；不要修改文件。")],
            [p("修复问题"), p("修复当前测试失败的问题。先说明原因，再修改并运行相关测试验证。")],
            [p("新增功能"), p("为用户列表增加按名称搜索。复用现有组件风格，并补充必要测试。")],
            [p("文档维护"), p("根据当前代码更新 README 的安装与启动说明，使用中文。")]]
t = Table(examples, colWidths=[3.0*cm, 12.0*cm], repeatRows=1)
t.setStyle(TableStyle([("BACKGROUND", (0,0), (-1,0), LIGHT_BLUE), ("GRID", (0,0), (-1,-1), .35, colors.HexColor("#D6E4F0")), ("VALIGN", (0,0), (-1,-1), "TOP"), ("LEFTPADDING", (0,0), (-1,-1), 9), ("RIGHTPADDING", (0,0), (-1,-1), 9), ("TOPPADDING", (0,0), (-1,-1), 7), ("BOTTOMPADDING", (0,0), (-1,-1), 7)]))
section3 += [t, Spacer(1, .4*cm), p("一个很实用的句式：<b>“完成 X；只改动 Y；不要改动 Z；最后运行 A 进行验证。”</b>", "BodyCN")]
story += [KeepTogether(section3)]

section4 = [p("4. 使用过程中的确认与命令", "H1CN"), p("Codex 会在需要运行命令、修改文件或执行可能影响系统的操作时展示计划或请求确认。请阅读命令、路径和改动范围后再批准。尤其注意删除文件、安装依赖、推送代码或访问外部服务等操作。", "BodyCN"), p("会话内常用命令", "H2CN")]
commands = [[p("命令", "H2CN"), p("用途", "H2CN")], [p("/help"), p("查看会话内命令与快捷键。")], [p("/exit"), p("退出当前会话。")], [p("/feedback"), p("提交体验或问题反馈。")]]
t = Table(commands, colWidths=[3.0*cm, 12.0*cm], repeatRows=1)
t.setStyle(TableStyle([("BACKGROUND", (0,0), (-1,0), LIGHT_BLUE), ("GRID", (0,0), (-1,-1), .35, colors.HexColor("#D6E4F0")), ("VALIGN", (0,0), (-1,-1), "TOP"), ("LEFTPADDING", (0,0), (-1,-1), 9), ("RIGHTPADDING", (0,0), (-1,-1), 9), ("TOPPADDING", (0,0), (-1,-1), 7), ("BOTTOMPADDING", (0,0), (-1,-1), 7)]))
section4 += [t]
story += [KeepTogether(section4)]

story += [KeepTogether([p("5. 推荐工作流程", "H1CN"), p("1. 进入仓库根目录，运行 <font name=\"Courier\">codex</font>。<br/>2. 说明目标和限制，并让它先分析项目。<br/>3. 审阅它给出的计划与可能的文件改动。<br/>4. 允许在可接受范围内的修改与测试。<br/>5. 最后查看改动内容、测试结果和未解决事项。", "BodyCN"), p("最佳实践", "H2CN"), p("• 在 Git 仓库中使用，便于随时查看差异和回退。<br/>• 小步提交任务；每次只聚焦一个清晰目标。<br/>• 不要在提示中粘贴 API 密钥、密码或其他敏感信息。<br/>• 对自动生成的代码仍应自行审阅，特别是安全、权限、数据库与支付逻辑。", "BodyCN")])]

story += [KeepTogether([p("6. 常见问题", "H1CN"), p("<b>找不到 codex 命令</b><br/>重新打开终端后再试；若仍无效，检查安装程序输出的 PATH 提示，或改用官方安装方式重新安装。", "BodyCN"), p("<b>npm 安装报权限错误</b><br/>不要用 sudo 强行全局安装。优先使用官方安装脚本，或修正 Node / npm 的全局目录权限配置。", "BodyCN"), p("<b>Codex 不该修改文件怎么办</b><br/>在任务开头明确写“只分析，不修改文件”；运行命令或编辑前仔细阅读确认信息。", "BodyCN"), p("<b>如何继续上一次工作</b><br/>回到同一项目目录并启动 Codex，根据界面提供的会话恢复选项继续；必要时在新会话中说明当前进度和已有改动。", "BodyCN"), Spacer(1, .3*cm), p("官方项目与最新安装说明：github.com/openai/codex", "SmallCN")])]

doc.build(story, onFirstPage=footer, onLaterPages=footer)
print(OUT)
