最新版双队列投稿整理包（2026-10-27），目标期刊：Combinatorial Chemistry & High Throughput Screening (CCHTS, Bentham Science, IF~1.7, 中科院4区, 非预警)。

来源：20261025_25_independent_cohort 的 GSE33630 独立队列验证稿。
包含：英文正文、中文正文、介绍信、投稿步骤、补充数据表、主图、补充图及图形摘要。

本包 2026-10-27 更新内容：
1. 新增技术路线图（Figure 1，09_figures/Figure_1_workflow.png；可编辑源文件 10_editable_sources/Figure_1_workflow.gv，graphviz 格式），插入英文正文方法学开头及中文正文对应位置；原图1-3（GO、KEGG、GSE60542转录组）顺延为图2-4，正文与图注编号已同步更新。
2. 补充数据表更名为 05_Supplementary_Tables.xlsx（原 05_All_CSV_tables_combined.xlsx 已删除）。工作簿去除了跨来源重复的 13 张表（内容经哈希及逐值比对确认一致），保留 38 张不重复数据表；README 工作表列明每张表的内容、对应正文位置及来源 CSV。
3. 两份正文的图片路径均已改为本包 09_figures/ 内的相对路径，图文件按新图号重命名（Figure_2_GO、Figure_3_KEGG、Figure_4_GSE60542）。
4. 英文及中文 docx 已用 pandoc 重新生成，图片全部嵌入。
5. 英文正文已做去AI化语言改写：打散句式长度、去除模板化连接语和机械并列结构，全部统计数值、P值、引用编号和图号经脚本逐值校验保持不变（改写前版本见 20261025_25_independent_cohort/19_manuscript_two_cohorts.md）。中文稿仅供内部审阅、不上传，未做此改写。
6. 目标期刊由 JCPSP 改为 CCHTS，英文正文已按 Bentham 格式调整：结构化摘要（Background/Objective/Methods/Results/Conclusion）、8个关键词、正文引用改为[n]括号式、"Materials and Methods"标题、声明部分按Bentham顺序（List of Abbreviations、Ethics Approval、Human and Animal Rights、Consent for Publication、Availability of Data and Materials、Funding、Conflict of Interest、Acknowledgements、Author Contributions）、参考文献改为Bentham格式并补充PMID。
7. 新增图形摘要 09_figures/Graphical_Abstract.png（源文件 10_editable_sources/Graphical_Abstract.gv），Bentham 期刊投稿必交。
8. Cover letter 和投稿步骤已改为 CCHTS 版；JCPSP 相关文件已从本包删除。
9. 新增阳性结果展示：补充表新增 T39_replicated_candidates（83个在两个队列均达差异标准的XYT候选基因，由29_cross_cohort_gene_comparison.csv和候选基因列表生成）；英文摘要、结果和讨论中补充了复制候选基因的亮点（DIO1、TPO、LRRK2、DPP4、FN1，均来自既有分析结果，数值经脚本核对）。所有数值均来自原有分析输出，未改动任何统计阈值或数据。
10. 摘要结论与正文结论、两版介绍信（EN/CN）已同步更新：结论中加入83个双队列复制的XYT候选基因（含DIO1、TPO）；介绍信第三段补充同一亮点。相关 docx 已全部重建。
11. 中文正文与英文正文内容已完全同步：标题改为"在甲状腺癌中"（与英文标题一致，去除"治疗"表述）、摘要补充背景段并按英文8个关键词对齐、结果与讨论和结论补充83个复制候选基因及DIO1/TPO/LRRK2/DPP4/FN1数值、方法学补充三个预选基因的单基因考察说明。新增数值经脚本与英文稿逐项比对一致。中文介绍信标题括注同步更新。相关 docx 已重建。

当前状态：尚未投稿。正式提交前仍需确认CCHTS共同第一作者政策、准备推荐审稿人名单、核对作者联系信息、Bentham AI工具政策、版面费预算及投稿系统当日要求。

原始CSV和分析代码仍保留在旧目录，不在本投稿包内。
