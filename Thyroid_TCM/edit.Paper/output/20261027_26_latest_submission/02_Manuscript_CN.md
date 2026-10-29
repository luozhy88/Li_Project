# 消瘿汤核心药组的网络药理学及甲状腺乳头状癌转录组评估

短标题：消瘿汤核心药组与甲状腺乳头状癌的基因关联

Tao Li^1^\*, Jiacai Li^2^\*, Yajie Ding^3^\*, Zhiyu Luo^4^, Xing Zhang^5^

^1^Department of Acupuncture and Moxibustion, Affiliated Hospital of Hangzhou Normal University, Hangzhou, Zhejiang, China.

^2^Practical Teaching Center, Shanghai University of Medicine & Health Sciences, Shanghai, China.

^3^Department of Oncology, Longhua Hospital Affiliated to Shanghai University of Traditional Chinese Medicine, Shanghai, China.

^4^Jiangsu Sanshu Biotechnology Co., Ltd., No. 259 Shibei Gaoxin Road, Tangzha Town Subdistrict, Chongchuan District, Nantong, Jiangsu, China.

^5^Department of Respiratory Disease, Shuguang Hospital Affiliated to Shanghai University of Traditional Chinese Medicine, Shanghai 201203, China.

\*Tao Li, Jiacai Li, and Yajie Ding contributed equally to this work and share first authorship.

Correspondence: Xing Zhang, Department of Respiratory Disease, Shuguang Hospital Affiliated to Shanghai University of Traditional Chinese Medicine, Shanghai 201203, China. Tel: 021-20256375; Email: zhangxing2012@shutcm.edu.cn

Additional author contact details: Yajie Ding, Email: dingyajie2018@163.com.

中文审阅说明：本稿与当前英文投稿主稿同步，包含GSE60542和GSE33630两个队列。供作者理解与审阅，不作为英文期刊投稿附件。作者姓名、单位英文信息及参考文献保留英文；补充图集中附于文后，图号与英文稿一致。

# 摘要

**目的：**评估所选消瘿汤核心药组的关联基因是否在甲状腺乳头状癌差异表达基因中富集。


**方法：**将消瘿汤核心药组（XYT）、浙贝母（ZBM）和夏枯草（XKC）的关联基因与GeneCards疾病基因列表取交集，开展功能富集分析。分别采用配对limma模型分析GSE60542的24例患者和GSE33630的44例患者。差异表达标准为Benjamini–Hochberg校正P<0.05且|log2倍数变化|≥1。分别以全部可测基因和其中的疾病相关基因为背景，检验候选集合的富集情况。

**结果：**XYT、ZBM和XKC的疾病交集分别包含1,595、78和842个基因。两个队列各有20,848个可测基因，分别有666和745个达到差异表达标准。XYT和XKC在两个队列中均显著富集：采用疾病限制背景时，GSE60542的比值比分别为1.80和2.08，GSE33630分别为1.59和1.80，校正P均<0.001。初始666个差异基因中，536个在GSE33630再次达到相同标准且方向一致，其中包括83个XYT候选基因。ZBM富集未达统计学显著性。POLB、CACNA1H和CYP3A4在两个队列中均未达到联合差异表达标准。

**结论：**所选XYT和XKC基因集与肿瘤表达特征的关联可重复，为功能实验提供了候选基因；本研究未评估治疗反应。

**关键词：**消瘿汤；网络药理学；甲状腺乳头状癌；转录组；分子对接；中医药；差异表达；β-谷甾醇

# 引言

甲状腺癌包含分子特征和临床过程不同的恶性肿瘤类型。[1]手术、适用于部分患者的放射性碘治疗及甲状腺激素管理是分化型甲状腺癌治疗的基础。[2]持续性病变和治疗抵抗推动了其他治疗策略的研究。不同亚型的分子驱动因素存在差异，因此，从广泛疾病数据库中提出的候选机制需要在明确的肿瘤类型中加以考察。[3]

网络药理学可将中药成分与基因及生物学过程联系起来，形成可检验的研究假设。数据库覆盖范围和预测关联会影响候选基因集的构成。外部肿瘤表达数据可用于判断这些基因是否在疾病组织中发生改变，而其对中药治疗的响应仍需进一步研究。

本研究围绕甲状腺疾病，从消瘿汤中选取核心草药进行分析，全文以XYT表示这一核心药组。另对浙贝母（ZBM；Fritillariae Thunbergii Bulbus）和夏枯草（XKC；Prunellae Spica）进行单独分析，以描述其在整体注释特征中的贡献。我们比较成分关联基因的覆盖范围和功能注释，再在两个独立队列的配对PTC与正常甲状腺样本中考察预先确定的候选集合，并通过β-谷甾醇的探索性对接分析考察候选化合物–蛋白相互作用。

# 方法

本研究整合了中药–成分及成分–基因数据、功能富集分析、配对肿瘤表达谱和探索性分子对接，总体技术路线见图1。

![](09_figures/Figure_1_workflow.png)

**图1：**研究技术路线。TCMSP中药成分记录映射至PubChem化合物标识符（CID），来自BindingDB、BATMAN-TCM和IUPHAR/BPS Guide to PHARMACOLOGY的成分关联基因与12,987个GeneCards甲状腺肿瘤相关基因取交集，形成ZBM、XKC和XYT候选集合。下游分析包括GO/KEGG富集分析、两个独立PTC队列（GSE60542和GSE33630）的配对limma差异表达分析及候选集合富集检验与跨队列验证，以及β-谷甾醇与CYP3A4（6UNI）和CaV3.2（9AYG）的探索性分子对接。

中药成分来自TCMSP，通过小写化合物名称与同义词表匹配，关联至PubChem化合物标识符（CID）。[4, 5]成分关联基因来自BindingDB、BATMAN-TCM的已知及预测相互作用，以及IUPHAR/BPS Guide to PHARMACOLOGY。[6, 7, 8]将BindingDB中的UniProt蛋白标识符映射为对应的基因符号。保留上述来源中的已记录关联，并纳入预测评分>0.80的BATMAN-TCM预测记录。合并CID–基因记录后，在每个CID内对基因符号去重。TCMSP提供中药成分信息，最终保留的基因关联来自上述三个靶点资源。

保留的XYT数据库药物清单包括莪术、白芍、当归、柴胡、青皮、黄药子、夏枯草、皂角刺、牡丹皮、白术、茯苓、白毛夏枯草、玄参和浙贝母（补充表T02）。历史脚本使用中文药名模式匹配，其中“夏枯草”同时匹配了另一种药物白毛夏枯草（Ajugae Decumbentis Herba）。因此，我们检查了排除后者对疾病交集基因的影响。ZBM和XKC采用精确中文药名提取，使用相同的靶点关联流程。各组内合并重复的中药–成分记录，分别统计全部已匹配成分、关联到保留基因的成分及关联到疾病交集基因的成分的唯一CID数量。未按口服生物利用度或类药性筛选，也未对药物剂量、制剂成分及浓度建模；XYT指数据库分析集合，不代表成分经过测定的制剂。

甲状腺肿瘤相关基因来自以“Thyroid tumor”为查询词的归档GeneCards导出数据，包含12,987个唯一基因符号。[9]使用完整导出列表，未额外施加相关性评分阈值。成分关联基因符号按分号拆分、去重后与疾病列表取交集，形成ZBM、XKC和XYT候选集合。分析基于归档数据库导出文件；现有分析文件未记录准确的原始检索日期及数据库发布版本。

表达分析使用两个Affymetrix GPL570队列：来自比利时和法国的GSE60542，[10]以及来自乌克兰切尔诺贝利组织库的GSE33630。[11]分别分析GSE60542的fRMA和GSE33630的RMA标准化log2表达矩阵，未再次进行log转换或跨队列标准化。纳入所有符合标准的配对样本，样本量由可用样本决定。仅纳入配对的原发PTC与非肿瘤甲状腺组织。GSE60542中同一患者同类组织的重复组织块取平均。GSE33630根据样本名称建立配对，并用临床组织库编号核实；排除5个不配对肿瘤、1个不配对正常样本及11个未分化癌。两个队列的招募来源不同，没有重复的GEO样本标识符。排除基因注释缺失或不明确的探针，同一基因的多个探针在log2尺度取平均，每个队列均获得20,848个基因。使用R 4.3.1、limma 3.58.1及稳健经验贝叶斯方法，拟合配对模型`~ patient_id + tissue_group`。[12]正log2倍数变化表示肿瘤中表达升高。差异标准为Benjamini–Hochberg校正P<0.05且|log2FC|≥1。沿用原XYT、XKC和ZBM集合，采用单侧Fisher精确检验，分别以全部可测基因和其中10,655个疾病相关基因为背景。在各背景内，候选集合归属与差异表达状态构成2 × 2列联表，对照为背景中不属于候选集合的基因。每个队列内的6项检验统一校正P值。条件比值比同时报告双侧精确95%置信区间。重复结果要求再次达到两项差异表达标准且方向一致；另以Spearman相关及单独的方向一致性作描述性评估。GSE33630按放射暴露状态开展敏感性分析，采用相同配对模型和阈值，12项富集检验统一校正P值。对原始结构筛选关注的POLB、CACNA1H和CYP3A4另行考察，但保留全基因组P值校正。补充表T39按两个队列绝对log2FC的均值降序排列；正文展示的候选基因实例在分析后选定。

使用`bitr`和org.Hs.eg.db将唯一基因符号映射为Entrez标识符，并将返回的全部标识符用于clusterProfiler富集分析。[13]GO分析涵盖生物学过程（BP）、分子功能（MF）和细胞组分（CC），KEGG分析使用人类注释。采用上尾超几何检验和Benjamini–Hochberg校正，参数为`pvalueCutoff = 0.05`、`qvalueCutoff = 0.2`。默认注释背景分别包含18,888个BP基因、18,522个MF基因、19,894个CC基因及8,859个KEGG基因。报告校正P<0.05的条目。

GO图取各药物集合在每个本体中校正P值最低的三个条目的并集。KEGG图取各组校正P值最低的五条通路的并集，并根据研究主题加入甲状腺癌通路（hsa05216）。GeneRatio的分母为对应分析的已注释输入基因数，因此可随本体而不同。缺失的柱形表示该条目未出现在相应显著结果表中。图示用于描述性比较注释特征。

配体使用Open Babel制备，受体结构来自RCSB PDB。[14, 15]初步筛选采用AutoDock Vina，β-谷甾醇的补充对接采用smina。[16, 17]结构分析聚焦人源CYP3A4（6UNI）及apo状态的CaV3.2/CACNA1H（9AYG）。[15, 18]4P5A为细菌UMP/dUMP甲基化酶PolB，并非人类DNA聚合酶β，因此排除在人源靶点解释之外。初始Vina方案采用20 × 20 × 20 Å搜索盒，以受体ATOM坐标均值四舍五入至整埃后作为中心，exhaustiveness为32，使用五个CPU线程。补充smina对接采用exhaustiveness 8、九个输出构象和一个CPU线程，基于提取的参考结构及`--autobox_add 4`定义搜索盒。结合受体身份、辅因子及搜索区域定义考察最低评分构象。这些未经验证的计算不用于确认结合亲和力或靶点占据。

终稿准备过程中使用OpenAI Codex（GPT-6）辅助文字修订、参考文献元数据核查、核验脚本及图形排版代码。其输出需由作者审阅，作者对提交的文字、分析和解释负责。

# 结果

ZBM、XKC和XYT成分表分别包含16、68和803个非缺失的唯一PubChem CID。其中8、57和501个具有保留的基因关联，8、57和492个与至少一个疾病交集基因关联。三个集合的唯一成分关联基因数分别为107、1,094和2,131，与疾病列表重叠的分别为78/107（72.9%）、842/1,094（77.0%）和1,595/2,131（74.8%）。仅属于疾病列表的基因数分别为12,909、12,145和11,392，仅属于成分关联列表的基因数为29、252和536（补充图S1）。ZBM的全部78个及XKC的全部842个疾病交集基因均包含于XYT的1,595个交集基因中。XYT的覆盖范围最大，三个集合与疾病列表重叠的比例相近。



GO分析中，ZBM、XKC和XYT分别有1,026、3,925和4,782个条目的校正P<0.05。按BP/MF/CC顺序，已注释输入基因数分别为ZBM 78/78/78、XKC 835/840/839及XYT 1,582/1,588/1,587。这些数量对应各本体的已注释输入集合，而非最初的疾病交集总数。氧化应激反应（GO:0006979）分别涉及ZBM的9/78个基因（校正P =5.98 × 10^-4^）、XKC的106/835个（6.35 × 10^-50^）及XYT的149/1,582个（7.06 × 10^-57^）。钙稳态和电压门控钙通道复合物相关注释也有富集。后者（GO:0005891）的GeneRatio分别为15/78、15/839和19/1,587，注释特征见图2。

![](09_figures/Figure_2_GO.png){height=8.4in}

**图2：**选取的GO条目，分为BP、MF和CC三个面板。展示各药物集合在每个本体中校正P最低的三个条目的并集。柱长表示GeneRatio，橙色、蓝绿色和紫色分别表示ZBM、XKC和XYT。缺失柱形表示条目未出现在相应显著结果表中。具体条目及统计值见配套数据表。

KEGG分析中，ZBM、XKC和XYT分别有83、221和239条通路的校正P<0.05，对应已注释输入基因数为74、724和1,346。三个集合均涉及钙信号、MAPK信号和细胞凋亡。甲状腺癌通路在XKC中包含16/724个基因（校正P =2.56 × 10^-8^），在XYT中包含20/1,346个（1.16 × 10^-7^），未出现在ZBM显著结果中。XKC和XYT校正P最低的KEGG条目均为脂质与动脉粥样硬化，ZBM则为尼古丁成瘾（图3）。

![](09_figures/Figure_3_KEGG.png)

**图3：**选取的KEGG通路。展示各组校正P最低的五条通路的并集，并加入hsa05216（甲状腺癌）。柱长表示GeneRatio，颜色及缺失柱形的含义同图2。完整富集结果和作图源数据随稿提供。

**配对PTC转录组分析。**最终纳入24例患者，共28个原发肿瘤组织块及27个正常甲状腺组织块，按患者和组织类型平均后形成48个观测。20,848个可测基因中，666个符合差异表达标准，其中347个在肿瘤中升高、319个降低（图4）。XYT的1,518个可测候选基因中有111个差异基因；以全部可测基因及疾病限制背景比较时，比值比分别为2.67和1.80，校正P分别为6.08 × 10^-16^和5.31 × 10^-7^。XKC有69/805个差异基因，对应比值比为3.05和2.08，校正P为2.44 × 10^-13^和5.31 × 10^-7^。ZBM有4/77个差异基因，两项检验均未显著（校正P =0.278和0.482）。POLB略有下降（log2FC =−0.070，校正P =0.00271），但幅度未达到预设标准。CACNA1H和CYP3A4变化较小且不显著，log2FC分别为−0.053和−0.031，校正P分别为0.273和0.338。因此，这三个基因均未同时满足两项差异标准。

![](09_figures/Figure_4_GSE60542.png)

**图4：**GSE60542中24例患者的配对PTC转录组分析。（A）20,848个可测基因的配对limma火山图，虚线标示校正P =0.05及|log2FC| =1。候选集合相互重叠，先绘制XYT（紫色），再绘制XKC（蓝绿色）和ZBM（橙色），可见颜色不代表互斥的集合归属。（B）候选集合在差异基因中过度代表的条件比值比及双侧精确95%置信区间，分别采用全部可测基因和疾病列表中的可测基因为背景。FDR为六项单侧Fisher检验经Benjamini–Hochberg校正后的P值。

GSE33630纳入44对匹配组织，共发现745个差异表达基因，其中407个在肿瘤中上调、338个下调。XYT的1,518个可测候选基因中有112个差异基因，XKC为68/805。以全部可测基因为背景时，比值比分别为2.35和2.64；以疾病相关可测基因为背景时，分别为1.59（95%CI 1.27–1.97）和1.80（1.36–2.35），后两项检验的校正P均为4.99 × 10^-5^。ZBM包含2/77个差异基因，未显著富集，两种背景的校正P均为0.913。20,848个基因在两个队列中的log2倍数变化相关，Spearman ρ=0.786。初始666个差异基因中，664个方向一致，536个在GSE33630再次达到相同标准，其中包括83/111个XYT基因和51/69个XKC基因（补充图S3）。可复制候选基因的示例包括DIO1（两个队列log2FC分别为−4.86和−3.24）和TPO（−4.61和−3.24）的下调，以及LRRK2（3.59和4.29）、DPP4（3.45和3.60）和FN1（3.02和2.51）的上调（完整83个候选基因见补充表T39，正文所列为示例）。POLB、CACNA1H和CYP3A4仍未达到联合标准。放射暴露及未暴露两层各22对样本，XYT和XKC在两层中均保持显著富集，校正P均<0.001，ZBM仍不显著。该分层分析未检验暴露与组织类型的交互作用。

排除白毛夏枯草后，成分关联基因从2,131个减少至2,128个；减少的EcR、MICAL1和usp均不在疾病列表中，因此1,595个XYT候选基因及基于该集合的下游结果不变。

β-谷甾醇（CID 222284）为三个成分集合共有，选用于结构分析。不同中药标签下重复出现的CID–受体组合代表相同的化合物–蛋白配对。

CYP3A4（6UNI）和CaV3.2（9AYG）的最低smina评分分别为−9.12和−9.68 kcal/mol（补充图S2）。CYP3A4受体制备时移除了血红素辅因子，而自动搜索盒参考同时包含血红素和晶体中的抑制剂。9AYG的参考包含多类非蛋白残基。这些受体制备和搜索区域问题使现有构象不能支持特异性结合的主张。



# 讨论

本研究将所选XYT核心药组关联到甲状腺疾病相关基因，并在配对PTC组织中考察其表达。功能注释涉及氧化应激、钙调节和信号转导。可复制候选基因包括在肿瘤中下调的DIO1和TPO。在两个队列中，XYT和XKC候选基因在肿瘤差异表达基因中富集，将背景限制为疾病列表中的可测基因后仍然显著。限制背景后的比值比较小，说明富集程度部分取决于比较基因的选择。持续存在的关联支持候选集合与PTC表达特征的联系，而草药治疗效应仍需进一步检验。

两个单味药候选集合均为XYT的子集，但彼此不构成包含关系。检验存在重叠，不能视为不同草药效应的独立证据。合并成分记录后XYT覆盖范围更大是可预期的，不能据此判断协同作用。ZBM仅有77个可测基因，在GSE60542和GSE33630中分别有4个和2个差异基因，富集估计的置信区间较宽。其富集不显著不等于它与XYT或XKC有显著差异，后者需要直接比较。交集比例及GeneRatio也不是治疗效果指标。

氧化应激注释与氧化还原调节在甲状腺癌中的作用相符。[19]实验研究已将TRPC1与甲状腺癌细胞迁移和增殖联系起来，因此钙相关注释也值得关注。[20]但TRPC1和CACNA1H编码不同类型的通道，共享钙相关注释并不能建立二者之间的机制联系。本次转录组分析评估的是整体候选集合，而非逐条GO或KEGG通路的独立活性。KEGG甲状腺癌通路的富集还部分取决于输入基因已先与甲状腺肿瘤基因列表取交集。

原始结构筛选关注的三个基因均未达到联合差异表达标准。POLB在两个队列中略有下降，log2FC分别为−0.070和−0.116；CACNA1H和CYP3A4未见显著差异。表达结果本身不能确认或排除蛋白结合。CYP3A4参与凡德他尼氧化代谢，[21]提示其关联也可能涉及药物代谢。现有对接无法解决这一问题：CYP3A4缺少血红素辅因子，自动搜索区域未经结合位点验证。评分不应解释为实测结合能，也不能用于比较两个受体的结合亲和力。

既往β-谷甾醇研究提供了参考，但实验体系有所不同。Cao等研究了非小细胞肺癌，[22]Chen等则在与桥本病相关的炎症性甲状腺细胞模型中研究夏枯草成分。[23]这些工作均不能确立本研究提出的PTC相互作用。Wnt和Keap1/Nrf2信号与甲状腺癌生物学相关，[24, 25]但本研究未测量所选草药对这些通路的调节。

DIO1和TPO均属于癌症基因组图谱PTC研究所定义的甲状腺分化评分基因。[26]它们在两个队列中的下调与甲状腺分化状态改变相符，但不构成药理靶点验证。正文突出展示的基因是在查看表达结果后选择的，属于探索性的实验候选。

靶点记录混合了证据强度不同的已记录关联和预测关联。BindingDB整合时未按物种或亲和力筛选；与人类疾病列表匹配到相同符号，不能证明对应的人类化合物–蛋白相互作用。较宽的疾病列表及全注释范围的GO/KEGG背景可能偏向研究较充分的基因。Fisher检验将基因视为可交换单位，未考虑表达相关性，其P值反映该模型下的重叠程度，不代表彼此独立的生物学效应。疾病限制背景解决了部分背景选择问题，但未消除这些依赖关系。不同患者队列之间仍可能共享数据库和平台偏倚。

两个队列均为PTC，非肿瘤组织来自癌症患者而非健康供者，组织成分差异可能影响表达结果。GSE33630患者年龄中位数为16.82岁，范围12.38–24.55岁，放射暴露与未暴露者人数相等，结果向其他年龄或甲状腺癌亚型推广仍需验证。归档药物清单未规定剂量、制备方法或可达到的成分暴露量，药名匹配也纳入了两种名称含“夏枯草”的不同药物。排除白毛夏枯草不改变疾病候选集合，但开展实验前仍需明确基原和制剂。数据库检索日期缺失还限制了从当前在线资源精确重建历史分析。

后续研究应结合浓度–反应实验、非恶性对照、靶点扰动及结合或功能测定。结构研究应保留必要辅因子，明确结合位点，并使用参考配体评价对接方案。

# 结论

XYT和XKC候选基因在两个独立队列的配对PTC组织差异表达基因中均富集，采用疾病限制背景后结论仍成立。83个XYT候选基因（包括甲状腺分化标志物DIO1和TPO）在两个队列均达到差异表达标准。ZBM的富集证据尚不明确，三个单独考察的基因均未达到联合差异表达标准。这些结果为进一步开展所选草药的功能研究提供了优先考察的疾病相关基因。

# 利益冲突

Zhiyu Luo为江苏三黍生物科技有限公司员工。其余作者声明不存在与本研究相关的利益冲突。

# 基金资助

本工作得到上海市科委上海中医肿瘤临床医学研究中心项目（21MC1930500）资助。另获 2025 年上海市老年教育理论研究中心课题（ZX0551-25）及 2025 年项目"数智化医疗转型背景下上海市复合型数字健康人才培养政策与路径"（2025JZ02）支持。所列资助来源均为本工作提供了支持。

# 患者知情同意

本研究对公开可用、去标识化的GEO数据开展二次分析。

# 伦理说明

本研究使用数据库记录、分子结构及公开可用的去标识化GEO表达数据，不涉及新招募参与者、干预或动物实验。

# 人体及动物研究

本研究未开展新的人体或动物实验，而是对公开、去标识化的人类组织表达数据进行二次分析。

# 作者贡献

Tao Li：数据收集；数据整理；统计分析。Jiacai Li：指导；稿件审阅与修改。Yajie Ding：数据整理；稿件审阅与修改。Zhiyu Luo：统计分析。Xing Zhang：数据整理；稿件审阅与修改。全体作者均已审阅并批准最终稿。

# 数据可用性

基因列表、样本和配对清单、差异表达结果、富集检验及质量控制图见补充材料。表达数据和平台注释可通过GEO编号GSE60542、GSE33630及GPL570获取。跨队列比较及放射暴露分层结果也见补充材料。底层网络药理学数据可向通讯作者提出合理申请获取。分析代码不公开。


# 补充材料

补充表T01–T39见05_Supplementary_Tables.xlsx；补充图S1–S3及英文图注见06_Supplementary_Figures.docx。图形摘要单独提供。

# 参考文献

[1] Cabanillas ME, McFadden DG, Durante C. Thyroid cancer. Lancet 2016; 388(10061): 2783-2795. doi: 10.1016/s0140-6736(16)30172-6. [PMID: 27240885]

[2] Haugen BR, Alexander EK, Bible KC, Doherty GM, Mandel SJ, Nikiforov YE, et al. 2015 American Thyroid Association Management Guidelines for Adult Patients with Thyroid Nodules and Differentiated Thyroid Cancer: The American Thyroid Association Guidelines Task Force on Thyroid Nodules and Differentiated Thyroid Cancer. Thyroid 2016; 26(1): 1-133. doi: 10.1089/thy.2015.0020. [PMID: 26462967]

[3] Xing M. Molecular pathogenesis and mechanisms of thyroid cancer. Nat Rev Cancer 2013; 13(3): 184-199. doi: 10.1038/nrc3431. [PMID: 23429735]

[4] Ru J, Li P, Wang J, Zhou W, Li B, Huang C, et al. TCMSP: a database of systems pharmacology for drug discovery from herbal medicines. J Cheminform 2014; 6: 13. doi: 10.1186/1758-2946-6-13. [PMID: 24735618]

[5] Kim S, Chen J, Cheng T, Gindulyte A, He J, He S, et al. PubChem 2023 update. Nucleic Acids Res 2023; 51(D1): D1373-D1380. doi: 10.1093/nar/gkac956. [PMID: 36305812]

[6] Gilson MK, Liu T, Baitaluk M, Nicola G, Hwang L, Chong J. BindingDB in 2015: A public database for medicinal chemistry, computational chemistry and systems pharmacology. Nucleic Acids Res 2016; 44(D1): D1045-53. doi: 10.1093/nar/gkv1072. [PMID: 26481362]

[7] Kong X, Liu C, Zhang Z, Cheng M, Mei Z, Li X, et al. BATMAN-TCM 2.0: an enhanced integrative database for known and predicted interactions between traditional Chinese medicine ingredients and target proteins. Nucleic Acids Res 2024; 52(D1): D1110-D1120. doi: 10.1093/nar/gkad926. [PMID: 37904598]

[8] Harding SD, Armstrong JF, Faccenda E, Southan C, Alexander SPH, Davenport AP, et al. The IUPHAR/BPS Guide to PHARMACOLOGY in 2024. Nucleic Acids Res 2024; 52(D1): D1438-D1449. doi: 10.1093/nar/gkad944. [PMID: 37897341]

[9] Safran M, Dalah I, Alexander J, Rosen N, Iny Stein T, Shmoish M, et al. GeneCards Version 3: the human gene integrator. Database (Oxford) 2010; 2010: baq020. doi: 10.1093/database/baq020. [PMID: 20689021]

[10] Tarabichi M, Saiselet M, Trésallet C, Hoang C, Larsimont D, Andry G, et al. Revisiting the transcriptional analysis of primary tumours and associated nodal metastases with enhanced biological and statistical controls: application to thyroid cancer. Br J Cancer 2015; 112(10): 1665-1674. doi: 10.1038/bjc.2014.665. [PMID: 25965298]

[11] Dom G, Tarabichi M, Unger K, Thomas G, Oczko-Wojciechowska M, Bogdanova T, et al. A gene expression signature distinguishes normal tissues of sporadic and radiation-induced papillary thyroid carcinomas. Br J Cancer 2012; 107(6): 994-1000. doi: 10.1038/bjc.2012.302. [PMID: 22828612]

[12] Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res 2015; 43(7): e47. doi: 10.1093/nar/gkv007. [PMID: 25605792]

[13] Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS 2012; 16(5): 284-287. doi: 10.1089/omi.2011.0118. [PMID: 22455463]

[14] O'Boyle NM, Banck M, James CA, Morley C, Vandermeersch T, Hutchison GR. Open Babel: An open chemical toolbox. J Cheminform 2011; 3: 33. doi: 10.1186/1758-2946-3-33. [PMID: 21982300]

[15] Burley SK, Bhikadiya C, Bi C, Bittrich S, Chen L, Crichlow GV, et al. RCSB Protein Data Bank: powerful new tools for exploring 3D structures of biological macromolecules for basic and applied research and education in fundamental biology, biomedicine, biotechnology, bioengineering and energy sciences. Nucleic Acids Res 2021; 49(D1): D437-D451. doi: 10.1093/nar/gkaa1038. [PMID: 33211854]

[16] Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. J Comput Chem 2010; 31(2): 455-461. doi: 10.1002/jcc.21334. [PMID: 19499576]

[17] Koes DR, Baumgartner MP, Camacho CJ. Lessons learned in empirical scoring with smina from the CSAR 2011 benchmarking exercise. J Chem Inf Model 2013; 53(8): 1893-1904. doi: 10.1021/ci300604z. [PMID: 23379370]

[18] Huang J, Fan X, Jin X, Lyu C, Guo Q, Liu T, et al. Structural basis for human Cav3.2 inhibition by selective antagonists. Cell Res 2024; 34(6): 440-450. doi: 10.1038/s41422-024-00959-8. [PMID: 38605177]

[19] Cazarin J, Dupuy C, Pires de Carvalho D. Redox Homeostasis in Thyroid Cancer: Implications in Na+/I− Symporter (NIS) Regulation. Int J Mol Sci 2022; 23(11): 6129. doi: 10.3390/ijms23116129. [PMID: 35682803]

[20] Asghar MY, Magnusson M, Kemppainen K, Sukumaran P, Löf C, Pulli I, et al. Transient Receptor Potential Canonical 1 (TRPC1) Channels as Regulators of Sphingolipid and VEGF Receptor Expression: IMPLICATIONS FOR THYROID CANCER CELL MIGRATION AND PROLIFERATION. J Biol Chem 2015; 290(26): 16116-16131. doi: 10.1074/jbc.m115.643668. [PMID: 25971967]

[21] Indra R, Pompach P, Martínek V, Takácsová P, Vavrová K, Heger Z, et al. Identification of Human Enzymes Oxidizing the Anti-Thyroid-Cancer Drug Vandetanib and Explanation of the High Efficiency of Cytochrome P450 3A4 in its Oxidation. Int J Mol Sci 2019; 20(14): E3392. doi: 10.3390/ijms20143392. [PMID: 31295928]

[22] Cao W, Yuan F, Liu T, Yin R. Network pharmacology analysis, molecular docking integrated experimental verification reveal β-sitosterol as the active anti-NSCLC ingredient of Polygonatum cyrtonema Hua by suppression of PI3K/Akt/HIF-1α signaling pathway. J Ethnopharmacol 2024; 328: 117900. doi: 10.1016/j.jep.2024.117900. [PMID: 38432577]

[23] Chen Y, Jiang B, Qu C, Jiang C, Zhang C, Wang Y, et al. Bioactive components in prunella vulgaris for treating Hashimoto's disease via regulation of innate immune response in human thyrocytes. Heliyon 2024; 10(16): e36103. doi: 10.1016/j.heliyon.2024.e36103. [PMID: 39253271]

[24] Sastre-Perona A, Santisteban P. Role of the wnt pathway in thyroid cancer. Front Endocrinol (Lausanne) 2012; 3: 31. doi: 10.3389/fendo.2012.00031. [PMID: 22645520]

[25] Renaud CO, Ziros PG, Chartoumpekis DV, Bongiovanni M, Sykiotis GP. Keap1/Nrf2 Signaling: A New Player in Thyroid Pathophysiology and Thyroid Cancer. Front Endocrinol (Lausanne) 2019; 10: 510. doi: 10.3389/fendo.2019.00510. [PMID: 31428048]

[26] Cancer Genome Atlas Research Network. Integrated genomic characterization of papillary thyroid carcinoma. Cell 2014; 159(3): 676-690. doi: 10.1016/j.cell.2014.09.050. [PMID: 25417114]

# 补充图（供对照阅读）

![](09_figures/Supplementary_Figure_S1_gene_overlap.png)

**补充图S1：**ZBM、XKC及XYT核心药组与疾病基因的交集，分别包含78、842和1,595个基因。三个比较采用相同的12,987个疾病基因。圆形为示意，面积不与集合大小成比例。XYT指针对甲状腺疾病选取的消瘿汤核心药组。

![](09_figures/Supplementary_Figure_S2_docking.png)

**补充图S2：**β-谷甾醇（CID 222284）的探索性对接构象。（A）人源CYP3A4（PDB 6UNI）；（B）人源CaV3.2/CACNA1H（PDB 9AYG）。左侧为配体绿色棒状结构及周围残基的局部视图，右侧为受体整体视图。构象对应最低smina评分−9.12及−9.68 kcal/mol。其解释应考虑正文所述受体制备和搜索区域的限制。

![](09_figures/Supplementary_Figure_S3_GSE33630.png)

**补充图S3：**GSE33630独立重复验证，共44对匹配组织。（A）配对limma差异表达结果，标准为BH校正P<0.05且|log2FC|≥1。（B）以全部20,848个可测基因或其中10,655个疾病相关基因为背景的候选集合富集。点为条件比值比，误差线为双侧精确95%置信区间，q为该队列6项单侧Fisher检验经BH校正后的P值。（C）两个队列的基因log2倍数变化；蓝点表示GSE60542中符合差异标准的666个基因，其中664个在GSE33630方向一致，536个方向一致且再次达到相同标准。相关系数依据全部20,848个基因计算。
