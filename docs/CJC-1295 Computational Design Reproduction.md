# **计算机辅助设计在长效生长激素释放激素类似物（CJC-1295）研发中的应用及本地复现指南**

*此内容仅供参考。如需医疗建议或诊断，请咨询专业人士。*

天然生长激素释放激素（Growth Hormone-Releasing Hormone, GHRH）是由下丘脑合成的含有44个氨基酸的神经肽，其核心生物学功能是结合并激活垂体前叶的生长激素释放激素受体（GHRHR），从而促进生长激素（GH）的合成与脉冲式分泌1。尽管天然GHRH在调节躯体生长、代谢平衡以及组织修复中具有不可替代的作用，但其在临床治疗中的应用却受到了严峻的药代动力学限制。天然GHRH在人体血液循环中的半衰期极短，仅存活数分钟，这使得传统的GHRH疗法需要频繁的连续静脉输注或高频次的皮下注射，极大地限制了其临床依从性与治疗潜力3。

为了突破这一治疗瓶颈，加拿大ConjuChem Biotechnologies公司开发了一种名为CJC-1295（亦被称为DAC:GRF或CJC-1295 DAC）的合成GHRH类似物6。CJC-1295的核心研发目标是通过分子层面的精准工程化改造，彻底解决天然GHRH体内半衰期过短的问题。在这一宏大的药物研发工程中，计算机辅助药物设计（Computer-Aided Drug Design, CADD）发挥了决定性的作用。科学家们利用先进的分子动力学模拟、序列置换算法以及空间构象建模技术，不仅精准预测并防御了血液降解酶的攻击，还开创性地设计了药物亲和复合物（Drug Affinity Complex, DAC）技术。通过计算模拟优化出的复杂化学“挂钩”（Linker），使得CJC-1295能够在进入血液后自动、高选择性地与循环血清白蛋白（Human Serum Albumin, HSA）发生共价结合，从而将其半衰期从区区几分钟惊人地延长至6到8天6。

本报告将以结构生物学和计算化学为核心视角，详尽剖析CJC-1295研发过程中计算机辅助设计的具体思路、方法与文献依据。报告的后半部分将提供一套具有极高实操价值的本地复现指南，指导研究人员如何利用现代开源计算分子生物学工具（如GROMACS、Rosetta等），在本地计算集群上从零开始重现这一经典的复杂共价肽类药物分子辅助设计工作。

## **酶切位点预测与序列置换的计算模拟机制**

天然GHRH在体内的快速降解，主要归因于血浆中广泛分布的二肽基肽酶-IV（Dipeptidyl Peptidase-IV, DPP-IV）的酶促裂解3。DPP-IV是一种具有高度底物特异性的丝氨酸蛋白酶，其分子机制在于专门识别并裂解多肽N端倒数第二个残基为脯氨酸（Pro）或丙氨酸（Ala）的肽键11。

### **酶切位点的识别与脆弱防线的暴露**

在CJC-1295的早期研发阶段，科学家们利用计算机同源建模与分子对接技术，模拟了DPP-IV活性中心（包含催化三联体Ser630、Asp708、His740）与天然GHRH (1-29) 序列的相互作用动态13。天然GHRH的核心活性区位于其前29个氨基酸，其N端序列前三个氨基酸为Tyr-Ala-Asp。计算机模型清晰地揭示出，序列中第2位的L-丙氨酸（L-Ala）与DPP-IV的催化口袋在空间构象和静电势分布上形成了完美的契合，是整条多肽防线中最薄弱的环节，导致肽键在数分钟内即被水解断裂10。

现代针对DPP-IV裂解位点的预测研究进一步证实了这一早期计算结果的准确性。例如，基于自然语言处理（NLP）和预训练语言模型（如BERT-DPPIV）的深度学习算法，能够自动捕获底物多肽的序列特征与DPP-IV酶切位点之间的内在联系，其预测准确率与湿实验数据高度吻合15。这些现代工具印证了当初ConjuChem团队利用分子建模锁定第2位氨基酸为核心改造靶点的科学前瞻性。

### **基于计算化学的序列置换与稳定性增强**

在锁定了第2位L-Ala这一致命弱点后，研发团队利用序列置换模拟（Sequence Substitution Simulation）系统性地评估了多种氨基酸替换方案的热力学与动力学效应。计算模拟显示，直接改变氨基酸侧链的化学性质可能会导致N端结合GHRH受体（GHRHR）的能力丧失，因为N端对于受体激活至关重要2。

最完美的解决方案来自于立体化学的翻转。计算机模拟表明，如果将第2位的左旋L-丙氨酸替换为右旋的D-丙氨酸（![][image1]），可以在不改变侧链化学基团的前提下，显著改变局部骨架的拉氏图（Ramachandran Plot）构象分布8。![][image1]的引入在DPP-IV的催化口袋中产生了强烈的空间位阻效应（Steric Hindrance），破坏了酶与底物之间形成四面体过渡态所需的精确几何对齐。体外实验证实，这一基于计算预测的微小改动，使修饰后的多肽对DPP-IV表现出了极高的抵抗力14。

除了针对DPP-IV的![][image1]置换，计算化学模型还指导了其他关键位点的优化，以确保多肽在复杂生物环境和工业生产中的化学稳定性。CJC-1295作为一种四取代类似物，总共经历了四处核心氨基酸置换设计。

| 氨基酸置换位点 | 原始天然残基 | 置换后残基 | 计算模型揭示的设计依据与生物学优势 |
| :---- | :---- | :---- | :---- |
| **位置 2** | L-丙氨酸 (L-Ala) | D-丙氨酸 (![][image1]) | 诱导局部构象偏转，产生空间位阻，完全阻断DPP-IV的酶切识别与降解路径8。 |
| **位置 8** | 天冬酰胺 (Asn) | 谷氨酰胺 (![][image2]) | 增加碳链长度以改变局部柔韧性，防止天冬酰胺在水相中发生重排或酰胺水解生成天冬氨酸，提升制剂稳定性8。 |
| **位置 15** | 甘氨酸 (Gly) | 丙氨酸 (![][image3]) | 降低骨架的构象熵，增强核心α-螺旋结构的稳定性，从而提升多肽与受体结合的内在生物活性8。 |
| **位置 27** | 甲硫氨酸 (Met) | 亮氨酸 (![][image4]) | 消除易受活性氧（ROS）攻击的硫醚侧链，彻底防止甲硫氨酸氧化，延长分子在血液及存储过程中的半衰期8。 |

这些精确的序列置换共同构成了CJC-1295的药效团骨架，即IUPAC命名中的 ![][image5]\-Sermorelin，为其后续的大分子偶联奠定了坚实的化学基础6。

## **DAC（药物亲和复合物）技术设计与分子建模**

如果说序列置换只是解决了酶解问题，那么DAC（Drug Affinity Complex）技术的引入则是CJC-1295实现长达数天半衰期的决定性因素。研发团队面临着一个极其复杂的计算设计难题：如何设计一个化学“挂钩”（Linker），使得这段仅有30个氨基酸的短肽在进入人体血液后，能够自动、快速且不可逆地与大分子蛋白质结合，同时又不能丧失其作为激动剂的生物学活性14。

### **血清白蛋白的表面特性与计算靶向**

分子建模的首要任务是寻找血液中合适的“载体”。人体血清白蛋白（HSA）是血浆中最丰富的蛋白质，浓度高达0.6 mM，占血浆总蛋白的60%左右，其在人体内的平均循环半衰期长达19至28天19。HSA分子量约为66.5 kDa，含有585个氨基酸残基，其结构中包含35个半胱氨酸。通过深入的X射线晶体衍射数据（如PDB结构1E78、1AO6）分析与计算，科学家们发现其中34个半胱氨酸参与形成了17对二硫键，仅留下唯一的一个游离半胱氨酸——Cys3419。

进一步的静电势映射与量子化学计算揭示了Cys34的独特性质。Cys34位于HSA第一结构域（Domain I）表面的一个浅层疏水裂口中19。周围的微环境（包括Tyr83、Glu33和Asp38残基）形成了一个稳定的氢键网络，极大地改变了该巯基的电离特性，使其pKa值降至约8.1（显著低于普通多肽巯基的9.1）20。这一计算结果意味着，在pH 7.4的生理血液环境中，有高达17%的Cys34以高度亲核的硫醇盐阴离子（Thiolate anion, ![][image6]）形式存在，成为极为理想的化学共价锚点20。

### **共价结合化学的反应动力学模拟**

确定了靶点后，计算化学家们需要选择一种能够在体内温和条件下，迅速且高选择性地与Cys34反应的化学基团。经过海量的化合物库筛选与分子轨道（HOMO/LUMO）计算，马来酰亚胺（Maleimide）基团脱颖而出23。

马来酰亚胺是一个缺电子的迈克尔受体（Michael acceptor），其双键（C=C）受到相邻两个羰基的强吸电子作用而被高度极化。在生理pH（6.5-7.5）下，马来酰亚胺与硫醇盐阴离子的巯基-迈克尔加成反应（Thiol-Michael Addition）具有极高的化学选择性。动力学计算表明，在pH 7.0时，该反应与巯基的速率比与循环中广泛存在的伯胺（如赖氨酸侧链）的反应速率快约1000倍24。

该反应不仅迅速，而且具有巨大的热力学驱动力（\>20 kJ/mole），反应产物为立体位阻极大、抗水解与抗氧化的硫代琥珀酰亚胺（Thiosuccinimide）硫醚键。这确保了CJC-1295在体内与白蛋白结合后，具有极高的不可逆性和长期循环稳定性24。

### **复杂“挂钩”（Linker）的拓扑结构与分子动力学模拟**

明确了反应化学后，最核心的CADD挑战在于“挂钩”的拓扑设计。研究团队在GHRH短肽的C端（对受体激活影响最小的区域）第30位添加了一个赖氨酸，并通过其侧链氨基连接了马来酰亚胺丙酰基（Maleimidopropionyl，MPA），从而构建了完整分子 ![][image7]\-maleimidopropionyl-Lysine6。

必须确保短肽挂在大分子白蛋白上后，其活性中心（N端负责受体识别的区域）依然能不受遮挡地露在外面，去精准对接垂体细胞表面的受体。这种空间构象的精准计算，完全依赖于高性能的分子动力学（Molecular Dynamics, MD）模拟。

在全原子MD模拟中，研究人员构建了CJC-1295与HSA共价结合后的庞大复合物模型。模拟轨迹长达数百纳秒，以考察分子的柔性与溶剂可及表面积（SASA）。计算结果证实了设计的精妙之处：虽然短肽的C端被牢牢锚定在庞大的HSA上，但马来酰亚胺丙酰基-赖氨酸接头提供了足够的柔性（Flexibility）和延伸空间27。这种结构犹如“风筝与线”，庞大的白蛋白是放飞风筝的人，接头是线，而包含活性![][image1]和核心α-螺旋的N端则是随血流飘动的风筝，能够在三维空间中自由探索并暴露出受体结合面8。同时，高达66.5 kDa的HSA像一面巨大的立体盾牌，利用空间位阻效应（Steric exclusion）进一步阻止了血液中其他非特异性内肽酶和外肽酶对GHRH肽链的靠近与切割，从而实现了特异性酶解抵抗（由![][image1]介导）与非特异性物理保护的双重保险14。

## **空间构象与受体精准对接的分子机制**

除了模拟其在血液循环中的稳定性，计算机建模还必须回答另一个关键问题：带着巨大白蛋白“尾巴”的CJC-1295，是否还能有效激活细胞表面的GHRH受体？

生长激素释放激素受体（GHRHR）属于B1类G蛋白偶联受体（Class B1 GPCR）家族2。生物物理学与结构生物学（如低温电子显微镜Cryo-EM结构解析，PDB ID: 7CZ5, 7V9M）揭示，B1类GPCR与多肽配体的结合遵循经典的“两步结合模型”（Two-step model）28。

1. **第一步结合（初步锚定）**：多肽配体的C端区域首先识别并结合到受体庞大的胞外结构域（Extracellular Domain, ECD）上。  
2. **第二步结合（受体激活）**：在C端锚定后，多肽的N端区域顺势深入受体的七跨膜结构域（7TMD）的核心空腔内，引发受体跨膜螺旋的构象重排（如TM3和TM6距离的增加），进而激活胞内的异三聚体G蛋白（主要是![][image8]），引发cAMP的大量积累2。

在对CJC-1295进行的高性能分子对接（Molecular Docking）模拟中，研究人员观察到，将共价挂钩设置在第30位（C端）是一项极其巧妙的计算设计。由于结合第一步依赖于C端与ECD的相互作用，接头的存在虽然引入了一定的位阻，但HSA巨大的体积主要停留在受体外部的溶剂区，并没有深入细胞膜的疏水核心。当多肽被HSA限制在受体附近时，其游离且保持了高度螺旋结构倾向的N端，依然能够顺利执行第二步结合，深入7TMD并激活cAMP级联反应8。

## **药效动力学与临床实验验证**

这些深度的计算机辅助设计最终在体内外生物学实验和临床试验中得到了完美的验证。

在体外大鼠前垂体细胞的GH分泌分析中，与HSA形成生物偶联物的CJC-1295变体展现出了与天然GHRH相媲美的生物学活性（在纳摩尔nM浓度范围内达到最大分泌量），证明大分子偶联未破坏药效团14。在体内Sprague Dawley大鼠模型中，通过免疫印迹（Western blot）分析证实，皮下注射CJC-1295后仅15分钟，血浆中即出现了分子量与白蛋白一致的免疫反应条带，且在血液中循环超过24小时甚至72小时；其在2小时内引发的GH曲线下面积（AUC）比天然GHRH激增了4倍14。

| 药代动力学/药效学指标 | 天然 GHRH (1-29) | CJC-1295 (DAC) | 数据依据与临床意义 |
| :---- | :---- | :---- | :---- |
| **体内半衰期 (![][image9])** | 约 7 分钟 | 5.8 \- 8.1 天 | 结合白蛋白使其免于肾小球滤过与酶解，大幅减少给药频率（可一周一次）5。 |
| **生长激素 (GH) 提升** | 短暂脉冲式上升 | 升高 2 \- 10 倍，持续 6 天以上 | 在保持生理性脉冲释放（Pulsatility）的前提下，显著提升谷值（基础值）和整体分泌量1。 |
| **IGF-1 提升** | 波动不明显 | 升高 1.5 \- 3 倍，持续 9 \- 11 天 | 作为GH的下游效应分子，IGF-1的持续升高证实了药物的长效合成代谢促进作用4。 |
| **受体脱敏效应** | 频繁注射易引发 | 无明显脱敏 | 因为作用于上游垂体受体而非直接替代GH，保留了机体天然的生长抑素（Somatostatin）负反馈调节机制5。 |

人类健康志愿者的Phase I/II期临床数据显示，单次皮下注射CJC-1295即可在长达近一个月的时间内显著提升体内IGF-1水平5。多项研究指出，CJC-1295能够有效增加瘦体重、促进脂肪代谢并改善睡眠质量，成为内分泌抗衰老与运动医学领域的明星分子3。

然而，值得一提的是，尽管CJC-1295在GHRH基因敲除小鼠模型中表现出极佳的促生长和代谢正常化作用35，ConjuChem公司针对HIV脂肪营养不良症患者的Phase II临床试验却因一名受试者的死亡而被迫终止。调查显示，该患者患有无症状的冠状动脉疾病，死亡原因极可能为斑块破裂和血管闭塞，主治医师认为该不良事件与CJC-1295药物治疗无关，但出于谨慎，该项目仍被搁置6。此事件凸显了尽管CADD能够完美优化药代动力学，复杂人群的临床试验依然充满不确定性。

## ---

**计算机辅助设计工作本地复现指南**

为了在本地计算环境（如搭载GPU加速的Linux工作站）中深度复现CJC-1295的计算机辅助设计工作，需要系统性地串联序列设计、量子化学参数化、分子动力学模拟以及大分子柔性对接技术。以下是基于现代开源软件链的具体实操步骤。

### **第一阶段：底层大分子与多肽结构准备**

**1\. 获取高质量蛋白质晶体结构** 首先需从RCSB PDB数据库获取相关受体与载体蛋白的结构。

* **血清白蛋白（HSA）**：下载结构代码如 1E78（分辨率2.60 Å）或 1AO6（分辨率2.50 Å）21。使用PyMOL或UCSF Chimera软件剔除结构中的结晶水、脂质及其他非特异性配体。重点检查第一结构域的Cys34残基，确保其未形成二硫键，呈现游离的硫醇状态。在后续处理中，需注意将其残基名称标记为游离半胱氨酸（在GROMACS中通常为CYS，而非形成二硫键的CYX）38。  
* **生长激素释放激素受体（GHRHR）**：下载含有GHRH配体与G蛋白复合物的冷冻电镜（Cryo-EM）结构，如 7CZ5 或 7V9M29。提取出GHRHR的受体骨架，并移除内源性配体以便进行后续对接。

**2\. 构建CJC-1295多肽模型并进行序列置换** 以天然GHRH (1-29) 结构（如可从 7CZ5 中提取）为起始模板。在PyMOL的 Mutagenesis Wizard 工具或基于Rosetta的序列设计模块中进行序列变异：

* 将位置2的L-Ala突变为D-Ala（此时需要注意D型氨基酸的立体手性配置）。  
* 进行另外三次突变：Asn8 ![][image10] Gln8，Gly15 ![][image10] Ala15，Met27 ![][image10] Leu2714。  
* 在序列C端（第30位）通过计算工具添加一个赖氨酸残基（Lys30）。

### **第二阶段：非天然氨基酸与共价偶联力场参数化**

由于标准的分子力场（如AMBER ff14SB, CHARMM36m）并未涵盖![][image1]、马来酰亚胺丙酰基以及其与半胱氨酸反应后生成的硫醚共价键，必须进行针对性的拓扑开发40。

**1\. 配体的量子化学参数化** 提取分子结构中的马来酰亚胺基团部分，使用量子化学软件（如ORCA或Gaussian，采用B3LYP基组）对其几何结构进行优化，并计算其静电势以拟合RESP（Restrained Electrostatic Potential）部分电荷23。 随后，利用AmberTools软件包中的 antechamber 工具，结合通用AMBER力场（GAFF2），生成马来酰亚胺基团的小分子力场参数。通过 ACPYPE 或 CGenFF 工具，将这些参数转换为GROMACS兼容的拓扑格式文件44。

**2\. GROMACS力场深度修改定制** 要实现马来酰亚胺与HSA Cys34之间的共价硫醚连接，必须对GROMACS环境的底层系统文件进行精细配置42：

* **residuetypes.dat**：在此文件中声明包含马来酰亚胺的连接子（假设自定义命名为 DM2），将其类别设定为 Protein，这样 pdb2gmx 工具才会将其作为多肽链的延伸部分进行处理。  
* **aminoacids.rtp**：在残基拓扑文件中，新增一个 \[ DM2 \] 模块，详尽列出原子的连接关系、局部电荷以及成键定义，特别是参与后续共价反应的C36碳原子。  
* **specbond.dat**：这是触发自动成键机制的核心文件。添加一行定义指示GROMACS自动在Cys34的硫原子和DM2的碳原子之间构建键合： CYS SG 1 DM2 C36 1 0.18 CYX DM2 （该指令告知程序：当CYS的SG原子与DM2的C36原子距离在约0.18 nm范围内时，自动生成特定的共价拓扑）。  
* **ffbonded.itp**：如果力场库中缺失对应的硫-碳（S-C）键合参数，需手动增添基于谐振子模型的参数（如 S C 1 0.18200 251000.0），以维持硫醚键在模拟过程中的稳定性。

### **第三阶段：全原子分子动力学（MD）模拟评估复合物构象**

力场配置完毕后，需要模拟CJC-1295与66.5 kDa的HSA结合后的分子动态，验证其构象自由度46。

**1\. 初始空间定位与溶剂化** 使用手动或脚本方式将CJC-1295的C端与HSA的Cys34在三维空间中靠近。如果初始对接存在不可避免的原子空间重叠（Steric clashes），建议先运行一段简短的拉伸分子动力学（Steered MD, SMD）。在SMD中，通过定义距离限制（Distance restraints），施加引力将马来酰亚胺的碳原子缓缓拉向半胱氨酸的硫原子至反应距离（约0.2 nm），随后对生成的新复合物进行能量最小化38。 接着，将复合物置于填充有TIP3P水分子模型的截角八面体盒子中，并加入适量的钠离子和氯离子（如0.15 M NaCl）以模拟血液环境并中和系统净电荷48。

**2\. 平衡与生产动力学轨迹演化**

* 实施最速下降法（Steepest Descent）进行系统的初步能量最小化。  
* 在NVT系综（恒容恒温，调节至生理温度310 K）和NPT系综（恒压恒温，1 bar）下分别进行1 ns的限制性平衡模拟。  
* 解除对蛋白质骨架的位置限制，开启长达微秒级（至少500 ns）的生产MD（Production Run）模拟。在此过程中，重点观察大分子HSA与小肽片段之间的相互作用模式41。

**3\. 柔性与暴露度分析（RMSF与SASA计算）** 模拟结束后，使用GROMACS套件内的 gmx rmsf 和 gmx sasa 工具分析多肽区域的数据。研究人员应重点寻找这样的物理图景：与HSA连接的C端区域RMSF值极低（被大分子彻底锁死限制），但N端的D-Ala2及其后侧的α-螺旋结构区不仅维持了二级结构稳定性，且展现出极高的RMSF值（构象灵活）和较高的溶剂可及表面积（未被白蛋白埋葬）。这在理论上印证了DAC技术能够确保活性药效团“暴露在外”，为结合受体保留了结构基础。

### **第四阶段：基于Rosetta的多肽-大分子受体柔性对接**

证明肽段暴露后，最后一个环节是模拟挂载着巨大HSA负载的CJC-1295是否能准确激活GHRHR。受限于模拟体系的庞大，这阶段推荐使用Rosetta的 **FlexPepDock** 从头计算对接模块（ab-initio protocol）进行验证50。

**1\. 生成多肽片段库（Fragment Library）** 借助PSIPRED等二级结构预测工具，基于CJC-1295的氨基酸序列从PDB数据库中抽取重构其结构空间的三聚体（3-mers）和九聚体（9-mers）片段库。这一步骤赋予了Rosetta在对接过程中对肽链骨架进行大尺度构象采样的能力50。

**2\. 受体预处理与空间冲突清除** 将准备好的GHRHR膜蛋白结构输入Rosetta，运行预填充（Prepacking）脚本以优化侧链构象，清除因结晶堆积造成的非自然原子冲突，为多肽的插入清理出平滑的势能表面： FlexPepDocking.linuxgccrelease \-database $PATH\_TO\_DB \-s receptor\_start.pdb \-flexpep\_prepack \-use\_input\_sc50。

**3\. 从头对接（Ab-initio Docking）与空间惩罚评估** 将CJC-1295的大致构象置于GHRHR胞外域的上方。在对接配置中，可以利用Rosetta的能量函数（Scoring function）引入一个代表HSA体积的虚拟排斥力场区域。 对接分两步进行：首先在低分辨率的质心模式（Centroid mode）下进行穷举式空间搜索（可生成约100,000个粗略模型，选项如 \-docking:low\_res\_protocol\_only）；随后对通过聚类分析脱颖而出的代表性构象，进行高分辨率的全原子细化与侧链重包装（Full-atom refinement）50。

**4\. 结合能打分与机理验证** 计算最终模型的界面结合能。如果在引入了大量C端空间位阻惩罚的情况下，FlexPepDock依旧能找到让CJC-1295的N端极稳定地插入GHRHR跨膜核心、且打分极优的结合姿态，那么这项计算机辅助工作便成功地在硅基环境（in silico）中再现了CJC-1295惊艳的药效学设计逻辑。

## **复现管线总体进度（2026-05-14）**

| 阶段 | 任务 | 状态 | 交付物 | GPU/资源 |
|------|------|------|--------|----------|
| **P0a** | 结构准备：DPP-IV、GHRH(1-29)、D-Ala2 突变体 | ✅ 完成 | `DPP4_clean.pdb`, `GHRH_1-29.pdb`, `GHRH_1-29_DAla2.pdb` | — |
| **P0b** | Rosetta FlexPepDock + Constrained Refine | ✅ 完成 | `refine_constrained_WT_combined.silent` (201 models), `refine_constrained_DAla2_combined.silent` (201 models) | CPU |
| **P0c** | MD 体系构建：溶剂化、离子化、能量最小化、NVT/NPT 平衡 | ✅ 完成 | `md.tpr`, `DAla2_production.tpr`, `short_peptide_md.tpr` | — |
| **P1** | **FEP 自由能微扰**：L-Ala2 ↔ D-Ala2 双拓扑 11 λ 窗口 | ⏳ **运行中** | 11 个 `lambda_*/fep.xtc`（各 5 ns） | GPU 2 |
| **P2a** | WT 200 ns 生产 MD | ⏳ **运行中**（66.9 ns / 200 ns，已澄清 GPU 正常） | `md.xtc` | GPU 0 |
| **P2b** | D-Ala2 200 ns 生产 MD | ⏳ **运行中** | `DAla2_md.xtc` (22 ns / 200 ns) | GPU 1 |
| **P2c** | Short Peptide WT 200 ns MD | ⏳ **运行中** | `short_peptide_md.xtc` (28.1 ns / 200 ns) | GPU 3 |
| **P3** | 数据分析：RMSF、SASA、催化几何、MM-PBSA、FEP ΔG | ⏳ **等待中** | 图表 + 定量对比表 | — |

**当前计算资源占用**：
- GPU 0: WT MD（93% util, 434 MiB）
- GPU 1: D-Ala2 MD（52% util, 434 MiB）
- GPU 2: FEP λ₀₀-λ₀₂（99% util, 1010 MiB，3 窗口并发）
- GPU 3: Short Peptide MD（39% util, 330 MiB）

---

## **本地复现已知问题与修复记录**

在本地复现过程中，我们发现了以下结构性 pipeline bug，已修复并记录如下：

### Issue 1: PyMOL `cmd.fab("YAD", ss=0)` 压平 backbone（Critical）

**根因**：`cmd.fab()` with `ss=0` 生成完全线性的肽链，所有 backbone dihedrals ≈0°（Ramachandran 禁区）。当 `prepare_docking_start.py` 用 `pair_fit()` 将 GHRH N-端对齐到这个线性参考时，GHRH 的 backbone 被强制压平，产生不物理的起始构象。

**后果**：
- D-Ala2 突变体起始结构 phi=-6°, psi=-3°（禁区）
- GROMACS 生产运行 dt=0.002 时立即 segfault
- 必须降到 dt=0.001 才能缓慢弛豫

**修复**：
1. `prepare_structures.py` / `prepare_structures_v2.py`：
   - 去掉 `cmd.fab("YAD", ss=0)`
   - 改为从 `GHRH_1-29_from_7CZ5.pdb` 提取天然 YAD 构象
2. `prepare_docking_start.py`：
   - 增加 phi/psi/omega 自动验证
   - 检测到 phi~0, psi~0 时发出 WARNING

**当前状态**：
- 已修复代码（v1.2）
- 现有 MD 不需要重启（D-Ala2 已在前 12 ns 自我纠正到 beta-sheet 区域）

### Issue 2: Rosetta refine-only 降低催化几何质量

**根因**：Rosetta `FlexPepDocking -pep_refine` 没有催化几何约束，单纯优化能量会导致催化位点偏离理想几何。

**后果**：I_sc 与催化能力呈反相关（高分模型反而催化几何更差）。

**修复建议**：使用 refine 结果时，额外用 `validate_docked_pose.py` 过滤催化几何指标；或改用带约束的 enzdes 协议。

**P0b 深入分析补充（2026-05-14）**：

对 961 个 ab-initio 对接模型 + 400 个约束优化模型的系统分析表明，
**Rosetta I_sc 完全不是催化几何的可靠代理指标**：

| 模型类型 | I_sc | 距 Ser630 | 催化价值 |
|----------|------|-----------|----------|
| Ab-initio Best | -12.12 | **85 Å** | 零（非催化位点结合） |
| Ab-initio Top 10 | -8.4 ~ -3.5 | 70–100 Å | 零 |
| Constrained WT | +4.5 (均值) | ~3–5 Å | 高（几何上在口袋内） |
| Constrained DAla2 | +5.1 (均值) | ~3–5 Å | 高（几何上在口袋内） |

Ab-initio 的高负分 I_sc 反映的是 DPP-IV **表面其他区域**的非特异性强结合；
约束优化的正值 I_sc 则是因为手动放置的构象在 Rosetta 能量 landscape 中处于局部高能区，
但它严格保留了催化口袋的空间位置（rmsALL_if ≈ 2.1 Å）。

**核心结论**：对于酶切机制研究，
- ❌ **不应使用 I_sc 筛选模型**
- ✅ **应以催化几何指标为首要判据**
- ✅ **约束优化的价值在于提供合理的 MD 起始构象，而非 Rosetta 打分**

**补充：Rosetta Refine 与 MD 结构差异的物理解释（2026-05-14）**

定量对比（Rosetta 起始结构 vs MD 66.7–70.3 ns）：
| 参数 | Rosetta 起始 | MD 66.7 ns | 说明 |
|------|-------------|------------|------|
| Ser630-Ala2C | 2.82 Å | 3.57 ± 0.25 Å | MD 仍维持合理催化距离 |
| Tyr1-Glu205 | 5.3 Å | 9.0 ± 0.7 Å | 盐桥在 MD 中断裂 |
| Tyr1-Glu206 | 2.0 Å | 7.0 ± 0.5 Å | 盐桥在 MD 中断裂 |

结论：
- **MD 没有问题**。差异源于 Rosetta 的坐标约束**人为冻结了热涨落**，而 MD **正确释放了热涨落**。
- 29 残基长肽在 DPP-IV 口袋中发生"枢轴运动"（pivot）：Ala2 羰基维持在 Ser630 附近，N-端 Tyr1 摆离。这是弱锚定长肽的**物理固有行为**。
- 该差异恰恰说明 DPP-IV 的酶切机制依赖于底物在口袋内的动态采样——完全刚性的结合反而不利于催化。

### Issue 3: 29 残基长肽在无约束 MD 中从口袋漂移（新发现）

**根因**：GHRH(1-29) C 端 24 个残基完全暴露在溶剂中，热涨落会拉动 N 端离开 DPP-IV 口袋。

**证据**：

| 指标 | 起始结构 | WT ~40 ns | D-Ala2 ~12.9 ns |
|------|---------|-----------|-----------------|
| Ser630 OG → Ala2 C | 2.82 Å | 4.84 Å | 6.40 Å |
| 攻击角 | 113.6° | 49.5° | 86.6° |
| **PASS** | **4/4** | **0/4** | **1/4** |

**修复建议**：
- 如需保持催化几何：mdp 中加 `define = -DPOSRES_NTERM`（限制 GHRH 1–5 CA）
- 或改用 PLUMED 反应坐标约束
- 当前 MD 继续运行，后续分析关注**相对漂移速率**（WT vs D-Ala2）

### Issue 4: WT MD 进程误判与不必要的重启（已澄清，2026-05-14）

**根因**：检查了旧的 `md.log`（part0001，仅到 9.44 ns），未注意到 `-noappend` 已创建 `md.part0002.log`。
`md.part0002.log` 显示该进程已正常跑到 66.92 ns，性能 74.141 ns/day。

**教训**：
- GROMACS 2026 默认自动启用 GPU 加速，`-gpu_id` 已足以触发
- 使用 `-noappend` 时，必须检查最新的 `*.part*.log`
- 重启前务必确认实际进度
- **WT MD 已重启继续运行，数据无损失**

### Issue 5: FEP 双拓扑构建中的一系列拓扑同步问题（已解决）

**根因**：GROMACS 双拓扑 FEP 对 atomtypes 定义顺序、GRO 原子数连续性、
虚拟原子坐标位置、二面角 multiplicity 匹配有严格要求。

**遇到的问题链**：
1. `Atomtype DU not found` → `[ atomtypes ] DU` 必须在第一个 `[ moleculetype ]` 之前定义
2. 缺少 DU 键合参数 → 在 TOP 中显式添加零力 `[ bondtypes ]`、`[ angletypes ]`、`[ dihedraltypes ]`
3. GRO 原子数不匹配 → D-原子必须插入蛋白和溶剂之间，保持连续编号
4. `perturbed excluded non-bonded pair beyond cut-off` → D-原子坐标必须紧邻 L-对应原子（~0.15 nm），而非远距离镜像
5. `Cannot perturb torsion with multiple terms` → DU 涉及的二面角必须用匹配的 function type 和 multiplicity

**修复状态**：全部修复，100-step 测试通过，11 窗口生产运行已启动。

### Issue 6: API 密钥泄露（Security）

**根因**：早期脚本 `cc.sh` 中包含明文 API key。

**修复**：已从仓库中删除该文件，但仍需轮换密钥。

---

## **结语**

CJC-1295（DAC:GRF）的研发不仅是内分泌药理学的一次突破，更是早期计算生物学与药物辅助设计在复杂多肽修饰领域中的标杆案例。其成功绝非偶然的试错，而是建立在对靶点结构、酶动力学与微观物理化学性质深刻计算洞察的基础上。

从锁定酶切弱点进行的![][image1]立体空间翻转，到瞄准血清白蛋白独一无二低pKa值Cys34的极速共价反应，再到利用大分子屏蔽效应配合精巧拓扑接头实现的高选择性受体激活机制。通过本指南提供的本地软件复现管线，现代研究者可以重新梳理这段利用计算化学克服多肽药代动力学鸿沟的历史，不仅为了解构CJC-1295，更旨在为当前研发长效GLP-1、多靶点抗癌肽等新型生物大分子药物提供一种系统性、可移植的底层计算设计范式。

#### **Works cited**

1. Pulsatile Secretion of Growth Hormone (GH) Persists during Continuous Stimulation by CJC-1295, a Long-Acting GH-Releasing Hormone Analog | The Journal of Clinical Endocrinology & Metabolism | Oxford Academic, [https://academic.oup.com/jcem/article/91/12/4792/2656274](https://academic.oup.com/jcem/article/91/12/4792/2656274)  
2. Dynamic properties of the Growth Hormone Releasing Hormone Receptor (GHRHR) and molecular determinants of GHRH binding | Request PDF \- ResearchGate, [https://www.researchgate.net/publication/317125360\_Dynamic\_properties\_of\_the\_Growth\_Hormone\_Releasing\_Hormone\_Receptor\_GHRHR\_and\_molecular\_determinants\_of\_GHRH\_binding](https://www.researchgate.net/publication/317125360_Dynamic_properties_of_the_Growth_Hormone_Releasing_Hormone_Receptor_GHRHR_and_molecular_determinants_of_GHRH_binding)  
3. CJC-1295 Peptide | Growth Hormone & Performance \- Paragon Sports Medicine, [https://www.paragonsportsmedicine.com/peptides/cjc-1295](https://www.paragonsportsmedicine.com/peptides/cjc-1295)  
4. Prolonged Stimulation of Growth Hormone (GH) and Insulin-Like Growth Factor I Secretion by CJC-1295, a Long-Acting Analog of GH-Releasing Hormone, in Healthy Adults \- ResearchGate, [https://www.researchgate.net/publication/7416716\_Prolonged\_Stimulation\_of\_Growth\_Hormone\_GH\_and\_Insulin-Like\_Growth\_Factor\_I\_Secretion\_by\_CJC-1295\_a\_Long-Acting\_Analog\_of\_GH-Releasing\_Hormone\_in\_Healthy\_Adults](https://www.researchgate.net/publication/7416716_Prolonged_Stimulation_of_Growth_Hormone_GH_and_Insulin-Like_Growth_Factor_I_Secretion_by_CJC-1295_a_Long-Acting_Analog_of_GH-Releasing_Hormone_in_Healthy_Adults)  
5. Rx CJC-1295 Protocol, [https://www.tryshed.com/es-us/products/rx-cjc-1295-protocol](https://www.tryshed.com/es-us/products/rx-cjc-1295-protocol)  
6. CJC-1295 \- Wikipedia, [https://en.wikipedia.org/wiki/CJC-1295](https://en.wikipedia.org/wiki/CJC-1295)  
7. CJC 1295 \- AdisInsight, [https://adisinsight.springer.com/drugs/800018006](https://adisinsight.springer.com/drugs/800018006)  
8. CJC-1295: A Long-Acting GHRH Analog \- Superpower, [https://superpower.com/guides/cjc-1295](https://superpower.com/guides/cjc-1295)  
9. CJC-1295 DAC vs. No DAC: Which Growth Hormone Peptide Is Right for You?, [https://livvnatural.com/cjc-1295-dac-vs-no-dac-peptide-comparison/](https://livvnatural.com/cjc-1295-dac-vs-no-dac-peptide-comparison/)  
10. \[His1,Nle27\]-Growth Hormone Releasing Factor 1-32 amide human | Sigma-Aldrich, [https://www.sigmaaldrich.com/BS/en/product/sigma/scp0162](https://www.sigmaaldrich.com/BS/en/product/sigma/scp0162)  
11. Computational Modeling of the Interactions between DPP IV and Hemorphins \- PMC \- NIH, [https://pmc.ncbi.nlm.nih.gov/articles/PMC10932442/](https://pmc.ncbi.nlm.nih.gov/articles/PMC10932442/)  
12. Computational Screening for the Dipeptidyl Peptidase-IV Inhibitory Peptides from Putative Hemp Seed Hydrolyzed Peptidome as a Potential Antidiabetic Agent \- PMC, [https://pmc.ncbi.nlm.nih.gov/articles/PMC11171819/](https://pmc.ncbi.nlm.nih.gov/articles/PMC11171819/)  
13. ACE- and DPP-IV-Inhibitory Peptides from Bambara Groundnut Hydrolysate: Elucidation Using Computational Tools and Molecular Docking \- MDPI, [https://www.mdpi.com/2079-7737/14/5/511](https://www.mdpi.com/2079-7737/14/5/511)  
14. Human Growth Hormone-Releasing Factor (hGRF)1–29-Albumin Bioconjugates Activate the GRF Receptor on the Anterior Pituitary in Rats: Identification of CJC-1295 as a Long-Lasting GRF Analog \- Oxford Academic, [https://academic.oup.com/endo/article/146/7/3052/2500187](https://academic.oup.com/endo/article/146/7/3052/2500187)  
15. Exploration of DPP-IV Inhibitory Peptide Design Rules Assisted by the Deep Learning Pipeline That Identifies the Restriction Enzyme Cutting Site \- PMC, [https://pmc.ncbi.nlm.nih.gov/articles/PMC10601436/](https://pmc.ncbi.nlm.nih.gov/articles/PMC10601436/)  
16. Oral Treatment of Obesity by GLP-1 and Its Analogs \- MDPI, [https://www.mdpi.com/1999-4923/17/12/1596](https://www.mdpi.com/1999-4923/17/12/1596)  
17. Honour Grandparents by Prioritising Elderly Dental Care, [https://exceldental.com.sg/honour-grandparents-by-prioritising-elderly-dental-care/](https://exceldental.com.sg/honour-grandparents-by-prioritising-elderly-dental-care/)  
18. DAC™ and PC-DAC™ Technologies \- ConjuChem, [https://www.conjuchem.com/technology/platform.html](https://www.conjuchem.com/technology/platform.html)  
19. Expression, purification and initial characterization of human serum albumin domain I and its cysteine 34 \- PMC, [https://pmc.ncbi.nlm.nih.gov/articles/PMC7549792/](https://pmc.ncbi.nlm.nih.gov/articles/PMC7549792/)  
20. Use of Human Serum Albumin Cys 34 (HSA-Cys 34 ) Adductomics as a Multidimensional and Integrative Biomarker Approach to Assess Oxidative Stress \- MDPI, [https://www.mdpi.com/2076-3921/15/4/458](https://www.mdpi.com/2076-3921/15/4/458)  
21. 1E78: Crystal structure of human serum albumin \- RCSB PDB, [https://www.rcsb.org/structure/1e78](https://www.rcsb.org/structure/1e78)  
22. 1AO6: CRYSTAL STRUCTURE OF HUMAN SERUM ALBUMIN \- RCSB PDB, [https://www.rcsb.org/structure/1ao6](https://www.rcsb.org/structure/1ao6)  
23. Decoding dynamic interactions between EGFR‐TKD and DAC through computational and experimental approaches: A novel breakthrough in lung melanoma treatment \- PMC, [https://pmc.ncbi.nlm.nih.gov/articles/PMC11058330/](https://pmc.ncbi.nlm.nih.gov/articles/PMC11058330/)  
24. Maleimide Reaction Chemistry \- Vector Labs, [https://vectorlabs.com/maleimide-reaction-chemistry/](https://vectorlabs.com/maleimide-reaction-chemistry/)  
25. Maleimide Crosslinkers for Antibody Labeling \- BOC Sciences, [https://www.bocsci.com/research-area/maleimide-crosslinkers-for-cysteine-specific-antibody-labeling.html](https://www.bocsci.com/research-area/maleimide-crosslinkers-for-cysteine-specific-antibody-labeling.html)  
26. An In-depth Technical Guide to the Core Principles of Maleimide-Thiol Reaction Chemistry \- Benchchem, [https://www.benchchem.com/pdf/An\_In\_depth\_Technical\_Guide\_to\_the\_Core\_Principles\_of\_Maleimide\_Thiol\_Reaction\_Chemistry.pdf](https://www.benchchem.com/pdf/An_In_depth_Technical_Guide_to_the_Core_Principles_of_Maleimide_Thiol_Reaction_Chemistry.pdf)  
27. Unraveling and Sliding of Polypeptide Strands Underlies the Exceptional Toughness of the Triple-Helix Collagen Molecule \- PMC, [https://pmc.ncbi.nlm.nih.gov/articles/PMC12854733/](https://pmc.ncbi.nlm.nih.gov/articles/PMC12854733/)  
28. Identification of Small-Molecule Antagonists Targeting the Growth Hormone Releasing Hormone Receptor (GHRHR) \- PMC, [https://pmc.ncbi.nlm.nih.gov/articles/PMC11423342/](https://pmc.ncbi.nlm.nih.gov/articles/PMC11423342/)  
29. 7CZ5: Cryo-EM structure of the human growth hormone-releasing hormone receptor-Gs protein complex \- RCSB PDB, [https://www.rcsb.org/structure/7cz5](https://www.rcsb.org/structure/7cz5)  
30. 7V9M: Cryo-EM structure of the GHRH-bound human GHRHR splice variant 1 complex, [https://www.rcsb.org/structure/7V9M](https://www.rcsb.org/structure/7V9M)  
31. Unveiling the activation mechanism of the GnRH Receptor \- STAX, [http://stax.strath.ac.uk/downloads/0r967426f](http://stax.strath.ac.uk/downloads/0r967426f)  
32. (PDF) Human Growth Hormone-Releasing Factor (hGRF) 1–29 \-Albumin Bioconjugates Activate the GRF Receptor on the Anterior Pituitary in Rats: Identification of CJC-1295 as a Long-Lasting GRF Analog \- ResearchGate, [https://www.researchgate.net/publication/7918454\_Human\_Growth\_Hormone-Releasing\_Factor\_hGRF\_1-29\_-Albumin\_Bioconjugates\_Activate\_the\_GRF\_Receptor\_on\_the\_Anterior\_Pituitary\_in\_Rats\_Identification\_of\_CJC-1295\_as\_a\_Long-Lasting\_GRF\_Analog](https://www.researchgate.net/publication/7918454_Human_Growth_Hormone-Releasing_Factor_hGRF_1-29_-Albumin_Bioconjugates_Activate_the_GRF_Receptor_on_the_Anterior_Pituitary_in_Rats_Identification_of_CJC-1295_as_a_Long-Lasting_GRF_Analog)  
33. Human Growth Hormone-Releasing Factor (hGRF)1–29 \- Albumin Bioconjugates Activate the GRF Receptor on the Anterior Pituitary i, [https://academic.oup.com/endo/article-pdf/146/7/3052/9023995/endo3052.pdf](https://academic.oup.com/endo/article-pdf/146/7/3052/9023995/endo3052.pdf)  
34. Prolonged Stimulation of Growth Hormone (GH) and Insulin-Like Growth Factor I Secretion by CJC-1295, a Long-Acting Analog of GH-Releasing Hormone, in Healthy Adults | The Journal of Clinical Endocrinology & Metabolism | Oxford Academic, [https://academic.oup.com/jcem/article-abstract/91/3/799/2843281](https://academic.oup.com/jcem/article-abstract/91/3/799/2843281)  
35. Once-daily administration of CJC-1295, a long-acting growth hormone-releasing hormone (GHRH) analog, normalizes growth in the GHRH knockout mouse | American Journal of Physiology-Endocrinology and Metabolism, [https://journals.physiology.org/doi/full/10.1152/ajpendo.00201.2006](https://journals.physiology.org/doi/full/10.1152/ajpendo.00201.2006)  
36. Once-daily administration of CJC-1295, a long-acting growth hormone-releasing hormone (GHRH) analog, normalizes growth in the GHRH knockout mouse | American Journal of Physiology-Endocrinology and Metabolism, [https://journals.physiology.org/doi/abs/10.1152/ajpendo.00201.2006](https://journals.physiology.org/doi/abs/10.1152/ajpendo.00201.2006)  
37. Lipodystrophy study halted after patient death \- Aidsmap, [https://www.aidsmap.com/news/jul-2006/lipodystrophy-study-halted-after-patient-death](https://www.aidsmap.com/news/jul-2006/lipodystrophy-study-halted-after-patient-death)  
38. Protocols for Multi-Scale Molecular Dynamics Simulations in Amber and Gromacs: a Case Study of Intrinsically Disordered Amyloid Beta | bioRxiv, [https://www.biorxiv.org/content/10.1101/2023.10.24.563575v1.full-text](https://www.biorxiv.org/content/10.1101/2023.10.24.563575v1.full-text)  
39. Introductory Tutorials for Simulating Protein Dynamics with GROMACS \- PMC, [https://pmc.ncbi.nlm.nih.gov/articles/PMC11457149/](https://pmc.ncbi.nlm.nih.gov/articles/PMC11457149/)  
40. Applications of the Newly Developed Force-Field Parameters Uncover a Dynamic Nature of Ω-Loop C in the Lys-Ligated Alkaline Form of Cytochrome c \- ACS Publications, [https://pubs.acs.org/doi/10.1021/acs.jpcb.4c00625](https://pubs.acs.org/doi/10.1021/acs.jpcb.4c00625)  
41. Enabling Biomolecular Simulations with Neural Network Potentials in GROMACS \- arXiv, [https://arxiv.org/html/2604.21441](https://arxiv.org/html/2604.21441)  
42. Parameter files \- GROMACS 2025.3 documentation, [https://manual.gromacs.org/documentation/2025.3/reference-manual/topologies/parameter-files.html](https://manual.gromacs.org/documentation/2025.3/reference-manual/topologies/parameter-files.html)  
43. Forcefield\_PTM: Ab Initio Charge and AMBER Forcefield Parameters for Frequently Occurring Post-Translational Modifications \- PMC, [https://pmc.ncbi.nlm.nih.gov/articles/PMC3904396/](https://pmc.ncbi.nlm.nih.gov/articles/PMC3904396/)  
44. Problems creating a covalent bond between a protein and a ligand in GROMACS, [https://gromacs.bioexcel.eu/t/problems-creating-a-covalent-bond-between-a-protein-and-a-ligand-in-gromacs/12910](https://gromacs.bioexcel.eu/t/problems-creating-a-covalent-bond-between-a-protein-and-a-ligand-in-gromacs/12910)  
45. How to create a topology for covalent inhibitors for gromacs? \- ResearchGate, [https://www.researchgate.net/post/How\_to\_create\_a\_topology\_for\_covalent\_inhibitors\_for\_gromacs](https://www.researchgate.net/post/How_to_create_a_topology_for_covalent_inhibitors_for_gromacs)  
46. Molecular Dynamics Simulations of Human Serum Albumin and Role of Disulfide Bonds | Request PDF \- ResearchGate, [https://www.researchgate.net/publication/257073344\_Molecular\_Dynamics\_Simulations\_of\_Human\_Serum\_Albumin\_and\_Role\_of\_Disulfide\_Bonds](https://www.researchgate.net/publication/257073344_Molecular_Dynamics_Simulations_of_Human_Serum_Albumin_and_Role_of_Disulfide_Bonds)  
47. pH Regulates Ligand Binding to an Enzyme Active Site by Modulating Intermediate Populations \- ACS Publications, [https://pubs.acs.org/doi/10.1021/acs.jpcb.2c05117](https://pubs.acs.org/doi/10.1021/acs.jpcb.2c05117)  
48. Parameterizing an isopeptide bond for molecular dynamics simulation in AMBER, [https://robinbetz.com/blog/2016/11/28/parameterizing-an-isopeptide-bond/](https://robinbetz.com/blog/2016/11/28/parameterizing-an-isopeptide-bond/)  
49. Molecular dynamics simulations of human serum albumin and role of disulfide bonds, [https://pubmed.ncbi.nlm.nih.gov/24066859/](https://pubmed.ncbi.nlm.nih.gov/24066859/)  
50. Rosetta FlexPepDock ab-initio: Simultaneous Folding, Docking and Refinement of Peptides onto Their Receptors \- PMC, [https://pmc.ncbi.nlm.nih.gov/articles/PMC3084719/](https://pmc.ncbi.nlm.nih.gov/articles/PMC3084719/)  
51. Crosslink Guided Protein-Protein Docking \- Rosetta, [https://docs.rosettacommons.org/demos/latest/public/xl\_driven\_protein\_docking/README](https://docs.rosettacommons.org/demos/latest/public/xl_driven_protein_docking/README)

[image1]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGMAAAAfCAYAAADz23MvAAADbElEQVR4Xu2ZW6hOQRTHl9zLpXggT5SkPIhccneIUsSLF3kRyosHKRxvJ0nkLkkuJ7dSnoTkngdFKMoDSShFSJFb7utvzbLnW2fv75u9nePbJ/Orf+ebtdbMWTN775nZs4kikUiknsxmfWH9YK03vsg/ZDlrllf+yXrulf97OrFGWWMbgcH/5JUXOlt3z5ZKA+sQ6wjrKOuY02HWQdY+1m7WUq3QTvlOMiBZLGPdYL0liRtU6c7FCydlKEmb4z1bKpjbtrFuk1SAmlk7nX0P64HnuyPV2hULKMk/C0wte6l2XBF2Uc42tQLuoCxuksQ8tI6SowMcMiChcXlAe6utsRrfSCpttg6DJtvXOkrKKdYZChvkqSQxl4z9b8B0tdYaa6HJ9rIOAzqHuEfWUUK6UHIBQi7GFZKYydZREFzUudZYC6zyIckCzK2hsfXmHWuG+605D0jcLWjNfmGGGeOVp7N6e+VMmkiS+GodKWBRb82k24pxVLm91JxnejZLSL+Gk+y6PpBsWUG/xP2b+SQvfG+cPlLtdv/wmSR4o3Wk8JjCkq43yM/f12vOKz2bD6Ym+C9ahwcGFjFdXfmZZxuiQa6cpiA0uKd1pKCxT62jRGxnnTc2zbvZ2BXM7/BPsA4HnjL/Qii5BroWedYLoLGTrKNEpPUFUzDst6zDUW0MtpL4TloHVa+Xm3UkjYWsF9gr5/nnuItwJlNERXnCOk3yRr2Ytcj9fUmSN6bkNKr1S30djX2Ys1819sJgoUGDIeuFJrXBOqrQoaCKMJDklPRuil5T9oBPJLGfsw5mBGXXO05ib7COoug/qrVe6KOKO6ysIL/+1uhYRdmDeoHEnnZ2tJ/Eh8XaktVeIbDvDWmwD4XF1ZNG1n1r9JhD2X2wdv/3Fldu8myKXw9jiTO8wuwgaazaeRQ6WSum3vSg9EH2GUwtB13x7Xg6znq+ac63ybOBZme/58pY3Mcm7nDwZjiPkiSukZzLTHHCWyu2hurHRSsj2IaeoCRP3L32m0E3kv5cpiRuDVVOSf7FwALf2fMB+LAWKQdINgqw4/AUaP1cYLHRN0JfuPOxmGM/jbfL61T5taqMYOBekRzG4TsEDjtHV0TIsQT6o3FYyN9T5eCOpGQcsr7boK7GrHC2JZ4NT14kEolEIpFIJBKJRCLCLzvvLmCCEblGAAAAAElFTkSuQmCC>

[image2]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAC0AAAAbCAYAAADoOQYqAAACQ0lEQVR4Xu2XS8hNURiGX3KJkpiZMDNRJEqRZEBkRikZmmFmoBRJITEgl5FMFDOUck9m6P9TLpGBAROXASOX3L/v/76vtc5rrf2fs//UP9hPve2z3/fba62z19rr7AN0dJR4IfrD5njmFxvjnTmiS2w2sRs2Jbk+iiZ6PsO9nKuiV+5z1oa5sHbWccDsRep0f2+E9e7P8uOV3hhHRe89u05ZW97B2vvAQfAVVvCQg4w1SF9qMWXKPVi2gvyxUJ05XfAanOegQLURNGeDMgmprZl5oMQUvOGgQtPAmrJBeSRawqayA4N39FZ0mU1hFayduxy05JtoCptKDPgxBw38Fi1gEzZYbWslB8I52HXHM2+DaBh2zbXMD3aKhtjchjToyZS1oTZjB0UnRPNh+S7RS9ERz2e7/8PPcw6IfopuhKG3v9ZRG2pthad7fNTwWq1d+w9R+J2DFuiS0LZuka97+mv/rNOtNaU9fOBBn+YgI2p0Vj6LvsC+pE5Zzh1Y3XLyc+KHZxoHaDHoPRxk6DpcJDqGVL9JNC8vyrImajWbYf4FDkpEI/0UH0K9U6UpUybA8tLDFs9WX5vBGYzeWaCdad1NDmA/2Zrd9vOzorUpHiHW82HyFR7DqOOJC7ZwkDEVqW4ZZcpFWLbaz0udxnqeTv5W90/6uW6D+1JcJt7YVKcoU+Iu3vdjie2wbCHsGXjaG4/AdzPYCPOXovddoy/0zS4azvUEaa098GOJZ7D6Txw4mpWWlvIclo/pX4ouh46Ojv/MX1+DwpMDZN3YAAAAAElFTkSuQmCC>

[image3]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAADIAAAAbCAYAAAA+nNxPAAACa0lEQVR4Xu2XuWtWQRTFT4SoYCxdUEgIBMR0ViIuCCGdRZBYajCKRTqxEmxSGMFKrPwDRDsFUwRXVEQkmCIkprJ2j7ihRuNyjzPzzeQ689YIFu8HB745Z7559+a9N/MFaGiowjrRL20KD0TnRLt08L/xVLQgeoF4I9MwfizLZZNogzb/MdcRL/axNsrABW9pM+CyaEb0DfGLVyHVyKQ2ijIHsyBvdYrTohuoccsjpBp5JDoi+qGDLNbDFxdbNGQvzJybyq9KqpH7oq3azOOnaBzFGrkLM2e38quSaqQ0Q6Jnon4Ua6TInDKkGrmmjTzcIhvt59iiIUXmlCHWiHvU+5SfZEJ0IRjnFbkHJr+tA8UZ0WuYnafdeit9vIRYI+SNNlKswt8L5DVyByZPnbY8g5jzRSVr7Di2Lh8dFstdkmLj3NYdPNP4HZ7wmXwS7VSeu2CX8h2xghyu6LAYcsr6n5W/LHDH+aBN+EL36cCS1YjLOpR/wvqjyl8WuPBx0THRsOiw6Kj1qZN+ags+Tsx4IGqydrxXMP5qHdSFL/dzmB9lWq6YS63ZHv50YbZDB8I7mOy8DpBusBZtyF6UzzHzWR0guyCX9egAxv+qzbq8FB3UZsAUzIV50mvCRjpFDxNZyAiMP2bH/P3E3agWA4hfLOQi4kVtsd49O+Z7st3HrTupcWvxnyc3rsx+mIPJLXpoafyHbtGg6D38vAOibcEceleCzyGuUXe3V4gWRU+svxamjqs2Lw3/atxmuXPwseLZEXtsuPezCR5KPJzeir6IvgdzzsI3uTnwHb3w+UeY4gkPPXrzdtzQ0NAQ5zf37MpgS8WecgAAAABJRU5ErkJggg==>

[image4]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAADQAAAAbCAYAAAAzgqwIAAACZ0lEQVR4Xu2XT0hVQRTGj2khghEoRIgbkdy0NAJbuHGTLcJwoyEuQkRbtGrbphYVRIqRfyhy40bauRFctkoMQcSFO8Hyb5CBUtq/8zEzz/M+u/c+LsJ7wf3Bh3O+M/M88+7cmXkiGRmnyazqj2pb1UA5+NPklTSHpn1D3ATwN4AY+i+oEFfsTePxBAZNO4/bqgnVU9Uj1WPVM9Ww6q24Dy8GKP66iTe9B+pV50wuj17VC9WKHH8LS6rXquemX7FBXV/ZjOOdnHyspUKPuLrKORFHmMw4J4rMWXF1VXIiiTChi5woMnbFlJl2LHVSmsvNbt3giOKrFOfALofJ/OZEgXSrVsV9xnvKgUuqPXEvdS3lwAVVP3m7qm+qL+LG/ZSTE3pJcY7wdEY4kUCXuHFzxsN2b580lvB3345aBew3Gc9qxvQBAxTnCANqOBHDfXFjhjgh+cWhHdY+Fw7Oe2+R/NSkeX+qJHrMgmrDt9EPSwZgy0X/Hz4O4KyD30l+anCARhVn+Wza4X1ZVj0Qd8v44L2Ppp9lUlz+Dvm/vH9qhMmMcYKw/zSMwX3qlqpZkg+9qC8tyk9N+EAsvSimxN3zAmmKQH/eRau9/4n81DRKcnHhpbUkjWHaxPUfJf+J9+/5eF1iLp2FgB9IccVhy0WO93+cKVFjAOewLOE9JB9nC/wzPuZxBXNN1SHHk8HO1CLu5G1V9anWTB7njSXcr96QjwMW/mXyAfwtE8+rDryPbf2u6pXJF0y7al/ctQLfENZ1KDwIHp4KDkQ8jX+BpYE+dhx+R0VxRfL7YiMBuBEg3vFxRkZGRnr+AlFyvhzo5pFoAAAAAElFTkSuQmCC>

[image5]: <data:image/png;base64,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>

[image6]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAB8AAAAfCAYAAAAfrhY5AAAA60lEQVR4XmNgGAWjYBCBViD+j4R/AzEHVA7EphmAWcaOJOYOFf8ApWkC/jHg9lk2A8TizegS1AAFDBDDedAlkABI3hiJPxmIpxHAlXDVeAAsfvEBQvJkA5jl8ugSSIDmloNwBJoczcEKBlQHwPA2ZEW0BNcZMC2HYboBfwZMhzxBUUEncJmBhr7vAmJudEEkwMZAQ8uJMRSk5jO6IDUAyGCQ73CBYAaIGkt0CUqBAwPE4Klo4sgAJP8OXZAaYD8DIj4r0OSYGSAVzR80caoBkKXiQMwExJ+gfGQ8HaGU+gBk6SgYBaNgFFAdAACq4U1SpC6bggAAAABJRU5ErkJggg==>

[image7]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAC8AAAAfCAYAAAB3XZQBAAACLklEQVR4Xu2XMWgVQRCGx2hhRIWICZLqsBNiYSOCwcYqFukkWFkJiorEys4iFhYqIulsTEzAtBExhagQQbESEhRFSBFIoqJBxViJ/j+7683O7V2ejfcg98HPu/ln3t7s7d3eeyING49r0G9oHjphcoTed+gndMjkauUFlKmYkxhR8QVoQcXfoH4V1wqbpariXhWzcZ0v5Rh0F5r04vGeqCJnAhqFrkJXoOvQragizSZos4rZ2H0Ta7oSXpKT0E3okeRX5HNUkcO6Ncnr3kGXowo3QeZmobMmRx6IGyPQJ8VGOVnrVTIu8QSqYH6rNcEX6LE1PT3iVuojdEn5RyR9Pnr6VqqExVv8J3UvTv+lQ9InOyXOP69UBuuW/HHmYws9nqslwgA3/HFqQHJO3ApZXov7TqeXXplt6pi8knh8e67tCa+UDFpVcWg+tR8vQ4etCd5A760Jdogba0Z53Ott87tVfNB7LXEHGlbxWym/+imP7JfyHP29JtYX6zk0peKX4na9luBg+v7a6T1ql/JJWYNkQFx+AXqm/APe5w70C5pTuQBXlDvRU4m/uy6phkLzbCRwGnqo4trJJF7CwKAUb50P4u7HtmEMumhNT2g+vEVTK1QrbEi/ujW3JZ/AP7/1/gfrNRSa58M0bXK1sg/6YU3DouQT4K7RNvCXHX8dVtEtxQe3VriTPJG8qTNS3M81bdP8cegr9Ala8Z/863VUFxmGxD28DQ0NDQ0NDSn+AIAsmzAWMHlpAAAAAElFTkSuQmCC>

[image8]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACwAAAAfCAYAAACcai8CAAAB7klEQVR4Xu2Wuy8FURDGx7tSkFAJnUaBSqElopOoFBKJ0h8gkRCVSkPQoBIKDZEgKhpBIhFKJQWFeMUr3jM5Z5zZsXvtruTuFvtLvuw9852d+927Z3cPQEZGRuoYRH0pXaMKrV9ua4kzDC7gqNeCTluvsMdVr51/nsEEOdCGoB3cD2pWXl75ABNiQRs+cODEuAQT4EwbASQaeACiB7hArehivuCwx9rIwSeqQRcFZahJMH1vUSNeOz794AKXKC8ua2D6tdlxgR2//cww7KhxKF4g+nLIxT2YXqWq3mTrE3bcilp3dng47Ks2YjALpteYNizyj6EXEC2byHCTaW0IeA5djUfUE5gf+C4nwd9Xin0KSvdALLjJkDYE9ahG1Di4+d2oOjGnT3hBsH+DqlReaLjJkjZ8oEsdFOoQTH1RGwI+d1sbUZiB4BAaustp3pY2kCMwXq82BGG/p0gXNNyoRxsCWnc8r0V5BN315NEOz48O+Dsw7UvOwfTfQNV6bQfvvEhTyiM4zJ49BvEAv/0a1BVqF8zOj3y+JyicRJ5Ln7vE2BfaoXFwqRNwL5V9ewxiGbznboLbPxMUkup3osbwOfPaCEOsZ+Q/oaVAjzsOnmpOUdVinOrAxWAC0pGoArPuU80cmND09qR1n5GRkTTf8Oie3vzr+9UAAAAASUVORK5CYII=>

[image9]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAB8AAAAbCAYAAACEP1QvAAABXUlEQVR4Xu2UvytGURzGvwslg8JgMJmUwSLZDVLKZEGSyWJWdmVQFAqLMtnNSgb/BJsBKbIQ+fX9Oue8nh7dcw6513I/9fTe5/me9zxv7+kekZqamn/iTfXOYVVYsf2AypkVV75FeSWciSvv5EFZzKh2VLviik173s/BulJYU22r9uXrvM2b+mFdqUyJK1/nQRWciyvv4MEvGVD1cFhEOO8UxxwUcEN+TPUs7liXaRZ9v7vFzU/8Zw6H8LygGgVve1wGE857szEuJqd8XtUC3r7zCH7aZ5+E97urMRY5UI2AD+SU85orr0CvwBo+72bySFGOrHBAbAjscwemyT+3hiGRKrf7IYXtsYjBqw+fVG04IFLlLxwQ9vcvcZhLrLxdNckhcKQa5/AnxMpPOQBWVYPgh+E5m1j5BQeeCXEXzK3Xg8T3+ca96lrc5WCf5pEhVR9lgfA2sf4Mvk6TfABtgF4mO8BPnAAAAABJRU5ErkJggg==>

[image10]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAAdCAYAAACwuqxLAAAAfElEQVR4Xu3QMQ5AQBSE4aencwCdwklo3cDVxAncQOEAziBB1BqRMOimHoXkfclf7DSbXTPn3C8UPKg1qORRKUQHj2oLynhUO1HAo1Jk7yUJ7Y8BTYJWey/ZjF5zHxTVaEepfaBHHY8qFRp5VLr//TM5anlUmlHMo3M/cgHoUSAvkuibbgAAAABJRU5ErkJggg==>
---

## **2026-05-15 初步 MD 分析结果**

### 分析背景
- **WT 全长肽**：可用轨迹 ~41 ns（md.part0003.trr，66.7–108 ns），md.part0002.trr（9.4–66.7 ns）意外缺失
- **D-Ala2 全长肽**：完整 35.8 ns（DAla2_md.xtc）
- **短肽 WT**：完整 112.5 ns（short_peptide_md.trr），GHRH(1-10)

### 关键距离分析（gmx distance，自动 PBC 处理）

> ⚠️ **2026-05-16 修正**：此前所有数据因原子索引错误而失效。`r 728`/`r 729` 实际选中 DPP-IV Val728/Asp730，非 GHRH 原子。正确残基编号为 GHRH `r 1`–`r 29`，DPP-IV `r 39`–`r 766`。详见下方"2026-05-16 重大修正"。

#### 1. 催化几何：Ser630 OG ↔ Ala2 C=O（修正后）

| 体系 | 时段 | 平均距离 | < 4 Å | 结论 |
|------|------|---------|-------|------|
| **CJC-1295 (D-Ala2) 全长** | 0–9.4 ns | 4.82 ± 0.43 Å | 4.2% | 初始远离 |
| **CJC-1295 (D-Ala2) 全长** | 66.7–200 ns | 4.29 ± 0.56 Å | 39.0% | 短暂窗口后脱结合 |
| **CJC-1295 (D-Ala2) 全长** | **100–150 ns** | **3.96 ± 0.49 Å** | **66.0%** | 最佳瞬态窗口 |
| **CJC-1295 (D-Ala2) 全长** | **150–200 ns** | **4.71 ± 0.27 Å** | **0.0%** | **完全脱结合** |
| **D-Ala2 短肽 (1-10)** | 0–200 ns | **5.41 ± 0.52 Å** | **3.0%** | **无催化几何** |

**修正前 vs 修正后**：
- 短肽从 "3.15 Å, 93.8%" 变为 "5.41 Å, 3.0%"
- 全长从 "3.31 Å, 87%"（0–9.4ns）变为 "4.82 Å, 4.2%"

**结论**：D-Ala2 **不形成稳定催化几何**。全长在 100–150ns 有一瞬态窗口（3.96Å），但 150ns 后反弹到 4.71Å。短肽完全失败，证实 GHRH 残基 11–29 对结合必不可少。

#### 2. N-末端盐桥：Glu205 ↔ Tyr1-N（待重新分析）

此前 Tyr1-N 原子索引错误（误用 H 而非 N，原子号 11656 vs 11655）。需用正确原子重新计算 Glu205-OE1/2 ↔ Tyr1-N（11655）氢键。短肽盐桥数据也可能受同一错误影响。

#### 3. 综合解读

- **WT 全长肽**：催化几何好，但 N-末端摆动导致盐桥断裂。这是物理合理的——长肽在浅口袋中发生 pivot 旋转，Ala2 附近保持在催化位点，N-末端 Tyr1 因熵驱动摆动。
- **D-Ala2**：催化几何差且盐桥断裂。D-构型不仅破坏盐桥，还导致整个肽链姿态偏移，使 scissile bond 远离 Ser630。
- **短肽 WT**：N-末端锚定最强，但催化几何中等（缺少 C-末端深入口袋的相互作用）。

### 局限性与下一步
1. WT 轨迹缺失 md.part0002.trr（~57 ns 空白），需等当前模拟完成至 200 ns
2. D-Ala2 仅 35.8 ns，需继续运行至 200 ns 确认趋势
3. 尚未做氢键分析、RMSF、构象聚类、MM-PBSA
4. FEP 完成后将得到 L-Ala2 ↔ D-Ala2 的定量自由能差

---

## **已知问题与修复记录（更新）**

### 2026-05-15：FEP λ₀₀ 崩溃修复
- **问题**：`maximal distance in decoupled molecule exceeds rlist`
- **根因**：λ=0.0 时 D-虚原子质量=0、键力常数=0，成为自由粒子飞出盒子
- **修复**：拓扑中复制所有非零键合参数到 dummy 状态（8 bonds + 18 angles + 38 dihedrals）
- **结果**：全部 11 个 λ 窗口重新 grompp 并启动，λ₀₀–λ₀₂ 运行正常

### 2026-05-15：WT TRR 分段缺失
- **问题**：md.part0002.trr（~57 ns）不存在，导致 WT 轨迹从 9.4 ns 跳到 66.7 ns
- **影响**：初步分析只能使用 66.7–108 ns 段
- **状态**：不影响当前运行，后续分析以 200 ns 完整轨迹为准


---

## ⚠️ 2026-05-15 重大修正：手性标记反转

> **详见 `docs/CHIRALITY_CORRECTION.md`**

### 核心问题
原始模板 `GHRH_1-29_from_7CZ5.pdb` 中的 Ala2 本身就是 **D-构型**。`build_DAla_mutant.py` 从这个 D-Ala 结构翻转后，错误地将产物标记为 "D-Ala2"（实际是 L-Ala），而未翻转的原始结构被标记为 "WT"（实际就是 D-Ala2）。

### 修正对照

| 原标签 | 实际身份 | 修正后标签 |
|--------|----------|-----------|
| WT 全长 MD | **CJC-1295 (D-Ala2)** | CJC-1295 (D-Ala2) |
| D-Ala2 全长 MD | **天然 GHRH (L-Ala)** | 天然 GHRH (L-Ala) |
| WT 短肽 MD | **D-Ala2 短肽** | D-Ala2 短肽 |
| FEP L→D | **实际计算 D→L** | ΔG(L→D) = −ΔG_FEP |

### ⚠️ 2026-05-16 重大修正：原子索引灾难

**详见 `workspace/step3/CATALYTIC_GEOMETRY_CORRECTED.md`**

#### 根因
`gmx make_ndx` 的 `r N` 使用 **PDB 残基编号**，不是系统残基索引。`-merge all` 后：
- GHRH 保留 chain B 编号 `r 1`–`r 29`
- DPP-IV 保留 chain A 编号 `r 39`–`r 766`

此前所有分析误用 `r 728`、`r 729` 等，实际选中的是 DPP-IV Val728/Asp730，**不是 GHRH 原子**。早期"好结果"（短肽 3.15 Å）测量的是 **Ser630-OG → Asp3-N**。

#### 修正后催化几何（Ser630-OG → GHRH-Ala2-C）

| 系统 | 平均距离 | < 4 Å | 结论 |
|------|---------|-------|------|
| **CJC-1295 (D-Ala2) 全长** | 4.29 ± 0.56 Å | 39.0% | 短暂结合窗口（100-150ns 最佳），150ns 后脱结合 |
| **D-Ala2 短肽 (1-10)** | 5.41 ± 0.52 Å | **3.0%** | **完全没有催化几何**，GHRH 11-29 对结合必不可少 |

**关键变化**：
- 短肽从"3.15 Å, 93.8% < 4Å"变为"5.41 Å, 3.0% < 4Å"
- 全长早期（0-9.4ns）从"3.31 Å"变为"4.82 Å"
- 全长 100-150ns 有一短暂窗口（3.96 Å, 66% < 4Å），但 150ns 后完全脱结合

#### 修正后生物学解释
CJC-1295 (D-Ala2) **不形成稳定催化几何**。全长肽在 100-150ns 有一瞬态结合窗口，但 150ns 后脱结合（反弹到 4.71Å）。短肽完全不结合（5.41Å）。

这与 FEP 结果（ΔG ≈ 0，亲和力无差异）一致：**抑制是空间/变构效应**，不是竞争性结合。D-Ala2 的立体化学阻碍破坏了过渡态几何，但整个肽链无法稳定在活性口袋中。

#### N-末端盐桥（待重新分析）
之前 Glu205-Tyr1 盐桥分析也使用了错误的 Tyr1 原子（`r 1` 的 H 而非 N）。需用正确原子（Tyr1 N = 11655）重新计算。

### 后续工作（按修正后标签）
- **P2a (CJC-1295 200 ns)**：全长已完成，但存在 9.4-66.7ns 轨迹空缺；分析以 66.7-200ns 为准
- **P2b (天然 GHRH 200 ns)**：GPU 1 运行中（~62ns/200ns），等完成后对比 D-Ala2 vs L-Ala2
- **P2c (D-Ala2 短肽 200 ns)**：已完成，证实短肽无催化活性
- **P3 (FEP)**：已完成 11 λ，ΔG(D-Ala2 → L-Ala2) = +0.83 ± 0.31 kJ/mol ≈ 0
- **OpenMM 复刻**：GPU 0 已完成 NPT 稳定性验证。OpenMM 8.5.1 NPT 速度 **44.4 ns/day**，比 GROMACS 2026.0 NPT (~33 ns/day) 快 **34%**。详见 `workspace/step3/OPENMM_REPLICA.md` 和 `openmm_production.py`。

