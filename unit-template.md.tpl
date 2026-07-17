
<div style="text-align: center; page-break-before: always; page-break-after: always; padding: 12% 5% 0 5%;">

  <h1 style="font-size: 2.8rem; margin-bottom: 1.8rem;">第一单元 统计信号处理理论基础</h1>
  <p style="font-size: 1.2rem; color: #555; margin-bottom: 0.2rem;">（讲义 01–09）</p>
  
  <p style="font-size: 1.2rem; max-width: 750px; margin: 2rem auto 2.5rem auto; line-height: 1.9; color: #333;">
    本单元从概率论与DSP基础出发，沿着<strong>“估计什么 → 多好算好 → 怎么估计 → 如何高效求解”</strong>的主线，系统构建统计信号处理的理论基石，覆盖经典估计、最优滤波与正则化核方法。
  </p>
  
 <div style="display: flex; align-items: center; justify-content: center; gap: 12px; margin: 1.8rem auto;">
  <span style="flex:1; max-width:80px; height:1px; background: linear-gradient(to right, transparent, #888);"></span>
  <span style="display:inline-block; width:6px; height:6px; background:#38bdf8; border-radius:50%;"></span>
  <span style="flex:1; max-width:80px; height:1px; background: linear-gradient(to left, transparent, #888);"></span>
</div>
  
  <div style="display: inline-block; text-align: left; font-size: 1.15rem; line-height: 2.6;">
    <b>01</b> 概率论与随机过程复习（测度、条件期望、方差分解）<br>
    <b>02</b> 确定性信号处理基础（最小相位、谱分解、匹配滤波）<br>
    <b>03</b> 最小方差无偏估计（MVUE）与充分统计量<br>
    <b>04</b> Cramér-Rao下界（CRLB）与线性估计<br>
    <b>05</b> 维纳滤波（Wiener Filtering）与正交投影<br>
    <b>06</b> 线性预测编码（LPC）与Levinson-Durbin递推<br>
    <b>07</b> 卡尔曼滤波（Kalman Filtering）及其扩展<br>
    <b>08</b> 正则化（岭回归、对角加载、偏差-方差权衡）<br>
    <b>09</b> 支持向量机（SVM）与核方法
  </div>

</div>

<div style="page-break-before: always; padding: 8% 8% 0 8%;">

  <h2 style="text-align: center; font-size: 2.2rem; margin-bottom: 2rem;">单元逻辑与内容梳理</h2>
  
  <p style="font-size: 1.1rem; line-height: 2; text-indent: 2em;">
    <strong>第一单元（01–09）</strong>是后续所有自适应滤波算法（第二单元）的“理论锚点”。整个单元按“基础→理论→方法→扩展”四个层次递进：
  </p>
  <p style="font-size: 1.1rem; line-height: 2; text-indent: 2em;">
    <strong>基础层（01–02）</strong>：回顾概率论与随机过程（均值、相关、功率谱、各态历经），以及数字信号处理的核心频域工具（Z变换、DTFT、最小相位与谱分解）。它们为后续所有估计与滤波提供了描述信号的语言。
  </p>
  <p style="font-size: 1.1rem; line-height: 2; text-indent: 2em;">
    <strong>理论层（03–04）</strong>：回答“什么是最优估计”以及“最优能有多好”。03从充分统计量、Rao-Blackwell到Lehmann-Scheffé定理，给出MVUE的存在性框架；04通过Fisher信息量与CRLB，给出无偏估计方差的理论下界，并延伸至线性估计与投影。两者共同构成经典估计理论的完整骨架。
  </p>
  <p style="font-size: 1.1rem; line-height: 2; text-indent: 2em;">
    <strong>方法层（05–07）</strong>：从平稳信号的最优线性估计（维纳滤波）出发，过渡到时变状态的最优递推（卡尔曼滤波），并在中间插入LPC与Levinson-Durbin高效算法，展现从“频域闭式解”到“时域递推估计”的演进脉络。
  </p>
  <p style="font-size: 1.1rem; line-height: 2; text-indent: 2em;">
    <strong>扩展层（08–09）</strong>：转向线性估计的约束与升维。08通过岭回归（Tikhonov正则化）揭示偏差-方差权衡与对角加载的工程本质；09借助SVM引入对偶与核技巧，完成从线性到非线性的跨越，并与正则化理论形成呼应。
  </p>
  <p style="font-size: 1.1rem; line-height: 2; text-indent: 2em;">
    整个第一单元为第二单元（自适应滤波算法：LMS/NLMS/RLS/QR-RLS/SVD/PCA/TLS等）提供了完整的理论预备——第二单元的所有算法，本质上都是在本单元所定义的最优准则、性能界与约束框架下，对“如何实时、高效地求解”这一工程问题的具体回应。
  </p>

</div>
<div style="page-break-before: always;"></div>