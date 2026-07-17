<div style="page-break-before: always; padding: 8% 8% 0 8%;">
 <h1 id="第六讲-时频分析基础-短时傅里叶变换与测不准原理" style="text-align: center; margin-bottom: 2rem; border-bottom: none; display: block;">第六讲 时频分析基础：短时傅里叶变换与测不准原理</h1> 
 <div style="display: flex; align-items: center; justify-content: center; gap: 12px; margin: 1.8rem auto;">
  <span style="flex:1; max-width:80px; height:1px; background: linear-gradient(to right, transparent, #888);"></span>
  <span style="display:inline-block; width:6px; height:6px; background:#38bdf8; border-radius:50%;"></span>
  <span style="flex:1; max-width:80px; height:1px; background: linear-gradient(to left, transparent, #888);"></span>
 </div>
</div>


<!-- # 第六讲 时频分析基础：短时傅里叶变换与测不准原理 -->

## 1. 背景

### 1.1 窄带假设的失效：从阵列信号到时频分析的必然过渡

在前面几讲中，已经系统研究了阵列信号处理的基本内容。从 MUSIC 算法到 ESPRIT 算法，再到相关信号处理与阵列处理的 CRLB，核心问题始终是同一类型：**空间中存在多个信号源，它们的频率已知（窄带），需要估计它们的方向。**

这些问题反复依赖一个核心假设——**窄带假设**：

\[
s(t) = a(t) \exp(j 2\pi f_c t), \qquad B \ll f_c
\tag{6.1}
\]

这个假设将"时延"转化为"相移"，从而把角度估计问题转化为相位差估计问题。窄带假设是 MUSIC、ESPRIT 等算法能够工作的前提——它把问题从"估计一个随时间变化的信号的到达方向"简化成了"估计一个恒定频率的信号的到达方向"。

**但现实中的信号并不总是窄带的。**

语音信号、音乐信号、雷达回波、生物医学信号、机械振动信号——这些信号的频率成分随时间变化。一段语音中，不同时刻有不同的音素；一段音乐中，不同时刻有不同的音符；一个移动目标的雷达回波，其多普勒频率随时间变化。

如果仍用窄带假设去处理这类信号，就会遇到根本性困难：**信号在频域中占据了较宽的带宽，无法用一个单一的载频来描述其方向信息。** 阵列信号处理的窄带方法无法直接应用。

更严重的是，实际问题往往不仅要求知道"信号从哪个方向来"，还要求知道"在某个特定时刻，信号包含了哪些频率成分"。这是一个更高维度的信息提取问题——需要同时刻画信号的时间演变和频率构成。

时频分析正是要回答这个问题。

### 1.2 时频分析的数学本质：信号在时间-频率平面上的能量分布

时频分析的核心任务：**在时间和频率两个维度上同时描述信号的特性。**

经典傅里叶分析能给出"信号中有哪些频率成分"，但无法判定"这些频率成分在什么时候出现"。经典时间域分析能描述"信号随时间如何变化"，但无法揭示"每个时刻的频率结构"。

时频分析将这两个视角合二为一。

**直观案例：音乐信号**

想象一段钢琴演奏的录音。音乐包含不同音符（频率）在不同时刻（时间）出现。如果你用傅里叶变换对整个录音做频谱分析，你会得到一个包含所有音符频率的频谱，但你无法知道哪个音符出现在哪个时刻。如果你只看时域波形，你能看到声音在哪些时刻响起和停止，但你无法分辨出是哪个音被弹奏。

![spectrogram](assets/06/01.png) 
时频分析提供了一种表示方式：**以时间为横轴、频率为纵轴、颜色表示能量强度**，形成一张"时频谱"（spectrogram）。在这张图上：

- 水平方向 → 时间演化
- 垂直方向 → 频率成分
- 颜色深浅 → 能量大小

这是一张二维的"信号地图"。沿着水平方向看，你能看到频率成分如何随时间变化；沿着垂直方向看，你能看到某个频率成分在不同时刻的能量分布。

从阵列信号处理到时频分析的视角转换：

在阵列信号处理中，面对的是空间维度——多个传感器在不同空间位置采集同一个信号的不同版本。在时频分析中，面对的是时间-频率维度——同一个传感器在不同时刻采集的信号，需要同时理解它的时间结构和频率结构。

两者的共同点在于：处理的都是二维结构。阵列信号处理是"空间 × 时间"，时频分析是"时间 × 频率"。

### 1.3 时频分辨率的下界：测不准原理的本质与来源

时频分析面临一个根本性的限制：**无法同时获得完美的时间分辨率和完美的频率分辨率。**

这就是信号处理中的**测不准原理**（Uncertainty Principle），来源于傅里叶变换的数学性质。

**直观案例：看乐谱与听音乐**

想象你正在听一首乐曲：

- 如果你想要**精确地知道某个时刻弹奏了哪个音符**（高时间分辨率），你就需要把注意力集中在一个极短的时间窗口内。但时间窗口越短，你听到的音符越不完整——你无法确定它的音高（频率）。一个极短的脉冲在频域中展开成一个很宽的频谱，你无法准确地判断它的频率。

- 如果你想要**精确地知道某个音符的音高**（高频率分辨率），你就需要听足够长的时间——听好几个周期才能判断它的频率。但时间越长，你越无法确定这个音符具体是在哪个时刻被弹奏的。

这就是测不准原理的本质：**时间分辨率与频率分辨率是不可兼得的。** 这是一个根源性的物理限制，与算法无关，与仪器精度无关。

**数学表达：**

对于任意信号 \( x(t) \) 及其傅里叶变换 \( X(f) \)，定义其时间宽度 \( \Delta t \) 和频率宽度 \( \Delta f \) 为：

\[
(\Delta t)^2 = \frac{\displaystyle \int_{-\infty}^{\infty} (t - t_0)^2 |x(t)|^2 dt}{\displaystyle \int_{-\infty}^{\infty} |x(t)|^2 dt}, \qquad
(\Delta f)^2 = \frac{\displaystyle \int_{-\infty}^{\infty} (f - f_0)^2 |X(f)|^2 df}{\displaystyle \int_{-\infty}^{\infty} |X(f)|^2 df}
\tag{6.2}
\]

其中 \( t_0 \) 和 \( f_0 \) 分别是信号在时域和频域的"中心"。

测不准原理说：

\[
\boxed{\Delta t \cdot \Delta f \ge \frac{1}{4\pi}}
\tag{6.3}
\]

等号成立当且仅当信号是高斯函数（Gabor 函数）。这意味着：信号的时宽和频宽的乘积存在一个下限——两者不可能同时任意小。

时频分析方法必须在时间分辨率和频率分辨率之间做出权衡。短时傅里叶变换使用窗函数截取信号——窗越短，时间分辨率越好，但频率分辨率越差；窗越长，频率分辨率越好，但时间分辨率越差。不存在一个"最佳"窗长，只有针对具体问题的最优选择。

本讲从傅里叶变换的基本性质出发，逐步引入短时傅里叶变换（STFT），并深入分析测不准原理对时频分析的深远影响。


## 2. 傅里叶变换及其重要的性质

短时傅里叶变换作为时频分析的一阶工具，需要用到之前学过的傅里叶分析工具。本节系统回顾傅里叶变换的定义与核心性质，重点推导泊松求和公式及其在采样定理中的关键应用，为后续时频分析打下数学基础。

---

### 2.1 傅里叶变换的定义

**连续时间傅里叶变换**（CTFT）将时域信号 \(x(t)\) 映射到频域 \(\tilde{X}(\omega)\)，其定义为

\[
\tilde{X}(\omega) = \mathcal{F}\{x(t)\} = \int_{-\infty}^{\infty} x(t) \, \exp(-j\omega t) \, dt , \tag{6.4}
\]

其中 \(\omega\) 为角频率（单位：弧度/秒）。对应的**傅里叶逆变换**为

\[
x(t) = \mathcal{F}^{-1}\{\tilde{X}(\omega)\} = \frac{1}{2\pi} \int_{-\infty}^{\infty} \tilde{X}(\omega) \, \exp(j\omega t) \, d\omega . \tag{6.5}
\]

**推导思路**：从周期信号的傅里叶级数出发，令周期 \(T \to \infty\)，求和变为积分，从而导出上述变换对。以下直接采用定义作为推导起点。

---

### 2.2 傅里叶变换的性质

#### 2.2.1 成对出现的对称性

傅里叶变换与逆变换具有高度对称性。若定义 \(\tilde{X}(\omega) = \mathcal{F}\{x(t)\}\)，则有以下配对关系（变换对）：

| 时域 \(x(t)\) | 频域 \(\tilde{X}(\omega)\) |
| :--- | :--- |
| \(x(t)\) | \(\tilde{X}(\omega)\) |
| \(\tilde{X}(t)\) | \(2\pi x(-\omega)\) （对偶性） |
| \(x(-t)\) | \(\tilde{X}(-\omega)\) |
| \(x^*(t)\) | \(\tilde{X}^*(-\omega)\) （共轭对称性，实信号时） |

特别地，当 \(x(t)\) 为实偶函数时，\(\tilde{X}(\omega)\) 也为实偶函数。这种"成对出现"的性质使得变换与逆变换在形式上几乎一致，便于相互转换。

#### 2.2.2 Parseval 定理

Parseval 定理说明傅里叶变换是正交变换（酉变换），正交变换保范数、保能量。两个函数 \(f, g\) 对应的傅里叶变换记作 \(\tilde{f}, \tilde{g}\)，则有内积保持关系

\[
\langle f, g \rangle = \frac{1}{2\pi} \langle \tilde{f}, \tilde{g} \rangle , \tag{6.6}
\]

其中内积定义为 \(\langle f, g \rangle = \int_{-\infty}^{\infty} f(t) \, g^*(t) \, dt\)，频域内积为 \(\langle \tilde{f}, \tilde{g} \rangle = \int_{-\infty}^{\infty} \tilde{f}(\omega) \, \tilde{g}^*(\omega) \, d\omega\)。

**详细计算**（推导）：

\[
\begin{aligned}
\langle f, g \rangle &= \int_{-\infty}^{\infty} f(t) \, g^*(t) \, dt \\
&= \int_{-\infty}^{\infty} \left[ \frac{1}{2\pi} \int_{-\infty}^{\infty} \tilde{f}(\omega) \exp(j\omega t) \, d\omega \right] g^*(t) \, dt \\
&= \frac{1}{2\pi} \int_{-\infty}^{\infty} \tilde{f}(\omega) \left[ \int_{-\infty}^{\infty} g^*(t) \exp(j\omega t) \, dt \right] d\omega \\
&= \frac{1}{2\pi} \int_{-\infty}^{\infty} \tilde{f}(\omega) \left[ \int_{-\infty}^{\infty} g(t) \exp(-j\omega t) \, dt \right]^* d\omega \\
&= \frac{1}{2\pi} \int_{-\infty}^{\infty} \tilde{f}(\omega) \, \tilde{g}^*(\omega) \, d\omega = \frac{1}{2\pi} \langle \tilde{f}, \tilde{g} \rangle .
\end{aligned}
\]

当 \(f = g\) 时，得到能量守恒：

\[
\int_{-\infty}^{\infty} |f(t)|^2 dt = \frac{1}{2\pi} \int_{-\infty}^{\infty} |\tilde{f}(\omega)|^2 d\omega .
\]

#### 2.2.3 Poisson Summation Formula（泊松求和公式）

从时域周期延拓出发，探究其与频域采样的关系。考虑信号 \(x(t)\)，构造周期为 \(T\) 的函数

\[
g(t) = \sum_{k=-\infty}^{\infty} x(t + kT) , \quad g(t) = g(t+T) .
\]

这个时域表示能否直接和某个频域表示建立关系呢？利用傅里叶级数对 \(g(t)\) 展开，得到

\[
g(t) = \sum_{k=-\infty}^{\infty} \alpha_k \exp\left(j \frac{2 k \pi}{T} t\right) .
\]

其中傅里叶系数 \(\alpha_k\) 的计算如下：

\[
\begin{aligned}
\alpha_k &= \frac{1}{T} \int_{0}^{T} g(t) \exp\left(-j \frac{2 k \pi}{T} t\right) dt \\
&= \frac{1}{T} \int_{0}^{T} \sum_{n=-\infty}^{\infty} x(t + nT) \exp\left(-j \frac{2 k \pi}{T} t\right) dt .
\end{aligned}
\]

交换求和与积分（在适当收敛条件下），令 \(u = t + nT\)，则 \(t = u - nT\)，\(dt = du\)，且当 \(t\) 从 \(0\) 到 \(T\) 时，\(u\) 从 \(nT\) 到 \((n+1)T\)。于是

\[
\begin{aligned}
\alpha_k &= \frac{1}{T} \sum_{n=-\infty}^{\infty} \int_{nT}^{(n+1)T} x(u) \exp\left(-j \frac{2 k \pi}{T} (u - nT)\right) du \\
&= \frac{1}{T} \sum_{n=-\infty}^{\infty} \int_{nT}^{(n+1)T} x(u) \exp\left(-j \frac{2 k \pi}{T} u\right) \exp(j 2 k \pi n) \, du .
\end{aligned}
\]

由于 \(\exp(j 2 k \pi n) = 1\)，求和与积分合并为整个实轴上的积分：

\[
\alpha_k = \frac{1}{T} \int_{-\infty}^{\infty} x(u) \exp\left(-j \frac{2 k \pi}{T} u\right) du = \frac{1}{T} \tilde{X}\left( \frac{2 k \pi}{T} \right) ,
\]

其中 \(\tilde{X}(\omega)\) 是 \(x(t)\) 的傅里叶变换。于是得到**泊松求和公式**：

\[
\boxed{ \sum_{k=-\infty}^{\infty} x(t + kT) = \frac{1}{T} \sum_{k=-\infty}^{\infty} \tilde{X}\left( \frac{2 k \pi}{T} \right) \exp\left( j \frac{2 k \pi}{T} t \right) } . \tag{6.7}
\]

特别地，取 \(t = 0\) 时，有常用形式

\[
\sum_{k=-\infty}^{\infty} x(kT) = \frac{1}{T} \sum_{k=-\infty}^{\infty} \tilde{X}\left( \frac{2 k \pi}{T} \right) .
\]

---

**应用举例：采样定理的推导**

以下展示泊松求和公式在采样定理中的关键作用。设 \(y(t)\) 是连续信号，需从 \(y(t)\) 采样得到离散序列 \(\{y_k\}\)。取 \(x(t) = \delta(t)\)（单位冲激），并假设采样周期为 \(T=1\)（简便起见）。此时

\[
\sum_{k=-\infty}^{\infty} x(t+k) = \sum_{k=-\infty}^{\infty} \delta(t+k) .
\]

理想采样过程可表示为连续信号乘以冲激串：

\[
\{ y_k \} \quad \text{（即采样值）} \quad \Longleftrightarrow \quad y_s(t) = y(t) \sum_{k=-\infty}^{\infty} \delta(t+k) .
\]

下图示意了一个连续函数 \(y(t)\)（光滑曲线），以及在其上等间隔（\(T=1\)）的采样竖线和对应的采样点（圆点）。竖线的高度代表该时刻的瞬时值，采样点为竖线与曲线的交点。

（图像描述：横轴为时间 \(t\)，纵轴为幅值。一条连续曲线 \(y(t)\) 上下波动；在整数点 \(t=0,\pm1,\pm2,\ldots\) 处，绘制竖直虚线，并在曲线与竖线交点处绘制实心圆点，表示采样值。）

为了书写简便，令 \((x(t))^{\hat{}} = \mathcal{F}\{x(t)\}\) 表示傅里叶变换。

\[
\left( y(t) \sum_{k=-\infty}^{\infty} \delta(t+k) \right)^{\hat{}} = \tilde{Y}(\omega) \ast \left( \sum_{k=-\infty}^{\infty} \delta(t+k) \right)^{\hat{}} ,
\]

其中 \(\ast\) 表示卷积。现在计算冲激串的傅里叶变换。直接积分

\[
\left( \sum_{k=-\infty}^{\infty} \delta(t+k) \right)^{\hat{}} = \int_{-\infty}^{\infty} \sum_{k=-\infty}^{\infty} \delta(t+k) \exp(-j\omega t) \, dt
\]

看似发散，但可利用泊松求和公式（或Parseval定理）来正则化。由泊松求和公式（6.7）取 \(x(t) = \delta(t)\)，此时 \(\tilde{X}(\omega) = 1\)，且 \(T=1\)，得到

\[
\sum_{k=-\infty}^{\infty} \delta(t+k) = \sum_{k=-\infty}^{\infty} \exp(j 2 k \pi t) .
\]

因此

\[
\left( \sum_{k=-\infty}^{\infty} \delta(t+k) \right)^{\hat{}} = \mathcal{F}\left\{ \sum_{k=-\infty}^{\infty} \exp(j 2 k \pi t) \right\} .
\]

而 \(\mathcal{F}\{\exp(j \omega_0 t)\} = 2\pi \delta(\omega - \omega_0)\)，所以

\[
\left( \sum_{k=-\infty}^{\infty} \delta(t+k) \right)^{\hat{}} = 2\pi \sum_{k=-\infty}^{\infty} \delta(\omega - 2k\pi) .
\]

代入卷积表达式：

\[
\tilde{y}_s(\omega) = \tilde{Y}(\omega) \ast \left[ 2\pi \sum_{k=-\infty}^{\infty} \delta(\omega - 2k\pi) \right] .
\]

利用卷积性质 \(f(\omega) \ast \delta(\omega - \omega_0) = f(\omega - \omega_0)\)（可由定义验证：\(\int_{-\infty}^{\infty} f(\omega - \omega') \delta(\omega' - \omega_0) d\omega' = f(\omega - \omega_0)\)），得到

\[
\tilde{y}_s(\omega) = 2\pi \sum_{k=-\infty}^{\infty} \tilde{Y}(\omega - 2k\pi) .
\]

忽略常数因子 \(2\pi\)（取决于变换定义），通常写作

\[
\boxed{ \tilde{y}_s(\omega) = \sum_{k=-\infty}^{\infty} \tilde{Y}(\omega - 2k\pi) } . \tag{6.8}
\]

上式表明：**理想采样信号的频谱是原信号频谱以 \(2\pi\) 为周期的周期延拓**。这就是**采样定理**的频域本质，其推导过程中关键的一步正是利用了泊松求和公式（将时域冲激串转换为频域冲激串）。若原信号带限且采样频率足够高（\(2\pi/T > 2\omega_{\max}\)），则周期延拓不发生混叠，从而可从采样信号完美重建原信号。



---

## 3. 短时傅里叶变换

短时傅里叶变换（STFT）是时频分析的一阶工具，其核心思想是在时域上对信号加窗，再对每一段进行傅里叶分析。为清晰描述这一过程，首先引入时域和频域的平移算子，它们是后续推导的基本工具。

---

### 3.1 时域、频域平移算子

在推导之前，先定义两个算子 \(\operatorname{T}\) 和 \(\operatorname{M}\)，它们将在本文及后续课程中大量使用。

**时域平移算子**（Time-shift operator）：

\[
\operatorname{T}_x \bigl( f(t) \bigr) \triangleq f(t - x) , \tag{6.9}
\]

其中 \(x \in \mathbb{R}\) 表示平移量。该算子将函数 \(f(t)\) 沿时间轴向右（\(x>0\)）或向左（\(x<0\)）移动。

**频域平移算子**（Frequency-shift operator，也称调制算子）：

\[
\operatorname{M}_{\omega} \bigl( f(t) \bigr) \triangleq \exp(j \omega t) \, f(t) , \tag{6.10}
\]

其中 \(\omega \in \mathbb{R}\) 表示调制角频率。该算子在时域上乘以复指数，相当于在频域上产生平移（由傅里叶变换的频移性质决定）。

---

#### 两个算子的非交换性

考察 \(\operatorname{T}_x\) 与 \(\operatorname{M}_{\omega}\) 的复合顺序是否可交换。先计算先频移后时移：

\[
\operatorname{T}_x \operatorname{M}_{\omega} \bigl( f(t) \bigr) = \operatorname{T}_x \Bigl( \operatorname{M}_{\omega} \bigl( f(t) \bigr) \Bigr) = \operatorname{T}_x \Bigl( \exp(j \omega t) f(t) \Bigr) = \exp\bigl(j \omega (t - x)\bigr) \, f(t - x) .
\]

再计算先时移后频移：

\[
\operatorname{M}_{\omega} \operatorname{T}_x \bigl( f(t) \bigr) = \operatorname{M}_{\omega} \Bigl( \operatorname{T}_x \bigl( f(t) \bigr) \Bigr) = \operatorname{M}_{\omega} \bigl( f(t - x) \bigr) = \exp(j \omega t) \, f(t - x) .
\]

比较两式，由于 \(\exp(j \omega (t - x)) \neq \exp(j \omega t)\) 一般成立（除非 \(x=0\) 或 \(\omega=0\)），因此

\[
\operatorname{T}_x \operatorname{M}_{\omega} \neq \operatorname{M}_{\omega} \operatorname{T}_x .
\]

这说明两个算子**没有交换性质**，其差异为

\[
\operatorname{T}_x \operatorname{M}_{\omega} = \exp(-j \omega x) \, \operatorname{M}_{\omega} \operatorname{T}_x .
\]

这个非交换关系在时频分析中至关重要，它体现了时间与频率观测的"不确定性"原理在代数层面的表现。

---

#### 平移后函数的傅里叶变换

下面计算复合算子 \(\operatorname{T}_x \operatorname{M}_{\omega}\) 作用在 \(f(t)\) 上所得新函数的傅里叶变换。记 \(\tilde{f}(\xi)\) 为 \(f(t)\) 的傅里叶变换（为避免与时间变量 \(t\) 混淆，这里频域变量用 \(\xi\) 表示，后续将统一使用 \(\omega\)），即

\[
\tilde{f}(\xi) = \int_{-\infty}^{\infty} f(t) \exp(-j \xi t) \, dt .
\]

要求

\[
\left( \operatorname{T}_x \operatorname{M}_{\omega} \bigl( f(t) \bigr) \right)^{\hat{}} \triangleq \mathcal{F}\left\{ \exp\bigl(j \omega (t - x)\bigr) f(t - x) \right\},
\]

注意这里的 \((\cdot)^{\hat{}}\) 表示傅里叶变换。按定义展开：

\[
\left( \operatorname{T}_x \operatorname{M}_{\omega} (f) \right)^{\hat{}}(\xi) 
= \int_{-\infty}^{\infty} \Bigl[ \exp\bigl(j \omega (t - x)\bigr) f(t - x) \Bigr] \exp(-j \xi t) \, dt .
\]

为计算这个积分，做变量代换。令

\[
u = t - x \quad \Longrightarrow \quad t = u + x, \quad dt = du .
\]

当 \(t \to -\infty\) 时 \(u \to -\infty\)，当 \(t \to +\infty\) 时 \(u \to +\infty\)，积分限不变。代入得

\[
\begin{aligned}
\left( \operatorname{T}_x \operatorname{M}_{\omega} (f) \right)^{\hat{}}(\xi)
&= \int_{-\infty}^{\infty} \exp\bigl(j \omega u\bigr) f(u) \exp\bigl(-j \xi (u + x)\bigr) \, du \\
&= \int_{-\infty}^{\infty} f(u) \exp\bigl(j \omega u\bigr) \exp(-j \xi u) \exp(-j \xi x) \, du \\
&= \exp(-j \xi x) \int_{-\infty}^{\infty} f(u) \exp\bigl(-j (\xi - \omega) u\bigr) \, du .
\end{aligned}
\]

注意到上述积分正是 \(f(u)\) 的傅里叶变换在频率点 \(\xi - \omega\) 处的取值，即

\[
\int_{-\infty}^{\infty} f(u) \exp\bigl(-j (\xi - \omega) u\bigr) \, du = \tilde{f}(\xi - \omega).
\]

因此

\[
\left( \operatorname{T}_x \operatorname{M}_{\omega} (f) \right)^{\hat{}}(\xi) = \exp(-j \xi x) \, \tilde{f}(\xi - \omega).
\]

为用算子形式表达这一结果，定义在频域上的平移和调制（由于频域变量记作 \(\omega\)，以 \(T_\omega\) 表示频域平移，\(M_{-x}\) 表示频域调制，注意符号约定）。上式右端可看作：先将 \(\tilde{f}\) 在频域平移（由 \(\tilde{f}(\xi - \omega)\) 体现），再乘以 \(\exp(-j \xi x)\)（即频域调制）。因此可以写成算子形式：

\[
\boxed{ \left( \operatorname{T}_x \operatorname{M}_{\omega} (f) \right)^{\hat{}} = \operatorname{M}_{-x} \operatorname{T}_{\omega} \bigl( \tilde{f} \bigr) , } \tag{6.11}
\]

其中 \(\operatorname{T}_{\omega}\) 作用于频域函数 \(\tilde{f}(\xi)\) 定义为 \(\operatorname{T}_{\omega}(\tilde{f})(\xi) = \tilde{f}(\xi - \omega)\)，而 \(\operatorname{M}_{-x}\) 作用于频域函数定义为 \(\operatorname{M}_{-x}(\tilde{f})(\xi) = \exp(-j \xi x) \tilde{f}(\xi)\)。

**详细中间步骤回顾**（展开书写以便理解）：

1. 写出傅里叶变换积分：
   \[
   \left( T_x M_\omega (f) \right)^{\hat{}}(\xi) = \int_{-\infty}^{\infty} \exp(j\omega(t-x)) f(t-x) \exp(-j\xi t) \, dt .
   \]
2. 变量替换 \(u = t-x\)，则 \(t = u+x\)，\(dt=du\)，得到
   \[
   = \int_{-\infty}^{\infty} \exp(j\omega u) f(u) \exp(-j\xi (u+x)) \, du .
   \]
3. 分离因子 \(\exp(-j\xi x)\)（与积分变量无关）：
   \[
   = \exp(-j\xi x) \int_{-\infty}^{\infty} f(u) \exp(j(\omega - \xi)u) \, du .
   \]
4. 改写指数为 \(\exp(-j(\xi-\omega)u)\)，即
   \[
   = \exp(-j\xi x) \int_{-\infty}^{\infty} f(u) \exp(-j(\xi-\omega)u) \, du = \exp(-j\xi x) \, \tilde{f}(\xi - \omega).
   \]
5. 最后解释右边等价于频域上先平移 \(\omega\) 再调制 \(-x\)，即式 (6.11)。

这个公式将时域的复合平移-调制操作映射为频域的相反操作（先频移，再乘以线性相位），是后续推导短时傅里叶变换时窗函数平移与调制特性的核心依据。

### 3.2 短时傅里叶变换的时-频双视角及其关系

从时间和频率两个维度分析信号的根本原因在于：经典傅里叶变换（FT）能够精确给出信号的全局频率成分，却完全丢失了时间信息——无法给出某个频率分量出现的时间。短时傅里叶变换（STFT）正是为解决这一不足而引入的一阶时频分析工具（后续还将学习二阶工具如Wigner-Ville变换，以及小波变换等）。STFT的核心思想是对信号加窗，在时间轴上逐段进行傅里叶分析，从而获得信号的时频局部化信息。

---

#### 3.2.1 从傅里叶变换到加窗傅里叶变换

经典傅里叶变换定义为

\[
\tilde{f}(\omega) = \int_{-\infty}^{\infty} f(t) \exp(-j\omega t) \, dt , \tag{6.12}
\]

它提供的是**频率局部化**（frequency localization），即 \(\tilde{f}(\omega)\) 反映了信号在整个时间轴上的频率成分，但无法获知这些频率成分发生在什么时刻。

为使傅里叶变换具备**时间局部化**（time localization）的能力，在被积函数中引入一个**窗函数** \(g(t - t')\)，使积分仅在窗函数支撑集的范围内有效。于是得到加窗傅里叶变换：

\[
\int_{-\infty}^{\infty} g(t - t') \, f(t') \exp(-j\omega t') \, dt' , \tag{6.13}
\]

其中 \(g(t - t')\) 是以时刻 \(t\) 为中心的窗函数。不同时刻 \(t\) 对应窗函数的不同位置，因此得到的傅里叶变换也随之变化，从而实现了"**随时间变化的频谱**"。

窗函数 \(g(\cdot)\) 通常应满足：在时域上具有紧支撑（或快速衰减），使得乘积 \(g(t - t') f(t')\) 仅在 \(t'\) 接近 \(t\) 时才有显著贡献。例如，取矩形窗

\[
g(\tau) = I_{[-a,a]}(\tau) = 
\begin{cases}
1, & |\tau| \le a, \\
0, & |\tau| > a,
\end{cases}
\]

代入式 (6.13) 得到

\[
\int_{t-a}^{t+a} f(t') \exp(-j\omega t') \, dt' . \tag{6.14}
\]

即加了一个窗，把原本从 \(-\infty\) 到 \(+\infty\) 的积分限改变为 \(t-a\) 到 \(t+a\)，即在时刻 \(t\) 附近进行局部傅里叶分析。时间 \(t\) 变化时，积分区间随之移动，从而得到随时间变化的频谱。

---

#### 3.2.2 短时傅里叶变换的两种等价定义

以下从两个角度来理解短时傅里叶变换：一个在时域加窗，另一个在频域加窗。这两种定义之间存在深刻的对偶关系。

**定义一（时域窗视角）**：以 \(g\) 为窗函数、\(f\) 为被分析信号的短时傅里叶变换定义为

\[
V_g f(t, \omega) \triangleq \int_{-\infty}^{\infty} f(t') \, \overline{g(t' - t)} \exp(-j\omega t') \, dt' . \tag{6.15}
\]

这里窗函数 \(g\) 在时域上以 \(t\) 为中心，\(\overline{g(\cdot)}\) 表示复共轭。该式直接刻画了"在时刻 \(t\) 附近，信号 \(f\) 的局部频谱"。

**定义二（频域窗视角）**：等价地，也可以在频域上对信号的傅里叶变换 \(\tilde{f}\) 加窗。设 \(\tilde{g}\) 为窗函数 \(g\) 的傅里叶变换，定义

\[
V_{\tilde{g}} \tilde{f}(\omega, t) \triangleq \int_{-\infty}^{\infty} \tilde{f}(\omega') \, \overline{\tilde{g}(\omega' - \omega)} \exp(-j t \omega') \, d\omega' . \tag{6.16}
\]

这里 \(\tilde{g}(\omega' - \omega)\) 是以频率 \(\omega\) 为中心的频域窗函数，\(\exp(-j t \omega')\) 实现了傅里叶变换中的相位因子。该式刻画了"在频率 \(\omega\) 附近，信号的频谱随时间的变化"。

---

#### 3.2.3 两种定义的等价关系推导

推导 \(V_{\tilde{g}} \tilde{f}(\omega, t)\) 与 \(V_g f(t, \omega)\) 之间的具体关系。这个推导展示了Parseval定理在时频分析中的深刻体现。

设傅里叶变换及其逆变换定义为

\[
\tilde{f}(\omega') = \int_{-\infty}^{\infty} f(t') \exp(-j\omega' t') \, dt', \qquad
f(t') = \frac{1}{2\pi} \int_{-\infty}^{\infty} \tilde{f}(\omega') \exp(j\omega' t') \, d\omega' .
\]

且 \(\tilde{g}(\cdot)\) 为 \(g(\cdot)\) 的傅里叶变换，即

\[
\tilde{g}(\omega' - \omega) = \int_{-\infty}^{\infty} g(\tau) \exp\bigl(-j(\omega' - \omega)\tau\bigr) \, d\tau .
\]

于是其共轭为

\[
\overline{\tilde{g}(\omega' - \omega)} = \int_{-\infty}^{\infty} \overline{g(\tau)} \exp\bigl(j(\omega' - \omega)\tau\bigr) \, d\tau . \tag{6.17}
\]

将式 (6.17) 代入定义二 (6.16)，得到

\[
\begin{aligned}
V_{\tilde{g}} \tilde{f}(\omega, t) 
&= \int_{-\infty}^{\infty} \tilde{f}(\omega') \, \exp(-j t \omega') \left[ \int_{-\infty}^{\infty} \overline{g(\tau)} \exp\bigl(j(\omega' - \omega)\tau\bigr) \, d\tau \right] d\omega' \\
&= \int_{-\infty}^{\infty} \int_{-\infty}^{\infty} \overline{g(\tau)} \, \tilde{f}(\omega') \exp\bigl(j(\omega' - \omega)\tau\bigr) \exp(-j t \omega') \, d\tau \, d\omega' .
\end{aligned}
\]

交换积分次序（在适当收敛条件下成立），整理指数项：

\[
V_{\tilde{g}} \tilde{f}(\omega, t) 
= \int_{-\infty}^{\infty} \int_{-\infty}^{\infty} \overline{g(\tau)} \, \tilde{f}(\omega') 
\exp\bigl(j\omega'(\tau - t) - j\omega \tau\bigr) \, d\omega' \, d\tau .
\]

将 \(\exp(-j\omega \tau)\) 因子移出内层积分：

\[
= \int_{-\infty}^{\infty} \overline{g(\tau)} \exp(-j\omega \tau) 
\left[ \int_{-\infty}^{\infty} \tilde{f}(\omega') \exp\bigl(j\omega'(\tau - t)\bigr) \, d\omega' \right] d\tau .
\]

注意到内层积分正是 \(\tilde{f}(\omega')\) 的傅里叶逆变换在点 \((\tau - t)\) 处的取值，即

\[
\int_{-\infty}^{\infty} \tilde{f}(\omega') \exp\bigl(j\omega'(\tau - t)\bigr) \, d\omega' = 2\pi \, f(\tau - t) .
\]

代入上式得

\[
V_{\tilde{g}} \tilde{f}(\omega, t) = 2\pi \int_{-\infty}^{\infty} f(\tau - t) \, \overline{g(\tau)} \exp(-j\omega \tau) \, d\tau .
\]

现在做变量代换：令 \(t' = \tau - t\)，则 \(\tau = t' + t\)，\(d\tau = dt'\)。当 \(\tau \to -\infty\) 时 \(t' \to -\infty\)，当 \(\tau \to +\infty\) 时 \(t' \to +\infty\)，积分限不变。代入得

\[
\begin{aligned}
V_{\tilde{g}} \tilde{f}(\omega, t) 
&= 2\pi \int_{-\infty}^{\infty} f(t') \, \overline{g(t' + t)} \exp\bigl(-j\omega (t' + t)\bigr) \, dt' \\
&= 2\pi \exp(-j\omega t) \int_{-\infty}^{\infty} f(t') \, \overline{g(t' + t)} \exp(-j\omega t') \, dt' .
\end{aligned}
\]

现在考察时域窗定义 (6.15) 在参数 \((-t, \omega)\) 处的取值：

\[
V_g f(-t, \omega) = \int_{-\infty}^{\infty} f(t') \, \overline{g(t' - (-t))} \exp(-j\omega t') \, dt' 
= \int_{-\infty}^{\infty} f(t') \, \overline{g(t' + t)} \exp(-j\omega t') \, dt' .
\]

上式右端正是前述积分中出现的那个积分。因此

\[
V_{\tilde{g}} \tilde{f}(\omega, t) = 2\pi \exp(-j\omega t) \, V_g f(-t, \omega) .
\]

在傅里叶变换采用"工程定义"（即正变换无 \(1/2\pi\)，逆变换有 \(1/2\pi\)）时，常数因子 \(2\pi\) 通常被吸收到逆变换定义中。若采用归一化酉变换（即正逆变换均含 \(1/\sqrt{2\pi}\)），则常数因子为 1。在本文的符号体系下（正变换无常数），保留 \(2\pi\)，但在实用中常写为

\[
\boxed{ V_{\tilde{g}} \tilde{f}(\omega, t) = \exp(-j\omega t) \, V_g f(-t, \omega) } , \tag{6.18}
\]

其含义是等式两边在相差一个常数因子（取决于傅里叶变换定义）下相等。

---

#### 3.2.4 结论与对偶性解读

最终得到如下重要关系：

\[
\boxed{ V_{\tilde{g}} \tilde{f}(\omega, t) = \exp(-j\omega t) \, V_g f(-t, \omega) } . \tag{6.18}
\]

这一结果与 Parseval 定理有非常相似的内核：**时域的内积（加窗）与频域的内积（加窗）通过傅里叶变换相互联系**，本质上都是酉变换保内积的体现。具体来看：

- 参数 \((\omega, t)\) 与 \((t, \omega)\) 交换了位置，体现了时间与频率的对偶性；
- 时间变量 \(t\) 变成了 \(-t\)，反映了逆变换中指数符号的反转（类似于傅里叶逆变换与正变换的关系）；
- 附加的相位因子 \(\exp(-j\omega t)\) 与先前时-频平移算子的非交换关系（见 §3.1）一脉相承，本质上来源于 \(\operatorname{T}_x\) 与 \(\operatorname{M}_{\omega}\) 的非交换性。

这一关系表明：**在时域加窗做短时傅里叶变换，与在频域加窗做短时傅里叶变换是等价的**，只是时间参数取了负号并多了一个相位因子。这为从频域视角理解和设计时频分析方法提供了重要依据。
### 3.3 短时傅里叶变换的卷积性质

上一节从时域和频域两个视角定义了短时傅里叶变换，并揭示了两者之间的对偶关系。本节进一步讨论STFT的另一个本质属性——**卷积性质**。这一性质将STFT与经典的线性时不变滤波理论联系起来，为理解STFT作为"滤波器组"的物理意义奠定基础。

从定义出发：

\[
V_g f(t, \omega) = \int_{-\infty}^{\infty} f(t') \, \overline{g(t' - t)} \exp(-j\omega t') \, dt' . \tag{6.19}
\]

这个积分形式上与卷积相似，但窗函数中的变量是 \(t' - t\) 而不是 \(t - t'\)，同时被积函数中还多了一个与频率 \(\omega\) 有关的复指数。目标是将它改写为标准卷积的形式，从而揭示其时域滤波的本质。

---

#### 3.3.1 卷积形式的数学推导

为了将式 (6.19) 转化为标准卷积，引入一个**依赖于频率 \(\omega\) 的核函数** \(h_\omega(\cdot)\)。观察积分中的窗函数部分 \(\overline{g(t' - t)}\)，它实际上是 \(\overline{g(-(t - t'))}\)。因此，定义

\[
h_\omega(\tau) \triangleq \overline{g(-\tau)} \exp(j\omega \tau) ,
\]

那么当 \(\tau = t - t'\) 时，有

\[
h_\omega(t - t') = \overline{g(t' - t)} \exp\bigl(j\omega (t - t')\bigr) .
\]

将上式代入式 (6.19)，可得

\[
\begin{aligned}
V_g f(t, \omega) 
&= \int_{-\infty}^{\infty} f(t') \, \overline{g(t' - t)} \exp(-j\omega t') \, dt' \\
&= \int_{-\infty}^{\infty} f(t') \, \Bigl[ \overline{g(t' - t)} \exp\bigl(j\omega (t - t')\bigr) \Bigr] \exp(-j\omega t) \, dt' \\
&= \exp(-j\omega t) \int_{-\infty}^{\infty} f(t') \, h_\omega(t - t') \, dt' .
\end{aligned}
\]

注意到最后的积分正是 \(f(t)\) 与 \(h_\omega(t)\) 的**卷积**在时刻 \(t\) 的取值，即

\[
(f * h_\omega)(t) = \int_{-\infty}^{\infty} f(t') \, h_\omega(t - t') \, dt' .
\]

因此得到短时傅里叶变换的**卷积表示**：

\[
\boxed{ V_g f(t, \omega) = \exp(-j\omega t) \, (f * h_\omega)(t) , } \tag{6.20}
\]

其中核函数

\[
h_\omega(t) = \overline{g(-t)} \exp(j\omega t) .
\]

![stft](assets/06/03.png)

---

#### 3.3.2 直观例释：时频平面上的"地形图"与滤波通道

直观地理解这个卷积性质：设想一个三维坐标系——水平面由两个轴构成，横轴为时间 \(t\)，纵轴为频率 \(\omega\)，这个水平面就是**时频平面**；竖直方向的坐标轴表示短时傅里叶变换的**幅度** \(|V_g f(t, \omega)|\)。整个三维空间中的曲面如同一幅**时频地形图**——山峰代表信号在某个时刻和某个频率上的能量集中区域，山谷代表能量稀疏区域。

固定某一个频率 \(\omega_0\)，即在地形图上沿着平行于时间轴的方向切一刀，得到一条"**等高线**"或"**剖面线**"：\(|V_g f(t, \omega_0)|\) 随 \(t\) 的变化。卷积性质 (6.20) 表明：**这条剖面线正是原始信号 \(f(t)\) 通过一个以 \(\omega_0\) 为中心的带通滤波器后的输出包络**。这个滤波器的冲激响应是 \(h_{\omega_0}(t) = \overline{g(-t)} \exp(j\omega_0 t)\)，它在频域上相当于把窗谱 \(\tilde{g}\) 平移到中心频率 \(\omega_0\) 处，并取共轭。

用一个更生活化的例子来说：想象你站在一个音乐厅里，手里拿着一个**可调谐的听筒**（就像老式收音机的调频旋钮）。你把听筒调到某个频率 \(\omega_0\)（比如中音C），然后你沿着时间轴行走（从乐曲的开头走到结尾），听筒里听到的响度变化就是 \(|V_g f(t, \omega_0)|\)。这个听筒的"带宽"和"形状"由你选择的窗函数 \(g\) 决定。你把旋钮调到不同的频率，就能得到整个时频平面上的响度分布——这就是一幅完整的地形图。而式 (6.20) 恰恰用数学语言精确描述了"调好频率的听筒"与"原始信号"之间的滤波关系：**先滤波（卷积），再解调（乘以 \(\exp(-j\omega t)\)）**，把高频成分搬回基带以便观察包络。

---

#### 3.3.3 专业的滤波器组理论阐释

从**滤波器组理论**（Filter Bank Theory）的角度来看，式 (6.20) 揭示了STFT的本质：对于每一个固定的频率 \(\omega\)，短时傅里叶变换 \(V_g f(t, \omega)\) 是信号 \(f(t)\) 经过一个**线性时不变（LTI）系统**后的输出，再乘以一个解调因子 \(\exp(-j\omega t)\)。该LTI系统的**冲激响应**为

\[
h_\omega(t) = \overline{g(-t)} \exp(j\omega t) ,
\]

其对应的**频率响应**为

\[
\tilde{h}_\omega(\xi) = \overline{\tilde{g}(\xi - \omega)} .
\]

这意味着，每个 \(\omega\) 都对应一个独立的带通滤波器，其通带中心位于 \(\omega\)，形状由窗函数 \(g\) 的频谱 \(\tilde{g}\) 决定。当 \(\omega\) 遍历整个实数轴时，这些滤波器构成了一个**均匀滤波器组**（Uniform Filter Bank）：所有子带滤波器具有完全相同的形状（都是 \(\tilde{g}\) 的平移），只是中心频率均匀地分布在频率轴上。

整个短时傅里叶变换 \(V_g f(t, \omega)\) 可以理解为：**将信号 \(f(t)\) 同时送入无穷多个中心频率为 \(\omega\) 的带通滤波器，每个滤波器的输出再乘以 \(\exp(-j\omega t)\) 进行解调，得到该频率通道的复包络**。这一视角在语音处理、雷达信号分析和多速率信号处理中具有根本意义——它把时频分析从"积分变换"的抽象数学拉回到了"滤波与解调"的工程直观。

核心结论：

\[
\boxed{ V_g f(t, \omega) = \exp(-j\omega t) \, (f * h_\omega)(t), \quad h_\omega(t) = \overline{g(-t)} \exp(j\omega t) } . \tag{6.20}
\]

该式的**后半部分** —— \((f * h_\omega)(t)\) —— 是在时域上完成的线性卷积（滤波），它决定了STFT的动态范围和时间分辨率；**前半部分** —— \(\exp(-j\omega t)\) —— 则是在时频域上的相位修正（解调），它将滤波后的高频信号搬回基带，使得 \(V_g f(t, \omega)\) 在复平面上具有明确的物理意义（即局部频谱的复振幅）。这两部分的结合，完整地刻画了STFT在时频平面上的卷积本质。


## 4. 测不准原理

### 4.1 时频分辨率的内在制约与能量归一化

从短时傅里叶变换（STFT）的分析中可以看出，时域上窗函数的划分可自由控制，但频域的有效范围却不可控。"时域变胖，频域变瘦"之间是否存在某种内在的、不可逾越的联系？海森堡测不准原理给出了确切的答案。

首先，为保证信号具有物理意义上的有限能量，并使后续的均值与方差定义有意义，假定信号能量有限：
$$
m_0 = \int_{-\infty}^{\infty} |f(t)|^2 dt < \infty \tag{6.21}
$$
> **注**：若信号已归一化，则 $m_0=1$。后续公式中的概率密度函数实际上隐含了除以 $m_0$ 的操作，为书写简洁，以下推导默认信号能量已归一化为1，即 $\int |f(t)|^2 dt = 1$。

### 4.2 时频域均值与方差的物理诠释

#### 4.2.1 时频中心：信号重心的度量

为量化信号在时域和频域的"集中程度"，需定义类似于统计学中期望值的概念。以下两个积分分别代表信号能量分布的时间中心和频率中心：
$$
\boxed{
m_t = \int_{-\infty}^{\infty} t |f(t)|^2 dt, \quad
m_{\omega} = \int_{-\infty}^{\infty} \omega |\hat{f}(\omega)|^2 d\omega
} \tag{6.22}
$$
*   **$m_t$（时间均值）**：表示信号能量在时间轴上的"重心"。若将 $|f(t)|^2$ 视为一个概率密度函数，$m_t$ 即为该随机变量的数学期望，反映信号主要出现的时刻。
*   **$m_{\omega}$（频率均值）**：表示信号能量在频率轴上的"重心"，反映信号的主载波频率或平均振荡速率。

进行时频分析时，关注的往往不是信号的绝对位置，而是信号相对于其中心的"展宽"程度。只有确定了中心 $m_t$ 和 $m_{\omega}$，才能进一步衡量信号偏离中心的离散程度。若不减去均值直接计算二阶矩，结果将包含中心位置的贡献，无法纯粹反映信号的"宽度"。

#### 4.2.2 时频展宽：不确定性的量化

定义了中心之后，用方差来描述信号能量围绕中心的分散程度，这正是"不确定性"或"分辨率"的数学表达：
$$
\boxed{
v_t = \int_{-\infty}^{\infty} (t - m_t)^2 |f(t)|^2 dt, \quad
v_{\omega} = \int_{-\infty}^{\infty} (\omega - m_{\omega})^2 |\hat{f}(\omega)|^2 d\omega
} \tag{6.23}
$$
*   **$v_t$（时间方差）**：衡量信号在时域上的有效持续时间。$v_t$ 越小，信号在时域上越集中，时间定位越精准。
*   **$v_{\omega}$（频率方差）**：衡量信号在频域上的有效带宽。$v_{\omega}$ 越小，频谱越窄，频率分辨力越高。

这两个方差正是测不准原理的主角，定量刻画了STFT中"窗宽"与"频谱宽"的矛盾。目标是探究 $v_t$ 和 $v_{\omega}$ 能否同时任意小。

### 4.3 海森堡测不准原理的结论

借鉴量子力学中海森堡关于粒子位置与动量的测不准关系，在信号处理领域得到完全同构的结论。对于任何有限能量的信号，其时域方差与频域方差的乘积存在一个下界：
$$
\boxed{
v_t v_{\omega} \geq C
} \tag{6.24}
$$
其中 $C$ 是一个正常数（具体数值取决于傅里叶变换的定义形式）。该不等式的物理含义是：**当一个域（如时域）的展宽缩小时，另一个域（如频域）的展宽必然变大。** 无法同时在时域和频域获得任意高的分辨率，这是信号本身的内在属性，而非测量仪器的缺陷。

### 4.4 测不准原理的严格证明

以下通过算子理论严格证明上述不等式。证明过程遵循线性代数与泛函分析的基本框架。

#### 4.4.1 伴随算子与自伴算子

设 $A$ 是定义在希尔伯特空间上的线性算子。其**伴随算子** $A^H$ 定义为满足以下内积关系的算子：
$$
\langle Af, g \rangle = \int_{-\infty}^{\infty} (Af)(t) \overline{g(t)} dt = \langle f, A^H g \rangle \tag{6.25}
$$
若算子满足 $A = A^H$，则称 $A$ 为**自伴算子**（Self-adjoint）。

对于熟悉矩阵运算的工科生而言，这一概念对应于共轭转置。在有限维空间中，内积可写为向量形式：
$$
\langle Af, g \rangle = (Af)^H g = f^H A^H g = \langle f, A^H g \rangle \tag{6.26}
$$
这表明伴随算子在无穷维函数空间中扮演了与矩阵共轭转置完全相同的角色。特别地，时间和频率算子在适当边界条件下均为自伴算子，这保证了它们的观测值（均值、方差）为实数。

#### 4.4.2 算子交换性与共同特征矢量

对于两个算子 $A$ 和 $B$，一个核心问题是它们是否可交换，即 $AB = BA$。
**命题**：若 $A$ 和 $B$ 拥有完全相同的特征矢量集，则 $AB = BA$。

**证明**：
设 $\{\phi_n\}$ 是 $A$ 和 $B$ 的共同完备特征矢量基，对应的特征值分别为 $\lambda_n$ 和 $\mu_n$，即：
$$
A \phi_n = \lambda_n \phi_n, \quad B \phi_n = \mu_n \phi_n
$$
对任意信号 $f$，将其展开为 $f = \sum_n c_n \phi_n$，则有：
$$
AB f = A \left( \sum_n c_n \mu_n \phi_n \right) = \sum_n c_n \mu_n \lambda_n \phi_n
$$
$$
BA f = B \left( \sum_n c_n \lambda_n \phi_n \right) = \sum_n c_n \lambda_n \mu_n \phi_n
$$
由于标量乘法可交换（$\mu_n \lambda_n = \lambda_n \mu_n$），故 $AB f = BA f$ 对任意 $f$ 成立，即 $AB = BA$。反之，若 $AB \neq BA$，则 $A$ 与 $B$ 不可能拥有完全相同的特征矢量集，这意味着它们无法被同时对角化，进而导致联合测量的不确定性。

#### 4.4.3 泊松括号（交换子）

为量化两个算子的非交换程度，定义**泊松括号**（Poisson Bracket，在量子力学中常称为交换子 Commutator）：
$$
[A, B] = AB - BA \tag{6.27}
$$
若 $[A, B] = 0$，则两算子可交换；若 $[A, B] \neq 0$，则存在内在的不确定性耦合。

#### 4.4.4 交换子期望值的虚部性质

计算交换子在状态 $f$ 下的期望值 $\langle [A, B]f, f \rangle$。假设 $A$ 和 $B$ 均为自伴算子（$A=A^H, B=B^H$），有如下关键恒等式：
$$
\boxed{
\langle [A, B]f, f \rangle = 2j \, \text{Im}\langle Bf, Af \rangle
} \tag{6.28}
$$
**详细推导**：
展开交换子：
$$
\langle [A, B]f, f \rangle = \langle ABf, f \rangle - \langle BAf, f \rangle
$$
利用伴随算子定义及 $A, B$ 的自伴性：
$$
\langle ABf, f \rangle = \langle Bf, A^H f \rangle = \langle Bf, Af \rangle
$$
$$
\langle BAf, f \rangle = \langle Af, B^H f \rangle = \langle Af, Bf \rangle = \overline{\langle Bf, Af \rangle}
$$
令复数 $z = \langle Bf, Af \rangle$，则：
$$
\langle [A, B]f, f \rangle = z - \bar{z} = 2j \, \text{Im}(z) = 2j \, \text{Im}\langle Bf, Af \rangle
$$
> **注**：该步骤揭示了交换子的期望值是纯虚数，且正比于 $\langle Bf, Af \rangle$ 的虚部。

#### 4.4.5 柯西-施瓦茨不等式的应用

对式 (6.28) 两边取模，并利用柯西-施瓦茨（Cauchy-Schwarz）不等式 $|\langle u, v \rangle| \leq \|u\| \|v\|$：
$$
|\langle [A, B]f, f \rangle| = |2j \, \text{Im}\langle Bf, Af \rangle| = 2 |\text{Im}\langle Bf, Af \rangle|
$$
由于复数的虚部绝对值不超过其模长，即 $|\text{Im}(z)| \leq |z|$，结合柯西-施瓦茨不等式得：
$$
\boxed{
|\langle [A, B]f, f \rangle| \leq 2 \|Bf\| \|Af\|
} \tag{6.29}
$$
这是测不准原理最一般的算子形式：两个算子交换子的期望值大小，受限于这两个算子作用于信号后的范数之积。

#### 4.4.6 代入时频算子导出测不准界

现在选取具体的时频算子。为使推导简洁，先假设信号已中心化，即 $m_t = 0, m_{\omega} = 0$（否则可用 $f(t)\exp(-j m_{\omega} t)$ 和时移预处理实现，不影响方差结果）。定义：
$$
(Af)(t) = t f(t), \quad (Bf)(t) = \frac{1}{2\pi j} \frac{d}{dt} f(t) \tag{6.30}
$$
> **说明**：此处 $B$ 的选择对应于以普通频率 $f$（Hz）为变量的傅里叶变换对 $\hat{f}(\xi) = \int f(t) \exp(-2\pi j \xi t) dt$。若使用角频率 $\omega$，则 $B = \frac{1}{j} \frac{d}{dt}$，常数 $C$ 会相应变化。以下推导基于式 (6.30) 的定义。

**第一步：计算交换子 $[A, B]$**
$$
[A, B]f = A(Bf) - B(Af) = t \cdot \frac{1}{2\pi j} f'(t) - \frac{1}{2\pi j} \frac{d}{dt}[t f(t)]
$$
$$
= \frac{t}{2\pi j} f'(t) - \frac{1}{2\pi j} [f(t) + t f'(t)] = -\frac{1}{2\pi j} f(t)
$$
即：
$$
[A, B] = -\frac{1}{2\pi j} I \tag{6.31}
$$
其中 $I$ 为恒等算子。

**第二步：计算交换子期望值**
$$
\langle [A, B]f, f \rangle = \left\langle -\frac{1}{2\pi j} f, f \right\rangle = -\frac{1}{2\pi j} \langle f, f \rangle = \frac{j}{2\pi} \|f\|^2 \tag{6.32}
$$
取模得：
$$
|\langle [A, B]f, f \rangle| = \frac{1}{2\pi} \|f\|^2 \tag{6.33}
$$

**第三步：关联方差与范数**
在 $m_t = 0, m_{\omega} = 0$ 的假设下：
$$
\|Af\|^2 = \int |t f(t)|^2 dt = v_t \tag{6.34}
$$
根据帕塞瓦尔定理及傅里叶变换的微分性质 $\mathcal{F}\{f'(t)\} = 2\pi j \xi \hat{f}(\xi)$：
$$
\|Bf\|^2 = \int \left| \frac{1}{2\pi j} f'(t) \right|^2 dt = \int |\xi \hat{f}(\xi)|^2 d\xi = v_{\omega} \tag{6.35}
$$
> **注**：此处 $v_{\omega}$ 是以 Hz 为单位的频率方差。若使用角频率 $\omega = 2\pi \xi$，则 $v_{\omega_{\text{rad}}} = (2\pi)^2 v_{\xi}$。

**第四步：得到最终界限**
将式 (6.33)、(6.34)、(6.35) 代入一般不等式 (6.29)：
$$
\frac{1}{2\pi} \|f\|^2 \leq 2 \sqrt{v_t} \sqrt{v_{\omega}}
$$
整理得：
$$
\boxed{
v_t v_{\omega} \geq \frac{1}{16\pi^2} \|f\|^4
} \tag{6.36}
$$
当信号能量归一化 $\|f\|^2 = 1$ 时，即为著名的海森堡测不准原理在信号处理中的标准形式：
$$
\boxed{
v_t v_{\omega} \geq \frac{1}{16\pi^2}
} \tag{6.37}
$$
若采用角频率 $\omega$ 定义傅里叶变换，对应的下界为 $v_t v_{\omega} \geq \frac{1}{4}$。无论哪种定义，核心结论不变：**时频方差的乘积存在一个严格的正下界**。等号成立当且仅当 $f(t)$ 为高斯函数 $\exp(-\alpha t^2)$，这也解释了为何高斯窗在STFT中具有最优的时频聚集性。
### 4.5 时频算子选取的物理与数学必然性

在式 (6.30) 中，特意选取了 $(Af)(t) = t f(t)$ 与 $(Bf)(t) = \frac{1}{2\pi j} \frac{d}{dt} f(t)$ 这一对算子代入通用不等式。这并非随意的数学凑配，而是由方差的定义、傅里叶变换的对偶性质以及量子力学与信号处理的同构性共同决定的。以下从三个层面解释这一选取的必然性。

#### 4.5.1 方差定义的直接对应（物理动机）

回顾式 (6.23) 中定义的方差：
$$
v_t = \int_{-\infty}^{\infty} (t - m_t)^2 |f(t)|^2 dt, \quad v_{\omega} = \int_{-\infty}^{\infty} (\omega - m_{\omega})^2 |\hat{f}(\omega)|^2 d\omega
$$
为了利用算子不等式 $|\langle [A, B]f, f \rangle| \leq 2 \|Af\| \|Bf\|$ 来约束 $v_t v_{\omega}$，必须找到两个算子，使得它们的范数平方 $\|Af\|^2$ 和 $\|Bf\|^2$ 恰好等于（或正比于）$v_t$ 和 $v_{\omega}$。

-   **对于 $v_t$**：被积函数包含因子 $t$（假设已中心化 $m_t=0$）。只有乘法算子 $(Af)(t) = t f(t)$ 能满足：
    $$
    \|Af\|^2 = \int |t f(t)|^2 dt = v_t
    $$
    任何其他形式的算子都无法直接生成 $t^2 |f(t)|^2$ 这一项。因此，$A$ 必须是时间乘法算子。

-   **对于 $v_{\omega}$**：被积函数包含因子 $\omega$（或 $\xi$），但积分是在频域进行的。需要一个作用在时域上的算子 $B$，使其能量等价于频域的加权能量。根据帕塞瓦尔定理和傅里叶变换的微分性质，频域的乘法对应时域的微分。具体地，若 $\mathcal{F}\{f(t)\} = \hat{f}(\xi)$，则：
    $$
    \mathcal{F}\left\{ \frac{1}{2\pi j} \frac{d}{dt} f(t) \right\} = \xi \hat{f}(\xi)
    $$
    因此，选取 $(Bf)(t) = \frac{1}{2\pi j} \frac{d}{dt} f(t)$ 可保证：
    $$
    \|Bf\|^2 = \int |\xi \hat{f}(\xi)|^2 d\xi = v_{\omega}
    $$
    这正是将频域方差"拉回"时域进行计算的唯一自然选择。

> **小结**：$A$ 和 $B$ 的选取是由 $v_t$ 和 $v_{\omega}$ 的数学结构唯一确定的。目标是界定这两个特定量的乘积，因此必须选用能精确生成这两个量的算子。

#### 4.5.2 傅里叶变换的对偶性与非交换根源（数学动机）

为什么偏偏是"乘法"和"微分"这一对？因为它们是傅里叶变换下的**对偶算子**。

傅里叶变换本质上是一个基变换，它将时域的平移不变性转换为频域的相位调制，将时域的尺度缩放转换为频域的反向缩放。在这种对偶关系下：
-   时域的**位置算子**（乘法）$\leftrightarrow$ 频域的**微分算子**
-   时域的**微分算子** $\leftrightarrow$ 频域的**位置算子**（乘法）

这种对偶性导致了它们天然不可交换，如式 (6.31) 所示：
$$
[A, B] = t \cdot \frac{1}{2\pi j} \frac{d}{dt} - \frac{1}{2\pi j} \frac{d}{dt} \cdot t = -\frac{1}{2\pi j} I
$$
这个非零常数结果（$-\frac{1}{2\pi j}$）正是测不准原理下界 $C$ 的来源。**如果 $A$ 和 $B$ 可交换，下界为零，测不准原理就不存在。** 而乘法与微分的非交换性，在数学上等价于"函数的自变量与其变化率不能同时被精确指定"，这正是时频分辨率矛盾的本质。

此外，只有当 $[A, B]$ 为标量算子（即常数乘以恒等算子）时，不等式右边才是一个与信号 $f$ 无关的普适常数。若选取其他算子组合，交换子可能依赖于 $f$ 本身，无法得到形如 $v_t v_{\omega} \geq C$ 的绝对下界。乘法与微分是唯一一对满足 $[A, B] = cI$（$c$ 为非零常数）的自伴算子对（在石川-冯·诺依曼唯一性定理的意义下）。

#### 4.5.3 与量子力学的同构映射（理论动机）

海森堡测不准原理最初诞生于量子力学，其表述为位置 $x$ 与动量 $p$ 的不确定性关系：
$$
\Delta x \Delta p \geq \frac{\hbar}{2}
$$
在量子力学中，位置算子为 $\hat{x} \psi(x) = x \psi(x)$，动量算子为 $\hat{p} \psi(x) = -j\hbar \frac{d}{dx} \psi(x)$。对比可知：

| 量子力学 | 信号处理 | 对应关系 |
| :--- | :--- | :--- |
| 波函数 $\psi(x)$ | 信号 $f(t)$ | 状态描述 |
| 位置算子 $\hat{x}$ | 时间算子 $A = t$ | 自变量 |
| 动量算子 $\hat{p} = -j\hbar \partial_x$ | 频率算子 $B = \frac{1}{2\pi j} \partial_t$ | 共轭变量 |
| 普朗克常数 $\hbar$ | $1/(2\pi)$ | 对偶尺度因子 |

信号处理中的"时间-频率"对，与量子力学中的"位置-动量"对在数学结构上完全同构。频率之所以对应微分，是因为在德布罗意关系中 $p = h/\lambda = h \xi$，动量正比于空间频率；而在信号中，瞬时频率正是相位的导数。

因此，选取 $A=t$ 和 $B=\frac{1}{2\pi j}\frac{d}{dt}$ 不仅是数学推导的需要，更是因为**时间和频率本身就是一对共轭物理量**，它们之间的关系由傅里叶变换这一数学结构所固化。任何试图绕过这对算子来证明时频测不准原理的尝试，都将失去物理根基。

> **核心结论**：选取这两个算子，是因为它们是唯一能够同时将"时域方差"和"频域方差"表达为算子范数、且其交换子为非零常数的自伴算子对。这一选取既由方差的定义所强制，也由傅里叶对偶性所保证，更由物理世界的共轭变量本质所决定。



## 5. 课后总结

### 5.1 核心逻辑链：从全局谱到时频谱

本讲从阵列信号处理的窄带假设出发，指出了经典傅里叶变换的局限性——有频率分辨率无时间分辨率，由此引入时频分析的根本动机。核心逻辑链条如下：

1. **窄带假设的失效**：阵列信号处理要求信号带宽远小于载频（$B \ll f_c$），才能将时延转化为相移。实际信号（语音、音乐、雷达回波）的频率随时间变化，窄带假设不再成立，需要同时刻画信号的时间演化和频率构成。

2. **时频分析的本质**：信号在时间-频率平面上的能量分布。以时间为横轴、频率为纵轴、颜色表示能量，形成时频谱（spectrogram）。

3. **短时傅里叶变换（STFT）**：对信号加窗，逐段做傅里叶分析。窗越短，时间分辨率越好，频率分辨率越差；窗越长，频率分辨率越好，时间分辨率越差。这是测不准原理的直接体现。

4. **STFT的两种等价定义**：时域加窗 $V_g f(t,\omega)$ 与频域加窗 $V_{\tilde{g}} \tilde{f}(\omega, t)$ 通过对偶性等价，相差一个相位因子和参数符号翻转。

5. **STFT的卷积性质**：$V_g f(t,\omega) = e^{-j\omega t}(f * h_\omega)(t)$，即STFT等价于一组中心频率均匀分布的带通滤波器组对信号滤波后再解调。这是时频分析联系工程实现（滤波器组）的关键桥梁。

6. **测不准原理**：$v_t \cdot v_\omega \geq 1/4$（角频率定义下），时宽和频宽的乘积有严格正下界。等号成立当且仅当信号是高斯函数。这是时频分析的根本限制，与算法无关。

### 5.2 关键公式总结

| 公式 | 含义 |
| :--- | :--- |
| $V_g f(t,\omega) = \int f(t') \overline{g(t'-t)} e^{-j\omega t'} dt'$ | STFT定义（时域窗） |
| $V_{\tilde{g}} \tilde{f}(\omega, t) = e^{-j\omega t} V_g f(-t, \omega)$ | 时频对偶关系 |
| $V_g f(t,\omega) = e^{-j\omega t} (f * h_\omega)(t)$ | STFT的卷积/滤波器组解释 |
| $v_t v_\omega \geq 1/4$ | 海森堡测不准原理（角频率） |
| $\Delta t \cdot \Delta f \geq 1/(4\pi)$ | 测不准原理（普通频率） |

### 5.3 重点概念

- **窗函数**：决定STFT的时间和频率分辨率。矩形窗频率泄漏严重；Hamming窗旁瓣低；高斯窗达到测不准下界。
- **时频分辨率权衡**：不存在最优窗长，只有针对具体问题的最优选择。
- **滤波器组解释**：每个频率 $\omega$ 对应一个带通滤波器，STFT就是并行滤波器组的输出。
- **测不准原理的物理意义**：信号不可能同时在时间和频率上任意集中，这是傅里叶变换的数学性质决定的根本限制。

---

## 6. 学习检查清单：自测核心知识点掌握情况

- [ ] 能解释为什么阵列信号处理需要窄带假设，以及窄带假设失效后怎么办
- [ ] 能写出STFT的定义式，并解释各符号的含义
- [ ] 能推导STFT的时频对偶关系 $V_{\tilde{g}} \tilde{f}(\omega, t) = e^{-j\omega t} V_g f(-t, \omega)$
- [ ] 能推导STFT的卷积性质，并解释"滤波器组"的物理意义
- [ ] 能陈述海森堡测不准原理，并解释其物理含义
- [ ] 能证明测不准原理（从算子交换子出发的完整证明）
- [ ] 能解释为什么高斯窗在STFT中具有最优时频聚集性
- [ ] 能区分时间分辨率与频率分辨率，并说明它们之间的权衡关系
- [ ] 能解释泊松求和公式，并说明它在采样定理中的作用
- [ ] 能写出Parseval定理，并说明其"保内积"的含义

---

## 7. 思考题：拓展与挑战

1. **窗函数的选择**：矩形窗、Hamming窗、Hanning窗、高斯窗，它们的时域和频域形状各有什么特点？的主瓣宽度和旁瓣水平如何权衡？如果用矩形窗做STFT，会出现什么问题？

2. **测不准原理的直观理解**：用一个正弦波乘以一个矩形窗，窗越长，频率估计越精确，但时间定位越模糊。用一个冲击信号（持续时间极短），时间定位精确，但频谱很宽。请用这个例子验证测不准原理 $v_t v_\omega \geq 1/4$。

3. **STFT vs 小波变换**：STFT使用固定长度的窗，小波变换使用长度随频率变化的窗（高频用短窗，低频用长窗）。为什么说小波变换在时频分析中具有"自适应"的分辨率？它是否突破了测不准原理？

4. **滤波器组解释的代码验证**：给定一个信号 $f(t) = \sin(2\pi f_1 t) + \sin(2\pi f_2 t)$，用STFT和用一组带通滤波器分别处理，比较结果是否一致。如果不一致，原因是什么？

5. **ESP32阵列的时频应用**：你有20多个ESP32组成的阵列。如果每个ESP32同时采集一段音频信号，你可以用STFT分析每个阵元接收到的信号的时频谱。请问：如何利用20个阵元的时频分析结果，估计声源的方向？（提示：结合MUSIC算法的思想，但用的是时频域的相位差而非频域的相位差。）

6. **2×2 MIMO SDR的时频实验**：SDR可以发射线性调频信号（chirp signal）$s(t) = \exp(j\pi \alpha t^2)$，其频率随时间线性变化。用STFT分析这个信号，时频谱应该是一条斜线。请推导这条斜线的斜率与调频参数 $\alpha$ 的关系。

---

## 8. 实验设计：基于ESP32阵列与2×2 MIMO SDR的时频分析实验

### 8.0 实验总览

本实验使用的硬件：
- **20+ ESP32组成的阵列**：每个ESP32可以通过I2S接口连接麦克风，实现同步采样，获取多通道音频信号。
- **2×2 MIMO SDR**：可以发射和接收任意波形的射频信号，用于验证STFT在雷达/chirp信号处理中的应用。

实验设计遵循"从软件仿真到硬件实现，从单通道到多通道"的递进逻辑。

| 实验编号 | 实验名称 | 使用硬件 | 预期结论 |
| :--- | :--- | :--- | :--- |
| 8.1 | STFT软件仿真：窗函数对比 | 纯软件（Python/MATLAB） | 不同窗函数的时频分辨率权衡 |
| 8.2 | ESP32阵列：多通道音频信号的STFT分析 | 20+ ESP32阵列 | 多通道时频分析，声源方向估计 |
| 8.3 | SDR chirp信号生成与STFT分析 | 2×2 MIMO SDR | 验证chirp信号的时频特性 |
| 8.4 | 测不准原理的硬件验证 | ESP32 + SDR | 定量验证 $v_t v_\omega \geq 1/4$ |

---

### 8.1 实验1：STFT软件仿真——窗函数对比与测不准原理可视化

**目的**：在软件中生成不同类型的信号，用不同窗函数做STFT，观察时频分辨率的权衡，验证测不准原理。

**步骤**：
1. 生成测试信号：
   - 信号1：两个正弦波叠加 $s(t) = \sin(2\pi f_1 t) + \sin(2\pi f_2 t)$，$f_1=100\text{ Hz}, f_2=120\text{ Hz}$
   - 信号2：线性调频信号 $s(t) = \exp(j\pi \alpha t^2)$，$\alpha=50\text{ Hz/s}$
   - 信号3：冲击信号（极短的高斯脉冲）
2. 分别用矩形窗、Hamming窗、高斯窗做STFT，窗长分别取 $N=64, 256, 1024$
3. 绘制时频谱图（时间×频率×能量），观察时间分辨率和频率分辨率的变化
4. 计算各信号的 $v_t$ 和 $v_\omega$，验证 $v_t v_\omega \geq 1/4$

**定量指标**：
- 两个正弦波的频率间隔 $\Delta f = 20\text{ Hz}$，多大窗长能区分这两个频率？
- 冲击信号的时宽 $v_t$ 和频宽 $v_\omega$，计算乘积并对比理论下界

**预期结果**：
- 窗长 $N$ 越大，频率分辨率越好（能区分更接近的频率），但时间分辨率越差
- 高斯窗的 $v_t v_\omega$ 最接近理论下界
- 矩形窗频率泄漏严重，时频谱中出现虚假频率成分

---

### 8.2 实验2：ESP32阵列多通道音频信号的STFT分析

**目的**：利用ESP32阵列采集多通道音频，用STFT分析各通道的时频特性，并尝试估计声源方向。

**硬件连接**：
- 20个ESP32，每个连接一个I2S麦克风（如INMP441）
- 通过WiFi或蓝牙将采样数据发送到主机（或其中一个ESP32作为主机）
- 所有ESP32共享同一个时钟信号（通过WiFi时间同步或硬件触发），实现同步采样

**步骤**：
1. 在阵列前方放置一个扬声器，播放一段音乐（频率随时间变化）
2. 所有ESP32同时采集音频信号，采样率 $f_s = 16\text{ kHz}$，采集时长 $T = 2\text{ s}$
3. 对每个通道的数据做STFT，得到20个时频谱 $V_g f_i(t,\omega), i=1,\dots,20$
4. 在每一个时间-频率点 $(t,\omega)$，比较20个通道的相位差 $\Delta \phi_i(t,\omega)$
5. 利用相位差估计声源的方向（类似MUSIC算法，但是时频域版本）

**定量指标**：
- 声源方向估计误差（与真实方向对比）
- 不同频率成分的方向估计精度（高频 vs 低频）

**预期结果**：
- 时频域的相位差可以用于估计声源方向
- 高频信号的方向估计精度高于低频信号（波长更短，相位差更明显）
- 如果声源是移动的，时频谱可以捕捉到频率和多普勒效应的变化

---

### 8.3 实验3：SDR chirp信号生成与STFT分析

**目的**：用SDR发射线性调频信号，用另一个SDR接收并做STFT分析，验证chirp信号的时频特性。

**硬件连接**：
- 一台SDR作为发射机（TX），一台SDR作为接收机（RX）
- 两台SDR通过电缆连接（加衰减器），避免电磁干扰

**步骤**：
1. 发射机生成线性调频信号 $s(t) = \exp(j\pi \alpha t^2)$，调频斜率 $\alpha = 1\text{ MHz/ms}$，中心频率 $f_c = 2.4\text{ GHz}$
2. 发射机通过SDR将信号上变频到 $f_c$ 并发射
3. 接收机接收信号，下变频到基带，得到复基带信号 $r(t)$
4. 对 $r(t)$ 做STFT，绘制时频谱
5. 从时频谱中提取频率随时间的变化曲线，与理论值 $\alpha t$ 对比

**定量指标**：
- 时频谱中频率斜线的斜率，与 $\alpha$ 的误差
- 不同窗长对chirp信号时频分析精度的影响

**预期结果**：
- 时频谱应显示一条清晰的斜线，斜率等于 $\alpha$
- 窗长越长，斜线越平滑，但瞬时频率的变化越模糊
- 如果发射机和接收机之间有相对运动，时频谱中会出现多普勒频移

---

### 8.4 实验4：测不准原理的硬件验证

**目的**：用硬件生成具有不同时宽的信号，通过STFT测量 $v_t$ 和 $v_\omega$，验证测不准原理。

**步骤**：
1. 用SDR生成三种信号：
   - 长脉冲：高斯包络，时宽 $v_t$ 大
   - 短脉冲：高斯包络，时宽 $v_t$ 小
   - 连续波：近似无限长，时宽 $v_t \to \infty$
2. 对每个信号做STFT，从时频谱中估计 $v_t$（时间方差）和 $v_\omega$（频率方差）
3. 计算 $v_t v_\omega$，与理论下界 $1/4$ 对比

**定量指标**：
- 三种信号的 $v_t, v_\omega, v_t v_\omega$
- $v_t v_\omega$ 是否始终大于等于 $1/4$

**预期结果**：
- 长脉冲：$v_t$ 大，$v_\omega$ 小，乘积接近 $1/4$
- 短脉冲：$v_t$ 小，$v_\omega$ 大，乘积接近 $1/4$
- 连续波：$v_t \to \infty$，$v_\omega \to 0$，乘积趋于 $1/4$（理论值）
- 高斯包络信号的 $v_t v_\omega$ 最接近理论下界

---

### 8.5 实验报告要求

1. 给出所有实验的时频谱图（时间×频率×能量），标注时间分辨率和频率分辨率
2. 对比不同窗函数的时频分析结果，分析其优缺点
3. 给出ESP32阵列多通道STFT的分析结果，包括各通道的时频谱和相位差
4. 给出SDR chirp信号的时频分析结果，提取频率斜线并验证调频斜率
5. 计算各信号的 $v_t v_\omega$，验证测不准原理
6. 讨论：基于你的硬件，STFT在实际应用中的主要限制是什么？（采样率、计算量、同步精度等）

---

<div style="page-break-before: always;"></div>
