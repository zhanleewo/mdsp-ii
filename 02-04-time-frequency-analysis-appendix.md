
<div style="page-break-before: always; padding: 8% 8% 0 8%;">
 <h1 id="附录1-自适应时频分析" style="text-align: center; margin-bottom: 2rem; border-bottom: none; display: block;">附录 1：自适应时频分析</h1> 
 <div style="display: flex; align-items: center; justify-content: center; gap: 12px; margin: 1.8rem auto;">
  <span style="flex:1; max-width:80px; height:1px; background: linear-gradient(to right, transparent, #888);"></span>
  <span style="display:inline-block; width:6px; height:6px; background:#38bdf8; border-radius:50%;"></span>
  <span style="flex:1; max-width:80px; height:1px; background: linear-gradient(to left, transparent, #888);"></span>
 </div>
</div>

<!-- # 附录：自适应时频分析 -->

## 1. 引言

前面三讲系统建立了时频分析的完整理论框架。然而，无论是STFT还是小波变换，其分析基函数**一旦选定便固定不变**，在整个时频平面上使用完全相同的"模板"去匹配所有位置的信号特征。这种"一刀切"的方式存在根本性局限：信号的非平稳性在时频平面上可能剧烈变化——某些区域需要高时间分辨率（瞬态脉冲），另一些区域需要高频率分辨率（稳态谐波），固定基函数无法适应这种局部变化。

自适应时频分析的核心思想是：**让分析基函数的参数根据信号的局部特征自适应地调整**，从而在时频平面的每个位置上获得最优的局部表示。本附录系统推导自适应时频分析的核心方法。

---

## 2. 自适应时频分析

自适应时频分析的关键在于引入**可调节的基函数参数**，并建立**优化准则**来确定每个时频位置上的最优参数。记可调节参数为 $\theta \in \Theta$（$\Theta$ 为参数空间）。

### 2.1 自适应连续小波变换

#### 2.1.1 参数化小波族

回顾CWT的标准形式。给定母小波 $\phi(t)$，基函数族为：

\[
\phi_{a,b}(t) = \frac{1}{\sqrt{a}} \, \phi\!\left(\frac{t - b}{a}\right)
\tag{A1.1}
\]

对应的小波变换系数：

\[
W_f(a, b) = \int_{-\infty}^{\infty} f(t) \, \overline{\phi_{a,b}(t)} \, dt
\tag{A1.2}
\]

在标准CWT中，母小波 $\phi(t)$ 是固定的。自适应CWT引入**参数化母小波族** $\{\phi_\theta(t) : \theta \in \Theta\}$，其中 $\theta$ 控制小波的形状特征。

最典型的参数化族是在Morlet小波基础上引入带宽参数 $\sigma$：

\[
\phi_\sigma(t) = \frac{1}{\sqrt[4]{\pi \sigma^2}} \exp\!\left(-\frac{t^2}{2\sigma^2}\right) \exp(j 2\pi \eta t)
\tag{A1.3}
\]

其中 $\sigma > 0$ 为高斯包络的宽度，$\eta > 0$ 为中心频率。更一般地，参数向量 $\theta = (\theta_1, \ldots, \theta_K)^T$ 描述小波的完整形状。

对于每个 $\theta$，母小波 $\phi_\theta(t)$ 必须满足可容许性条件：

\[
C_{\phi_\theta} = \int_{-\infty}^{\infty} \frac{|\hat{\phi}_\theta(\omega)|^2}{|\omega|} \, d\omega < \infty
\tag{A1.4}
\]

参数化基函数的尺度-平移生成规则与标准CWT一致：

\[
\phi_{\theta, a, b}(t) = \frac{1}{\sqrt{a}} \, \phi_\theta\!\left(\frac{t - b}{a}\right)
\tag{A1.5}
\]

#### 2.1.2 自适应CWT定义

对每个尺度-平移坐标 $(a, b)$ 独立选择最优参数 $\theta^\star(a, b)$。

**定义1（自适应CWT系数）**：

\[
\boxed{
W_f^\star(a, b) = \int_{-\infty}^{\infty} f(t) \, \overline{\phi_{\theta^\star(a,b), a, b}(t)} \, dt
}
\tag{A1.6}
\]

其中最优参数由局部幅度最大化准则确定：

\[
\theta^\star(a, b) = \arg\max_{\theta \in \Theta} \; |W_f^\theta(a, b)|
\tag{A1.7}
\]

$W_f^\theta(a, b)$ 为使用参数 $\theta$ 的小波变换系数：

\[
W_f^\theta(a, b) = \int_{-\infty}^{\infty} f(t) \, \overline{\phi_{\theta, a, b}(t)} \, dt
\tag{A1.8}
\]

上述准则的直观含义：在尺度 $a$ 和平移 $b$ 处，遍历所有可能的 $\theta$，选择使内积绝对值最大的参数——即选择与信号局部片段"最相似"的小波形状。

**全局优化准则**：对于整个时频平面统一选择参数 $\theta^\star$，可采用Renyi熵最小化。首先计算归一化分布：

\[
P_f^\theta(a, b) = \frac{|W_f^\theta(a, b)|^2}{\iint |W_f^\theta(a, b)|^2 \, \frac{da \, db}{a^2}}
\tag{A1.9}
\]

$\alpha$ 阶Renyi熵定义为：

\[
H_\alpha(W_f^\theta) = \frac{1}{1 - \alpha} \log_2 \iint \left[P_f^\theta(a, b)\right]^\alpha \, \frac{da \, db}{a^2}, \quad \alpha > 0, \alpha \neq 1
\tag{A1.10}
\]

Renyi熵越小，时频能量越集中，表示质量越好。全局最优参数为：

\[
\theta^\star = \arg\min_{\theta \in \Theta} \; H_\alpha(W_f^\theta)
\tag{A1.11}
\]

#### 2.1.3 解析信号的还原

对于解析信号 $f_a(t)$（$\hat{f}_a(\omega)$ 仅在 $\omega \ge 0$ 处非零），标准CWT的重建公式为：

\[
f_a(t) = \frac{1}{C_\phi} \int_{0}^{\infty} \int_{-\infty}^{\infty} W_f(a, b) \, \phi_{a,b}(t) \, \frac{da \, db}{a^2}
\tag{A1.12}
\]

**逐点自适应重建**：由于每个 $(a, b)$ 处使用的参数 $\theta^\star(a, b)$ 不同，重建公式为：

\[
\boxed{
f_a(t) = \int_{0}^{\infty} \int_{-\infty}^{\infty} \frac{1}{C_{\phi_{\theta^\star(a,b)}}} \,
W_f^\star(a, b) \, \phi_{\theta^\star(a,b), a, b}(t) \, \frac{da \, db}{a^2}
}
\tag{A1.13}
\]

其中 $C_{\phi_\theta}$ 是参数为 $\theta$ 的母小波的可容许性常数：

\[
C_{\phi_\theta} = \int_{0}^{\infty} \frac{|\hat{\phi}_\theta(\omega)|^2}{\omega} \, d\omega
\tag{A1.14}
\]

**推导思路**：对任意固定的 $\theta$，可容许性条件保证了：

\[
f_a(t) = \frac{1}{C_{\phi_\theta}} \int_{0}^{\infty} \int_{-\infty}^{\infty} \langle f_a, \phi_{\theta, a, b} \rangle \, \phi_{\theta, a, b}(t) \, \frac{da \, db}{a^2}
\tag{A1.15}
\]

在逐点自适应下，对每个 $(a, b)$ 选取使 $|\langle f_a, \phi_{\theta, a, b} \rangle|$ 最大的 $\theta$，并用对应的 $C_{\phi_\theta}$ 归一化。由于 $C_{\phi_\theta}$ 依赖于 $\theta$ 而不能提取到积分号外，(9.13) 中的 $C_{\phi_{\theta^\star(a,b)}}$ 被保留在积分内部。

**全局统一参数重建**：若采用全局优化确定统一参数 $\theta^\star$，则：

\[
\boxed{
f_a(t) = \frac{1}{C_{\phi_{\theta^\star}}} \int_{0}^{\infty} \int_{-\infty}^{\infty}
W_f^{\theta^\star}(a, b) \, \phi_{\theta^\star, a, b}(t) \, \frac{da \, db}{a^2}
}
\tag{A1.16}
\]

#### 2.1.4 实信号的还原

对于实信号 $f(t)$，利用频谱的Hermitian对称性 $\hat{f}(-\omega) = \overline{\hat{f}(\omega)}$，标准CWT的重建为：

\[
f(t) = \frac{2}{C_\phi} \, \Re \left\{ \int_{0}^{\infty} \int_{-\infty}^{\infty} W_f(a, b) \, \phi_{a,b}(t) \, \frac{da \, db}{a^2} \right\}
\tag{A1.17}
\]

自适应CWT对实信号的还原公式为：

\[
\boxed{
f(t) = 2 \, \Re \left\{ \int_{0}^{\infty} \int_{-\infty}^{\infty} \frac{1}{C_{\phi_{\theta^\star(a,b)}}} \,
W_f^\star(a, b) \, \phi_{\theta^\star(a,b), a, b}(t) \, \frac{da \, db}{a^2} \right\}
}
\tag{A1.18}
\]

取实部是因为实信号的负频率部分是正频率部分的共轭，仅需计算正频率部分的积分并取实部的两倍。

---

### 2.2 自适应短时傅里叶变换

#### 2.2.1 参数化窗函数族

标准STFT以 $g(t)$ 为窗函数：

\[
V_g f(t, \omega) = \int_{-\infty}^{\infty} f(\tau) \, \overline{g(\tau - t)} \, \exp(-j\omega \tau) \, d\tau
\tag{A1.19}
\]

自适应STFT引入**参数化窗函数族** $\{g_\theta(t) : \theta \in \Theta\}$。

**高斯窗族**是最常见的构造，引入宽度参数 $\sigma$：

\[
g_\sigma(t) = \frac{1}{\sqrt[4]{2\pi \sigma^2}} \exp\!\left(-\frac{t^2}{4\sigma^2}\right)
\tag{A1.20}
\]

归一化因子 $1/\sqrt[4]{2\pi\sigma^2}$ 保证 $\|g_\sigma\| = 1$。对应的频域窗为：

\[
\tilde{g}_\sigma(\omega) = \sqrt[4]{8\pi \sigma^2} \, \exp(-\sigma^2 \omega^2)
\tag{A1.21}
\]

（推导：对 (9.20) 做傅里叶变换，利用 $\int \exp(-at^2 - j\omega t) dt = \sqrt{\pi/a} \, \exp(-\omega^2/(4a))$，代入 $a = 1/(4\sigma^2)$ 即得。）

高斯窗族满足测不准原理的下界：

\[
\Delta t_g \cdot \Delta \omega_g = \frac{1}{2}
\tag{A1.22}
\]

参数 $\sigma$ 的优化选择相当于在时间分辨率与频率分辨率之间寻找最优折中。

#### 2.2.2 自适应STFT定义

**定义2（自适应STFT系数）**：

\[
\boxed{
V_f^\star(t, \omega) = \int_{-\infty}^{\infty} f(\tau) \, \overline{g_{\theta^\star(t,\omega)}(\tau - t)} \, \exp(-j\omega \tau) \, d\tau
}
\tag{A1.23}
\]

其中最优参数 $\theta^\star(t, \omega)$ 由局部幅度最大化确定：

\[
\theta^\star(t, \omega) = \arg\max_{\theta \in \Theta} \; |V_f^\theta(t, \omega)|
\tag{A1.24}
\]

其中 $V_f^\theta(t, \omega)$ 是使用窗函数 $g_\theta$ 的STFT系数：

\[
V_f^\theta(t, \omega) = \int_{-\infty}^{\infty} f(\tau) \, \overline{g_\theta(\tau - t)} \, \exp(-j\omega \tau) \, d\tau
\tag{A1.25}
\]

#### 2.2.3 解析信号的还原

对于解析信号 $f_a(t)$，标准STFT的逆变换（由Moyal公式导出）为：

\[
f_a(t) = \frac{1}{2\pi \, \|g\|^2} \int_{-\infty}^{\infty} \int_{-\infty}^{\infty}
V_g f(t, \omega) \, g(\tau - t) \, \exp(j\omega t) \, d\omega \, dt
\tag{A1.26}
\]

对于自适应STFT，逐点自适应重建公式为：

\[
\boxed{
f_a(\tau) = \frac{1}{2\pi} \int_{-\infty}^{\infty} \int_{-\infty}^{\infty}
\frac{1}{\|g_{\theta^\star(t,\omega)}\|^2} \,
V_f^\star(t, \omega) \, g_{\theta^\star(t,\omega)}(\tau - t) \, \exp(j\omega \tau) \, d\omega \, dt
}
\tag{A1.27}
\]

当窗函数族归一化为单位能量 $\|g_\theta\| = 1$ 时，简化为：

\[
f_a(\tau) = \frac{1}{2\pi} \int_{-\infty}^{\infty} \int_{-\infty}^{\infty}
V_f^\star(t, \omega) \, g_{\theta^\star(t,\omega)}(\tau - t) \, \exp(j\omega \tau) \, d\omega \, dt
\tag{A1.28}
\]

#### 2.2.4 实信号的还原

对于实信号 $f(t)$，取实部重建：

\[
\boxed{
f(\tau) = \frac{1}{\pi} \, \Re \left\{ \int_{0}^{\infty} \int_{-\infty}^{\infty}
\frac{1}{\|g_{\theta^\star(t,\omega)}\|^2} \,
V_f^\star(t, \omega) \, g_{\theta^\star(t,\omega)}(\tau - t) \, \exp(j\omega \tau) \, d\omega \, dt \right\}
}
\tag{A1.29}
\]

当窗函数族归一化后：

\[
f(\tau) = \frac{1}{\pi} \, \Re \left\{ \int_{0}^{\infty} \int_{-\infty}^{\infty}
V_f^\star(t, \omega) \, g_{\theta^\star(t,\omega)}(\tau - t) \, \exp(j\omega \tau) \, d\omega \, dt \right\}
\tag{A1.30}
\]

---

## 3. 自适应同步压缩变换

同步压缩变换（Synchrosqueezing Transform, SST）是一种后处理技术：在计算得到时频表示后，将系数沿频率方向"挤压"到瞬时频率曲线上，从而获得更锐利的时频表示。将自适应参数选择与同步压缩结合，可以在参数优化的基础上进一步提升集中度。

### 3.1 自适应同步压缩小波变换

#### 3.1.1 同步压缩小波变换回顾

**瞬时频率估计**：对于CWT系数 $W_f(a, b) \neq 0$，定义瞬时频率估计：

\[
\boxed{
\omega_f(a, b) = \frac{\partial_b W_f(a, b)}{j \, W_f(a, b)}
}
\tag{A1.31}
\]

**推导验证**（纯谐波信号 $f(t) = A \exp(j\omega_0 t)$）：其CWT为：

\[
\begin{aligned}
W_f(a, b) &= \int A \exp(j\omega_0 t) \, \frac{1}{\sqrt{a}} \, \overline{\phi\!\left(\frac{t - b}{a}\right)} \, dt \\
&= \frac{A}{\sqrt{a}} \int \exp(j\omega_0(au + b)) \, \overline{\phi(u)} \, a \, du \qquad (u = (t-b)/a) \\
&= A \sqrt{a} \, \exp(j\omega_0 b) \int \overline{\phi(u)} \, \exp(j a \omega_0 u) \, du \\
&= A \sqrt{a} \, \exp(j\omega_0 b) \, \overline{\hat{\phi}(a\omega_0)}
\end{aligned}
\tag{A1.32}
\]

对 $b$ 求偏导：

\[
\partial_b W_f(a, b) = j\omega_0 \cdot A \sqrt{a} \, \exp(j\omega_0 b) \, \overline{\hat{\phi}(a\omega_0)} = j\omega_0 \, W_f(a, b)
\tag{A1.33}
\]

因此 $\omega_f(a, b) = j\omega_0 W_f(a, b) / (j W_f(a, b)) = \omega_0$，精确还原信号频率。

**同步压缩操作**：将CWT系数从尺度坐标 $a$ 重新分配（reassign）到频率坐标 $\omega = \omega_f(a, b)$：

\[
\boxed{
S_f(\omega, b) = \int_{\{a:\, W_f(a,b) \neq 0\}} W_f(a, b) \, a^{-3/2} \, \delta(\omega - \omega_f(a, b)) \, da
}
\tag{A1.34}
\]

权重 $a^{-3/2}$ 来自重建公式中 $da/a^2$ 的测度与 $1/\sqrt{a}$ 归一化的组合。

同步压缩的重建公式保持了信号的可重建性：

\[
f(t) = \frac{1}{C_\phi} \, \Re \left\{ \frac{1}{2\pi} \int_{0}^{\infty} S_f(\omega, t) \, d\omega \right\}
\tag{A1.35}
\]

#### 3.1.2 自适应同步压缩CWT

将自适应参数选择融入同步压缩框架：

**第一步**：计算自适应CWT系数 $W_f^\star(a, b)$（见 (9.6)）。

**第二步**：基于自适应系数估计瞬时频率：

\[
\omega_f^\star(a, b) = \frac{\partial_b W_f^\star(a, b)}{j \, W_f^\star(a, b)}, \qquad W_f^\star(a, b) \neq 0
\tag{A1.36}
\]

这里 $\partial_b W_f^\star(a, b)$ 的计算需要注意：由于自适应系数 $W_f^\star(a, b)$ 在相邻的 $(a, b)$ 网格点可能使用不同的母小波参数 $\theta^\star$，直接数值差分可能引入不连续性。实践中，对给定的 $(a, b)$，使用 $\theta^\star(a, b)$ 对应的母小波计算CWT在 $b$ 附近的值，再进行数值微分。

**第三步**：执行同步压缩操作：

\[
\boxed{
S_f^\star(\omega, b) = \int_{\{a:\, W_f^\star(a,b) \neq 0\}}
W_f^\star(a, b) \, a^{-3/2} \, \delta(\omega - \omega_f^\star(a, b)) \, da
}
\tag{A1.37}
\]

**数值离散化**：将频率轴离散化为 $\{\omega_\ell\}_{\ell=1}^{L}$，$\delta$ 函数用有限宽度近似：

\[
S_f^\star(\omega_\ell, b) = \sum_{a_k:\, |\omega_f^\star(a_k, b) - \omega_\ell| \le \Delta\omega/2}
W_f^\star(a_k, b) \, a_k^{-3/2} \, \Delta a_k
\tag{A1.38}
\]

其中 $\Delta\omega$ 为频率分辨率，$\Delta a_k$ 为尺度步长。

**重建公式**：

\[
f(t) = \Re \left\{ \frac{1}{2\pi} \sum_{\ell} S_f^\star(\omega_\ell, t) \, \Delta\omega \right\}
\tag{A1.39}
\]

---

### 3.2 自适应同步压缩短时傅里叶变换

#### 3.2.1 同步压缩STFT回顾

**局部瞬时频率估计**：对于STFT系数 $V_g f(t, \omega) \neq 0$，定义：

\[
\boxed{
\hat{\omega}_f(t, \omega) = \omega - \Im\left\{ \frac{V_{g'} f(t, \omega)}{V_g f(t, \omega)} \right\}
}
\tag{A1.40}
\]

其中 $g'(t) = dg(t)/dt$ 为窗函数的导数，$V_{g'} f(t, \omega)$ 为使用导数窗的STFT：

\[
V_{g'} f(t, \omega) = \int_{-\infty}^{\infty} f(\tau) \, \overline{g'(\tau - t)} \, \exp(-j\omega \tau) \, d\tau
\tag{A1.41}
\]

**推导验证**（纯谐波信号 $f(t) = A \exp(j\omega_0 t)$）：

第一步：计算 $V_g f(t, \omega)$。利用傅里叶变换的性质：

\[
\begin{aligned}
V_g f(t, \omega) &= \int A \exp(j\omega_0 \tau) \, \overline{g(\tau - t)} \, \exp(-j\omega \tau) \, d\tau \\
&= A \exp(-j(\omega - \omega_0)t) \int \overline{g(u)} \, \exp(-j(\omega - \omega_0)u) \, du \qquad (u = \tau - t) \\
&= A \exp(-j(\omega - \omega_0)t) \, \overline{\tilde{g}(\omega - \omega_0)}
\end{aligned}
\tag{A1.42}
\]

其中 $\tilde{g}(\omega) = \int g(t) \exp(-j\omega t) dt$。

第二步：计算 $V_{g'} f(t, \omega)$。对 $g'$ 的傅里叶变换，利用 $\int g'(t) \exp(-j\omega t) dt = j\omega \, \tilde{g}(\omega)$（分部积分，边界项为零）。因此：

\[
\begin{aligned}
V_{g'} f(t, \omega) &= A \exp(-j(\omega - \omega_0)t) \, \overline{j(\omega - \omega_0) \, \tilde{g}(\omega - \omega_0)} \\
&= -j(\omega - \omega_0) \, A \exp(-j(\omega - \omega_0)t) \, \overline{\tilde{g}(\omega - \omega_0)} \\
&= -j(\omega - \omega_0) \, V_g f(t, \omega)
\end{aligned}
\tag{A1.43}
\]

第三步：代入 (9.40) 验证：

\[
\hat{\omega}_f(t, \omega) = \omega - \Im\left\{ \frac{-j(\omega - \omega_0) V_g f(t, \omega)}{V_g f(t, \omega)} \right\}
= \omega - \Im\{-j(\omega - \omega_0)\} = \omega - (-(\omega - \omega_0)) = \omega_0
\tag{A1.44}
\]

即对于纯谐波信号，局部瞬时频率估计精确还原了 $\omega_0$。

**同步压缩STFT的定义**：

\[
\boxed{
S_f(t, \eta) = \int_{\{\omega:\, V_g f(t,\omega) \neq 0\}} V_g f(t, \omega) \, \delta(\eta - \hat{\omega}_f(t, \omega)) \, d\omega
}
\tag{A1.45}
\]

**重建公式**：同步压缩STFT保持了信号的可重建性：

\[
f(\tau) = \frac{1}{2\pi \, \overline{g(0)}} \int_{-\infty}^{\infty} S_f(t, \eta) \, \exp(j\eta \tau) \, d\eta
\tag{A1.46}
\]

其中 $g(0)$ 是窗函数在原点的取值（非零假设）。

#### 3.2.2 自适应同步压缩STFT

将自适应窗参数选择与STFT同步压缩结合：

**第一步**：计算自适应STFT系数 $V_f^\star(t, \omega)$（见 (9.23)）。

**第二步**：基于自适应系数估计局部瞬时频率。由于每个 $(t, \omega)$ 处使用不同的窗函数 $g_{\theta^\star(t,\omega)}$，瞬时频率估计需要对应使用该位置的导数窗：

\[
\hat{\omega}_f^\star(t, \omega) = \omega - \Im\left\{ \frac{V_{g'_{\theta^\star}} f(t, \omega)}{V_f^\star(t, \omega)} \right\}
\tag{A1.47}
\]

其中 $g'_\theta(t) = dg_\theta(t)/dt$，且：

\[
V_{g'_{\theta^\star}} f(t, \omega) = \int_{-\infty}^{\infty} f(\tau) \, \overline{g'_{\theta^\star(t,\omega)}(\tau - t)} \, \exp(-j\omega \tau) \, d\tau
\tag{A1.48}
\]

**第三步**：执行同步压缩：

\[
\boxed{
S_f^\star(t, \eta) = \int_{\{\omega:\, V_f^\star(t,\omega) \neq 0\}}
V_f^\star(t, \omega) \, \delta(\eta - \hat{\omega}_f^\star(t, \omega)) \, d\omega
}
\tag{A1.49}
\]

**数值离散化**：

\[
S_f^\star(t, \eta_\ell) = \sum_{\omega_k:\, |\hat{\omega}_f^\star(t, \omega_k) - \eta_\ell| \le \Delta\eta/2}
V_f^\star(t, \omega_k) \, \Delta\omega_k
\tag{A1.50}
\]

**重建公式**（利用归一化窗函数族 $\|g_\theta\| = 1$）：

\[
f(\tau) = \frac{1}{2\pi \, \overline{g_{\theta^\star}(0)}} \int_{-\infty}^{\infty}
S_f^\star(t, \eta) \, \exp(j\eta \tau) \, d\eta
\tag{A1.51}
\]

---

## 4. 自适应二阶同步压缩变换

一阶同步压缩变换依赖于纯谐波信号的假设——它假定在每个局部区域内信号可以近似为 $A \exp(j\omega_0 t)$，即频率恒定。当信号分量具有**线性频率调制**（chirp信号 $f(t) = A \exp(j(\omega_0 t + \frac{1}{2}\alpha t^2))$）时，一阶瞬时频率估计存在偏差，导致压缩后的能量仍然沿频率轴扩散。二阶同步压缩变换通过引入**调频率估计**来修正这一偏差。

### 4.1 二阶同步压缩变换

#### 4.1.1 二阶瞬时频率估计

**CWT框架下的二阶估计**：

对于CWT系数 $W_f(a, b) \neq 0$，定义**群延迟估计**（group delay）：

\[
\hat{t}_f(a, b) = b + \Re\left\{ \frac{a \cdot W_f^{t\phi}(a, b)}{W_f(a, b)} \right\}
\tag{A1.52}
\]

其中 $W_f^{t\phi}(a, b)$ 是使用 $t \cdot \phi(t)$（小波乘以时间变量）作为母小波的CWT：

\[
W_f^{t\phi}(a, b) = \int_{-\infty}^{\infty} f(t) \, \frac{1}{\sqrt{a}} \, \overline{\frac{t-b}{a} \, \phi\!\left(\frac{t-b}{a}\right)} \, dt
\tag{A1.53}
\]

进一步，定义**调频率估计**（chirp rate estimate）：

\[
\boxed{
\hat{\alpha}_f(a, b) = \frac{\partial_b \omega_f(a, b)}{1 - \partial_b \hat{t}_f(a, b)}
}
\tag{A1.54}
\]

其中 $\omega_f(a, b)$ 是一阶瞬时频率估计 (9.31)，$\partial_b$ 表示对平移参数 $b$ 的偏导数。

**STFT框架下的二阶估计**：

对于STFT系数 $V_g f(t, \omega) \neq 0$，定义群延迟估计：

\[
\hat{t}_f(t, \omega) = t + \Re\left\{ \frac{V_{tg} f(t, \omega)}{V_g f(t, \omega)} \right\}
\tag{A1.55}
\]

其中 $V_{tg} f(t, \omega)$ 是使用 $t \cdot g(t)$ 作为窗函数的STFT。

调频率估计：

\[
\boxed{
\hat{\alpha}_f(t, \omega) = \frac{\partial_t \hat{\omega}_f(t, \omega)}{1 - \partial_t \hat{t}_f(t, \omega)}
}
\tag{A1.56}
\]

其中 $\hat{\omega}_f(t, \omega)$ 是一阶瞬时频率估计 (9.40)，$\partial_t$ 表示对时间参数 $t$ 的偏导数。

**二阶瞬时频率**：在一阶估计的基础上加入调频率修正：

对于CWT：

\[
\boxed{
\tilde{\omega}_f(a, b) = \omega_f(a, b) + \hat{\alpha}_f(a, b) \cdot (b - \hat{t}_f(a, b))
}
\tag{A1.57}
\]

对于STFT：

\[
\boxed{
\tilde{\omega}_f(t, \omega) = \hat{\omega}_f(t, \omega) + \hat{\alpha}_f(t, \omega) \cdot (t - \hat{t}_f(t, \omega))
}
\tag{A1.58}
\]

#### 4.1.2 二阶同步压缩CWT

将一阶瞬时频率 $\omega_f(a, b)$ 替换为二阶瞬时频率 $\tilde{\omega}_f(a, b)$：

\[
\boxed{
S_f^{(2)}(\omega, b) = \int_{\{a:\, W_f(a,b) \neq 0\}}
W_f(a, b) \, a^{-3/2} \, \delta(\omega - \tilde{\omega}_f(a, b)) \, da
}
\tag{A1.59}
\]

**对于chirp信号的改进**：设 $f(t) = A \exp(j(\omega_0 t + \frac{1}{2}\alpha t^2))$。一阶瞬时频率估计 $\omega_f(a, b)$ 存在与 $\alpha$ 成比例的偏差，导致同步压缩后的能量分布在 $\omega_0 + \alpha t$ 附近而非精确集中。二阶修正 $\tilde{\omega}_f(a, b)$ 通过减去 $\hat{\alpha}_f(a, b) \cdot \hat{t}_f(a, b)$ 项来补偿这一偏差，使得压缩后的能量精确沿 $\omega_0 + \alpha t$ 分布。

#### 4.1.3 二阶同步压缩STFT

类似地，对STFT执行二阶同步压缩：

\[
\boxed{
S_f^{(2)}(t, \eta) = \int_{\{\omega:\, V_g f(t,\omega) \neq 0\}}
V_g f(t, \omega) \, \delta(\eta - \tilde{\omega}_f(t, \omega)) \, d\omega
}
\tag{A1.60}
\]

---

### 4.2 自适应二阶同步压缩

将自适应参数选择与二阶同步压缩结合。核心思路：在参数优化的基础上，进一步利用二阶调频率信息来修正瞬时频率估计。

#### 4.2.1 自适应二阶同步压缩CWT

**第一步**：计算自适应CWT系数 $W_f^\star(a, b)$ 和最优参数 $\theta^\star(a, b)$（见 (9.6)-(9.7)）。

**第二步**：基于自适应系数估计一阶瞬时频率 $\omega_f^\star(a, b)$（见 (9.36)）。

**第三步**：基于自适应系数估计群延迟和调频率。注意由于每个 $(a, b)$ 处使用的母小波不同，需要统一使用对应参数的基函数：

\[
\hat{t}_f^\star(a, b) = b + \Re\left\{ \frac{a \cdot W_f^{\star, t\phi}(a, b)}{W_f^\star(a, b)} \right\}
\tag{A1.61}
\]

其中 $W_f^{\star, t\phi}(a, b)$ 是使用 $t \cdot \phi_{\theta^\star(a,b)}(t)$ 作为母小波的CWT：

\[
W_f^{\star, t\phi}(a, b) = \int_{-\infty}^{\infty} f(t) \, \frac{1}{\sqrt{a}} \, \overline{\frac{t-b}{a} \, \phi_{\theta^\star(a,b)}\!\left(\frac{t-b}{a}\right)} \, dt
\tag{A1.62}
\]

调频率估计：

\[
\hat{\alpha}_f^\star(a, b) = \frac{\partial_b \omega_f^\star(a, b)}{1 - \partial_b \hat{t}_f^\star(a, b)}
\tag{A1.63}
\]

**第四步**：计算自适应二阶瞬时频率：

\[
\boxed{
\tilde{\omega}_f^\star(a, b) = \omega_f^\star(a, b) + \hat{\alpha}_f^\star(a, b) \cdot (b - \hat{t}_f^\star(a, b))
}
\tag{A1.64}
\]

**第五步**：执行自适应二阶同步压缩：

\[
\boxed{
S_f^{\star(2)}(\omega, b) = \int_{\{a:\, W_f^\star(a,b) \neq 0\}}
W_f^\star(a, b) \, a^{-3/2} \, \delta(\omega - \tilde{\omega}_f^\star(a, b)) \, da
}
\tag{A1.65}
\]

#### 4.2.2 自适应二阶同步压缩STFT

**第一步**：计算自适应STFT系数 $V_f^\star(t, \omega)$ 和最优参数 $\theta^\star(t, \omega)$。

**第二步**：计算一阶瞬时频率 $\hat{\omega}_f^\star(t, \omega)$（见 (9.47)）。

**第三步**：估计群延迟和调频率：

\[
\hat{t}_f^\star(t, \omega) = t + \Re\left\{ \frac{V_{tg_{\theta^\star}} f(t, \omega)}{V_f^\star(t, \omega)} \right\}
\tag{A1.66}
\]

其中 $V_{tg_{\theta^\star}} f(t, \omega)$ 是使用 $t \cdot g_{\theta^\star(t,\omega)}(t)$ 作为窗函数的STFT。

\[
\hat{\alpha}_f^\star(t, \omega) = \frac{\partial_t \hat{\omega}_f^\star(t, \omega)}{1 - \partial_t \hat{t}_f^\star(t, \omega)}
\tag{A1.67}
\]

**第四步**：计算自适应二阶瞬时频率：

\[
\boxed{
\tilde{\omega}_f^\star(t, \omega) = \hat{\omega}_f^\star(t, \omega) + \hat{\alpha}_f^\star(t, \omega) \cdot (t - \hat{t}_f^\star(t, \omega))
}
\tag{A1.68}
\]

**第五步**：执行自适应二阶同步压缩：

\[
\boxed{
S_f^{\star(2)}(t, \eta) = \int_{\{\omega:\, V_f^\star(t,\omega) \neq 0\}}
V_f^\star(t, \omega) \, \delta(\eta - \tilde{\omega}_f^\star(t, \omega)) \, d\omega
}
\tag{A1.69}
\]

---

## 5. 自适应参数选择

前面的推导中，参数优化准则 $\theta^\star = \arg\max_\theta |\langle f, \phi_\theta \rangle|$（局部幅度最大化）是最直接的思路，但并非唯一选择。本节系统讨论自适应参数选择的多种准则。

### 5.1 浓度测度最大化

**定义（浓度测度）**：对于时频表示 $T_f(\xi)$（$\xi$ 为时频坐标），定义 $\gamma$ 阶浓度测度：

\[
\mathcal{C}_\gamma(T_f) = \frac{\iint |T_f(\xi)|^{2\gamma} \, d\xi}{\left(\iint |T_f(\xi)|^2 \, d\xi\right)^\gamma}, \qquad \gamma > 1
\tag{A1.70}
\]

该测度衡量的是时频能量分布的"峰值度"：当能量高度集中于少数点时，分子的 $\gamma$ 次幂放大了峰值效应，使得 $\mathcal{C}_\gamma$ 较大（$\gamma = 2$ 时即为常用的kurtosis型测度）。

对于自适应参数选择，定义局部浓度测度（在 $(a, b)$ 附近的小邻域内）：

\[
\mathcal{C}_\gamma^{\text{local}}(a, b; \theta) = \frac{\iint_{\mathcal{N}(a,b)} |W_f^\theta(a', b')|^{2\gamma} \, \frac{da' db'}{a'^2}}{\left(\iint_{\mathcal{N}(a,b)} |W_f^\theta(a', b')|^2 \, \frac{da' db'}{a'^2}\right)^\gamma}
\tag{A1.71}
\]

其中 $\mathcal{N}(a, b)$ 为以 $(a, b)$ 为中心的邻域。参数选择为：

\[
\theta^\star(a, b) = \arg\max_{\theta \in \Theta} \; \mathcal{C}_\gamma^{\text{local}}(a, b; \theta)
\tag{A1.72}
\]

### 5.2 熵最小化准则

**定义（$\alpha$ 阶Renyi熵）**：见 (9.10)。对于时频表示 $T_f$，其归一化分布为：

\[
P_f(\xi) = \frac{|T_f(\xi)|^2}{\iint |T_f(\xi')|^2 \, d\xi'}
\tag{A1.73}
\]

$\alpha$ 阶Renyi熵为：

\[
H_\alpha(T_f) = \frac{1}{1 - \alpha} \log_2 \iint [P_f(\xi)]^\alpha \, d\xi'
\tag{A1.74}
\]

常用的 $\alpha = 3$ 赋予峰值更高权重，有利于选择使时频表示最集中的参数。

**局部熵最小化**：在时频平面的每个点附近选择使局部熵最小的参数：

\[
\theta^\star(a, b) = \arg\min_{\theta \in \Theta} \; H_\alpha^{\text{local}}(a, b; \theta)
\tag{A1.75}
\]

其中 $H_\alpha^{\text{local}}$ 是在邻域 $\mathcal{N}(a, b)$ 内计算的局部Renyi熵。

### 5.3 局部自适应策略

实际计算中，遍历整个参数空间 $\Theta$ 对所有 $(a, b)$ 进行穷举搜索的计算代价为 $O(M \cdot |\Theta|)$ 倍于标准变换，这在很多应用中难以承受。以下是几种降低计算复杂度的策略。

**策略一：参数空间粗-精细搜索**

先在粗糙的参数网格上进行搜索，确定最优参数的大致范围，再在该范围内进行精细搜索：

\[
\Theta^{\text{coarse}} = \{\theta_1, \theta_2, \ldots, \theta_{K_c}\}, \quad \theta^\star_{\text{coarse}} = \arg\max_{\theta \in \Theta^{\text{coarse}}} \mathcal{C}(\theta)
\tag{A1.76}
\]

\[
\Theta^{\text{fine}} = \{\theta : |\theta - \theta^\star_{\text{coarse}}| < \delta\}, \quad \theta^\star = \arg\max_{\theta \in \Theta^{\text{fine}}} \mathcal{C}(\theta)
\tag{A1.77}
\]

**策略二：分段常数参数**

假设信号在时频平面上的特征变化是缓慢的，将时频平面划分为较大的"超像素"块，在每块内使用统一的参数：

\[
\theta^\star(\mathcal{R}_m) = \arg\max_{\theta \in \Theta} \; \iint_{\mathcal{R}_m} |W_f^\theta(a, b)|^2 \, \frac{da \, db}{a^2}
\tag{A1.78}
\]

其中 $\{\mathcal{R}_m\}$ 是时频平面的一个分划。这既减少了参数选择次数，也避免了相邻点参数突变带来的不连续性。

**策略三：基于信号特征的直接推断**

对于某些参数化族，最优参数可以直接从信号的局部特征推断，无需穷举搜索。例如，对于高斯窗族 $g_\sigma$，在时频点 $(t, \omega)$ 附近的最优窗宽 $\sigma$ 与信号的局部调频率 $\alpha$ 和频率 $\omega$ 有关：

\[
\sigma^\star(t, \omega) \propto \frac{1}{\sqrt{|\alpha(t, \omega)|}}
\tag{A1.79}
\]

即调频率越大（频率变化越快），需要的窗越窄（时间分辨率要求越高）。这一关系可以从测不准原理和 chirp 信号的局部表示特性导出。

**策略四：迭代自适应**

先使用初始参数 $\theta^{(0)}$ 计算时频表示，从中估计信号的局部特征（如瞬时频率、调频率），再根据这些特征更新参数 $\theta^{(1)}$，重复至收敛：

\[
\theta^{(n+1)} = \mathcal{F}(\theta^{(n)}, \hat{\omega}_f^{(n)}, \hat{\alpha}_f^{(n)})
\tag{A1.80}
\]

其中 $\hat{\omega}_f^{(n)}$ 和 $\hat{\alpha}_f^{(n)}$ 是从第 $n$ 轮时频表示中提取的特征。这个迭代过程结合了自适应参数选择和信号特征提取，可以逐步逼近最优解。

### 5.4 参数化族的可容许性检查

无论采用何种参数选择策略，最终选定的参数 $\theta^\star$ 必须确保对应的小波满足可容许性条件 (9.4)。对于不直接满足可容许性的参数，需要引入修正因子或用最接近的可容许参数近似替代。

典型的情况是：当 $\sigma$ 过小（时域极窄的Morlet小波）时，频域可能包含显著的低频成分，导致可容许性积分发散。实践中通过设置参数下界 $\theta_{\min}$ 来避免这一问题：

\[
\Theta = \{\theta : C_{\phi_\theta} < \infty\} \cap \{\theta : \theta \ge \theta_{\min}\}
\tag{A1.81}
\]

---

**本附录的核心结论**：

1. **自适应CWT/STFT**：通过局部参数优化，使基函数形状与信号的局部时频特征相匹配，获得比固定基函数更集中的时频表示。重建公式在逐点自适应和全局自适应两种情形下具有不同的形式。

2. **自适应同步压缩**：在自适应参数选择的基础上施加同步压缩操作，进一步将时频能量挤压到瞬时频率曲线上，是参数优化和频率重分配两种增强策略的叠加。重建公式保持了信号的可重建性。

3. **自适应二阶同步压缩**：引入调频率的二阶修正，克服了一阶同步压缩对chirp类信号的频率估计偏差，在自适应框架下同时实现了参数自适应和高阶瞬时频率修正。

4. **参数选择策略**：从穷举搜索到基于特征的推断，计算代价与表示精度之间存在折中。分段常数策略和迭代自适应策略在实际应用中最为实用。



<div style="page-break-before: always;"></div>