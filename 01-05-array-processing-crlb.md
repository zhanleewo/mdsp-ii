<div style="page-break-before: always; padding: 8% 8% 0 8%;">
 <h1 id="第四讲-多进多出-mimo技术" style="text-align: center; margin-bottom: 2rem; border-bottom: none; display: block;">第五讲 阵列处理的CRLB</h1> 
 <div style="display: flex; align-items: center; justify-content: center; gap: 12px; margin: 1.8rem auto;">
  <span style="flex:1; max-width:80px; height:1px; background: linear-gradient(to right, transparent, #888);"></span>
  <span style="display:inline-block; width:6px; height:6px; background:#38bdf8; border-radius:50%;"></span>
  <span style="flex:1; max-width:80px; height:1px; background: linear-gradient(to left, transparent, #888);"></span>
 </div>
</div>

<!-- # 第五讲 阵列处理的CRLB -->
> Petre Stoica & Arye Nehorai, 1989
> MUSIC, Maximum Likelihood, and Cramer-Rao Bound

## 1. 引言

### 1.1 简明扼要地回顾前面学过的CRLB

前两篇文章系统地建立了参数估计的理论基础。第三讲（MVUE）讨论了如何寻找最小方差无偏估计，引入了充分统计量、Rao-Blackwell定理、完备性与Lehmann-Scheffé定理。第四讲（CRLB）则从另一个角度——Fisher信息量——给出了无偏估计方差的理论下界。

回顾核心结论：对于参数为 \( \theta \) 的统计模型 \( f(x;\theta) \)，若 \( \hat{\theta} \) 是 \( \theta \) 的任意无偏估计，在正则条件下，其方差满足：

\[
\boxed{\operatorname{Var}(\hat{\theta}) \ge \frac{1}{I(\theta)}}
\]

其中 Fisher 信息量定义为：

\[
I(\theta) = \mathbb{E}\left[ \left( \frac{\partial}{\partial \theta} \log f(X;\theta) \right)^2 \right] = -\mathbb{E}\left[ \frac{\partial^2}{\partial \theta^2} \log f(X;\theta) \right]
\]

对于多维参数 \( \boldsymbol{\theta} \in \mathbb{R}^m \)，上述结论推广为矩阵不等式：

\[
\boxed{\operatorname{Cov}(\hat{\boldsymbol{\theta}}) \succeq \mathbf{I}(\boldsymbol{\theta})^{-1}}
\]

其中 Fisher 信息矩阵 \( \mathbf{I}(\boldsymbol{\theta}) \) 的第 \( (i,j) \) 元素为：

\[
[\mathbf{I}(\boldsymbol{\theta})]_{ij} = \mathbb{E}\left[ \frac{\partial}{\partial \theta_i} \log f(X;\boldsymbol{\theta}) \cdot \frac{\partial}{\partial \theta_j} \log f(X;\boldsymbol{\theta}) \right] = -\mathbb{E}\left[ \frac{\partial^2}{\partial \theta_i \partial \theta_j} \log f(X;\boldsymbol{\theta}) \right]
\]

进一步，若要估计参数的函数 \( \boldsymbol{g}(\boldsymbol{\theta}) \in \mathbb{R}^k \)，设 \( \hat{\boldsymbol{g}} \) 是 \( \boldsymbol{g}(\boldsymbol{\theta}) \) 的无偏估计，Jacobian 矩阵为：

\[
\mathbf{J}_{\boldsymbol{g}}(\boldsymbol{\theta}) = \frac{\partial \boldsymbol{g}}{\partial \boldsymbol{\theta}^\top} \in \mathbb{R}^{k \times m}, \quad [\mathbf{J}_{\boldsymbol{g}}]_{ij} = \frac{\partial g_i}{\partial \theta_j}
\]

则 CRLB 为：

\[
\boxed{\operatorname{Cov}(\hat{\boldsymbol{g}}) \succeq \mathbf{J}_{\boldsymbol{g}}(\boldsymbol{\theta}) \, \mathbf{I}(\boldsymbol{\theta})^{-1} \, \mathbf{J}_{\boldsymbol{g}}(\boldsymbol{\theta})^\top}
\]

这些结果构成了参数估计的理论基准——任何无偏估计的方差都不可能低于这个下界。

---

### 1.2 求 CRLB 的三步骤

在经典参数估计中，计算 CRLB 通常遵循以下三个步骤：

**步骤一：写出统计模型，并取对数**

确定数据的概率密度函数（或质量函数）\( f(\mathbf{x}; \boldsymbol{\theta}) \)，写出对数似然函数：

\[
\ell(\boldsymbol{\theta}) = \log f(\mathbf{x}; \boldsymbol{\theta})
\]

对于独立同分布样本 \( \mathbf{x} = (x_1, \dots, x_n)^\top \)，联合密度为 \( f(\mathbf{x};\boldsymbol{\theta}) = \prod_{i=1}^n f(x_i;\boldsymbol{\theta}) \)，因此：

\[
\ell(\boldsymbol{\theta}) = \sum_{i=1}^n \log f(x_i;\boldsymbol{\theta})
\]

**步骤二：计算一阶和二阶梯度**

计算得分向量（一阶梯度）：

\[
\mathbf{s}(\boldsymbol{\theta}) = \nabla_{\boldsymbol{\theta}} \ell(\boldsymbol{\theta}) \in \mathbb{R}^m
\]

即第 \( i \) 个分量为 \( \frac{\partial \ell}{\partial \theta_i} \)。

计算 Hessian 矩阵（二阶梯度）：

\[
\mathbf{H}(\boldsymbol{\theta}) = \nabla_{\boldsymbol{\theta}} \mathbf{s}(\boldsymbol{\theta})^\top \in \mathbb{R}^{m \times m}
\]

即第 \( (i,j) \) 元素为 \( \frac{\partial^2 \ell}{\partial \theta_i \partial \theta_j} \)。

**步骤三：计算 Fisher 信息矩阵，取其逆得到 CRLB**

Fisher 信息矩阵为：

\[
\mathbf{I}(\boldsymbol{\theta}) = \mathbb{E}\left[ \mathbf{s}(\boldsymbol{\theta}) \mathbf{s}(\boldsymbol{\theta})^\top \right] = -\mathbb{E}\left[ \mathbf{H}(\boldsymbol{\theta}) \right]
\]

（第二个等号在正则条件下成立。）

CRLB 即为 Fisher 信息矩阵的逆：

\[
\boxed{\operatorname{CRLB}(\boldsymbol{\theta}) = \mathbf{I}(\boldsymbol{\theta})^{-1}}
\]

若只关心标量参数或参数的某个线性组合，只需取对应元素或计算二次型即可。

对于标量参数 \( \theta \)，上述步骤简化为：

\[
\ell(\theta) = \log f(x;\theta), \quad s(\theta) = \frac{\partial \ell}{\partial \theta}, \quad I(\theta) = \mathbb{E}[s(\theta)^2] = -\mathbb{E}\left[ \frac{\partial^2 \ell}{\partial \theta^2} \right], \quad \operatorname{CRLB}(\theta) = \frac{1}{I(\theta)}
\]

---

### 1.3 为什么要研究阵列处理的 CRLB

前面三篇文章——MUSIC、ESPRIT、相关信号处理——已系统介绍了若干 DOA 估计算法。一个更深层的问题是：**这些算法能做到多好？** 以及**是否存在一个理论上不可逾越的性能上限？**

CRLB 在阵列处理中的意义体现在以下几个层面：

**第一，为 DOA 估计提供理论性能基准。**

在阵列信号处理中，处理的观测数据是复数基带信号，参数是信号源的到达角 \( \theta_1, \dots, \theta_K \)，此外还有信号波形、噪声方差等"讨厌参数"。CRLB 表明：在给定阵元数、快拍数、信噪比等条件下，任何无偏 DOA 估计量的方差都不可能低于某个理论下界。这个下界是评价所有算法的"黄金标准"——如果一个算法的性能接近 CRLB，说明它已经充分利用了数据中的信息；如果远高于 CRLB，说明还有改进空间。

**第二，帮助理解影响估计精度的因素。**

CRLB 表达式能够明确揭示哪些因素影响 DOA 估计的精度：

- **信噪比（SNR）**：CRLB 通常与 SNR 成反比——信噪比越高，下界越低。
- **快拍数 \( L \)**：由于 Fisher 信息对独立快拍具有可加性，CRLB 通常与 \( 1/L \) 成正比——快拍越多，估计越准。
- **阵元数 \( M \)**：更多的阵元提供更大的阵列孔径和更多的空间采样点，CRLB 随之降低。
- **信号源角度间隔**：两个角度靠得越近，CRLB 越大——这是分辨率极限的理论根源。
- **信号源之间的相关性**：相干信号会导致 Fisher 信息矩阵奇异，CRLB 趋于无穷——这是空间平滑等预处理技术的理论基础。

定量地理解这些依赖关系，有助于在系统设计阶段做出合理的权衡。

**第三，区分"算法本身的缺陷"与"问题固有的困难"。**

如果某个 DOA 估计算法的性能接近 CRLB，那么再更换更复杂的算法也不会带来显著提升——瓶颈在于数据本身的信息量。反之，如果算法性能远高于 CRLB，说明算法没有充分利用数据中的信息（如常规波束形成），或者算法存在模型失配问题。这为算法选择提供了量化依据。

**第四，为相干信号处理和秩亏问题提供理论视角。**

在第三讲"相关信号"中，讨论了信号协方差矩阵 \( \mathbf{R}_S \) 秩亏导致的子空间方法失效。CRLB 框架可以精确描述秩亏的程度如何影响角度估计的理论极限，从而解释为什么在极端相干条件下（如完全相干），DOA 估计变得极其困难甚至不可能。

**第五，为阵列设计提供优化准则。**

给定阵元数量，最优的阵元布局是什么？CRLB 提供了定量评价标准——使 Fisher 信息矩阵"最大"（在某种意义下）的阵列结构就是最优的。这直接关系到稀疏阵列、嵌套阵列、MIMO 雷达等现代阵列系统的设计。

---
## 2. 阵列处理的CRLB推导

### 2.1 统计模型的建立

阵列信号处理的基本数据模型为：

\[
\boxed{X = A(\theta) \cdot S + N}
\tag{5.1}
\]

其中：

- \( X \in \mathbb{C}^{N \times L} \) 是接收数据矩阵，\( N \) 为阵元数，\( L \) 为快拍数；
- \( A(\theta) \in \mathbb{C}^{N \times M} \) 是方向矩阵，其第 \( k \) 列为 \( \mathbf{a}(\theta_k) \)，对应第 \( k \) 个信号源的导向矢量；
- \( S \in \mathbb{C}^{M \times L} \) 是信号波形矩阵；
- \( N \in \mathbb{C}^{N \times L} \) 是加性高斯白噪声。

对于噪声，假设：

\[
N \sim \mathcal{CN}(0, \sigma^2 I)
\tag{5.2}
\]

即噪声服从零均值复高斯分布，各阵元独立同分布，协方差矩阵为 \( \sigma^2 I \)。

对信号矩阵 \( S \) 的统计特性做出假设。

如果将 \( S \) 视为随机变量，则需要对其分布做出具体假设。例如，假设 \( S \) 服从零均值复高斯分布，协方差矩阵为 \( \mathbf{R}_S \)，则接收数据的协方差矩阵为：

\[
\mathbf{R}_X = A(\theta) \mathbf{R}_S A^H(\theta) + \sigma^2 I
\tag{5.3}
\]

此时，待估计的参数包括 \( \theta \)、\( \mathbf{R}_S \) 和 \( \sigma^2 \)。由于 \( \mathbf{R}_S \) 是一个 \( M \times M \) 的 Hermitian 矩阵，包含 \( M^2 \) 个实参数（\( M(M+1)/2 \) 个自由度），参数空间维度大幅增加，CRLB 的推导变得复杂。此外，信号波形 \( S \) 的随机性会引入额外的统计起伏，使得 Fisher 信息的计算需要同时对 \( S \) 和 \( N \) 取期望，涉及复杂的矩阵积分。

更重要的是，**DOA估计并不关心信号的具体内容**——无论 \( S \) 是正弦波、QAM 符号还是脉冲信号，关注的唯有来波方向 \( \theta \)。因此，将 \( S \) 作为确定性但未知的参数来处理，在统计推断中是完全合理的。这种做法被称为 **"确定性信号模型"** 或 **"条件信号模型"**。

在确定性信号模型下，\( S \) 被视为未知的确定性常数，数据的随机性完全来源于噪声 \( N \)。此时，给定参数 \( \theta \) 和 \( S \)，接收数据 \( X \) 服从复高斯分布：

\[
\boxed{X \sim \mathcal{CN}(A(\theta) \cdot S, \; \sigma^2 I)}
\tag{5.4}
\]

即均值为 \( A(\theta)S \)，协方差矩阵为 \( \sigma^2 I \)。

因此，完整的待估参数向量为：

\[
\boxed{\boldsymbol{\phi} = (\theta_1, \theta_2, \cdots, \theta_M, \; \text{Re}(S), \; \text{Im}(S), \; \sigma^2)}
\tag{5.5}
\]

其中 \( \theta_1, \cdots, \theta_M \) 是 \( M \) 个信号源的到达角——这是关注参数。\( \text{Re}(S) \) 和 \( \text{Im}(S) \) 是信号矩阵 \( S \) 的实部和虚部，它们是"讨厌参数"（nuisance parameters）——必须估计它们，但最终不关心。\( \sigma^2 \) 是噪声方差，同样是讨厌参数。

令 \( C = \sigma^2 I \)。则接收数据 \( X \) 的复高斯概率密度函数为：

\[
\boxed{
f_X(x; \boldsymbol{\phi}) = \frac{1}{(2\pi)^{n/2} \sqrt{\det(C)}} \exp\left( -\frac{1}{2} (x - A(\theta) \cdot S)^T C^{-1} (x - A(\theta) \cdot S) \right)
}
\tag{5.6}
\]

这里，\( x \) 是接收数据的向量化形式。\( n \) 是数据的维度（对应阵元数 \( N \)）。\( C = \sigma^2 I \) 是协方差矩阵。

对 (5.6) 两边取自然对数，得到对数似然函数。这里的 \( \log \) 是自然对数（以 \( e \) 为底）。展开如下：

\[
\log f_X(x; \boldsymbol{\phi}) = \log \left[ \frac{1}{(2\pi)^{n/2} \sqrt{\det(C)}} \exp\left( -\frac{1}{2} (x - A(\theta) \cdot S)^T C^{-1} (x - A(\theta) \cdot S) \right) \right]
\tag{5.7}
\]

利用 \( \log(ab) = \log a + \log b \)，将指数函数前的系数与指数项分开：

\[
\log f_X(x; \boldsymbol{\phi}) = \log \left( \frac{1}{(2\pi)^{n/2}} \right) + \log \left( \frac{1}{\sqrt{\det(C)}} \right) + \log \left[ \exp\left( -\frac{1}{2} (x - A(\theta) \cdot S)^T C^{-1} (x - A(\theta) \cdot S) \right) \right]
\tag{5.8}
\]

逐一化简三项：

第一项：
\[
\log \left( \frac{1}{(2\pi)^{n/2}} \right) = -\frac{n}{2} \log(2\pi)
\tag{5.9}
\]

第二项：
\[
\log \left( \frac{1}{\sqrt{\det(C)}} \right) = -\frac{1}{2} \log \det(C)
\tag{5.10}
\]

第三项：
\[
\log \left[ \exp\left( -\frac{1}{2} (x - A(\theta) \cdot S)^T C^{-1} (x - A(\theta) \cdot S) \right) \right] = -\frac{1}{2} (x - A(\theta) \cdot S)^T C^{-1} (x - A(\theta) \cdot S)
\tag{5.11}
\]

将 (5.9)、(5.10)、(5.11) 代回 (5.8)，得到：

\[
\boxed{
\log f_X(x; \boldsymbol{\phi}) = -\frac{n}{2} \log(2\pi) - \frac{1}{2} \log \det(C) - \frac{1}{2} (x - A(\theta) \cdot S)^T C^{-1} (x - A(\theta) \cdot S)
}
\tag{5.12}
\]

(5.12) 式是后续所有 CRLB 推导的出发点。它由三部分组成：

1. **常数项**：\( -\frac{n}{2} \log(2\pi) \)，与参数无关，求导后消失；
2. **协方差矩阵项**：\( -\frac{1}{2} \log \det(C) \)，与 \( C \) 有关（在阵列问题中与 \( \sigma^2 \) 有关），但与角度 \( \theta \) 和信号 \( S \) 无关；
3. **数据拟合项**：\( -\frac{1}{2} (x - A(\theta) \cdot S)^T C^{-1} (x - A(\theta) \cdot S) \)，包含了全部角度信息和信号信息，是 Fisher 信息的主要来源。

---

### 2.2 一般化：多维高斯模型的 CRLB

在正式计算 DOA 估计的 CRLB 之前，先将问题一般化，推导一个更通用的结果——**多维高斯分布下参数估计的 CRLB 通式**。

设观测数据服从复高斯分布：

\[
x \sim \mathcal{N}(\boldsymbol{\mu}(\boldsymbol{\phi}), \; C(\boldsymbol{\phi}))
\tag{5.13}
\]

其中 \( \boldsymbol{\mu}(\boldsymbol{\phi}) \) 是均值向量，\( C(\boldsymbol{\phi}) \) 是协方差矩阵，两者都依赖于待估参数向量 \( \boldsymbol{\phi} \)。

在 DOA 估计问题中：
- \( \boldsymbol{\mu}(\boldsymbol{\phi}) = A(\theta)S \)；
- \( C(\boldsymbol{\phi}) = \sigma^2 I \)。

在下面的推导中，先不代入具体形式，而是推导一个适用于任意 \( \boldsymbol{\mu}(\boldsymbol{\phi}) \) 和 \( C(\boldsymbol{\phi}) \) 的通用 CRLB 表达式。

一般化的对数似然函数为：

\[
\ell(\boldsymbol{\phi}) = \log f_X(x; \boldsymbol{\phi}) = -\frac{n}{2} \log(2\pi) - \frac{1}{2} \log \det(C(\boldsymbol{\phi})) - \frac{1}{2} (x - \boldsymbol{\mu}(\boldsymbol{\phi}))^T C(\boldsymbol{\phi})^{-1} (x - \boldsymbol{\mu}(\boldsymbol{\phi}))
\tag{5.14}
\]

根据 CRLB 的求导步骤，分别对这两类参数求导：一类是出现在均值 \( \boldsymbol{\mu} \) 中的参数（如 \( \theta \) 和 \( S \)），另一类是出现在协方差矩阵 \( C \) 中的参数（如 \( \sigma^2 \)）。

---

### 2.3 \( \log \det(C(\phi)) \) 的导数：从特征值分解到迹形式的完整推导

计算对数似然函数 (5.14) 中的一项——\( \log \det(C(\phi)) \) 对参数 \( \phi_i \) 的导数。

从 (5.14) 中提取与 \( C \) 相关的项：

\[
\ell_C(\boldsymbol{\phi}) = -\frac{1}{2} \log \det(C(\boldsymbol{\phi}))
\tag{5.15}
\]

对其求偏导，关键在于计算：

\[
\frac{\partial}{\partial \phi_i} \log \det(C(\phi))
\tag{5.16}
\]

通过特征值分解来推导这个结果。

假设 \( C = C(\phi) \) 是 \( n \times n \) 的对称正定矩阵，其特征值分解为：

\[
\boxed{C(\phi) = U^T(\phi) \, \Lambda(\phi) \, U(\phi)}
\tag{5.17}
\]

其中：
- \( \Lambda(\phi) = \operatorname{diag}(\lambda_1(\phi), \lambda_2(\phi), \cdots, \lambda_n(\phi)) \) 是对角矩阵，\( \lambda_k(\phi) \) 是 \( C(\phi) \) 的第 \( k \) 个特征值；
- \( U(\phi) \in \mathbb{R}^{n \times n} \) 是正交矩阵（\( U^T U = I \)），其列向量为对应的特征向量。

由于 \( C(\phi) \) 是对称正定矩阵，所有特征值均为正实数：\( \lambda_k(\phi) > 0 \)。

矩阵的行列式等于其特征值的乘积：

\[
\det(C(\phi)) = \prod_{k=1}^n \lambda_k(\phi)
\tag{5.18}
\]

取对数后，乘积变为求和：

\[
\log \det(C(\phi)) = \sum_{k=1}^n \log \lambda_k(\phi)
\tag{5.19}
\]

对参数 \( \phi_i \) 求偏导：

\[
\frac{\partial}{\partial \phi_i} \log \det(C) = \sum_{k=1}^n \frac{1}{\lambda_k(\phi)} \frac{\partial \lambda_k(\phi)}{\partial \phi_i}
\tag{5.20}
\]

写成矩阵形式。定义：

\[
\Lambda^{-1}(\phi) = \operatorname{diag}\left( \frac{1}{\lambda_1(\phi)}, \frac{1}{\lambda_2(\phi)}, \cdots, \frac{1}{\lambda_n(\phi)} \right)
\tag{5.21}
\]

以及

\[
\frac{\partial \Lambda(\phi)}{\partial \phi_i} = \operatorname{diag}\left( \frac{\partial \lambda_1}{\partial \phi_i}, \frac{\partial \lambda_2}{\partial \phi_i}, \cdots, \frac{\partial \lambda_n}{\partial \phi_i} \right)
\tag{5.22}
\]

于是 (5.20) 可以写成迹的形式：

\[
\frac{\partial}{\partial \phi_i} \log \det(C) = \operatorname{Tr}\left( \Lambda^{-1}(\phi) \frac{\partial \Lambda(\phi)}{\partial \phi_i} \right)
\tag{5.23}
\]

或者显式地写出矩阵乘法：

\[
\operatorname{Tr}\left[
\begin{pmatrix}
\lambda_1^{-1}(\phi) & & & \\
& \lambda_2^{-1}(\phi) & & \\
& & \ddots & \\
& & & \lambda_n^{-1}(\phi)
\end{pmatrix}
\begin{pmatrix}
\frac{\partial \lambda_1(\phi)}{\partial \phi_i} & & & \\
& \frac{\partial \lambda_2(\phi)}{\partial \phi_i} & & \\
& & \ddots & \\
& & & \frac{\partial \lambda_n(\phi)}{\partial \phi_i}
\end{pmatrix}
\right]
= \operatorname{Tr}\left( \Lambda^{-1}(\phi) \frac{\partial \Lambda(\phi)}{\partial \phi_i} \right)
\tag{5.24}
\]

以上是从特征值角度得到的初步结果——这个形式仍然与 \( \Lambda \) 有关，希望在最终结果中直接看到与 \( C \) 的关系。

利用 \( C = U^T \Lambda U \)，尝试将 (5.23) 改写为与 \( C \) 相关的形式。在两个矩阵之间插入 \( U^T U = I \)，不改变迹的值：

\[
\operatorname{Tr}\left( \Lambda^{-1} \frac{\partial \Lambda}{\partial \phi_i} \right)
= \operatorname{Tr}\left( \Lambda^{-1} \, U^T U \, \frac{\partial \Lambda}{\partial \phi_i} \, U^T U \right)
\tag{5.25}
\]

利用迹的循环性质 \( \operatorname{Tr}(ABC) = \operatorname{Tr}(BCA) = \operatorname{Tr}(CAB) \)，将矩阵循环排列：

\[
\operatorname{Tr}\left( \Lambda^{-1} U^T U \frac{\partial \Lambda}{\partial \phi_i} U^T U \right)
= \operatorname{Tr}\left( U \Lambda^{-1} U^T U \frac{\partial \Lambda}{\partial \phi_i} U^T \right)
\tag{5.26}
\]

这里 \( U \Lambda^{-1} U^T = C^{-1} \)。因此：

\[
\operatorname{Tr}\left( \Lambda^{-1} \frac{\partial \Lambda}{\partial \phi_i} \right)
= \operatorname{Tr}\left( C^{-1} \cdot U \frac{\partial \Lambda}{\partial \phi_i} U^T \right)
\tag{5.27}
\]

要把 \( U \frac{\partial \Lambda}{\partial \phi_i} U^T \) 替换成 \( \frac{\partial C}{\partial \phi_i} \)。问题是：

\[
U \frac{\partial \Lambda}{\partial \phi_i} U^T \xrightarrow{?} \frac{\partial}{\partial \phi_i} \left( U \Lambda U^T \right) = \frac{\partial C}{\partial \phi_i}
\tag{5.28}
\]

这里存在一个关键问题：\( U \) 是 \( C \) 的特征向量矩阵，而 \( C \) 依赖于参数 \( \phi \)，因此 \( U = U(\phi) \) 本身也是 \( \phi \) 的函数。当对 \( U \Lambda U^T \) 求导时，需要考虑 \( U \) 随 \( \phi \) 的变化：

\[
\frac{\partial C}{\partial \phi_i} = \frac{\partial}{\partial \phi_i} \left( U \Lambda U^T \right)
= \frac{\partial U}{\partial \phi_i} \Lambda U^T + U \frac{\partial \Lambda}{\partial \phi_i} U^T + U \Lambda \frac{\partial U^T}{\partial \phi_i}
\tag{5.29}
\]

(5.29) 包含三项，而所需要的只是中间那一项 \( U \frac{\partial \Lambda}{\partial \phi_i} U^T \)。多出来的两项是：

\[
\frac{\partial U}{\partial \phi_i} \Lambda U^T + U \Lambda \frac{\partial U^T}{\partial \phi_i}
\tag{5.30}
\]

问题是：这两项是否可以忽略？答案是否定的——它们不能忽略。但当取迹并乘以 \( C^{-1} \) 时，这两项会相互抵消。

需验证：

\[
\operatorname{Tr}\left( C^{-1} \frac{\partial C}{\partial \phi_i} \right)
= \operatorname{Tr}\left( \Lambda^{-1} \frac{\partial \Lambda}{\partial \phi_i} \right)
\tag{5.31}
\]

这个等式的左边是：

\[
\operatorname{Tr}\left( C^{-1} \frac{\partial C}{\partial \phi_i} \right)
= \operatorname{Tr}\left( C^{-1} \frac{\partial U}{\partial \phi_i} \Lambda U^T \right)
+ \operatorname{Tr}\left( C^{-1} U \frac{\partial \Lambda}{\partial \phi_i} U^T \right)
+ \operatorname{Tr}\left( C^{-1} U \Lambda \frac{\partial U^T}{\partial \phi_i} \right)
\tag{5.32}
\]

证明 (5.32) 等于 \( \operatorname{Tr}\left( \Lambda^{-1} \frac{\partial \Lambda}{\partial \phi_i} \right) \)，即第一项和第三项之和为零。

将 \( C^{-1} = U \Lambda^{-1} U^T \) 代入第一项和第三项。

第一项：

\[
\operatorname{Tr}\left( C^{-1} \frac{\partial U}{\partial \phi_i} \Lambda U^T \right)
= \operatorname{Tr}\left( U \Lambda^{-1} U^T \frac{\partial U}{\partial \phi_i} \Lambda U^T \right)
\tag{5.33}
\]

利用迹的循环性质，将 \( U \Lambda^{-1} \) 移到最右边：

\[
= \operatorname{Tr}\left( \Lambda^{-1} U^T \frac{\partial U}{\partial \phi_i} \Lambda U^T U \right)
= \operatorname{Tr}\left( \Lambda^{-1} U^T \frac{\partial U}{\partial \phi_i} \Lambda \right)
\tag{5.34}
\]

第三项：

\[
\operatorname{Tr}\left( C^{-1} U \Lambda \frac{\partial U^T}{\partial \phi_i} \right)
= \operatorname{Tr}\left( U \Lambda^{-1} U^T U \Lambda \frac{\partial U^T}{\partial \phi_i} \right)
= \operatorname{Tr}\left( U \Lambda^{-1} \Lambda \frac{\partial U^T}{\partial \phi_i} \right)
= \operatorname{Tr}\left( U \frac{\partial U^T}{\partial \phi_i} \right)
\tag{5.35}
\]

利用 \( U^T U = I \) 的导数关系。对 \( U^T U = I \) 两边求导：

\[
\frac{\partial}{\partial \phi_i} (U^T U) = \frac{\partial U^T}{\partial \phi_i} U + U^T \frac{\partial U}{\partial \phi_i} = 0
\tag{5.36}
\]

因此：

\[
\frac{\partial U^T}{\partial \phi_i} U = - U^T \frac{\partial U}{\partial \phi_i}
\tag{5.37}
\]

利用这个关系，重新审视 (5.34) 和 (5.35) 的和。

第一项 + 第三项：

\[
\operatorname{Tr}\left( \Lambda^{-1} U^T \frac{\partial U}{\partial \phi_i} \Lambda \right) + \operatorname{Tr}\left( U \frac{\partial U^T}{\partial \phi_i} \right)
\tag{5.38}
\]

利用迹的循环性质，将第一项的 \( \Lambda \) 移到最前面：

\[
\operatorname{Tr}\left( \Lambda^{-1} U^T \frac{\partial U}{\partial \phi_i} \Lambda \right)
= \operatorname{Tr}\left( \Lambda \Lambda^{-1} U^T \frac{\partial U}{\partial \phi_i} \right)
= \operatorname{Tr}\left( U^T \frac{\partial U}{\partial \phi_i} \right)
\tag{5.39}
\]

于是第一项 + 第三项变为：

\[
\operatorname{Tr}\left( U^T \frac{\partial U}{\partial \phi_i} \right) + \operatorname{Tr}\left( U \frac{\partial U^T}{\partial \phi_i} \right)
= \operatorname{Tr}\left( U^T \frac{\partial U}{\partial \phi_i} + U \frac{\partial U^T}{\partial \phi_i} \right)
= \operatorname{Tr}\left( \frac{\partial}{\partial \phi_i} (U^T U) \right)
\tag{5.40}
\]

由于 \( U^T U = I \)，而 \( I \) 与 \( \phi_i \) 无关，所以：

\[
\operatorname{Tr}\left( \frac{\partial}{\partial \phi_i} (U^T U) \right) = \operatorname{Tr}\left( \frac{\partial}{\partial \phi_i} I \right) = \operatorname{Tr}(0) = 0
\tag{5.41}
\]

因此，第一项和第三项之和为零。

于是，从 (5.32) 到 (5.41) 的完整推导证明了：

\[
\boxed{
\operatorname{Tr}\left( C^{-1} \frac{\partial C}{\partial \phi_i} \right) = \operatorname{Tr}\left( \Lambda^{-1} \frac{\partial \Lambda}{\partial \phi_i} \right) = \frac{\partial}{\partial \phi_i} \log \det(C)
}
\tag{5.42}
\]

关键在于：虽然 \( U(\phi) \) 与参数有关，不能简单地直接与求导交换顺序，但在取迹并乘以 \( C^{-1} \) 后，\( U \) 的导数项完全抵消了。因此，从特征值表示可以安全地回到与 \( C \) 直接相关的表示。

因此，最终恒等式为：

\[
\boxed{
\nabla_{\phi} \log \det(C(\phi)) = \operatorname{Tr}\left( C^{-1}(\phi) \frac{\partial C(\phi)}{\partial \phi} \right)
}
\tag{5.43}
\]

回到对数似然函数 (5.14)。由 (5.43)：

\[
\frac{\partial}{\partial \phi_i} \left( -\frac{1}{2} \log \det(C) \right) = -\frac{1}{2} \operatorname{Tr}\left( C^{-1} \frac{\partial C}{\partial \phi_i} \right)
\tag{5.44}
\]

这个结果将在下一节中用于构建 Fisher 信息矩阵。

### 2.4 计算二次型的导数

计算对数似然函数 (5.14) 的第三部分——二次型——对参数 \( \phi_i \) 的导数。这是 CRLB 推导中的关键计算步骤。

从 (5.14) 中提取第三部分：

\[
Q(\boldsymbol{\phi}) = -\frac{1}{2} (x - \boldsymbol{\mu}(\boldsymbol{\phi}))^T C(\boldsymbol{\phi})^{-1} (x - \boldsymbol{\mu}(\boldsymbol{\phi}))
\tag{5.45}
\]

的目标是计算：

\[
\frac{\partial}{\partial \phi_i} Q(\boldsymbol{\phi}) = -\frac{1}{2} \frac{\partial}{\partial \phi_i} \left[ (x - \boldsymbol{\mu}(\boldsymbol{\phi}))^T C(\boldsymbol{\phi})^{-1} (x - \boldsymbol{\mu}(\boldsymbol{\phi})) \right]
\tag{5.46}
\]

先关注括号内的导数上。设：

\[
R(\boldsymbol{\phi}) = (x - \boldsymbol{\mu}(\boldsymbol{\phi}))^T C(\boldsymbol{\phi})^{-1} (x - \boldsymbol{\mu}(\boldsymbol{\phi}))
\tag{5.47}
\]

\( R(\boldsymbol{\phi}) \) 是一个标量，它依赖于 \( \boldsymbol{\mu} \) 和 \( C \)，而 \( \boldsymbol{\mu} \) 和 \( C \) 都依赖于参数 \( \boldsymbol{\phi} \)。

由于 \( R(\boldsymbol{\phi}) \) 是三个因子的乘积——两个 \( (x - \boldsymbol{\mu}) \) 因子夹着中间的 \( C^{-1} \)——使用乘积法则 \( (uv)' = u'v + uv' \) 来计算其导数。但是这里有三个因子，需要仔细展开。

将 \( R(\boldsymbol{\phi}) \) 写为：

\[
R = (x - \boldsymbol{\mu})^T \, C^{-1} \, (x - \boldsymbol{\mu})
\tag{5.48}
\]

现在对 \( \phi_i \) 求偏导。根据乘积法则，求导会对三个因子分别作用：

- 对第一个 \( (x - \boldsymbol{\mu}) \) 求导，得到 \( -\frac{\partial \boldsymbol{\mu}}{\partial \phi_i} \)；
- 对中间的 \( C^{-1} \) 求导，得到 \( \frac{\partial C^{-1}}{\partial \phi_i} \)；
- 对第二个 \( (x - \boldsymbol{\mu}) \) 求导，得到 \( -\frac{\partial \boldsymbol{\mu}}{\partial \phi_i} \)。

展开：

\[
\frac{\partial R}{\partial \phi_i} = 
\left( \frac{\partial}{\partial \phi_i} (x - \boldsymbol{\mu})^T \right) C^{-1} (x - \boldsymbol{\mu})
+ (x - \boldsymbol{\mu})^T \left( \frac{\partial}{\partial \phi_i} C^{-1} \right) (x - \boldsymbol{\mu})
+ (x - \boldsymbol{\mu})^T C^{-1} \left( \frac{\partial}{\partial \phi_i} (x - \boldsymbol{\mu}) \right)
\tag{5.49}
\]

由于 \( \frac{\partial}{\partial \phi_i} (x - \boldsymbol{\mu}) = -\frac{\partial \boldsymbol{\mu}}{\partial \phi_i} \)，而 \( \left( \frac{\partial}{\partial \phi_i} (x - \boldsymbol{\mu})^T \right) = -\left( \frac{\partial \boldsymbol{\mu}}{\partial \phi_i} \right)^T \)，代入 (5.49) 得到：

\[
\frac{\partial R}{\partial \phi_i} = 
-\left( \frac{\partial \boldsymbol{\mu}}{\partial \phi_i} \right)^T C^{-1} (x - \boldsymbol{\mu})
+ (x - \boldsymbol{\mu})^T \left( \frac{\partial C^{-1}}{\partial \phi_i} \right) (x - \boldsymbol{\mu})
- (x - \boldsymbol{\mu})^T C^{-1} \left( \frac{\partial \boldsymbol{\mu}}{\partial \phi_i} \right)
\tag{5.50}
\]

观察第一项和第三项。第一项是一个标量，第三项也是一个标量。由于 \( \frac{\partial \boldsymbol{\mu}}{\partial \phi_i} \) 是列向量，\( (x - \boldsymbol{\mu}) \) 也是列向量，\( C^{-1} \) 是对称矩阵，第一项和第三项实际上是相等的（互为转置）：

第一项：\( -\left( \frac{\partial \boldsymbol{\mu}}{\partial \phi_i} \right)^T C^{-1} (x - \boldsymbol{\mu}) \)

第三项：\( -(x - \boldsymbol{\mu})^T C^{-1} \left( \frac{\partial \boldsymbol{\mu}}{\partial \phi_i} \right) \)

因为 \( C^{-1} \) 是对称矩阵（\( C^{-1} \) 与 \( C \) 有相同的特征向量，且 \( C \) 对称，所以 \( C^{-1} \) 也对称），所以：

\[
\left( \frac{\partial \boldsymbol{\mu}}{\partial \phi_i} \right)^T C^{-1} (x - \boldsymbol{\mu})
= (x - \boldsymbol{\mu})^T C^{-1} \left( \frac{\partial \boldsymbol{\mu}}{\partial \phi_i} \right)
\tag{5.51}
\]

因此，第一项和第三项可以合并为：

\[
-2 \left( \frac{\partial \boldsymbol{\mu}}{\partial \phi_i} \right)^T C^{-1} (x - \boldsymbol{\mu})
\tag{5.52}
\]

于是 (5.50) 化简为：

\[
\frac{\partial R}{\partial \phi_i} = 
-2 \left( \frac{\partial \boldsymbol{\mu}}{\partial \phi_i} \right)^T C^{-1} (x - \boldsymbol{\mu})
+ (x - \boldsymbol{\mu})^T \left( \frac{\partial C^{-1}}{\partial \phi_i} \right) (x - \boldsymbol{\mu})
\tag{5.53}
\]

计算 \( \frac{\partial C^{-1}}{\partial \phi_i} \)。

**计算 \( \frac{\partial C^{-1}}{\partial \phi_i} \)**：利用恒等式 \( C^{-1} C = I \)。两边对 \( \phi_i \) 求导：

\[
\frac{\partial}{\partial \phi_i} (C^{-1} C) = \frac{\partial}{\partial \phi_i} I = 0
\tag{5.54}
\]

左边用乘积法则展开：

\[
\frac{\partial C^{-1}}{\partial \phi_i} \cdot C + C^{-1} \cdot \frac{\partial C}{\partial \phi_i} = 0
\tag{5.55}
\]

将第一项移到等号右边：

\[
\frac{\partial C^{-1}}{\partial \phi_i} \cdot C = - C^{-1} \frac{\partial C}{\partial \phi_i}
\tag{5.56}
\]

两边右乘 \( C^{-1} \)：

\[
\boxed{
\frac{\partial C^{-1}}{\partial \phi_i} = - C^{-1} \frac{\partial C}{\partial \phi_i} C^{-1}
}
\tag{5.57}
\]

将 (5.57) 代入 (5.53)：

\[
\frac{\partial R}{\partial \phi_i} = 
-2 \left( \frac{\partial \boldsymbol{\mu}}{\partial \phi_i} \right)^T C^{-1} (x - \boldsymbol{\mu})
- (x - \boldsymbol{\mu})^T C^{-1} \left( \frac{\partial C}{\partial \phi_i} \right) C^{-1} (x - \boldsymbol{\mu})
\tag{5.58}
\]

---

为了简化符号，引入更紧凑的记法，因为后续的推导会越来越复杂，这样写起来更清楚一些。

令：

\[
(\boldsymbol{\mu})_i' \triangleq \frac{\partial \boldsymbol{\mu}}{\partial \phi_i}, \qquad (C)_i' \triangleq \frac{\partial C}{\partial \phi_i}
\tag{5.59}
\]

这个记法表示：下标 \( i \) 表示对第 \( i \) 个参数 \( \phi_i \) 求偏导，撇号 \( ' \) 表示这是导数。

在这些记法下，(5.58) 可以写为：

\[
\frac{\partial R}{\partial \phi_i} = 
-2 \left( (\boldsymbol{\mu})_i' \right)^T C^{-1} (x - \boldsymbol{\mu})
- (x - \boldsymbol{\mu})^T C^{-1} (C)_i' C^{-1} (x - \boldsymbol{\mu})
\tag{5.60}
\]

回到对数似然函数的导数 (5.46)。目标是 \( -\frac{1}{2} R \)，因此：

\[
\frac{\partial \ell}{\partial \phi_i} = -\frac{1}{2} \frac{\partial R}{\partial \phi_i}
= -\frac{1}{2} \left[ -2 \left( (\boldsymbol{\mu})_i' \right)^T C^{-1} (x - \boldsymbol{\mu}) - (x - \boldsymbol{\mu})^T C^{-1} (C)_i' C^{-1} (x - \boldsymbol{\mu}) \right]
\tag{5.61}
\]

将 \( -\frac{1}{2} \) 分配到每一项：

第一项：\( -\frac{1}{2} \times (-2) = +1 \)，所以第一项变为：

\[
\left( (\boldsymbol{\mu})_i' \right)^T C^{-1} (x - \boldsymbol{\mu})
\tag{5.62}
\]

第二项：\( -\frac{1}{2} \times (-1) = +\frac{1}{2} \)，所以第二项变为：

\[
\frac{1}{2} (x - \boldsymbol{\mu})^T C^{-1} (C)_i' C^{-1} (x - \boldsymbol{\mu})
\tag{5.63}
\]

因此：

\[
\boxed{
\frac{\partial \ell}{\partial \phi_i} = 
\left( (\boldsymbol{\mu})_i' \right)^T C^{-1} (x - \boldsymbol{\mu})
+ \frac{1}{2} (x - \boldsymbol{\mu})^T C^{-1} (C)_i' C^{-1} (x - \boldsymbol{\mu})
}
\tag{5.64}
\]

---

### 2.5 合并两项：\( \ell(\boldsymbol{\phi}) \) 的完整导数

将 2.3 节和 2.4 节得到的两个结果合并，得到对数似然函数的完整导数。

从 (5.14) 出发：

\[
\ell(\boldsymbol{\phi}) = -\frac{n}{2} \log(2\pi) - \frac{1}{2} \log \det(C) - \frac{1}{2} (x - \boldsymbol{\mu})^T C^{-1} (x - \boldsymbol{\mu})
\tag{5.65}
\]

对 \( \phi_i \) 求偏导：

- 第一项 \( -\frac{n}{2} \log(2\pi) \) 是常数，导数为 0。
- 第二项 \( -\frac{1}{2} \log \det(C) \) 的导数，由 (5.44) 得：\( -\frac{1}{2} \operatorname{Tr}\left( C^{-1} (C)_i' \right) \)。
- 第三项 \( -\frac{1}{2} (x - \boldsymbol{\mu})^T C^{-1} (x - \boldsymbol{\mu}) \) 的导数，由 (5.61)-(5.63) 得：\( \left( (\boldsymbol{\mu})_i' \right)^T C^{-1} (x - \boldsymbol{\mu}) + \frac{1}{2} (x - \boldsymbol{\mu})^T C^{-1} (C)_i' C^{-1} (x - \boldsymbol{\mu}) \)。

将这三项相加：

\[
\boxed{
\frac{\partial \ell}{\partial \phi_i} = 
\left( (\boldsymbol{\mu})_i' \right)^T C^{-1} (x - \boldsymbol{\mu})
-\frac{1}{2} \operatorname{Tr}\left( C^{-1} (C)_i' \right)
+ \frac{1}{2} (x - \boldsymbol{\mu})^T C^{-1} (C)_i' C^{-1} (x - \boldsymbol{\mu})
}
\tag{5.66}
\]

或者等价地写成：

\[
\boxed{
\ell(\boldsymbol{\phi})' = 
\left( (\boldsymbol{\mu})_i' \right)^T C^{-1} (x - \boldsymbol{\mu})
-\frac{1}{2} \operatorname{Tr}\left( C^{-1} (C)_i' \right)
+ \frac{1}{2} (x - \boldsymbol{\mu})^T C^{-1} (C)_i' C^{-1} (x - \boldsymbol{\mu})
}
\tag{5.67}
\]

(5.66) 或 (5.67) 是 CRLB 推导中的关键中间结果。它给出了任意多维高斯分布的对数似然函数对任意参数 \( \phi_i \) 的导数的通用表达式。

这个表达式由三项组成：

1. **均值项**：\( \left( (\boldsymbol{\mu})_i' \right)^T C^{-1} (x - \boldsymbol{\mu}) \)，它来自均值 \( \boldsymbol{\mu} \) 对参数的依赖。这一项是线性的，因为 \( (x - \boldsymbol{\mu}) \) 的指数是 1。

2. **协方差矩阵的迹项**：\( -\frac{1}{2} \operatorname{Tr}\left( C^{-1} (C)_i' \right) \)，它来自协方差矩阵的行列式 \( \log \det(C) \) 的导数。这一项与数据 \( x \) 无关，只与参数 \( \phi \) 有关，它是 Fisher 信息中关于协方差矩阵参数的贡献。

3. **协方差矩阵的二次型项**：\( \frac{1}{2} (x - \boldsymbol{\mu})^T C^{-1} (C)_i' C^{-1} (x - \boldsymbol{\mu}) \)，它也来自协方差矩阵 \( C \) 对参数的依赖。这一项是二次型的，包含了数据的随机波动。它的期望（在 Fisher 信息矩阵的计算中）将与第二项相抵消。

### 2.6 计算二阶导数（Hessian 矩阵）

为计算 Fisher 信息矩阵，需要对一阶导数 (5.66) 再次求导，得到 Hessian 矩阵的元素。记 Hessian 矩阵的第 \( (i,j) \) 个元素为：

\[
\mathbf{H}_{ij} = \frac{\partial^2 \ell}{\partial \phi_i \partial \phi_j}
\tag{5.68}
\]

从 (5.66) 出发：

\[
\frac{\partial \ell}{\partial \phi_i} = 
\left( (\boldsymbol{\mu})_i' \right)^T C^{-1} (x - \boldsymbol{\mu})
-\frac{1}{2} \operatorname{Tr}\left( C^{-1} (C)_i' \right)
+ \frac{1}{2} (x - \boldsymbol{\mu})^T C^{-1} (C)_i' C^{-1} (x - \boldsymbol{\mu})
\tag{5.69}
\]

对 \( \phi_j \) 再次求偏导，得到一个很长的表达式：

\[
\frac{\partial^2 \ell}{\partial \phi_i \partial \phi_j}
=
\left(
\left( (\boldsymbol{\mu})_i' \right)^T C^{-1} (x - \boldsymbol{\mu}) 
-\frac{1}{2} \operatorname{Tr}\left( C^{-1} (C)_i' \right)
+ \frac{1}{2} (x - \boldsymbol{\mu})^T C^{-1} (C)_i' C^{-1} (x - \boldsymbol{\mu})
\right)'
\tag{5.70}
\]

逐项展开。为简化取期望操作（\( \mathbb{E}[x - \boldsymbol{\mu}] = 0 \)），在展开过程中把项分为两类：

- **不含 \( (x - \boldsymbol{\mu}) \) 的项**：这些项在取期望时保持不变；
- **含 \( (x - \boldsymbol{\mu}) \) 的项**：这些项在取期望时会消失，因为 \( \mathbb{E}[x - \boldsymbol{\mu}] = 0 \) 以及 \( \mathbb{E}[(x - \boldsymbol{\mu})(x - \boldsymbol{\mu})^T] = C \)，它们会贡献到 Fisher 信息矩阵中。

展开过程非常繁长，按顺序逐项处理。

**第一项的导数：**

第一项为：

\[
T_1 = \left( (\boldsymbol{\mu})_i' \right)^T C^{-1} (x - \boldsymbol{\mu})
\tag{5.71}
\]

对 \( \phi_j \) 求导（用乘积法则）：

\[
\frac{\partial T_1}{\partial \phi_j}
=
\left( (\boldsymbol{\mu})_{ij}'' \right)^T C^{-1} (x - \boldsymbol{\mu})
+ \left( (\boldsymbol{\mu})_i' \right)^T \left( -C^{-1} (C)_j' C^{-1} \right) (x - \boldsymbol{\mu})
- \left( (\boldsymbol{\mu})_i' \right)^T C^{-1} (\boldsymbol{\mu})_j'
\tag{5.72}
\]

这里 \( (\boldsymbol{\mu})_{ij}'' \) 表示 \( \boldsymbol{\mu} \) 对 \( \phi_i \) 和 \( \phi_j \) 的二阶混合偏导。前两项含有 \( (x - \boldsymbol{\mu}) \)，第三项不含。

**第二项的导数：**

第二项为：

\[
T_2 = -\frac{1}{2} \operatorname{Tr}\left( C^{-1} (C)_i' \right)
\tag{5.73}
\]

对 \( \phi_j \) 求导：

\[
\frac{\partial T_2}{\partial \phi_j}
=
-\frac{1}{2} \operatorname{Tr}\left( \frac{\partial}{\partial \phi_j} \left( C^{-1} (C)_i' \right) \right)
\tag{5.74}
\]

利用 \( \frac{\partial C^{-1}}{\partial \phi_j} = -C^{-1} (C)_j' C^{-1} \)，展开：

\[
\frac{\partial T_2}{\partial \phi_j}
=
-\frac{1}{2} \operatorname{Tr}\left( -C^{-1} (C)_j' C^{-1} (C)_i' + C^{-1} (C)_{ij}'' \right)
\tag{5.75}
\]

化简：

\[
\frac{\partial T_2}{\partial \phi_j}
=
\frac{1}{2} \operatorname{Tr}\left( C^{-1} (C)_j' C^{-1} (C)_i' \right)
- \frac{1}{2} \operatorname{Tr}\left( C^{-1} (C)_{ij}'' \right)
\tag{5.76}
\]

**第三项的导数：**

第三项为：

\[
T_3 = \frac{1}{2} (x - \boldsymbol{\mu})^T C^{-1} (C)_i' C^{-1} (x - \boldsymbol{\mu})
\tag{5.77}
\]

对 \( \phi_j \) 求导。令 \( A = C^{-1} (C)_i' C^{-1} \)，则 \( T_3 = \frac{1}{2} (x - \boldsymbol{\mu})^T A (x - \boldsymbol{\mu}) \)。对 \( \phi_j \) 求导，需要同时对 \( (x - \boldsymbol{\mu}) \)、\( A \) 和 \( (x - \boldsymbol{\mu}) \) 三个因子求导。

\[
\frac{\partial T_3}{\partial \phi_j}
=
-\left( (\boldsymbol{\mu})_j' \right)^T A (x - \boldsymbol{\mu})
+ \frac{1}{2} (x - \boldsymbol{\mu})^T \left( \frac{\partial A}{\partial \phi_j} \right) (x - \boldsymbol{\mu})
- \frac{1}{2} (x - \boldsymbol{\mu})^T A (\boldsymbol{\mu})_j'
\tag{5.78}
\]

其中 \( \frac{\partial A}{\partial \phi_j} = \frac{\partial}{\partial \phi_j} \left( C^{-1} (C)_i' C^{-1} \right) \)。展开 \( \frac{\partial A}{\partial \phi_j} \)：

\[
\frac{\partial A}{\partial \phi_j}
=
-C^{-1} (C)_j' C^{-1} (C)_i' C^{-1}
+ C^{-1} (C)_{ij}'' C^{-1}
- C^{-1} (C)_i' C^{-1} (C)_j' C^{-1}
\tag{5.79}
\]

将 (5.79) 代入 (5.78)，得到：

\[
\frac{\partial T_3}{\partial \phi_j}
=
-\left( (\boldsymbol{\mu})_j' \right)^T C^{-1} (C)_i' C^{-1} (x - \boldsymbol{\mu})
- \frac{1}{2} (x - \boldsymbol{\mu})^T C^{-1} (C)_j' C^{-1} (C)_i' C^{-1} (x - \boldsymbol{\mu})
+ \frac{1}{2} (x - \boldsymbol{\mu})^T C^{-1} (C)_{ij}'' C^{-1} (x - \boldsymbol{\mu})
- \frac{1}{2} (x - \boldsymbol{\mu})^T C^{-1} (C)_i' C^{-1} (C)_j' C^{-1} (x - \boldsymbol{\mu})
- \frac{1}{2} (x - \boldsymbol{\mu})^T A (\boldsymbol{\mu})_j'
\tag{5.80}
\]

注意 (5.72) 中也有一个类似的项 \( -\left( (\boldsymbol{\mu})_i' \right)^T C^{-1} (\boldsymbol{\mu})_j' \)。把这些所有项合并，表达式将非常冗长。

---

**Fisher 信息矩阵可以用 \( -\mathbb{E}[\mathbf{H}] \) 或 \( \mathbb{E}[\mathbf{s} \mathbf{s}^T] \) 计算，其中 \( \mathbf{s} \) 是得分向量。第二种方法通常更简单，因为它避免了二阶导数的复杂展开，而且 (5.66) 已给出得分向量的表达式。**

然而，为实现简洁表达，从得分向量的二阶矩出发。

从 (5.66) 开始，得分向量的第 \( i \) 个分量为：

\[
s_i = \frac{\partial \ell}{\partial \phi_i}
= \left( (\boldsymbol{\mu})_i' \right)^T C^{-1} (x - \boldsymbol{\mu})
-\frac{1}{2} \operatorname{Tr}\left( C^{-1} (C)_i' \right)
+ \frac{1}{2} (x - \boldsymbol{\mu})^T C^{-1} (C)_i' C^{-1} (x - \boldsymbol{\mu})
\tag{5.81}
\]

Fisher 信息矩阵的元素为：

\[
[\mathbf{I}(\boldsymbol{\phi})]_{ij} = \mathbb{E}[s_i s_j]
\tag{5.82}
\]

其中期望是对数据 \( x \) 取的。在计算这个乘积的期望时，用到以下标准结果（针对零均值高斯向量 \( z = x - \boldsymbol{\mu} \)）：

\[
\mathbb{E}[z] = 0, \quad \mathbb{E}[z z^T] = C
\tag{5.83}
\]

以及对于高斯随机向量的四阶矩（Isserlis 定理 / Wick 定理）：

\[
\mathbb{E}[z_i z_j z_k z_l] = C_{ij} C_{kl} + C_{ik} C_{jl} + C_{il} C_{jk}
\tag{5.84}
\]

利用这些结果，可以直接计算 (5.82)。

注意，在 (5.81) 中，得分向量由三项组成：
- 第一项是线性的；
- 第二项是常数（与 \( z \) 无关）；
- 第三项是二次型的。

计算 \( \mathbb{E}[s_i s_j] \) 时，发现：
- 线性和常数项的乘积在取期望时为零（因为 \( \mathbb{E}[z] = 0 \)）；
- 线性项与线性项的乘积给出 \( \left( (\boldsymbol{\mu})_i' \right)^T C^{-1} (\boldsymbol{\mu})_j' \)；
- 常数项与常数项的乘积在取期望时与协方差矩阵的二阶导数项相结合，最终与二次型项的期望中的相应部分相抵消；
- 二次型项与常数项的乘积也涉及 \( \mathbb{E}[z z^T] \)，这些项会与常数项的导数项抵消。

经过这些计算（推导过程较长但每一步都是标准的），Fisher 信息矩阵的最终表达式为：

\[
\boxed{
[\mathbf{I}(\boldsymbol{\phi})]_{ij} =
\left( (\boldsymbol{\mu})_i' \right)^T C^{-1} (\boldsymbol{\mu})_j'
+ \frac{1}{2} \operatorname{Tr}\left( (C)_i' \, C^{-1} \, (C)_j' \, C^{-1} \right)
}
\tag{5.85}
\]

其中 \( (\boldsymbol{\mu})_i' = \frac{\partial \boldsymbol{\mu}}{\partial \phi_i} \)，\( (C)_i' = \frac{\partial C}{\partial \phi_i} \)。

---
### 2.7 Fisher 信息矩阵的期望计算：从得分函数到 Slepian-Bangs 公式

2.6 节中已得到得分函数（一阶导数）的表达式，下面严格计算 Fisher 信息矩阵：

\[
[\mathbf{I}(\boldsymbol{\phi})]_{ij} = \mathbb{E}\left[ \frac{\partial \ell}{\partial \phi_i} \frac{\partial \ell}{\partial \phi_j} \right]
\tag{5.86}
\]

这个期望是对数据 \( x \) 取的，在确定性信号模型下，数据的随机性完全来源于 \( z = x - \boldsymbol{\mu} \)，其中 \( z \sim \mathcal{N}(0, C) \)。

为了计算这个期望，采用**得分向量外积法**，而不是 Hessian 的负期望。这能让更清晰地看到各项如何贡献于最终结果。

---

#### 2.7.1 将得分函数分解为线性部分和中心化二次型

从 (5.66) 出发，得分函数的第 \( i \) 个分量为：

\[
s_i = \frac{\partial \ell}{\partial \phi_i} = 
\underbrace{\left( (\boldsymbol{\mu})_i' \right)^T C^{-1} z}_{L_i}
-\frac{1}{2} \operatorname{Tr}\left( C^{-1} (C)_i' \right)
+ \frac{1}{2} z^T C^{-1} (C)_i' C^{-1} z
\tag{5.87}
\]

其中 \( z = x - \boldsymbol{\mu} \)，\( \mathbb{E}[z] = 0 \)，\( \mathbb{E}[z z^T] = C \)。

先计算二次型 \( z^T C^{-1} (C)_i' C^{-1} z \) 的期望。令 \( A_i = C^{-1} (C)_i' C^{-1} \)，则

\[
\mathbb{E}\left[ z^T A_i z \right] = \operatorname{Tr}(A_i C) = \operatorname{Tr}\left( C^{-1} (C)_i' \right)
\tag{5.88}
\]

这个结果在推导中心化二次型时将会用到。

将 \( s_i \) 重新组合，定义一个新的变量 \( U_i \)，使得 \( s_i \) 分解为一个线性项 \( L_i \) 和一个**零均值的中心化二次型** \( U_i \)：

\[
s_i = L_i + U_i
\tag{5.89}
\]

其中：

\[
L_i = \left( (\boldsymbol{\mu})_i' \right)^T C^{-1} z
\tag{5.90}
\]

\[
U_i = \frac{1}{2} z^T C^{-1} (C)_i' C^{-1} z - \frac{1}{2} \operatorname{Tr}\left( C^{-1} (C)_i' \right)
= \frac{1}{2} \left[ z^T C^{-1} (C)_i' C^{-1} z - \operatorname{Tr}\left( C^{-1} (C)_i' \right) \right]
\tag{5.91}
\]

注意 \( U_i \) 的定义：它是在二次型 \( z^T C^{-1} (C)_i' C^{-1} z \) 的基础上，减去了它的期望 \( \operatorname{Tr}(C^{-1} (C)_i') \)。因此：

\[
\boxed{\mathbb{E}[U_i] = 0}
\tag{5.92}
\]

这一构造消除了 \( s_i \) 中的常数项，使得 \( s_i \) 成为零均值项 \( L_i \) 和零均值项 \( U_i \) 的和。

---

#### 2.7.2 交叉项 \( \mathbb{E}[L_i U_j] \) 的计算

Fisher 信息矩阵的元素为：

\[
[\mathbf{I}]_{ij} = \mathbb{E}[s_i s_j] = \mathbb{E}[(L_i + U_i)(L_j + U_j)]
\tag{5.93}
\]

展开：

\[
\mathbb{E}[s_i s_j] = \mathbb{E}[L_i L_j] + \mathbb{E}[L_i U_j] + \mathbb{E}[U_i L_j] + \mathbb{E}[U_i U_j]
\tag{5.94}
\]

这里 \( L_i \) 是 \( z \) 的线性函数，而 \( U_j \) 是 \( z \) 的中心化二次型。由于高斯分布的奇数阶矩为零（特别地，线性项与二次型项的乘积的期望为零）：

\[
\mathbb{E}[L_i U_j] = 0, \qquad \mathbb{E}[U_i L_j] = 0
\tag{5.95}
\]

因此，(5.94) 简化为：

\[
\boxed{[\mathbf{I}]_{ij} = \mathbb{E}[L_i L_j] + \mathbb{E}[U_i U_j]}
\tag{5.96}
\]

---

#### 2.7.3 计算 \( \mathbb{E}[L_i L_j] \)

由 (5.90)：

\[
L_i = \left( (\boldsymbol{\mu})_i' \right)^T C^{-1} z
\tag{5.97}
\]

根据二次型的期望公式 \( \mathbb{E}[(a^T z)(b^T z)] = a^T \mathbb{E}[z z^T] b = a^T C b \)，有：

\[
\mathbb{E}[L_i L_j] = \left( (\boldsymbol{\mu})_i' \right)^T C^{-1} \cdot C \cdot C^{-1} (\boldsymbol{\mu})_j'
= \left( (\boldsymbol{\mu})_i' \right)^T C^{-1} (\boldsymbol{\mu})_j'
\tag{5.98}
\]

---

#### 2.7.4 计算 \( \mathbb{E}[U_i U_j] \)

由 (5.91)：

\[
U_i = \frac{1}{2} \left[ z^T A_i z - \operatorname{Tr}(A_i C) \right], \qquad A_i = C^{-1} (C)_i' C^{-1}
\tag{5.99}
\]

注意 \( \operatorname{Tr}(A_i C) = \operatorname{Tr}(C^{-1} (C)_i') \)，即 (5.88) 中的结果。

下面计算两个中心化二次型乘积的期望。经典结果如下：对于零均值高斯随机向量 \( z \sim \mathcal{N}(0, C) \)，若 \( A \) 和 \( B \) 是对称矩阵，则

\[
\mathbb{E}\left[ \left( z^T A z - \operatorname{Tr}(A C) \right) \left( z^T B z - \operatorname{Tr}(B C) \right) \right] = 2 \operatorname{Tr}(A C B C)
\tag{5.100}
\]

这个公式来自 Isserlis 定理（Wick 定理）的四阶矩展开，它直接给出了两个中心化二次型之间的协方差。

在 \( U_i U_j \) 中：

\[
A = A_i = C^{-1} (C)_i' C^{-1}, \qquad B = A_j = C^{-1} (C)_j' C^{-1}
\tag{5.101}
\]

代入 (5.100)：

\[
\mathbb{E}[U_i U_j]
= \frac{1}{4} \mathbb{E}\left[ \left( z^T A_i z - \operatorname{Tr}(A_i C) \right) \left( z^T A_j z - \operatorname{Tr}(A_j C) \right) \right]
= \frac{1}{4} \cdot 2 \operatorname{Tr}(A_i C A_j C)
\tag{5.102}
\]

化简系数：

\[
\mathbb{E}[U_i U_j] = \frac{1}{2} \operatorname{Tr}(A_i C A_j C)
\tag{5.103}
\]

现在代入 \( A_i = C^{-1} (C)_i' C^{-1} \)：

\[
A_i C = C^{-1} (C)_i'
\tag{5.104}
\]

\[
A_i C A_j C = \left( C^{-1} (C)_i' \right) \left( C^{-1} (C)_j' \right) = C^{-1} (C)_i' C^{-1} (C)_j'
\tag{5.105}
\]

因此：

\[
\boxed{\mathbb{E}[U_i U_j] = \frac{1}{2} \operatorname{Tr}\left( C^{-1} (C)_i' C^{-1} (C)_j' \right)}
\tag{5.106}
\]

---

#### 2.7.5 合并得到 Fisher 信息矩阵的最终形式

将 (5.98) 和 (5.106) 代入 (5.96)：

\[
\boxed{
[\mathbf{I}(\boldsymbol{\phi})]_{ij} =
\left( (\boldsymbol{\mu})_i' \right)^T C^{-1} (\boldsymbol{\mu})_j'
+ \frac{1}{2} \operatorname{Tr}\left( (C)_i' \, C^{-1} \, (C)_j' \, C^{-1} \right)
}
\tag{5.107}
\]

即 **Slepian-Bangs 公式**，它给出了多维高斯分布参数估计的 Fisher 信息矩阵的标准形式。

---

#### 2.7.6 推导中的关键要点

1. **中心化二次型的构造**：在 (5.91) 中，通过减去 \( \frac{1}{2} \operatorname{Tr}(C^{-1} (C)_i') \)，将 \( s_i \) 分解为线性项 \( L_i \) 和零均值项 \( U_i \) 的和。这一步消除了常数项，简化了后续的期望计算。

2. **交叉项的消失**：线性项与中心化二次型的乘积的期望为零（\( \mathbb{E}[L_i U_j] = 0 \)），因为高斯分布的奇数阶矩为零。这使得 Fisher 信息矩阵简化为线性项和二次型项的贡献之和。

3. **中心化二次型的协方差公式**：(5.100) 是高斯分布最重要的矩公式之一。它直接来自四阶矩的 Isserlis 定理。如果不用这个公式，直接展开 \( \mathbb{E}[z_i z_j z_k z_l] \) 并利用 \( \mathbb{E}[z_i z_j z_k z_l] = C_{ij} C_{kl} + C_{ik} C_{jl} + C_{il} C_{jk} \)，也可以得到同样的结果，但 (5.100) 更简洁。

4. **最终形式**：(5.107) 虽然推导过程较长，但最终表达式非常简洁：它由"均值贡献"和"协方差贡献"两部分组成。这一形式是后续计算 DOA 估计 CRLB 的出发点。

---
### 2.8 DOA 估计问题的 Fisher 信息矩阵

将一般性的 Slepian-Bangs 公式 (5.107) 具体应用到 DOA 估计问题中。

对于 DOA 估计，有：

\[
\boldsymbol{\mu} = A(\theta) S, \qquad C = \sigma^2 I
\tag{5.108}
\]

其中 \( A(\theta) \in \mathbb{R}^{N \times M} \) 是方向矩阵，\( S \in \mathbb{R}^{M \times 1} \) 是信号向量（为简化推导，先考虑单快拍情况，多快拍情况可以在后面通过 Kronecker 积扩展）。

完整的参数向量为：

\[
\boldsymbol{\phi} = (\theta_1, \cdots, \theta_M, \; s_1, \cdots, s_M, \; \sigma^2)^T
\tag{5.109}
\]

其中前 \( M \) 个是角度参数（关心的），中间 \( M \) 个是信号参数（讨厌参数），最后一个是噪声方差参数（讨厌参数）。

---

#### 2.8.1 Fisher 信息矩阵的四块结构

由于 \( C = \sigma^2 I \) 只依赖于 \( \sigma^2 \)，而 \( \boldsymbol{\mu} = A(\theta) S \) 只依赖于 \( \theta \) 和 \( S \)，Slepian-Bangs 公式 (5.107) 中的两项贡献如下：

**第一项（来自均值的贡献）：**

\[
[\mathbf{I}_{\mu\mu}]_{ij} = \left( \frac{\partial \boldsymbol{\mu}}{\partial \phi_i} \right)^T C^{-1} \left( \frac{\partial \boldsymbol{\mu}}{\partial \phi_j} \right)
\tag{5.110}
\]

**第二项（来自协方差矩阵的贡献）：**

由于 \( C = \sigma^2 I \)，只有当参数是 \( \sigma^2 \) 时，\( \frac{\partial C}{\partial \phi_i} \neq 0 \)。当参数是角度 \( \theta \) 或信号 \( S \) 时，\( \frac{\partial C}{\partial \phi_i} = 0 \)。

\[
[\mathbf{I}_{CC}]_{ij} = \frac{1}{2} \operatorname{Tr}\left( \frac{\partial C}{\partial \phi_i} C^{-1} \frac{\partial C}{\partial \phi_j} C^{-1} \right)
\tag{5.111}
\]

当 \( \phi_i = \phi_j = \sigma^2 \) 时：

\[
\frac{\partial C}{\partial \sigma^2} = I, \qquad C^{-1} = \frac{1}{\sigma^2} I
\tag{5.112}
\]

代入 (5.112)：

\[
[\mathbf{I}_{CC}]_{\sigma^2\sigma^2} = \frac{1}{2} \operatorname{Tr}\left( I \cdot \frac{1}{\sigma^2} I \cdot I \cdot \frac{1}{\sigma^2} I \right)
= \frac{1}{2} \operatorname{Tr}\left( \frac{1}{\sigma^4} I \right)
= \frac{1}{2} \cdot \frac{N}{\sigma^4}
= \frac{N}{2\sigma^4}
\tag{5.113}
\]

对于所有其他情况（即参数不是 \( \sigma^2 \)），该项为零。

因此，Fisher 信息矩阵呈现分块结构：

\[
\mathbf{I}(\boldsymbol{\phi}) = 
\begin{pmatrix}
\mathbf{I}_{\theta\theta} & \mathbf{I}_{\theta S} & \mathbf{0} \\
\mathbf{I}_{S\theta} & \mathbf{I}_{SS} & \mathbf{0} \\
\mathbf{0} & \mathbf{0} & \frac{N}{2\sigma^4}
\end{pmatrix}
\tag{5.114}
\]

其中各块定义如下：

\[
\mathbf{I}_{\theta\theta} = \frac{1}{\sigma^2} \left( \frac{\partial (A(\theta)S)}{\partial \theta} \right)^T \frac{\partial (A(\theta)S)}{\partial \theta}
\in \mathbb{R}^{M \times M}
\tag{5.115}
\]

\[
\mathbf{I}_{\theta S} = \frac{1}{\sigma^2} \left( \frac{\partial (A(\theta)S)}{\partial \theta} \right)^T \frac{\partial (A(\theta)S)}{\partial S}
\in \mathbb{R}^{M \times M}
\tag{5.116}
\]

\[
\mathbf{I}_{S\theta} = \mathbf{I}_{\theta S}^T
\tag{5.117}
\]

\[
\mathbf{I}_{SS} = \frac{1}{\sigma^2} \left( \frac{\partial (A(\theta)S)}{\partial S} \right)^T \frac{\partial (A(\theta)S)}{\partial S}
\in \mathbb{R}^{M \times M}
\tag{5.118}
\]

---

#### 2.8.2 用雅可比矩阵表达

为简化记号，定义以下雅可比矩阵：

\[
F \triangleq \frac{\partial}{\partial \theta} (A(\theta)S) \in \mathbb{R}^{N \times M}
\tag{5.119}
\]

\[
G \triangleq \frac{\partial}{\partial S} (A(\theta)S) \in \mathbb{R}^{N \times M}
\tag{5.120}
\]

注意 \( G = A(\theta) \)，因为 \( A(\theta)S \) 对 \( S \) 的偏导数就是 \( A(\theta) \) 本身。而 \( F \) 的每一列是 \( A(\theta)S \) 对 \( \theta_i \) 的偏导数：

\[
F_{:,i} = \frac{\partial A(\theta)}{\partial \theta_i} S
\tag{5.121}
\]

则 Fisher 信息矩阵的分块可以简洁地写为：

\[
\mathbf{I}(\boldsymbol{\phi}) = 
\begin{pmatrix}
\frac{1}{\sigma^2} F^T F & \frac{1}{\sigma^2} F^T G & \mathbf{0} \\
\frac{1}{\sigma^2} G^T F & \frac{1}{\sigma^2} G^T G & \mathbf{0} \\
\mathbf{0} & \mathbf{0} & \frac{N}{2\sigma^4}
\end{pmatrix}
\tag{5.122}
\]

将 \( G = A(\theta) \) 代入：

\[
\boxed{
\mathbf{I}(\boldsymbol{\phi}) = 
\begin{pmatrix}
\frac{1}{\sigma^2} F^T F & \frac{1}{\sigma^2} F^T A & \mathbf{0} \\
\frac{1}{\sigma^2} A^T F & \frac{1}{\sigma^2} A^T A & \mathbf{0} \\
\mathbf{0} & \mathbf{0} & \frac{N}{2\sigma^4}
\end{pmatrix}
}
\tag{5.123}
\]

---
#### 2.8.3 角度参数的 CRLB：用分块矩阵求逆公式推导

在 (5.123) 中，得到了 Fisher 信息矩阵的分块结构：

\[
\mathbf{I}(\boldsymbol{\phi}) = 
\begin{pmatrix}
\mathbf{I}_{\theta\theta} & \mathbf{I}_{\theta S} & \mathbf{0} \\
\mathbf{I}_{S\theta} & \mathbf{I}_{SS} & \mathbf{0} \\
\mathbf{0} & \mathbf{0} & d
\end{pmatrix}
\tag{5.124}
\]

由于第三块与 $\theta$ 和 $S$ 没有耦合，只需要对左上角的 2×2 分块求逆：

\[
\mathbf{I}_{\text{red}} = 
\begin{pmatrix}
A & B \\
C & D
\end{pmatrix}
=
\begin{pmatrix}
\mathbf{I}_{\theta\theta} & \mathbf{I}_{\theta S} \\
\mathbf{I}_{S\theta} & \mathbf{I}_{SS}
\end{pmatrix}
\tag{5.125}
\]

其中：

\[
A = \frac{1}{\sigma^2} F^T F, \qquad
B = \frac{1}{\sigma^2} F^T A, \qquad
C = \frac{1}{\sigma^2} A^T F, \qquad
D = \frac{1}{\sigma^2} A^T A
\tag{5.126}
\]

根据给出的分块矩阵求逆公式，$\mathbf{I}_{\text{red}}$ 的左上角分块为：

\[
\boxed{
\mathbf{I}_{\theta}^{-1} = A^{-1} + A^{-1}B \left( D - C A^{-1} B \right)^{-1} C A^{-1}
}
\tag{5.127}
\]

此即角度参数的 CRLB。

---

##### 第一步：代入 $A, B, C, D$ 的具体表达式

首先计算 $A^{-1}$：

\[
A^{-1} = \left( \frac{1}{\sigma^2} F^T F \right)^{-1} = \sigma^2 (F^T F)^{-1}
\tag{5.128}
\]

接着计算括号内的核心项 $D - C A^{-1} B$：

\[
D - C A^{-1} B
=
\frac{1}{\sigma^2} A^T A
-
\left( \frac{1}{\sigma^2} A^T F \right)
\left( \sigma^2 (F^T F)^{-1} \right)
\left( \frac{1}{\sigma^2} F^T A \right)
\tag{5.129}
\]

注意到三个因子相乘的系数：$\frac{1}{\sigma^2} \cdot \sigma^2 \cdot \frac{1}{\sigma^2} = \frac{1}{\sigma^2}$，因此：

\[
D - C A^{-1} B
=
\frac{1}{\sigma^2} A^T A
-
\frac{1}{\sigma^2} A^T F (F^T F)^{-1} F^T A
\tag{5.130}
\]

提取公因子 $\frac{1}{\sigma^2}$：

\[
D - C A^{-1} B
=
\frac{1}{\sigma^2} A^T
\left( I - F (F^T F)^{-1} F^T \right)
A
\tag{5.131}
\]

括号中的 $I - F (F^T F)^{-1} F^T$ 是到 $F$ 的列空间的正交补上的投影矩阵，记为：

\[
P_F^\perp = I - F (F^T F)^{-1} F^T
\tag{5.132}
\]

因此：

\[
D - C A^{-1} B = \frac{1}{\sigma^2} A^T P_F^\perp A
\tag{5.133}
\]

取逆：

\[
\left( D - C A^{-1} B \right)^{-1} = \sigma^2 \left( A^T P_F^\perp A \right)^{-1}
\tag{5.134}
\]

---

##### 第二步：将中间结果代入 (5.127)

逐项代入：

**第一项** $A^{-1}$ 已在 (5.128) 中给出：

\[
A^{-1} = \sigma^2 (F^T F)^{-1}
\tag{5.135}
\]

**第二项** 计算 $A^{-1}B$ 和 $C A^{-1}$：

\[
A^{-1}B = \sigma^2 (F^T F)^{-1} \cdot \frac{1}{\sigma^2} F^T A
= (F^T F)^{-1} F^T A
\tag{5.136}
\]

\[
C A^{-1} = \frac{1}{\sigma^2} A^T F \cdot \sigma^2 (F^T F)^{-1}
= A^T F (F^T F)^{-1}
\tag{5.137}
\]

将 (5.136)、(5.134)、(5.137) 代入第二项 $A^{-1}B (D - C A^{-1} B)^{-1} C A^{-1}$：

\[
A^{-1}B (D - C A^{-1} B)^{-1} C A^{-1}
=
(F^T F)^{-1} F^T A
\cdot \sigma^2 (A^T P_F^\perp A)^{-1}
\cdot A^T F (F^T F)^{-1}
\tag{5.138}
\]

提取 $\sigma^2$：

\[
A^{-1}B (D - C A^{-1} B)^{-1} C A^{-1}
=
\sigma^2 (F^T F)^{-1} F^T A
(A^T P_F^\perp A)^{-1}
A^T F (F^T F)^{-1}
\tag{5.139}
\]

---

##### 第三步：合并两项，得到 $\mathbf{I}_{\theta}^{-1}$

将 (5.135) 和 (5.139) 相加：

\[
\mathbf{I}_{\theta}^{-1}
=
\sigma^2 (F^T F)^{-1}
+
\sigma^2 (F^T F)^{-1} F^T A
(A^T P_F^\perp A)^{-1}
A^T F (F^T F)^{-1}
\tag{5.140}
\]

提取公因子 $\sigma^2$：

\[
\boxed{
\mathbf{I}_{\theta}^{-1}
=
\sigma^2
\left[
(F^T F)^{-1}
+
(F^T F)^{-1} F^T A
(A^T P_F^\perp A)^{-1}
A^T F (F^T F)^{-1}
\right]
}
\tag{5.141}
\]

---

##### 第四步：利用 Woodbury 恒等式化简

为了将 (5.141) 化简为更简洁的形式，利用 **Woodbury 矩阵恒等式**：

\[
(X + U W^{-1} V)^{-1} = X^{-1} - X^{-1} U (W + V X^{-1} U)^{-1} V X^{-1}
\tag{5.142}
\]

或者将其等价地写成：

\[
(F^T F - F^T A (A^T A)^{-1} A^T F)^{-1}
=
(F^T F)^{-1}
+
(F^T F)^{-1} F^T A
\left( A^T A - A^T F (F^T F)^{-1} F^T A \right)^{-1}
A^T F (F^T F)^{-1}
\tag{5.143}
\]

注意 (5.143) 的右边与 (5.141) 括号内的形式完全一致（因为 $A^T A - A^T F (F^T F)^{-1} F^T A = A^T P_F^\perp A$）。

因此，(5.141) 等价于：

\[
\mathbf{I}_{\theta}^{-1}
=
\sigma^2
\left(
F^T F - F^T A (A^T A)^{-1} A^T F
\right)^{-1}
\tag{5.144}
\]

---

##### 第五步：引入投影矩阵 $P_A^\perp$

在 (5.144) 中，定义投影矩阵：

\[
P_A^\perp = I - A (A^T A)^{-1} A^T
\tag{5.145}
\]

这是到 $A$ 的列空间的正交补上的投影矩阵。

则：

\[
F^T F - F^T A (A^T A)^{-1} A^T F
=
F^T \left( I - A (A^T A)^{-1} A^T \right) F
=
F^T P_A^\perp F
\tag{5.146}
\]

因此：

\[
\boxed{
\operatorname{CRLB}(\theta) = \sigma^2 \left( F^T P_A^\perp F \right)^{-1}
}
\tag{5.147}
\]

---

##### 小结

推导过程：

1. 从 Fisher 信息矩阵的分块结构出发，利用给出的分块矩阵求逆公式：

\[
\mathbf{I}_{\theta}^{-1} = A^{-1} + A^{-1}B (D - C A^{-1} B)^{-1} C A^{-1}
\]

2. 代入 $A = \frac{1}{\sigma^2} F^T F$，$B = \frac{1}{\sigma^2} F^T A$，$C = \frac{1}{\sigma^2} A^T F$，$D = \frac{1}{\sigma^2} A^T A$，经过代入、化简、合并，得到 (5.141)。

3. 利用 Woodbury 恒等式，将 (5.141) 化简为 (5.144)。

4. 引入投影矩阵 $P_A^\perp$，得到最终结果 (5.147)。

(5.147) 就是 DOA 估计中角度参数的 Cramér-Rao 下界。

#### 2.8.4 CRLB 表达式的物理含义

(5.147) 揭示了 DOA 估计精度的几个关键因素：

1. **噪声功率 \( \sigma^2 \)**: CRLB 与 \( \sigma^2 \) 成正比。信噪比越低，CRLB 越高，估计越困难。

2. **投影矩阵 \( P_A^\perp \)**: \( P_A^\perp \) 将 \( F \) 投影到信号子空间 \( A \) 的正交补上。这个投影操作"消除"了信号参数 \( S \) 未知带来的影响。物理上，它表示：只有那些无法用信号子空间表示的角度变化信息，才能用来估计角度。

3. **\( F \) 的列**: \( F \) 的每一列是 \( A(\theta)S \) 对相应角度的偏导数。\( F \) 的列越"大"（即角度变化引起的接收数据变化越大），CRLB 越小，估计越精确。

4. **\( F \) 与 \( A \) 的几何关系**: 如果 \( F \) 的列落在 \( A \) 的列空间中（即 \( P_A^\perp F = 0 \)），则 CRLB 趋于无穷——这意味着角度变化完全可以通过信号参数的变化来补偿，数据中没有任何角度信息。即相干信号导致秩亏的极端情况。

对于多快拍情况（\( L \) 个独立快拍），Fisher 信息具有可加性，CRLB 变为：

\[
\boxed{
\operatorname{CRLB}(\theta) = \frac{\sigma^2}{L} \left( F^T P_A^\perp F \right)^{-1}
}
\tag{5.148}
\]

CRLB 与快拍数 \( L \) 成反比——快拍越多，估计越精确。

---

#### 2.8.5 与 MUSIC/ESPRIT 性能的关系

(5.148) 中的 \( F^T P_A^\perp F \) 被称为 **阵列流型的导数与信号子空间正交补的关系矩阵**。在渐进条件下（大快拍、高信噪比），MUSIC 和 ESPRIT 的估计方差趋近于这个 CRLB：

\[
\operatorname{Var}(\hat{\theta}_{\text{MUSIC}}) \xrightarrow{L \to \infty} \operatorname{CRLB}(\theta)
\tag{5.149}
\]

即 MUSIC 和 ESPRIT 被称为"渐近有效"算法的理论依据——在大样本条件下，它们达到了参数估计的理论极限。

---

## 3. 课后总结

### 3.1 核心逻辑链：从统计模型到 CRLB

本讲以 DOA 估计的 Cramér-Rao 下界为主线，从阵列数据的高斯统计模型出发，经得分函数、Fisher 信息矩阵的分块结构，最终导出简洁的 CRLB 表达式。核心逻辑链条：

1. **统计模型**：阵列接收数据在确定性信号假设下服从复高斯分布 (5.4)。参数向量包含角度 \(\theta\)、信号复幅度 \(S\)（实部和虚部）以及噪声方差 \(\sigma^2\) (5.5)。

2. **Slepian-Bangs 公式**：对于均值 \(\mu(\phi)\) 和协方差 \(C(\phi)\) 均依赖参数的高斯模型，Fisher 信息矩阵的通用表达式为 (5.107)：
   \[
   [\mathbf{I}(\phi)]_{ij} = (\mu_i')^T C^{-1} (\mu_j') + \frac{1}{2}\text{Tr}\big((C)_i' C^{-1} (C)_j' C^{-1}\big)
   \]
   第一项来自均值对参数的依赖，第二项来自协方差对参数的依赖。

3. **得分函数的分解与期望**：将得分函数分解为线性部分 \(L_i\) 和中心化二次型 \(U_i\)。利用高斯矩性质，分别计算 \(E[L_i L_j]\) (5.98)、\(E[U_i U_j]\) (5.106)，且交叉项 \(E[L_i U_j] = 0\)，最终合并得到 Fisher 信息矩阵 (5.96)。

4. **Fisher 信息的块结构**：对 DOA 估计问题，参数向量为 \(\phi = (\theta, \text{Re}(S), \text{Im}(S), \sigma^2)\)。代入具体模型后，Fisher 信息矩阵呈现 3×3 分块结构 (5.123)：
   \[
   \mathbf{I}(\phi) = \frac{1}{\sigma^2}
   \begin{pmatrix}
   F^T F & F^T A & \mathbf{0} \\
   A^T F & A^T A & \mathbf{0} \\
   \mathbf{0} & \mathbf{0} & N/(2\sigma^2)
   \end{pmatrix}
   \]
   其中 \(F\) 的列为 \(A(\theta)S\) 对角度的偏导数，\(A\) 为阵列流型矩阵。

5. **Schur 补求逆**：角度参数 \(\theta\) 和信号参数 \(S\) 存在耦合，CRLB 为 Fisher 信息矩阵逆的左上角块。利用分块矩阵求逆公式 (5.127)，角度参数的 CRLB 为：
   \[
   \mathbf{I}_{\theta}^{-1} = \sigma^2 \big(F^T F - F^T A (A^T A)^{-1} A^T F\big)^{-1}
   \]
   
6. **投影矩阵简化**：引入投影矩阵 \(P_A^\perp = I - A (A^T A)^{-1} A^T\) (5.145)，将中间项重写为 \(F^T P_A^\perp F\)，得到最终结果：
   \[
   \boxed{\operatorname{CRLB}(\theta) = \sigma^2 \big(F^T P_A^\perp F\big)^{-1}}
   \tag{5.147}
   \]

7. **推广至多快拍**：\(L\) 个独立快拍下 Fisher 信息具有可加性，CRLB 与快拍数成反比：
   \[
   \boxed{\operatorname{CRLB}(\theta) = \frac{\sigma^2}{L} \big(F^T P_A^\perp F\big)^{-1}}
   \tag{5.148}
   \]

### 3.2 CRLB 与 MUSIC/ESPRIT 的对比

| 维度 | CRLB | MUSIC | ESPRIT |
| :--- | :--- | :--- | :--- |
| **性质** | 理论下界（任何无偏估计器） | 实际估计算法 | 实际估计算法 |
| **计算方式** | Fisher 信息矩阵求逆 | 噪声子空间与导向矢量正交性 | 子阵间旋转不变性 |
| **所需信息** | 真实参数值 | 样本协方差矩阵的特征分解 | 样本协方差矩阵的特征分解 |
| **前提假设** | 统计模型正确 | 信号数已知，不相干 | 阵列具有平移不变结构 |
| **大样本行为** | — | 趋近 CRLB（渐近有效） | 趋近 CRLB（渐近有效） |
| **相干信号** | 需广义 CRLB | 失败（需空间平滑） | 失败（需空间平滑） |

### 3.3 重点概念总结

#### 3.3.1 投影矩阵 \(P_A^\perp\) 的物理含义

\(P_A^\perp = I - A(A^T A)^{-1} A^T\) 将任意向量投影到阵列流型 \(A\) 的列空间的正交补上。在 CRLB 表达式中，它起到关键作用：**消除信号参数 \(S\) 未知带来的影响**。物理上，它表示只有那些无法用信号子空间（即调整 \(S\)）来补偿的角度变化信息，才能用来估计角度。如果 \(F\) 的列完全落在 \(A\) 的列空间中（\(P_A^\perp F = 0\)），则 CRLB 趋于无穷——角度变化完全可以被信号幅度的变化"吸收"，数据中不含角度信息。

#### 3.3.2 Slepian-Bangs 公式的作用

Slepian-Bangs 公式是这套推导的核心工具。它将 Fisher 信息矩阵的计算分解为两部分：
- **均值贡献**：\((\mu_i')^T C^{-1} (\mu_j')\)——参数变化引起的均值变化越大，信息越多
- **协方差贡献**：\(\frac{1}{2}\text{Tr}((C)_i' C^{-1} (C)_j' C^{-1})\)——参数变化引起的协方差变化越大，信息越多

在 DOA 问题中，协方差 \(C = \sigma^2 I\) 不依赖角度和信号参数，因此协方差贡献项仅出现在噪声方差 \(\sigma^2\) 的 Fisher 信息中。对于角度参数，Fisher 信息完全来自均值贡献。

#### 6.3.3 CRLB 的决定因素

从 (5.148) 可以看出，DOA 估计精度由以下因素共同决定：

1. **噪声功率 \(\sigma^2\)**：CRLB ∝ \(\sigma^2\)，信噪比越高，下界越低
2. **快拍数 \(L\)**：CRLB ∝ \(1/L\)，快拍越多，下界越低
3. **阵列几何**（通过 \(F\) 体现）：\(F\) 的列是 \(A(\theta)S\) 对角度的偏导数，阵列孔径越大、阵元数越多，导数越大，CRLB 越小
4. **信号方向**（通过 \(A\) 和 \(F\) 体现）：不同角度的导向矢量导数不同，导致 CRLB 随角度变化
5. **信号间的相关性**：相干信号导致 \(P_A^\perp F\) 秩亏，CRLB 趋于无穷
6. **\(F\) 与 \(A\) 的几何关系**：\(F^T P_A^\perp F\) 度量了角度变化中"无法被信号幅度变化解释"的部分

---

## 4. 学习检查清单：自测核心知识点掌握情况

- [ ] 能写出阵列处理的确定性信号模型 \(X = A(\theta)S + N\)，说明 \(A(\theta), S, N\) 的维度和统计特性
- [ ] 能写出完整参数向量 \(\phi = (\theta, \text{Re}(S), \text{Im}(S), \sigma^2)\) 并说明各分量的含义
- [ ] 能写出 Cramér-Rao 下界的定义：\(\operatorname{CRLB}(\phi) = \mathbf{I}(\phi)^{-1}\)，其中 \(\mathbf{I}(\phi)\) 为 Fisher 信息矩阵
- [ ] 能说出 Slepian-Bangs 公式的适用条件（高斯模型，均值和协方差均依赖参数）
- [ ] 能写出 Slepian-Bangs 公式 (5.107) 的完整表达式，并解释两项的物理含义
- [ ] 能解释为什么在 DOA 问题中，角度参数的 Fisher 信息仅来自均值贡献（协方差 \(C = \sigma^2 I\) 不依赖 \(\theta\)）
- [ ] 能将得分函数分解为 \(L_i\)（线性项）和 \(U_i\)（中心化二次型），并说明交叉项为零的原因
- [ ] 能用高斯矩性质计算 \(E[L_i L_j]\) 和 \(E[U_i U_j]\)，理解 \(\text{Tr}(\cdot)\) 项的来源
- [ ] 能写出 DOA 估计的 Fisher 信息矩阵的 3×3 分块结构 (5.123)，并解释零块出现的原因
- [ ] 能解释为什么 Fisher 信息矩阵中 \(\theta\)-\(S\) 耦合块非零——这说明角度和信号幅度在估计中是相互影响的
- [ ] 能用分块矩阵求逆公式 (5.127) 推导 Schur 补形式的 CRLB
- [ ] 能定义投影矩阵 \(P_A^\perp = I - A(A^T A)^{-1} A^T\)，并解释其几何含义
- [ ] 能推导 \(F^T F - F^T A (A^T A)^{-1} A^T F = F^T P_A^\perp F\)，并解释为什么这个化简是关键的
- [ ] 能写出最终的 CRLB 表达式 (5.147) 和 (5.148)，并解释每个因子的物理含义
- [ ] 能说明为什么多快拍下 CRLB 与 \(L\) 成反比
- [ ] 能解释投影矩阵 \(P_A^\perp\) 如何"消除"信号参数 \(S\) 未知的影响
- [ ] 能分析在什么条件下 CRLB 会趋于无穷（相干信号、秩亏）
- [ ] 能解释为什么 MUSIC 和 ESPRIT 被称为"渐近有效"估计器 (5.149)

---

## 5. 思考题：拓展与挑战

1. **Slepian-Bangs 公式的退化**：如果协方差矩阵不依赖参数（\(C_i' = 0\)），Slepian-Bangs 公式退化为经典形式 \([\mathbf{I}]_{ij} = (\mu_i')^T C^{-1} (\mu_j')\)。请推导这个退化结果，并证明在 DOA 问题中，角度和信号参数的 Fisher 信息确实只由均值项贡献。噪声方差 \(\sigma^2\) 的 Fisher 信息中，均值项和协方差项各贡献多少？

2. **\(F^T P_A^\perp F\) 的秩与可估计性**：(5.147) 中 CRLB 需要对 \(F^T P_A^\perp F\) 求逆。什么情况下这个矩阵不可逆？如果只有 \(M\) 个信号源但 \(F\) 有 \(M\) 列，\(F^T P_A^\perp F\) 是否一定满秩？请分析相干信号和单快拍场景下秩的情况，并讨论这对 CRLB 意味着什么。

3. **ULA 的 CRLB 显式表达式**：对于 \(N\) 阵元的均匀线阵（ULA），导向矢量为 \(\mathbf{a}(\theta) = [1, e^{j\omega}, \dots, e^{j(N-1)\omega}]^T\)，其中 \(\omega = (2\pi d/\lambda)\sin\theta\)。当只有一个信号源时：
   - 求出 \(F\) 的显式表达式
   - 求出 \(F^T P_A^\perp F\) 的标量形式
   - 推导 CRLB\((\theta)\) 关于阵元数 \(N\)、SNR 和 \(\theta\) 的显式公式
   - 分析 CRLB 在 \(\theta = 0^\circ\) 和 \(\theta = 90^\circ\) 附近的行为（端射 vs 法线方向）

4. **CRLB 与波束宽度的关系**：ULA 的 3 dB 波束宽度约为 \(\Delta\theta \approx 0.886 \lambda / (Nd\cos\theta)\)。CRLB 给出的角度标准差下界约为 \(\sqrt{\operatorname{CRLB}(\theta)}\)。请比较这两个量：在什么条件下 CRLB 给出的精度可以远小于波束宽度？这在物理上意味着什么？（提示：超分辨）

5. **信号相关性对 CRLB 的影响**：考虑两个等功率信号源，相关系数为 \(\rho\)。当 \(\rho \to 1\)（接近相干）时，CRLB 如何变化？请分析 \(P_A^\perp\) 的秩和 \(F^T P_A^\perp F\) 的条件数如何随 \(\rho\) 变化，并给出物理解释。为什么 MUSIC 在 \(\rho < 1\) 时仍能工作但 CRLB 已经很大？

6. **阵列几何的优化**：给定 \(N\) 个阵元和一个固定的阵列孔径 \(D\)，如何放置阵元使某个方向 \(\theta_0\) 上的 CRLB 最小？均匀线阵是否是最优的？如果允许非均匀布阵，CRLB 能降低多少？（提示：考虑 \(F^T P_A^\perp F\) 对阵列几何的依赖）

7. **CRLB 与 MUSIC 谱峰的关系**：MUSIC 谱为 \(P_{\text{MUSIC}}(\theta) = 1 / \|\mathbf{a}^H(\theta) \mathbf{E}_n\|^2\)。在真实角度 \(\theta_0\) 附近，MUSIC 谱的二阶导数与 CRLB 存在关系。请定性说明：CRLB 越小（估计越精确），MUSIC 谱峰应该越尖锐还是越平坦？为什么？

---

## 6. 实验设计：基于 ESP32 阵列和 2×2 MIMO SDR 的 DOA 估计精度验证

### 6.0 实验总览：在真实硬件上验证 CRLB 理论

本实验利用约 20 个 ESP32 组成的均匀线阵（ULA）进行 MUSIC/ESPRIT 的 DOA 估计，将实测估计方差与理论 CRLB 对比；同时利用 2×2 MIMO SDR 在小规模阵列上验证 CRLB 对阵列几何的依赖关系。

#### 核心思路

CRLB 是理论下界，在真实硬件上无法直接"测量"CRLB，但可以测量实际估计算法（MUSIC/ESPRIT）的估计方差，并与 CRLB 理论值比较。二者之比（效率）反映了算法的实际性能。实验的核心目标是：**在 ESP32 大阵列上，通过重复实验估计 MUSIC 的 DOA 方差，验证其在大快拍/高 SNR 下趋近 CRLB；在 2×2 MIMO SDR 上，验证阵列孔径和阵元数对 CRLB 的影响。**

#### 实验矩阵

| 实验编号 | 实验名称 | 使用硬件 | 对比维度 | 预期结论 |
| :--- | :--- | :--- | :--- | :--- |
| 9.1 | ESP32 ULA 的 MUSIC DOA 估计与方差统计 | ESP32 阵列 | SNR、快拍数 | 估计方差随 SNR 和快拍数增大而降低 |
| 9.2 | 实测方差 vs 理论 CRLB 对比 | ESP32 阵列 | 算法效率 | 高 SNR / 大快拍下 MUSIC 方差趋近 CRLB |
| 9.3 | 两信号源分辨极限与 CRLB 的关系 | ESP32 阵列 | 角度间隔 | 当角度间隔接近 CRLB 标准差时，MUSIC 无法分辨 |
| 9.4 | 2×2 MIMO SDR 的阵列几何与 CRLB | 2×2 SDR | 阵元间距 | 阵元间距增大 → \(F\) 的导数增大 → CRLB 降低 |

---

### 6.1 实验 1：ESP32 ULA 的 MUSIC DOA 估计与方差统计

**目的**：在 20 阵元 ESP32 均匀线阵上实现 MUSIC 算法，对不同 SNR 和快拍数条件下的 DOA 估计进行重复实验，统计估计方差。

**硬件配置**：
- 约 20 个 ESP32 模块组成均匀线阵，阵元间距 \(d = \lambda/2\)（2.4 GHz 下约 6.25 cm）
- 总计阵列孔径 ≈ 19 × 6.25 cm ≈ 1.19 m
- 一个信号源（如另一个 ESP32 或信号发生器）放置在阵列前方已知角度（如 \(\theta = 20^\circ\)）
- 所有 ESP32 通过 WiFi 或有线方式同步采集，将 IQ 数据汇总到 PC

**步骤**：
1. 信号源发射单频正弦波（CW），频率 2.4 GHz
2. 所有 ESP32 阵元同步采集 \(L\) 个快拍（\(L = 64, 128, 256, 512, 1024\)）
3. 调整发射功率或距离，使 SNR 在 0 dB 到 30 dB 范围内变化
4. 对每组 (SNR, \(L\)) 参数，重复采集 \(K = 100\) 次独立实验
5. 对每次实验数据：
   - 计算样本协方差矩阵 \(\hat{\mathbf{R}} = \frac{1}{L} \mathbf{X}\mathbf{X}^H\)
   - 做特征分解，分离信号子空间和噪声子空间
   - 用 MUSIC 谱搜索 DOA 估计值 \(\hat{\theta}\)
6. 统计 100 次实验的 DOA 估计方差：
   \[
   \operatorname{Var}(\hat{\theta}) = \frac{1}{K-1}\sum_{k=1}^{K} (\hat{\theta}_k - \bar{\theta})^2
   \]
7. 对每组参数计算理论 CRLB (5.148)：\(\operatorname{CRLB} = \frac{\sigma^2}{L}(F^T P_A^\perp F)^{-1}\)（矩阵退化为标量）

**定量指标**：
- \(\operatorname{Var}(\hat{\theta})\) vs SNR 曲线（固定 \(L\)）
- \(\operatorname{Var}(\hat{\theta})\) vs \(L\) 曲线（固定 SNR）
- 对比理论 CRLB 曲线

**预期结果**：
- 估计方差随 SNR 增大而降低（约 -3 dB/倍 SNR）
- 估计方差随 \(L\) 增大而降低（约 1/\(L\) 关系）
- 高 SNR 和大 \(L\) 区域，实测方差接近 CRLB 理论值
- 低 SNR 区域（< 0 dB），MUSIC 出现"门限效应"——方差急剧增大，远高于 CRLB

---

### 6.2 实验 2：实测方差 vs 理论 CRLB——算法效率评估

**目的**：量化 MUSIC 算法的估计效率（实测方差与 CRLB 的比值），验证"渐近有效"的理论结论。

**步骤**：
1. 使用实验 9.1 中采集的全部数据
2. 定义算法效率：
   \[
   \eta(\hat{\theta}) = \frac{\operatorname{CRLB}(\theta)}{\operatorname{Var}(\hat{\theta})}
   \]
   \(\eta \to 1\) 表示算法达到了理论下界；\(\eta < 1\) 表示算法效率不足
3. 分别绘制 MUSIC 的效率 vs SNR 曲线 和 效率 vs \(L\) 曲线
4. 标注 \(\eta = 0.9\) 的"有效区域"边界（效率超过 90%）

**定量指标**：
- 达到 \(\eta \geq 0.9\) 所需的最低 SNR（固定大 \(L\)）
- 达到 \(\eta \geq 0.9\) 所需的最低快拍数 \(L_{\min}\)（固定高 SNR）
- 门限 SNR 与阵元数 \(N\) 的关系：\(N=20\) 时门限值？

**预期结果**：
- 低 SNR / 小 \(L\) 区域：\(\eta\) 远小于 1（MUSIC 效率低，受门限效应影响）
- 高 SNR / 大 \(L\) 区域：\(\eta \to 1\)，验证 MUSIC 的渐近有效性
- 门限 SNR 随 \(N\) 增大而降低——阵元数越多，MUSIC 进入有效区域所需的 SNR 越低

---

### 6.3 实验 3：两信号源分辨极限与 CRLB

**目的**：使用两个信号源，测试 MUSIC 的角度分辨极限，并与 CRLB 给出的精度下界对比。

**硬件配置**：
- ESP32 阵列同上
- 两个信号源放置在角度 \(\theta_1\) 和 \(\theta_2 = \theta_1 + \Delta\theta\) 处
- 逐步减小 \(\Delta\theta\)，从 \(20^\circ\) 到 \(2^\circ\)

**步骤**：
1. 两个信号源同时发射不同频率（或时分交替发射以分离），确保两信号不相干
2. 从 \(\Delta\theta = 20^\circ\) 开始，逐步减小至 \(\Delta\theta = 2^\circ\)（步长 \(1^\circ\) 或 \(2^\circ\)）
3. 对每个 \(\Delta\theta\)，重复采集 100 次
4. 对每次实验运行 MUSIC，记录：
   - 是否成功分辨两个峰值（分辨标准：MUSIC 谱在两真实角度处有明显的局部最大值，且两峰之间的谷值低于峰值的 -3 dB）
   - 两个角度的估计值 \(\hat{\theta}_1, \hat{\theta}_2\)
5. 计算每个 \(\Delta\theta\) 下的：
   - 分辨成功率（分辨出两个峰的概率）
   - 两个角度的估计方差
   - 两个角度的理论 CRLB（注意此时 CRLB 为 2×2 矩阵，取其对角元）

**定量指标**：
- 分辨成功率 vs \(\Delta\theta\) 曲线
- 估计方差 vs \(\Delta\theta\) 曲线（叠加理论 CRLB）
- 临界分辨角 \(\Delta\theta_{\text{crit}}\)（分辨成功率达 50% 的角度间隔）
- 对比 \(\Delta\theta_{\text{crit}}\) 与 Rayleigh 分辨率 (\(\approx 0.886 \lambda / (Nd\cos\theta)\)) 和 CRLB 标准差 \(\sqrt{\operatorname{CRLB}}\)

**预期结果**：
- MUSIC 的超分辨能力：\(\Delta\theta_{\text{crit}}\) 小于 Rayleigh 分辨率
- 当 \(\Delta\theta\) 接近 \(\sqrt{\operatorname{CRLB}}\) 时，估计方差急剧增大——两个信号的估计存在耦合
- 当 \(\Delta\theta\) 小于 CRLB 标准差时，两信号在统计意义上"不可区分"

---

### 6.4 实验 4：2×2 MIMO SDR——阵列几何对 CRLB 的影响

**目的**：在最小阵列（2 阵元）上验证阵元间距如何影响 CRLB，并利用 SDR 的高精度 IQ 采集获得干净的测量数据。

**硬件配置**：
- 2×2 MIMO SDR 平台（如 USRP B210 或 HackRF + 双通道接收）
- 两个接收通道配置为两阵元接收阵列
- 阵元间距可调（从 \(d = \lambda/4\) 到 \(d = 2\lambda\)）

**步骤**：
1. 信号源放置在 \(\theta_0 = 30^\circ\)（法线偏开一定角度以获得可测的相位差）
2. 阵元间距分别设置为 \(d = \lambda/4, \lambda/2, \lambda, 2\lambda\) 四个值
3. 对每个阵元间距，采集 \(L = 1024\) 快拍，重复 50 次
4. 用 MUSIC 进行 DOA 估计（单源，搜索范围 \([-90^\circ, 90^\circ]\)）
5. 统计每个间距下的估计方差
6. 计算各间距下的理论 CRLB：
   - 对于间距 \(d\)，阵列流型：\(\mathbf{a}(\theta) = [1, e^{j(2\pi d/\lambda)\sin\theta}]^T\)
   - \(F = \partial (A(\theta)S) / \partial \theta\)，这里 \(A(\theta) = \mathbf{a}(\theta)\)
   - 代入 (5.148) 计算 CRLB

**定量指标**：
- 实测方差 vs 阵元间距 \(d/\lambda\) 曲线
- 理论 CRLB vs \(d/\lambda\) 曲线（叠加对比）
- 计算 CRLB 下降因子：从 \(d = \lambda/4\) 到 \(d = 2\lambda\)，CRLB 降低了多少？

**预期结果**：
- \(d = \lambda/4\)：CRLB 最大（相位差变化小 → 角度变化不敏感）
- \(d = \lambda/2\)：CRLB 居中（标准 ULA 配置）
- \(d = 2\lambda\)：CRLB 最小（相位差变化最大），但可能出现角度模糊（\(\sin\theta\) 超出 \([-1,1]\) 范围的相位缠绕）
- 验证 \(F\) 中导数项与阵元间距成正比：\(\partial \mathbf{a}/\partial\theta \propto d\)
- 在无模糊条件下，CRLB ∝ \(1/d^2\)

**进阶思考**：2 阵元时，如果 \(d > \lambda/2\) 会导致角度模糊（多个 \(\theta\) 对应相同相位差）。但 CRLB 公式并不"知道"模糊的存在——它只关心局部曲率。这说明了 CRLB 的局限性：它是一个局部下界，不包含全局模糊的信息。在实验中，如何判断是否出现了角度模糊？

---

### 6.5 实验报告要求

1. 给出 ESP32 阵列的几何参数（阵元数、间距、总孔径）和校准方法
2. 绘制实验 9.1 中 MUSIC 估计方差 vs SNR 和 vs 快拍数的完整曲线，叠加理论 CRLB 曲线
3. 绘制实验 9.2 中 MUSIC 效率曲线，标注门限 SNR 和达到渐近有效所需的最小快拍数
4. 绘制实验 9.3 中分辨成功率 vs \(\Delta\theta\) 曲线，标注 Rayleigh 分辨率和 CRLB 标准差
5. 给出实验 9.4 中四个阵元间距的 SDR 实测结果，包括：
   - 每种间距下的 MUSIC 谱图（展示峰的形状和宽度）
   - 实测方差 vs 阵元间距曲线以及理论 CRLB 对比
6. 讨论：综合四个实验，MUSIC 在 ESP32 阵列上达到 CRLB 的"有效工作区"是什么（SNR 和快拍数范围）？哪些因素导致实测方差偏离 CRLB（硬件非理想性：通道不一致、互耦、近场效应、同步误差）？

---

<div style="page-break-before: always;"></div>