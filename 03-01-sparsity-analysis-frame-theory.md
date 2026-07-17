
<div style="page-break-before: always; padding: 8% 8% 0 8%;">
 <h1 id="第九讲-框架理论" style="text-align: center; margin-bottom: 2rem; border-bottom: none; display: block;">第九讲 框架理论</h1> 
 <div style="display: flex; align-items: center; justify-content: center; gap: 12px; margin: 1.8rem auto;">
  <span style="flex:1; max-width:80px; height:1px; background: linear-gradient(to right, transparent, #888);"></span>
  <span style="display:inline-block; width:6px; height:6px; background:#38bdf8; border-radius:50%;"></span>
  <span style="flex:1; max-width:80px; height:1px; background: linear-gradient(to left, transparent, #888);"></span>
 </div>
</div>


<!-- # 第九讲 框架理论 -->

从本文开始，后续五篇文章将系统讲解稀疏信号处理，包括框架定理、稀疏表示、稀疏恢复、稀疏编码等内容。前面的文章已涵盖统计信号处理基础、自适应滤波、谱分析等内容，这些理论几乎都建立在一个核心假设之上：**基是正交的**。例如，傅里叶变换中的复指数基、PCA中的特征向量基、小波变换中的正交小波基，都是正交基的典型代表。

然而，稀疏信号处理面对的核心问题是**非正交**的情况。在实际应用中，常常需要处理冗余的、非正交的表示。本章要回答的核心问题是：**能否找到一组非正交但有冗余的基（即框架），它们仍然保留标准正交基的一些优良性质？**

答案就是**框架理论**。

---

## 1. 标准正交基

在引入框架之前，先回顾标准正交基的性质。这些性质是框架理论试图推广的"黄金标准"。

设 \(\{x_1, x_2, \cdots, x_n\}\) 是 \(n\) 维向量空间中的一组基，它们是**线性无关**的，即不存在非零系数组合使得：

\[
\sum_{k=1}^n \alpha_k x_k = 0 \quad \Longleftrightarrow \quad \alpha_1 = \alpha_2 = \cdots = \alpha_n = 0
\tag{9.1}
\]

线性无关保证了基的**唯一表示性**——任意向量 \(X\) 都可以由这组基唯一地表示为：

\[
\boxed{
X = \sum_{k=1}^n \alpha_k x_k
}
\tag{9.2}
\]

这种一一对应关系 \(X \leftrightarrow (\alpha_1, \alpha_2, \cdots, \alpha_n)\) 建立了向量与其系数之间的坐标映射。这是标准正交基的第一个重要性质：**唯一表示**。

---

### 1.1 标准正交基的定义

如果上述基还满足**正交性**条件，则称为标准正交基（Orthonormal Basis, ONB）：

\[
\boxed{
\langle x_i, x_j \rangle = 
\begin{cases}
1, & i = j \\
0, & i \neq j
\end{cases}
}
\tag{9.3}
\]

其中 \(\langle \cdot, \cdot \rangle\) 表示内积，对于实向量空间为 \(\langle x, y \rangle = x^T y\)，对于复向量空间为 \(\langle x, y \rangle = x^H y\)。

正交性意味着这组基向量在空间中"各走各的路"，互不干扰；归一化则保证了每个基向量的长度为 1。

---

### 1.2 标准正交基下系数的求解

对于标准正交基，求系数 \(\alpha_j\) 的过程极其简洁。对 (9.2) 两边同时与 \(x_j\) 做内积：

\[
\langle X, x_j \rangle = \left\langle \sum_{k=1}^n \alpha_k x_k, \; x_j \right\rangle
\tag{9.4}
\]

利用内积的线性性质，将求和提出内积外：

\[
\langle X, x_j \rangle = \sum_{k=1}^n \alpha_k \langle x_k, x_j \rangle
\tag{9.5}
\]

根据正交性条件 (9.3)，当 \(k \neq j\) 时，\(\langle x_k, x_j \rangle = 0\)；当 \(k = j\) 时，\(\langle x_j, x_j \rangle = 1\)。因此，求和式中只有 \(k = j\) 这一项保留：

\[
\langle X, x_j \rangle = \alpha_j \langle x_j, x_j \rangle = \alpha_j \cdot 1 = \alpha_j
\tag{9.6}
\]

于是：

\[
\boxed{
\alpha_j = \langle X, x_j \rangle
}
\tag{9.7}
\]

这个公式说明：**在标准正交基下，信号的展开系数就是信号与对应基向量的内积。** 这正是匹配滤波、傅里叶变换、PCA 投影等操作背后的统一原理。

将 (9.7) 代入 (9.2)，得到**正交展开公式**：

\[
\boxed{
X = \sum_{k=1}^n \langle X, x_k \rangle \, x_k
}
\tag{9.8}
\]

这就是信号在标准正交基下的"重构公式"——它由所有投影分量叠加而成。

---

### 1.3 标准正交基的几何意义

从几何角度看，标准正交基的系数求解就是**正交投影**。向量 \(X\) 在 \(x_j\) 方向上的投影长度就是 \(\langle X, x_j \rangle\)，(9.8) 说明 \(X\) 是它在所有基方向上投影的向量和。这是 Parseval 定理的直接体现——信号的能量等于其所有系数能量之和：

\[
\|X\|^2 = \sum_{k=1}^n |\alpha_k|^2
\tag{9.9}
\]

---

### 1.4 标准正交基的局限：缺乏冗余

标准正交基的这些优点——唯一表示、系数求解简便、能量保持——背后隐藏着一个根本性的缺陷：基向量个数恰好等于空间的维度，每个向量只有一个唯一的表示，这意味着**信息没有任何冗余**。

这个缺点的后果是：**一旦丢失一路信息，就无法恢复原始信号。**

设想一个通信或存储系统：信号 \(X\) 在标准正交基 \(\{x_1, \cdots, x_n\}\) 下被分解为系数 \(\{\alpha_1, \cdots, \alpha_n\}\)。如果某个系数 \(\alpha_j\) 在传输或存储过程中丢失或损坏了，那么重构时：

\[
\hat{X} = \sum_{k \neq j} \alpha_k x_k
\tag{9.10}
\]

永远无法精确恢复 \(X\)，因为我们没有关于 \(\alpha_j\) 的任何额外信息。正交基提供了"最紧凑"的表示，但也意味着"最脆弱"——每一个分量都是必需的。

**对比冗余表示：**

如果在一个更高维的空间中表示信号，使得基向量的数量 \(m\) 大于空间的维度 \(n\)（即过完备表示），那么同一个信号可以有多种表示方式。即使丢失了某些分量，仍然可能从剩余分量中恢复出原始信号。这种**冗余**提供了对数据丢失、噪声和损坏的**鲁棒性**。

总结标准正交基的优缺点：

| 优点 | 缺点 |
| :--- | :--- |
| 唯一表示 | 无冗余 |
| 系数求解简便（内积） | 不鲁棒（丢失即无法恢复） |
| 能量保持（Parseval） | 每个分量都必不可少 |
| 计算高效 | 无法容忍数据损坏 |

**标准正交基的"脆弱性"促使我们寻找一种能够保留正交基的优点、同时引入冗余的新工具——这就是框架。** 下一节将从标准正交基出发，逐步推广到框架的定义和性质。

## 2. 框架

在上一节中，指出了标准正交基的根本局限：**无冗余，不鲁棒。** 只要丢失一个系数，信号就无法恢复。与此同时，介绍了由三个均匀分布的单位向量（0°、120°、240°）组成的冗余表示系统——向量的数量（3）大于空间维度（2），这天然具备了冗余性。

一个自然的问题浮现出来：**能否找到一组基，它既能提供冗余带来的稳定性，又能像标准正交基一样保持简洁的表示和重建公式？**

框架（Frame）理论正是要回答这个问题。

---

### 2.1 二维示例：三向量冗余表示

考虑二维空间 \( \mathbb{R}^2 \)，标准正交基为：

\[
a = (1, 0), \qquad b = (0, 1)
\tag{9.11}
\]

定义一组新的向量 \( \{x_1, x_2, x_3\} \)，它们是 \( a \) 和 \( b \) 的线性组合：

\[
x_1 = a = (1, 0)
\tag{9.12}
\]

\[
x_2 = -\frac{1}{2}a + \frac{\sqrt{3}}{2}b = \left(-\frac{1}{2}, \frac{\sqrt{3}}{2}\right)
\tag{9.13}
\]

\[
x_3 = \frac{1}{2}a - \frac{\sqrt{3}}{2}b = \left(\frac{1}{2}, -\frac{\sqrt{3}}{2}\right)
\tag{9.14}
\]

这三个向量长度均为 1，彼此之间的夹角为 120°，均匀分布在单位圆上。它们的关系可以用一张图直观地表示：

![框架示例](assets/09/02.png)

上图展示了这三个向量在 \( \mathbb{R}^2 \) 平面上的分布。这三个向量的数量（3 个）大于空间的维度（2 维），因此它们构成一个**过完备**（overcomplete）的表示系统。这正是框架的典型特征——基向量个数多于空间维度。

目标是：**对于任意向量 \( X = (x, y) \in \mathbb{R}^2 \)，能否像标准正交基一样，通过内积得到一组系数，再通过这些系数重建 \( X \)？**

---

### 2.2 投影系数的计算

按照标准正交基的做法，首先计算 \( X \) 与每个框架向量的内积，得到三个系数。

**与 \( x_1 \) 的内积：**

\[
\langle X, x_1 \rangle = \langle (x, y), (1, 0) \rangle = x \cdot 1 + y \cdot 0 = x
\tag{9.15}
\]

**与 \( x_2 \) 的内积：**

\[
\langle X, x_2 \rangle = \left\langle (x, y), \left(-\frac{1}{2}, \frac{\sqrt{3}}{2}\right) \right\rangle
= x \cdot \left(-\frac{1}{2}\right) + y \cdot \frac{\sqrt{3}}{2}
= -\frac{1}{2}x + \frac{\sqrt{3}}{2}y
\tag{9.16}
\]

**与 \( x_3 \) 的内积：**

\[
\langle X, x_3 \rangle = \left\langle (x, y), \left(\frac{1}{2}, -\frac{\sqrt{3}}{2}\right) \right\rangle
= x \cdot \frac{1}{2} + y \cdot \left(-\frac{\sqrt{3}}{2}\right)
= \frac{1}{2}x - \frac{\sqrt{3}}{2}y
\tag{9.17}
\]

这三个系数 \( \langle X, x_1 \rangle, \langle X, x_2 \rangle, \langle X, x_3 \rangle \) 构成了 \( X \) 在框架 \( \{x_1, x_2, x_3\} \) 下的**分析系数**（analysis coefficients）。

---

### 2.3 由投影系数重建信号

在标准正交基中，重建信号只需将每个系数乘以对应的基向量再求和（即 (9.8) 式）。尝试同样的操作——将三个投影向量相加：

\[
\sum_{k=1}^{3} \langle X, x_k \rangle \, x_k
= \langle X, x_1 \rangle x_1 + \langle X, x_2 \rangle x_2 + \langle X, x_3 \rangle x_3
\tag{9.18}
\]

逐项代入：

第一项：

\[
\langle X, x_1 \rangle x_1 = x \cdot (1, 0) = (x, 0)
\tag{9.19}
\]

第二项：

\[
\langle X, x_2 \rangle x_2 = \left(-\frac{1}{2}x + \frac{\sqrt{3}}{2}y\right) \left(-\frac{1}{2}, \frac{\sqrt{3}}{2}\right)
\tag{9.20}
\]

将其写成分量形式：

\[
= \left( -\frac{1}{2} \left(-\frac{1}{2}x + \frac{\sqrt{3}}{2}y\right), \; \frac{\sqrt{3}}{2} \left(-\frac{1}{2}x + \frac{\sqrt{3}}{2}y\right) \right)
\tag{9.21}
\]

化简：

\[
= \left( \frac{1}{4}x - \frac{\sqrt{3}}{4}y, \; -\frac{\sqrt{3}}{4}x + \frac{3}{4}y \right)
\tag{9.22}
\]

第三项：

\[
\langle X, x_3 \rangle x_3 = \left(\frac{1}{2}x - \frac{\sqrt{3}}{2}y\right) \left(\frac{1}{2}, -\frac{\sqrt{3}}{2}\right)
\tag{9.23}
\]

写成分量形式：

\[
= \left( \frac{1}{2} \left(\frac{1}{2}x - \frac{\sqrt{3}}{2}y\right), \; -\frac{\sqrt{3}}{2} \left(\frac{1}{2}x - \frac{\sqrt{3}}{2}y\right) \right)
\tag{9.24}
\]

化简：

\[
= \left( \frac{1}{4}x - \frac{\sqrt{3}}{4}y, \; -\frac{\sqrt{3}}{4}x + \frac{3}{4}y \right)
\tag{9.25}
\]

将三项相加，分别合并 \( x \) 分量和 \( y \) 分量。

**\( x \) 分量：**

第一项的 \( x \) 分量为 \( x \)。

第二项的 \( x \) 分量：\( \frac{1}{4}x - \frac{\sqrt{3}}{4}y \)。

第三项的 \( x \) 分量：\( \frac{1}{4}x - \frac{\sqrt{3}}{4}y \)。

相加：

\[
x + \left(\frac{1}{4}x - \frac{\sqrt{3}}{4}y\right) + \left(\frac{1}{4}x - \frac{\sqrt{3}}{4}y\right)
= x + \frac{1}{4}x + \frac{1}{4}x - \frac{\sqrt{3}}{4}y - \frac{\sqrt{3}}{4}y
= \frac{3}{2}x - \frac{\sqrt{3}}{2}y
\tag{9.26}
\]

**\( y \) 分量：**

第一项的 \( y \) 分量为 \( 0 \)。

第二项的 \( y \) 分量：\( -\frac{\sqrt{3}}{4}x + \frac{3}{4}y \)。

第三项的 \( y \) 分量：\( -\frac{\sqrt{3}}{4}x + \frac{3}{4}y \)。

相加：

\[
0 + \left(-\frac{\sqrt{3}}{4}x + \frac{3}{4}y\right) + \left(-\frac{\sqrt{3}}{4}x + \frac{3}{4}y\right)
= -\frac{\sqrt{3}}{2}x + \frac{3}{2}y
\tag{9.27}
\]

合并 (9.26) 和 (9.27)，得到：

\[
\sum_{k=1}^{3} \langle X, x_k \rangle x_k
= \left( \frac{3}{2}x - \frac{\sqrt{3}}{2}y, \; -\frac{\sqrt{3}}{2}x + \frac{3}{2}y \right)
\tag{9.28}
\]

这个结果并不是 \( X = (x, y) \)，而是某种线性变换作用于 \( X \)。将 (9.28) 写成矩阵形式：

\[
\sum_{k=1}^{3} \langle X, x_k \rangle x_k
=
\begin{pmatrix}
\frac{3}{2} & -\frac{\sqrt{3}}{2} \\
-\frac{\sqrt{3}}{2} & \frac{3}{2}
\end{pmatrix}
\begin{pmatrix}
x \\
y
\end{pmatrix}
\tag{9.29}
\]

如果进一步将这个矩阵作用于 \( X \)，能得到一个更干净的结果。事实上，这正是下一节要引入的**框架算子**——在这个例子中，框架算子恰好是 \( \frac{3}{2} \) 乘以恒等算子。关键结论是：无论 \( X \) 取什么值，这个线性组合总是给出 \( X \) 的一个固定倍数。

---

### 2.4 重建公式与关键发现

上面已经计算出，对于这个三向量框架，投影和与原始向量 \( X \) 之间的关系是：

\[
\boxed{
\sum_{k=1}^{3} \langle X, x_k \rangle x_k = \frac{3}{2} X
}
\tag{9.30}
\]

这个结果的含义是：**虽然这三个向量不是正交基（它们之间有冗余），但通过内积和求和的操作，仍然可以稳定地重建原始信号——只需要乘上一个常数因子 \( 2/3 \)。**

因此，重建公式为：

\[
\boxed{
X = \frac{2}{3} \sum_{k=1}^{3} \langle X, x_k \rangle x_k
}
\tag{9.31}
\]

这正是框架理论的精髓所在。将此例与标准正交基的重建公式对比：

| | 标准正交基 | 三向量框架 |
| :--- | :--- | :--- |
| 基/框架向量个数 | 2（等于维度） | 3（大于维度） |
| 重建公式 | \( X = \sum \langle X, x_k \rangle x_k \) | \( X = \frac{2}{3} \sum \langle X, x_k \rangle x_k \) |
| 冗余 | 无 | 有 |
| 稳定性 | 丢失一个系数即失败 | 有冗余，有鲁棒性 |
| 框架界 | \( A = B = 1 \) | \( A = B = \frac{3}{2} \) |

这个例子中，三个向量构成一个**紧框架**（Tight Frame）——框架界的下界 \( A \) 和上界 \( B \) 相等，均为 \( 3/2 \)。特别地，当 \( A = B = 1 \) 时，框架退化为标准正交基；当 \( A = B > 1 \) 时，框架是冗余的，但仍然保持 \( X = \frac{1}{A} \sum \langle X, x_k \rangle x_k \) 的简洁重建形式。这类框架称为 **Parseval 框架**（Parseval Frame）。

紧框架保留了标准正交基最核心的性质——简洁的重建公式——同时提供了冗余带来的鲁棒性。

---

### 2.5 框架的一般定义

将上面的例子推广到一般情况。设 \( \mathcal{H} \) 是一个有限维或无限维希尔伯特空间（在通常语境中，\( \mathcal{H} \) 指 \( \mathbb{R}^n \) 或 \( \mathbb{C}^n \)）。一族向量 \( \{x_k\}_{k \in \mathcal{I}} \)（\( \mathcal{I} \) 为指标集）被称为一个 **框架**，如果存在两个常数 \( 0 < A \leq B < \infty \)，使得对于任意 \( X \in \mathcal{H} \)，都有：

\[
\boxed{
A \|X\|^2 \leq \sum_{k \in \mathcal{I}} |\langle X, x_k \rangle|^2 \leq B \|X\|^2
}
\tag{9.32}
\]

这个条件称为 **框架条件**。它保证了表示系统的两个关键性质：

- **稳定性（下界 \( A > 0 \)）**：如果 \( X \neq 0 \)，则至少有一些系数 \( \langle X, x_k \rangle \) 是非零的。即使某些系数丢失，只要框架足够密集，仍能稳定重建。
- **能量有界（上界 \( B < \infty \)）**：系数能量不会无限放大，系统不会产生不稳定或数值发散。

当 \( A = B \) 时，框架称为 **紧框架**（Tight Frame）。当 \( A = B = 1 \) 时，框架就是标准正交基；当 \( A = B > 1 \) 时，框架是冗余的紧框架——重建公式为：

\[
\boxed{
X = \frac{1}{A} \sum_{k} \langle X, x_k \rangle x_k
}
\tag{9.33}
\]

这恰好是标准正交基重建公式的直接推广——唯一的区别是在求和前乘上一个归一化因子 \( 1/A \)。

在下一节中，将引入框架算子（Frame Operator），它是连接框架系数与原始信号的桥梁，也是理解框架性质的核心数学工具。

### 2.6 框架的基本性质

#### 2.6.1 Parseval

在标准正交基中，Parseval等式是一个核心性质：

\[
\|X\|^2 = \sum_{k=1}^{n} |\langle X, e_k \rangle|^2
\tag{9.34}
\]

它与重建公式 \( X = \sum_{k=1}^{n} \langle X, e_k \rangle e_k \) 是等价的。

问题是：**对于框架 \( \{x_k\} \)，是否也有类似的性质？** 换句话说，下面的两个命题是否等价？

\[
X = \sum_{k} \langle X, x_k \rangle x_k
\tag{9.35}
\]

\[
\|X\|^2 = \sum_{k} |\langle X, x_k \rangle|^2
\tag{9.36}
\]

答案是**肯定的**——对于紧框架，这两个命题完全等价。下面从两个方向分别证明。

---

#### 2.6.1a 充要性：重建公式 ⇒ 能量恒等式

**已知**：对于任意 \( X \)，都有

\[
X = \sum_{k} \langle X, x_k \rangle x_k
\tag{9.37}
\]

**要证明**：

\[
\|X\|^2 = \sum_{k} |\langle X, x_k \rangle|^2
\tag{9.38}
\]

**证明**：

对 (9.37) 两边同时与 \( X \) 做内积：

\[
\langle X, X \rangle = \left\langle X, \sum_{k} \langle X, x_k \rangle x_k \right\rangle
\]

左边就是 \( \|X\|^2 \)。右边利用内积的双线性性质——内积对第二个变量是线性的，因此求和可以提到内积外面：

\[
\|X\|^2 = \sum_{k} \langle X, \langle X, x_k \rangle x_k \rangle
\]

这里 \( \langle X, x_k \rangle \) 是一个标量（系数），它可以提到内积外面（内积对第二个变量是线性的，常数因子可以自由进出）。对于实内积，直接提取即可；对复内积需要取共轭，但最终结果都是模平方。写成实内积的形式：

\[
\|X\|^2 = \sum_{k} \langle X, x_k \rangle \langle X, x_k \rangle
= \sum_{k} |\langle X, x_k \rangle|^2
\]

因此：

\[
\boxed{
\|X\|^2 = \sum_{k} |\langle X, x_k \rangle|^2
}
\tag{9.38}
\]

这样就证明了：**如果重建公式成立，那么能量恒等式必然成立。**

---

#### 2.6.1b 充要性：能量恒等式 ⇒ 重建公式

**已知**：对于任意 \( X \)，都有

\[
\|X\|^2 = \sum_{k} |\langle X, x_k \rangle|^2
\tag{9.39}
\]

**要证明**：

\[
X = \sum_{k} \langle X, x_k \rangle x_k
\tag{9.40}
\]

**证明思路**：能量恒等式对所有 \( X \) 成立，意味着映射 \( X \mapsto \sum_{k} \langle X, x_k \rangle x_k \) 是保范数的。保范数的线性算子一定是等距同构，从而一定是恒等算子。

为此，定义**框架算子**：

\[
S(X) = \sum_{k} \langle X, x_k \rangle x_k
\tag{9.41}
\]

这个算子将 \( X \) 映射为它在框架 \( \{x_k\} \) 下的"投影和"。(9.39) 表明：

\[
\|X\|^2 = \sum_{k} |\langle X, x_k \rangle|^2 = \langle S(X), X \rangle
\tag{9.42}
\]

验证最后一个等号：

\[
\langle S(X), X \rangle = \left\langle \sum_{k} \langle X, x_k \rangle x_k, \; X \right\rangle
= \sum_{k} \langle X, x_k \rangle \langle x_k, X \rangle
= \sum_{k} \langle X, x_k \rangle \overline{\langle X, x_k \rangle}
= \sum_{k} |\langle X, x_k \rangle|^2
\]

因此 (9.42) 可以写为：

\[
\boxed{
\langle S(X), X \rangle = \|X\|^2
}
\tag{9.43}
\]

现在，需要证明 \( S = I \)（即 \( S(X) = X \) 对所有 \( X \) 成立）。由于 (9.43) 对所有 \( X \) 成立，可以利用**极化恒等式**来将内积用范数表示：

\[
\langle X, Y \rangle = \frac{1}{4} \left( \|X+Y\|^2 - \|X-Y\|^2 \right)
\tag{9.44}
\]

对任意 \( X, Y \)，将 (9.43) 应用于向量 \( X+Y \)：

\[
\langle S(X+Y), X+Y \rangle = \|X+Y\|^2
\]

由于 \( S \) 是线性的，展开左边：

\[
\langle S(X) + S(Y), X+Y \rangle = \|X+Y\|^2
\]

利用内积的双线性性质展开：

\[
\langle S(X), X \rangle + \langle S(X), Y \rangle + \langle S(Y), X \rangle + \langle S(Y), Y \rangle = \|X+Y\|^2
\]

同样地，对 \( X-Y \) 有：

\[
\langle S(X), X \rangle - \langle S(X), Y \rangle - \langle S(Y), X \rangle + \langle S(Y), Y \rangle = \|X-Y\|^2
\]

两式相减：

\[
2\langle S(X), Y \rangle + 2\langle S(Y), X \rangle = \|X+Y\|^2 - \|X-Y\|^2
\]

由极化恒等式 (9.44)，右边等于 \( 4\langle X, Y \rangle \)。因此：

\[
2\langle S(X), Y \rangle + 2\langle S(Y), X \rangle = 4\langle X, Y \rangle
\]

两边除以 2：

\[
\langle S(X), Y \rangle + \langle S(Y), X \rangle = 2\langle X, Y \rangle
\tag{9.45}
\]

又因为 \( S \) 是自伴算子，所以 \( \langle S(Y), X \rangle = \langle S(X), Y \rangle \)。代入 (9.45)：

\[
2\langle S(X), Y \rangle = 2\langle X, Y \rangle
\quad \Longrightarrow \quad \langle S(X), Y \rangle = \langle X, Y \rangle
\]

由于这对任意 \( Y \) 成立，因此 \( S(X) = X \) 对所有 \( X \) 成立。于是：

\[
\boxed{
X = \sum_{k} \langle X, x_k \rangle x_k
}
\tag{9.40}
\]

---

#### 2.6.1c 等价性结论

综合两个方向的证明，得到以下结论：

\[
\boxed{
X = \sum_{k} \langle X, x_k \rangle x_k \quad \Longleftrightarrow \quad \|X\|^2 = \sum_{k} |\langle X, x_k \rangle|^2
}
\tag{9.46}
\]

这个等价关系说明：**对于一个紧框架，重建公式和能量恒等式是一体两面——知道一个，就能推出另一个。**

这与标准正交基的性质完全一致。标准正交基是 \( A = B = 1 \) 的紧框架，而一般的紧框架是 \( A = B > 1 \) 的过完备表示。两者都保持 Parseval 等式，唯一的区别是框架引入了冗余，从而获得了稳定性和鲁棒性。

---

#### 2.6.2 维度

**规范正交基的维度表示**

设 \(\{e_k\}_{k=1}^n\) 是希尔伯特空间 \(H\) 的一组规范正交基，其中 \(n = \dim(H)\)。由于每个基矢量的范数满足 \(\|e_k\| = 1\)，因此

\[
\dim(H) - \sum_{k=1}^n \|e_k\|^2 = n - \sum_{k=1}^n 1 = n - n = 0,
\]

即

\[
\dim(H) = \sum_{k=1}^n \|e_k\|^2. \tag{9.47}
\]

这表明，在规范正交基下，空间的维度等于所有基矢量范数平方之和。这个结果看似平凡，但它揭示了一个重要事实：维数 \(n\) 可以通过一组基向量的"能量"（范数平方）来度量。

---

**框架下的推广问题**

自然要问：对于一般的框架 \(\{x_l\}_{l=1}^N\)，上述等式是否还成立？即，是否仍有

\[
\dim(H) = \sum_{l=1}^N \|x_l\|^2 ? \tag{9.48}
\]

为回答这个问题，借助规范正交基 \(\{e_k\}_{k=1}^n\) 作为分析工具，将框架向量 \(\{x_l\}\) 在这些基上展开，然后逐步推导。下面的推导将揭示框架情形与规范正交基情形的本质区别。

---

**基于规范正交基的展开与推导**

首先，根据帕塞瓦尔恒等式（Parseval's identity），对于每一个框架向量 \(x_l\)，它在规范正交基 \(\{e_k\}\) 下的展开系数满足

\[
\|x_l\|^2 = \sum_{k=1}^n |\langle x_l, e_k \rangle|^2. \tag{9.49}
\]

这个恒等式是对每一个 \(x_l\) 都严格成立的，因为 \(\{e_k\}\) 是完备的规范正交基。

将上式对所有 \(l\) 求和，得到

\[
\sum_{l=1}^N \|x_l\|^2
= \sum_{l=1}^N \sum_{k=1}^n |\langle x_l, e_k \rangle|^2.
\]

交换求和顺序，并利用内积的共轭对称性（模长平方不受共轭影响，即 \(|\langle x_l, e_k \rangle|^2 = |\langle e_k, x_l \rangle|^2\)），有

\[
\sum_{l=1}^N \|x_l\|^2
= \sum_{k=1}^n \sum_{l=1}^N |\langle e_k, x_l \rangle|^2. \tag{9.50}
\]

式 (9.50) 是一个恒等式，它在任何框架下都成立。然而，这个双重和 \(\sum_{k=1}^n \sum_{l=1}^N |\langle e_k, x_l \rangle|^2\) 是否等于 \(\dim(H)\)，则需要进一步考察。

---

**框架算子的引入与一般结论**

为分析式 (9.50) 中的双重和，引入框架算子（Frame Operator）\(S: H \to H\)，其定义为

\[
Sx = \sum_{l=1}^N \langle x, x_l \rangle x_l, \quad \forall x \in H.
\]

将 \(x = e_k\) 代入框架算子的定义，并与 \(e_k\) 做内积，得到

\[
\langle S e_k, e_k \rangle
= \left\langle \sum_{l=1}^N \langle e_k, x_l \rangle x_l, \; e_k \right\rangle
= \sum_{l=1}^N \langle e_k, x_l \rangle \langle x_l, e_k \rangle
= \sum_{l=1}^N |\langle e_k, x_l \rangle|^2.
\]

因此，式 (9.50) 中的双重和可以写成

\[
\sum_{k=1}^n \sum_{l=1}^N |\langle e_k, x_l \rangle|^2
= \sum_{k=1}^n \langle S e_k, e_k \rangle.
\]

由于 \(\{e_k\}_{k=1}^n\) 是规范正交基，上式右端正是框架算子 \(S\) 的迹（trace），即

\[
\sum_{k=1}^n \langle S e_k, e_k \rangle = \text{Tr}(S).
\]

结合式 (9.50)，得到一般框架下的重要结论：

\[
\boxed{\sum_{l=1}^N \|x_l\|^2 = \text{Tr}(S)}. \tag{9.51}
\]

式 (9.51) 表明：在一般框架下，框架向量的范数平方和并不直接等于空间的维度 \(\dim(H)\)，而是等于框架算子 \(S\) 的迹。这意味着式 (9.48) 所提出的等式在一般情况下不成立。

---

**帕塞瓦尔框架下的维度公式**

式 (9.48) 能够成立的特殊情况是：\(\dim(H) = \text{Tr}(I)\)，其中 \(I\) 是单位算子。因此，若要使

\[
\sum_{l=1}^N \|x_l\|^2 = \dim(H),
\]

结合式 (9.51)，就需要 \(\text{Tr}(S) = \text{Tr}(I)\)。对于正定的框架算子 \(S\)，迹相等并不自动推出 \(S = I\)。但是，如果框架是**紧框架**（Tight Frame），即 \(S = A I\)（其中 \(A\) 是框架边界），则 \(\text{Tr}(S) = A \dim(H)\)。此时，\(\sum_{l=1}^N \|x_l\|^2 = A \dim(H)\)。若进一步要求 \(A = 1\)，则框架称为**帕塞瓦尔框架**（Parseval Frame）。

对于帕塞瓦尔框架，有 \(S = I\)。将其代入上式，得到

\[
\sum_{k=1}^n \sum_{l=1}^N |\langle e_k, x_l \rangle|^2
= \sum_{k=1}^n \langle e_k, e_k \rangle
= \sum_{k=1}^n 1
= n
= \dim(H). \tag{9.52}
\]

结合式 (9.50) 与式 (9.52)，得到帕塞瓦尔框架下的完整维度公式：

\[
\boxed{
\dim(H)
= \sum_{k=1}^n \sum_{l=1}^N |\langle e_k, x_l \rangle|^2
= \sum_{l=1}^N \sum_{k=1}^n |\langle x_l, e_k \rangle|^2
= \sum_{l=1}^N \|x_l\|^2
}. \tag{9.53}
\]

式 (9.53) 完整地给出了帕塞瓦尔框架下维度与框架向量范数平方和之间的等价关系。它表明，尽管帕塞瓦尔框架可能是过完备的（即框架向量个数 \(N\) 可能大于维度 \(n\)），但其所有向量的范数平方之和恰好等于空间的维度。这个结果在形式上与规范正交基的式 (9.47) 完全一致，因此，帕塞瓦尔框架是规范正交基在过完备情形下的自然推广。

---

**一般框架与帕塞瓦尔框架的对比总结**

一般框架与帕塞瓦尔框架的结论对比如下：

| 框架类型 | 维度与范数平方和的关系 | 成立条件 |
| :--- | :--- | :--- |
| 规范正交基 | \(\dim(H) = \sum_{k=1}^n \|e_k\|^2\) | 恒成立（\(\|e_k\|=1\)） |
| 一般框架 | \(\sum_{l=1}^N \|x_l\|^2 = \text{Tr}(S)\) | 恒成立，但一般不等于 \(\dim(H)\) |
| 帕塞瓦尔框架 | \(\dim(H) = \sum_{l=1}^N \|x_l\|^2\) | 紧框架且框架边界 \(A = 1\)（即 \(S = I\)） |

特别地，对于一般框架，还可以利用框架边界 \(A\) 和 \(B\) 给出 \(\sum_{l=1}^N \|x_l\|^2\) 与 \(\dim(H)\) 之间的上下界关系。由于 \(A I \le S \le B I\)，两边取迹，得到

\[
A \dim(H) \le \sum_{l=1}^N \|x_l\|^2 \le B \dim(H). \tag{9.54}
\]

式 (9.54) 表明，即使是一般框架，其向量的范数平方和也与维度保持在一个由框架边界确定的线性范围内。当 \(A = B = 1\) 时，上下界同时等于 \(\dim(H)\)，式 (9.54) 便退化为式 (9.53) 的等号情形。

---

**回到最初的推导**：原本试图证明

\[
\dim(H)
= \sum_{k=1}^n \sum_{l=1}^N |\langle e_k, x_l \rangle|^2
= \sum_{l=1}^N \sum_{k=1}^n |\langle x_l, e_k \rangle|^2
= \sum_{l=1}^N \|x_l\|^2.
\]

经过上面的分析，可以看到：这个证明链条中的**第一个等号**（即 \(\dim(H) = \sum_{k=1}^n \sum_{l=1}^N |\langle e_k, x_l \rangle|^2\)）在一般框架下并不成立。它仅在帕塞瓦尔框架（\(S = I\)）下才成立。而后面的两个等号（交换求和顺序和帕塞瓦尔恒等式）则是恒成立的。因此，整个证明链条是**帕塞瓦尔框架下的特殊结论**，而非一般框架下的普适结论。对于一般框架，正确的结论是式 (9.51) 和式 (9.54)。
#### 2.6.3 线性变换迹的框架表示

在规范正交基下已经熟知，任意线性变换 \(B: H \to H\) 的迹（trace）可以表示为

\[
\operatorname{Tr}(B) = \sum_{k=1}^{n} \langle B e_k, e_k \rangle, \tag{9.55}
\]

其中 \(\{e_k\}_{k=1}^n\) 是 \(H\) 的任意一组规范正交基。这个定义是基无关的，但形式上依赖于基的选择。自然要问：如果空间中的向量用框架 \(\{x_l\}_{l=1}^N\) 来表示，那么迹是否也能用这些框架向量直接表达？本节首先推导 Parseval 框架下的简洁表示，然后推广至一般框架。

---

**Parseval 框架下的迹公式及其详细推导**

设 \(\{x_l\}_{l=1}^N\) 是 \(H\) 的一个 Parseval 框架，即其框架算子 \(S = I\)。根据 Parseval 框架的重构公式，对任意 \(x \in H\)，有

\[
x = \sum_{l=1}^{N} \langle x, x_l \rangle \, x_l, \tag{9.56}
\]

等价地，用外积表示即为 \(\sum_{l=1}^N |x_l\rangle\langle x_l| = I\)。期望证明：

\[
\operatorname{Tr}(B) = \sum_{l=1}^{N} \langle B x_l, x_l \rangle, \tag{9.57}
\]

但在复数域上，更准确的形式是 \(\operatorname{Tr}(B) = \sum_{l=1}^N \langle x_l, B x_l \rangle\)。这是因为迹的循环性质：\(\operatorname{Tr}(B) = \operatorname{Tr}(B I) = \sum_l \operatorname{Tr}(B |x_l\rangle\langle x_l|) = \sum_l \langle x_l, B x_l \rangle\)。若 \(B\) 是自伴算子（或工作于实内积空间），则 \(\langle B x_l, x_l \rangle = \langle x_l, B x_l \rangle\)，两种写法等价。下面以 \(\langle B x_l, x_l \rangle\) 形式给出结论，并采用规范正交基逐步展开计算，以验证这一关系。

**第一步：用规范正交基表示迹**

由式 (9.55) 可知

\[
\operatorname{Tr}(B) = \sum_{k=1}^{n} \langle B e_k, e_k \rangle. \tag{9.58}
\]

**第二步：利用 Parseval 框架展开基矢量**

对每一个基矢量 \(e_k\)，应用式 (9.56) 的重构公式，得到

\[
e_k = \sum_{l=1}^{N} \langle e_k, x_l \rangle \, x_l. \tag{9.59}
\]

将式 (9.59) 代入式 (9.58) 中的 \(B e_k\)（注意线性性），有

\[
\langle B e_k, e_k \rangle
= \left\langle B \left( \sum_{l=1}^{N} \langle e_k, x_l \rangle \, x_l \right), \; e_k \right\rangle
= \sum_{l=1}^{N} \langle e_k, x_l \rangle \, \langle B x_l, e_k \rangle. \tag{9.60}
\]

**第三步：对所有 \(k\) 求和并交换求和顺序**

将式 (9.60) 左右两边对 \(k\) 从 1 到 \(n\) 求和：

\[
\operatorname{Tr}(B)
= \sum_{k=1}^{n} \sum_{l=1}^{N} \langle e_k, x_l \rangle \, \langle B x_l, e_k \rangle
= \sum_{l=1}^{N} \sum_{k=1}^{n} \langle e_k, x_l \rangle \, \langle B x_l, e_k \rangle. \tag{9.61}
\]

**第四步：化简内积乘积**

对于固定的 \(l\)，有

\[
\sum_{k=1}^{n} \langle e_k, x_l \rangle \, \langle B x_l, e_k \rangle.
\]

利用内积的共轭对称性，\(\langle e_k, x_l \rangle = \overline{\langle x_l, e_k \rangle}\)。若期望得到 \(\langle B x_l, x_l \rangle\)，则需将上述和化为

\[
\sum_{k=1}^{n} \langle B x_l, e_k \rangle \, \langle e_k, x_l \rangle.
\]

而 \(\langle B x_l, x_l \rangle = \left\langle B x_l, \sum_{k=1}^{n} \langle x_l, e_k \rangle e_k \right\rangle = \sum_{k=1}^{n} \langle x_l, e_k \rangle \, \langle B x_l, e_k \rangle\)。比较可知，如果 \(\langle e_k, x_l \rangle = \langle x_l, e_k \rangle\)（即内积为实内积，或内积是对称的），则两者一致。在复数域中，正确的无歧义表达式应为

\[
\operatorname{Tr}(B) = \sum_{l=1}^{N} \langle x_l, B x_l \rangle. \tag{9.62}
\]

为得到式 (9.62)，可从式 (9.61) 出发，将 \(\langle B x_l, e_k \rangle\) 写成 \(\langle e_k, B^* x_l \rangle^*\)（其中 \(B^*\) 是 \(B\) 的伴随），然后利用 Parseval 恒等式。但更简洁的推导是利用外积的迹循环性质：

\[
I = \sum_{l=1}^{N} |x_l\rangle\langle x_l|,
\]

于是

\[
\operatorname{Tr}(B) = \operatorname{Tr}(B I) = \sum_{l=1}^{N} \operatorname{Tr}\bigl(B |x_l\rangle\langle x_l|\bigr).
\]

对于任意向量 \(u,v\)，有 \(\operatorname{Tr}(B |u\rangle\langle v|) = \langle v, B u \rangle\)。因此

\[
\operatorname{Tr}(B) = \sum_{l=1}^{N} \langle x_l, B x_l \rangle. \tag{9.63}
\]

若 \(B\) 是自伴的，则 \(\langle x_l, B x_l \rangle = \langle B x_l, x_l \rangle\)，从而得到前述形式。因此，在实空间或自伴算子的情形下，Parseval 框架下的迹公式可写为

\[
\boxed{
\operatorname{Tr}(B) = \sum_{l=1}^{N} \langle B x_l, x_l \rangle
} \quad (\text{自伴或实内积}). \tag{9.64}
\]

对于一般复数域，使用式 (9.63) 的规范形式。

---

**一般框架下的迹表示**

将上述结论推广至一般框架 \(\{x_l\}_{l=1}^N\)，其框架算子为 \(S\)，边界为 \(0 < A \le B < \infty\)。此时 Parseval 的重构公式不再成立，需要引入**对偶框架** \(\{\tilde{x}_l\}_{l=1}^N\)。对偶框架满足重构公式

\[
x = \sum_{l=1}^{N} \langle x, \tilde{x}_l \rangle x_l = \sum_{l=1}^{N} \langle x, x_l \rangle \tilde{x}_l, \quad \forall x \in H. \tag{9.65}
\]

在算子语言中，这意味着

\[
\sum_{l=1}^{N} |x_l\rangle\langle \tilde{x}_l| = I, \quad \text{以及} \quad \sum_{l=1}^{N} |\tilde{x}_l\rangle\langle x_l| = I. \tag{9.66}
\]

利用与 Parseval 框架完全相同的迹循环论证，将 \(I\) 表示为 \(\sum_l |x_l\rangle\langle \tilde{x}_l|\)，有

\[
\operatorname{Tr}(B) = \operatorname{Tr}\left(B \sum_{l=1}^{N} |x_l\rangle\langle \tilde{x}_l| \right)
= \sum_{l=1}^{N} \operatorname{Tr}\bigl(B |x_l\rangle\langle \tilde{x}_l|\bigr)
= \sum_{l=1}^{N} \langle \tilde{x}_l, B x_l \rangle. \tag{9.67}
\]

同样，利用另一个分解 \(I = \sum_l |\tilde{x}_l\rangle\langle x_l|\)，可得

\[
\operatorname{Tr}(B) = \sum_{l=1}^{N} \langle x_l, B \tilde{x}_l \rangle. \tag{9.68}
\]

式 (9.67) 和 (9.68) 是**一般框架下的迹公式**。它们表明，迹可以用框架向量及其对偶框架的内积组合来表示。当 \(\{x_l\}\) 本身就是 Parseval 框架时，对偶框架等于自身（\(\tilde{x}_l = x_l\)），式 (9.67) 退化为 \(\operatorname{Tr}(B) = \sum_l \langle x_l, B x_l \rangle\)，与式 (9.63) 一致。

---

**紧框架（非 Parseval）的特殊情形**

若框架是紧框架，即 \(S = A I\)（\(A > 0\)），则其对偶框架为 \(\tilde{x}_l = \frac{1}{A} x_l\)。代入式 (9.67)，得到

\[
\operatorname{Tr}(B) = \sum_{l=1}^{N} \left\langle \frac{1}{A} x_l, B x_l \right\rangle
= \frac{1}{A} \sum_{l=1}^{N} \langle x_l, B x_l \rangle. \tag{9.69}
\]

等价地（若 \(B\) 自伴），

\[
\operatorname{Tr}(B) = \frac{1}{A} \sum_{l=1}^{N} \langle B x_l, x_l \rangle. \tag{9.70}
\]

当 \(A = 1\) 时，即 Parseval 框架，式 (9.70) 还原为式 (9.64)。

---

**总结与对比**

三种情况下的迹公式汇总如下：

| 框架类型 | 迹的表达式 | 条件 |
| :--- | :--- | :--- |
| 规范正交基 | \(\operatorname{Tr}(B) = \sum_{k=1}^{n} \langle B e_k, e_k \rangle\) | 恒成立（任意 \(B\)） |
| Parseval 框架 | \(\operatorname{Tr}(B) = \sum_{l=1}^{N} \langle x_l, B x_l \rangle\) | 复数域通用；自伴时 \(\langle B x_l, x_l \rangle\) |
| 紧框架 (\(S = A I\)) | \(\operatorname{Tr}(B) = \frac{1}{A} \sum_{l=1}^{N} \langle x_l, B x_l \rangle\) | 同上 |
| 一般框架 | \(\operatorname{Tr}(B) = \sum_{l=1}^{N} \langle \tilde{x}_l, B x_l \rangle\) | \(\{\tilde{x}_l\}\) 是对偶框架 |

---

**关键结论**

对于一般框架，迹不能简单地用框架向量自身的内积和来表示，必须引入对偶框架。只有当框架为紧框架时，才可以通过一个常数因子 \(1/A\) 进行修正。特别地，当 \(A = 1\)（Parseval 框架）时，修正因子消失，得到最简洁的形式。这一结果与上一节中维度公式的讨论一脉相承：**Parseval 框架在保持迹的简单代数形式方面具有与规范正交基相同的地位**。


#### 2.6.4 三种基本性质（重构、维度、迹）的框架对比

在前三小节中，分别研究了**框架算子与重构（2.6.1）**、**维度公式（2.6.2）** 以及**线性变换的迹（2.6.3）**。为更清晰地揭示从**规范正交基**到**Parseval框架**，再到**紧框架**和**一般框架**的推广脉络，本节将上述三个核心性质并列对比。

下表汇总了四种表示系统在重构公式、框架算子、对偶框架、维度计算、迹计算以及冗余性方面的完整对比：

| 对比性质 | **规范正交基 (ONB)** | **Parseval 框架** | **紧框架 (\(S = A I\))** | **一般框架** |
| :--- | :--- | :--- | :--- | :--- |
| **框架算子 \(S\)** | \(S = I\) | \(S = I\) | \(S = A I \quad (A>0)\) | \(S = \displaystyle\sum_{l=1}^{N} \langle \cdot, x_l \rangle x_l\) |
| **对偶框架 \(\tilde{x}_l\)** | \(\tilde{x}_l = e_l\)（对偶即自身） | \(\tilde{x}_l = x_l\)（对偶即自身） | \(\tilde{x}_l = \dfrac{1}{A} x_l\)（对偶为自身的缩放） | \(\tilde{x}_l = S^{-1} x_l\)（需通过逆算子求解） |
| **2.6.1 重构公式** | \(x = \displaystyle\sum_{k=1}^{n} \langle x, e_k \rangle e_k\) | \(x = \displaystyle\sum_{l=1}^{N} \langle x, x_l \rangle x_l\) | \(x = \dfrac{1}{A} \displaystyle\sum_{l=1}^{N} \langle x, x_l \rangle x_l\) | \(x = \displaystyle\sum_{l=1}^{N} \langle x, \tilde{x}_l \rangle x_l = \displaystyle\sum_{l=1}^{N} \langle x, x_l \rangle \tilde{x}_l\) |
| **2.6.2 维度公式**（范数平方和） | \(\displaystyle\sum_{k=1}^{n} \|e_k\|^2 = n = \dim(H)\) | \(\displaystyle\sum_{l=1}^{N} \|x_l\|^2 = n = \dim(H)\) | \(\displaystyle\sum_{l=1}^{N} \|x_l\|^2 = A \cdot \dim(H)\) | \(\displaystyle\sum_{l=1}^{N} \|x_l\|^2 = \operatorname{Tr}(S) \quad (\text{一般 } \neq \dim(H))\) |
| **2.6.3 迹公式**（线性算子 \(B\)） | \(\operatorname{Tr}(B) = \displaystyle\sum_{k=1}^{n} \langle B e_k, e_k \rangle\) | \(\operatorname{Tr}(B) = \displaystyle\sum_{l=1}^{N} \langle x_l, B x_l \rangle\)<br>（复数域通用形式） | \(\operatorname{Tr}(B) = \dfrac{1}{A} \displaystyle\sum_{l=1}^{N} \langle x_l, B x_l \rangle\)<br>（复数域通用形式） | \(\operatorname{Tr}(B) = \displaystyle\sum_{l=1}^{N} \langle \tilde{x}_l, B x_l \rangle\)<br>（需引入对偶框架） |
| **框架边界** | \(A = B = 1\) | \(A = B = 1\) | \(A = B\)（上界等于下界） | \(0 < A \le B < \infty\)（上下界可以不等） |
| **冗余度（过完备性）** | 无冗余，\(N = n\) | 允许过完备，\(N \ge n\) | 允许过完备，\(N \ge n\) | 允许过完备，\(N \ge n\) |

---

**主要观察与结论**

1.  **Parseval框架的特殊地位**：从表中可以清楚地看到，Parseval框架在**重构公式**、**维度公式**和**迹公式**三个核心性质上，与规范正交基保持了完全一致的形式。唯一的区别在于它允许过完备（\(N \ge n\)）。因此，Parseval框架可以视为规范正交基在过完备情形下的**最自然推广**。

2.  **紧框架的缩放修正**：当框架为紧框架（\(S = A I\)）但 \(A \neq 1\) 时，重构公式和迹公式中均多出一个因子 \(1/A\)，而维度公式中则多出因子 \(A\)。这表明紧框架在代数结构上与规范正交基只差一个**全局缩放因子**，其所有性质都可以通过该因子统一修正。

3.  **一般框架的本质区别**：对于一般框架，由于其框架算子 \(S\) 不是单位算子的倍数，**重构**和**迹**的计算不再能仅用原框架向量本身完成，必须引入**对偶框架** \(\{\tilde{x}_l\}\)。特别是迹公式，从 \(\sum \langle x_l, B x_l \rangle\) 变成了 \(\sum \langle \tilde{x}_l, B x_l \rangle\)，体现了对偶框架在一般线性运算中的核心作用。此外，维度公式中的范数平方和等于 \(\operatorname{Tr}(S)\) 而非 \(\dim(H)\)，说明此时"框架向量的总能量"偏离了空间维度。

4.  **推广脉络总结**：
    - **规范正交基**：\(S = I\)，对偶 = 自身，形式最简洁。
    - **Parseval框架**：\(S = I\)，对偶 = 自身，允许过完备，简洁性保持不变。
    - **紧框架**：\(S = A I\)，对偶 = 自身缩放，简洁性略有损失（引入常数因子）。
    - **一般框架**：\(S \neq c I\)，对偶 \(\neq\) 自身，必须借助对偶框架才能表达重构与迹。


## 3. 过完备字典

在前两节中，从规范正交基出发，逐步推广到 Parseval 框架和一般框架。Parseval 框架虽然允许过完备（即向量个数 \(N\) 可以大于空间维度 \(n\)），但在重构、维度和迹等核心性质上仍然保持了与规范正交基完全一致的形式。这一节进一步放宽约束，进入**过完备字典**的一般框架。此时，框架算子 \(S\) 不再等于单位算子，甚至不再是单位算子的倍数。本节揭示这种一般性带来的**优势**与**代价**，并引入**对偶框架**（或称互补基）的概念。

---

### 3.1 推广路线：标准正交基 → Parseval 框架 → 过完备字典

这条推广路线的逻辑如下：

- **规范正交基**：\(\{e_k\}_{k=1}^n\)，\(N = n\)，\(S = I\)，重构 \(x = \sum_{k=1}^n \langle x, e_k \rangle e_k\)。
- **Parseval 框架**：\(\{x_l\}_{l=1}^N\)，\(N \ge n\)，但仍有 \(S = I\)，重构 \(x = \sum_{l=1}^N \langle x, x_l \rangle x_l\)。这是过完备但保持简单性的特例。
- **一般框架（过完备字典）**：\(\{x_l\}_{l=1}^N\)，\(N \gg n\)，框架算子 \(S = \sum_{l=1}^N |x_l\rangle\langle x_l| \neq c I\)。此时重构不能直接用自身完成，而需要借助对偶框架。

进一步推广到过完备字典，主要有以下原因：

1. **冗余带来健壮性**：当系统受到噪声、干扰或数据丢失时，过完备表示能够利用冗余信息进行恢复。例如，在无线通信中，如果信号通过两个独立的信道传输，每个信道的误差会叠加；而如果只使用一组过完备的表示进行传输，误差仅影响一次，误差不会被放大。
   
2. **高精度**：过完备字典可以更灵活地逼近信号的稀疏结构或非线性特征，从而获得更高的近似精度。

3. **牺牲的是计算方便性**：因为重构需要求解对偶框架或进行矩阵求逆，运算量远大于正交基下的简单内积。

这些权衡让人想起经济学中的"蒙代尔不可能三角"——稳定性、精度和计算方便性三者难以同时满足。规范正交基在计算上最方便，但在某些应用场景下稳定性或精度不足；而过完备字典则在稳定性和精度上表现出色，却以计算复杂度为代价。

---

### 3.2 一般框架下的分析算子、合成算子及其复合

为推导一般框架下的重构公式，引入三个基本算子，它们是框架理论中的核心工具。

设 \(\{x_k\}_{k=1}^N\) 是希尔伯特空间 \(H\)（\(\dim H = n\)）中的一组向量（构成框架，但一般不满足 Parseval 条件）。定义一个**分析算子**（Analysis Operator）\(\Theta: H \to \mathbb{C}^N\)，它将向量 \(X \in H\) 映射为其在所有框架向量上的内积序列：

\[
\Theta X = \bigl( \langle X, x_1 \rangle, \langle X, x_2 \rangle, \ldots, \langle X, x_N \rangle \bigr)^T. \tag{9.71}
\]

分析算子的作用相当于对信号 \(X\) 进行"特征提取"，将其转换为一个 \(N\) 维的系数向量。

与之对应的，定义**合成算子**（Synthesis Operator）\(\Theta^*: \mathbb{C}^N \to H\)，它是 \(\Theta\) 的伴随算子。对于任意系数向量 \((c_1, \ldots, c_N)^T\)，合成算子将其组合为

\[
\Theta^* (c_1, \ldots, c_N)^T = \sum_{k=1}^N c_k \, x_k. \tag{9.72}
\]

特别地，当系数向量是标准基向量 \(\tilde{e}_k\)（即第 \(k\) 个分量为 1，其余为 0）时，有

\[
\Theta^* \tilde{e}_k = x_k. \tag{9.73}
\]

其中 \(\{\tilde{e}_k\}_{k=1}^N\) 是 \(\mathbb{C}^N\) 中的标准正交基。这个关系表明，合成算子将单位向量映射为对应的框架向量。

考虑复合算子 \(\Theta^* \Theta: H \to H\)。对于任意 \(X \in H\)，首先由 \(\Theta\) 得到系数序列，然后用合成算子重建：

\[
\Theta^* \Theta X = \Theta^* \bigl( \langle X, x_1 \rangle, \ldots, \langle X, x_N \rangle \bigr)^T
= \sum_{k=1}^N \langle X, x_k \rangle \, x_k. \tag{9.74}
\]

这正是框架算子 \(S\) 的定义，因此有

\[
S = \Theta^* \Theta. \tag{9.75}
\]

Parseval 框架对应于 \(S = I\)，此时重构即为 \(X = \sum_{k=1}^N \langle X, x_k \rangle x_k\)。如果 Parseval 条件不成立，则 \(S \neq I\)，此时式 (9.74) 的重建结果不等于 \(X\)，而是 \(S X\)。为恢复 \(X\)，只需将上述结果左乘 \(S^{-1}\)（假设框架是紧的或一般的可逆框架）：

\[
X = S^{-1} (S X) = S^{-1} \left( \sum_{k=1}^N \langle X, x_k \rangle \, x_k \right)
= \sum_{k=1}^N \langle X, x_k \rangle \, S^{-1} x_k. \tag{9.76}
\]

令

\[
y_k = S^{-1} x_k, \quad k = 1, 2, \ldots, N, \tag{9.77}
\]

则式 (9.76) 可以写成

\[
\boxed{X = \sum_{k=1}^N \langle X, x_k \rangle \, y_k}. \tag{9.78}
\]

式 (9.78) 给出了过完备字典下的重构公式：**用 \(x_k\) 计算系数（内积），用 \(y_k\) 作为重建原子**。这里 \(\{y_k\}_{k=1}^N\) 就是 \(\{x_k\}\) 的**对偶框架**（也称为正则对偶）。对偶框架 \(\{y_k\}\) 本身也是 \(H\) 的一个框架。

---

### 3.3 对偶重构公式的推导

通过 \(y_k = S^{-1} x_k\) 定义了对偶框架，并得到了重构公式 \(X = \sum_{k=1}^N \langle X, x_k \rangle y_k\)（见式 (9.78)）。以下证明**对偶框架反过来也成立**，即信号也可以用 \(y_k\) 计算系数、用 \(x_k\) 来重建：

\[
X = \sum_{k=1}^{N} \langle X, y_k \rangle x_k. \tag{9.79}
\]

这个关系并不是显然的，因为它涉及到系数计算原子与重建原子的互换。为得到这一结果，从一个恒等式出发：

\[
X = S \left( S^{-1} X \right). \tag{9.80}
\]

式 (9.80) 恒成立，因为 \(S^{-1} S = I\)，无论 \(X\) 取什么值都成立。将等号右边的 \(S\) 算子展开成显式的求和形式。

---

**第一步：将框架算子 \(S\) 作用于向量 \(S^{-1}X\)**

根据框架算子 \(S\) 的定义，对任意向量 \(Z \in H\)，有

\[
S Z = \sum_{k=1}^{N} \langle Z, x_k \rangle x_k. \tag{9.81}
\]

令式 (9.81) 中的 \(Z\) 取特定的向量 \(S^{-1}X\)。代入后得到

\[
S \left( S^{-1} X \right) = \sum_{k=1}^{N} \left\langle S^{-1}X, \; x_k \right\rangle \, x_k. \tag{9.82}
\]

结合式 (9.80) 和式 (9.82)，立即得到

\[
X = \sum_{k=1}^{N} \left\langle S^{-1}X, \; x_k \right\rangle \, x_k. \tag{9.83}
\]

目前只是把恒等式换了一种写法，但系数中仍然含有 \(S^{-1}\)，形式还不够简洁。

---

**第二步：利用 \(S^{-1}\) 的自伴性将算子从左边挪到右边**

由于框架算子 \(S = \Theta^*\Theta\) 是正定且自伴的（即 \(S = S^*\)），其逆算子 \(S^{-1}\) 也必然是自伴的，即

\[
(S^{-1})^* = S^{-1}. \tag{9.84}
\]

自伴性允许将内积中的线性算子从第一个位置"挪"到第二个位置，且保持内积结果不变。具体地，对于任意向量 \(X\) 和 \(x_k\)，有

\[
\left\langle S^{-1}X, \; x_k \right\rangle
= \left\langle X, \; S^{-1} x_k \right\rangle. \tag{9.85}
\]

这是一个标准的内积性质：\(\langle A^* u, v \rangle = \langle u, A v \rangle\)，当 \(A^* = A\) 时即为此形式。

---

**第三步：代入对偶框架的定义 \(y_k = S^{-1} x_k\)**

根据式 (9.77) 中对对偶框架原子的定义：

\[
y_k = S^{-1} x_k. \tag{9.86}
\]

将式 (9.86) 代入式 (9.85) 的右端，得到

\[
\left\langle S^{-1}X, \; x_k \right\rangle = \langle X, \; y_k \rangle. \tag{9.87}
\]

---

**第四步：将式 (9.87) 回代到式 (9.83) 中**

将式 (9.87) 的结果替换式 (9.83) 中的内积项，得到

\[
X = \sum_{k=1}^{N} \langle X, y_k \rangle \, x_k.
\]

这正是要证明的公式。

---

**完整推导过程的连贯展示**

将上述四步合并成一组连贯的等式链：

\[
\boxed{
\begin{aligned}
X &= S \left( S^{-1} X \right) && \text{(平凡的恒等变形)} \\
&= \sum_{k=1}^{N} \left\langle S^{-1}X, \; x_k \right\rangle \, x_k && \text{(展开框架算子的定义)} \\
&= \sum_{k=1}^{N} \left\langle X, \; S^{-1} x_k \right\rangle \, x_k && \text{(利用 }S^{-1}\text{ 的自伴性)} \\
&= \sum_{k=1}^{N} \langle X, \; y_k \rangle \, x_k. && \text{(代入对偶定义 } y_k = S^{-1}x_k \text{)}
\end{aligned}
} \tag{9.88}
\]

式 (9.88) 给出了互换系数计算和重建原子的重构公式。它和式 (9.78) 共同构成了互补基的完整图像：

- 用 \(x_k\) 算系数，用 \(y_k\) 重建：\(X = \sum \langle X, x_k \rangle y_k\)；
- 用 \(y_k\) 算系数，用 \(x_k\) 重建：\(X = \sum \langle X, y_k \rangle x_k\)。

因此，\(\{x_k\}\) 和 \(\{y_k\}\) 互为对偶。在 OFDM 通信中，这就像是用一组已知的导频波形（对应 \(x_k\)）去接收信号并估计信道系数，然后在接收端利用信道响应的逆（对应 \(S^{-1}\)）计算出另一组等效波形（对应 \(y_k\)），从而重构出发送的原始符号，实现稳健的信号恢复。

---

### 3.4 互补基的概念

从式 (9.78) 和 (9.88) 可以看到，两个框架 \(\{x_k\}\) 和 \(\{y_k\}\) 互为对偶，它们配合起来实现了完全重构。具体来说：

- **系数计算**：用 \(x_k\) 与信号做内积得到系数 \(\langle X, x_k \rangle\)；
- **信号重建**：用这些系数乘以 \(y_k\) 并求和。

或者反过来：

- 用 \(y_k\) 计算系数 \(\langle X, y_k \rangle\)；
- 用 \(x_k\) 重建。

这种关系被称为 **"互补基"**（Complementary Bases）或对偶框架。虽然 \(\{x_k\}\) 和 \(\{y_k\}\) 各自可能都不是正交基，甚至不是紧框架，但它们之间满足双正交关系：

\[
\langle x_j, y_k \rangle = \delta_{jk} \quad \text{?}
\]

实际上，对于过完备框架，这个双正交关系并不一定以简单的克罗内克形式成立，因为框架是冗余的，\(N > n\)。更准确地说，有两组向量，它们满足重构恒等式，但并不要求两两正交。互补基的含义是：**它们彼此充当对方的对偶，共同构成一个完备的表示系统，允许信号在这两个系统之间进行无损的系数-信号转换**。

在信号处理中，一个典型的应用是**正交频分复用（OFDM）**。在 OFDM 系统中，使用一组子载波（通常为正弦波）作为基函数来传输数据。然而，实际信道会造成频率偏移和符号间干扰。为对抗这些干扰，OFDM 常常采用循环前缀和冗余的导频结构。可以将发送端使用的调制波形视为一个过完备字典，而接收端的均衡和信道估计则对应于寻找对偶框架，以从畸变的接收信号中恢复原始数据。例如，在存在多径效应的信道中，接收信号是发送波形经过线性系统后的叠加，此时若将信道响应纳入框架设计，则发送的符号系数（如 QAM 符号）可以看作是在某个过完备字典上的展开系数，接收端则通过已知的导频估计对偶波形，从而实现稳健的重构。虽然 OFDM 通常不直接使用显式的框架公式，但其核心思想——**利用冗余信息（循环前缀、导频）提高鲁棒性，同时需要额外的计算（信道估计、均衡）**——与过完备字典的权衡完全一致。

---

### 3.5 一般框架下的重构步骤总结

过完备字典下的重构过程总结为以下步骤：

1. **给定**：一组过完备字典 \(\{x_k\}_{k=1}^N\)（可能不是 Parseval 框架），以及信号 \(X\)。
2. **计算框架算子**：\(S = \sum_{k=1}^N |x_k\rangle\langle x_k|\)（或等价地，\(S X = \sum_{k=1}^N \langle X, x_k \rangle x_k\)）。
3. **求逆**：若框架是紧的（\(S = A I\)），则 \(S^{-1} = \frac{1}{A} I\)；否则需计算 \(S^{-1}\)（可能通过矩阵求逆或迭代方法）。
4. **构造对偶框架**：\(y_k = S^{-1} x_k\)。
5. **重构**：用 \(X = \sum_{k=1}^N \langle X, x_k \rangle \, y_k\) 或 \(X = \sum_{k=1}^N \langle X, y_k \rangle \, x_k\)。

如果框架是 Parseval 的，则 \(S = I\)，\(y_k = x_k\)，重构退化为最简单的形式。

---

### 3.6 过完备字典：优势与代价

以下表格总结三种系统的特点：

| 系统类型 | 框架算子 \(S\) | 重构公式 | 冗余度 | 计算复杂度 | 稳健性/精度 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| 规范正交基 | \(I\) | \(X = \sum_k \langle X, e_k \rangle e_k\) | \(N = n\) | 最低 | 一般 |
| Parseval 框架 | \(I\) | \(X = \sum_l \langle X, x_l \rangle x_l\) | \(N \ge n\) | 低 | 较好 |
| 紧框架 (\(S = A I\)) | \(A I\) | \(X = \frac{1}{A} \sum_l \langle X, x_l \rangle x_l\) | \(N \ge n\) | 中（乘常数） | 较好 |
| 一般过完备字典 | \(S \neq cI\) | \(X = \sum_l \langle X, x_l \rangle (S^{-1} x_l)\) | \(N \gg n\) | 高（求逆） | 最高 |

通过这个对比，可以直观地理解为何过完备字典被广泛用于压缩感知、雷达信号处理、无线通信等对稳定性和精度要求较高的领域，尽管其代价是计算上的挑战。这就是"蒙代尔不可能三角"在信号处理中的体现：**稳健性、精度和计算方便性三者不可兼得**。

---


## 4. 课后总结

### 4.1 核心逻辑链：从正交基到过完备表示的泛化

本讲以"基的泛化"为主线，从最熟悉的标准正交基出发，逐步放宽条件，最终到达过完备字典这一最通用的表示框架。核心逻辑如下：

1. **标准正交基**：$\{e_k\}_{k=1}^n$ 满足 $\langle e_i, e_j\rangle = \delta_{ij}$，重构公式为 $X = \sum_k \langle X, e_k\rangle e_k$。这是所有信号表示理论的起点。

2. **框架**：将标准正交基的"正交性"和"完备性"两个要求放宽为一组不等式 $A\|X\|^2 \le \sum_k |\langle X, x_k\rangle|^2 \le B\|X\|^2$。框架允许 $N \ge n$（冗余）且不要求原子正交。

3. **Parseval 框架**：当 $A = B = 1$ 时，框架算子 $S = I$，重构公式与标准正交基完全一致——$X = \sum_k \langle X, x_k\rangle x_k$。这是"最接近正交基"的框架。

4. **对偶框架与重构**：一般框架下 $S \neq cI$，需要通过对偶框架 $\{y_k = S^{-1}x_k\}$ 实现重构：$X = \sum_k \langle X, x_k\rangle y_k = \sum_k \langle X, y_k\rangle x_k$。$\{x_k\}$ 与 $\{y_k\}$ 互为对偶。

5. **过完备字典的权衡**：冗余度越高，表示越稳健，但计算复杂度也越高。这是信号处理中的"不可能三角"——正交基计算方便但缺乏冗余；过完备字典稳健性强但需要求 $S^{-1}$。

### 4.2 三种表示系统对比

| 系统类型 | 框架算子 $S$ | 重构公式 | 冗余度 | 计算复杂度 | 稳健性 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| 标准正交基 | $I$ | $X = \sum \langle X, e_k\rangle e_k$ | $N = n$ | 最低 | 一般 |
| Parseval 框架 | $I$ | $X = \sum \langle X, x_k\rangle x_k$ | $N \ge n$ | 低 | 较好 |
| 紧框架 | $AI$ | $X = \frac{1}{A}\sum \langle X, x_k\rangle x_k$ | $N \ge n$ | 中 | 较好 |
| 一般过完备字典 | $S \neq cI$ | $X = \sum \langle X, x_k\rangle S^{-1}x_k$ | $N \gg n$ | 高 | 最高 |

### 4.3 重点公式

**框架条件：**
$$
\boxed{A\|X\|^2 \le \sum_{k=1}^N |\langle X, x_k\rangle|^2 \le B\|X\|^2}
$$

**框架算子与对偶：**
$$
\boxed{S = \sum_{k=1}^N |x_k\rangle\langle x_k|, \quad y_k = S^{-1}x_k}
$$

**重构公式（对偶形式）：**
$$
\boxed{X = \sum_{k=1}^N \langle X, x_k\rangle y_k = \sum_{k=1}^N \langle X, y_k\rangle x_k}
$$


## 5. 学习检查清单：自测核心知识点掌握情况

- [ ] 能写出框架的定义（上下界不等式），解释 $A, B$ 的物理含义
- [ ] 能区分标准正交基、Parseval 框架、紧框架、一般框架之间的包含关系
- [ ] 能写出框架算子 $S$ 的定义，并解释 $S$ 为什么是正定算子
- [ ] 能说明 Parseval 框架的充要条件：$S = I \iff \forall X, \|X\|^2 = \sum_k |\langle X, x_k\rangle|^2$
- [ ] 能推导对偶框架的构造：$y_k = S^{-1}x_k$，并理解 $S^{-1}$ 的自伴性在推导中的作用
- [ ] 能写出两种对偶重构公式，并说明系数计算原子与重建原子可以互换
- [ ] 能通过二维三向量示例，具体计算框架算子、框架上下界、对偶原子和重构过程
- [ ] 能解释"互补基"的概念，以及它与双正交基的区别与联系
- [ ] 能阐述过完备字典的三元权衡：正交基、Parseval 框架、一般过完备字典在计算复杂度与稳健性之间的递变关系
- [ ] 能类比 OFDM 中的对偶框架思想：导频波形 ↔ $x_k$，信道均衡 ↔ $S^{-1}$


## 6. 思考题：拓展与挑战

1. **框架界的几何意义**：框架界 $A$ 和 $B$ 分别度量了最坏方向和最好方向上的"能量保存比"。试证明：当 $A = B$ 时，框架是紧框架，且重构可简化为 $X = \frac{1}{A}\sum_k \langle X, x_k\rangle x_k$。为什么条件数 $B/A$ 大的框架在数值上不稳定？

2. **对偶框架的唯一性**：在过完备框架中（$N > n$），对偶框架 $\{y_k = S^{-1}x_k\}$（正则对偶）并不是唯一的——存在其他对偶框架也能实现重构。试给出一个 $N > n$ 框架存在非正则对偶的充要条件，并举一个简单的 $\mathbb{R}^2$ 三向量例子，构造出不同于正则对偶的另一组对偶。

3. **框架与压缩感知的联系**：压缩感知的测量过程 $\mathbf{y} = \mathbf{\Phi}\mathbf{x}$ 中，如果 $\mathbf{x}$ 在某个框架 $\mathbf{D}$ 下是稀疏的，即 $\mathbf{x} = \mathbf{D}\mathbf{s}$（$\mathbf{s}$ 稀疏），那么 $\mathbf{y} = \mathbf{\Phi}\mathbf{D}\mathbf{s}$。这里的 $\mathbf{\Phi}\mathbf{D}$ 是否也构成一个框架？这个框架的性质（特别是相干性）如何影响稀疏恢复的性能？

4. **SVD 与框架算子**：框架算子 $S = \Theta^*\Theta$（其中 $\Theta$ 是分析算子）是一个正定 Hermitian 矩阵，其特征分解与框架的 SVD 之间存在直接关系。试说明 $\Theta$ 的奇异值与 $S$ 的特征值之间的关系，并解释为何奇异值决定了框架的"紧密性"。

5. **Welch Bound 与框架设计**：Welch Bound 给出了任意框架相干性的下界。如果要在 $\mathbb{R}^n$ 中构造一个 $N$ 个原子的框架，且希望相干性尽可能小，能达到的理论极限是什么？等角紧框架（ETF）为什么是最优的？ETF 在什么维度和原子数组合下存在（或不存在）？





<div style="page-break-before: always;"></div>