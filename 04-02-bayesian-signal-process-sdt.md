<div style="page-break-before: always; padding: 8% 8% 0 8%;">
 <h1 id="第十五讲-贝叶斯决策与采样基础" style="text-align: center; margin-bottom: 2rem; border-bottom: none; display: block;">第十五讲：贝叶斯决策与采样基础</h1> 
 <div style="display: flex; align-items: center; justify-content: center; gap: 12px; margin: 1.8rem auto;">
  <span style="flex:1; max-width:80px; height:1px; background: linear-gradient(to right, transparent, #888);"></span>
  <span style="display:inline-block; width:6px; height:6px; background:#38bdf8; border-radius:50%;"></span>
  <span style="flex:1; max-width:80px; height:1px; background: linear-gradient(to left, transparent, #888);"></span>
 </div>
</div>

<!-- # 第十五讲：贝叶斯决策与采样基础 -->

## 1. 导言

### 1.1 贝叶斯思想回顾

上一讲确立了贝叶斯方法的哲学根基：参数不是未知的常数，而是对世界认知的一种概率表达。先验分布 \(P(\theta)\) 编码了观测数据之前的信念，似然函数 \(P(X|\theta)\) 描述了数据与参数的关联机制，而后验分布 \(P(\theta|X)\) 则是数据到来之后更新了的认知状态。三者通过贝叶斯公式构成一个完整的认知更新循环：

\[
\boxed{
P(\theta|X) = \frac{P(X|\theta)P(\theta)}{P(X)}
}
\]

这个公式本身是简洁的，但它的含义是深远的——它把"学习"这件事翻译成了概率演算。参数 \(\theta\) 的随机性不是物理意义上的随机，而是认知意义上的随机：它度量的是"不知道什么"以及"知道多少"。

然而，贝叶斯公式给出了认知更新的形式，却没有回答两个工程上最紧迫的问题。第一，当需要输出一个具体的决策或估计值时，应该从后验分布中取哪个点？第二，当后验分布无法解析计算时，该怎么办？本讲正是要回答这两个问题。

### 1.2 内容概述

本讲分为两大模块：**统计决策理论**（Statistical Decision Theory，也称风险理论）和 **Monte Carlo 方法**（采样方法）。前者回答"从后验中取哪个点"的问题，后者回答"后验无法解析计算时怎么办"的问题。

**第一个模块：统计决策理论。**

统计决策理论的核心框架是：选择一个损失函数 \(L(\hat{\theta}, \theta)\) 来衡量估计误差的代价，然后寻找使期望损失最小化的估计量。不同的损失函数对应不同的最优决策：

- **均方误差（L2损失）**：\(L(\hat{\theta}, \theta) = (\hat{\theta} - \theta)^2\)。它的最优估计是后验均值，即条件期望 \(\mathbb{E}[\theta|X]\)。这是上一卷中非常熟悉的结论——最小均方误差估计就是条件期望。

- **平均绝对误差（L1损失）**：\(L(\hat{\theta}, \theta) = |\hat{\theta} - \theta|\)。它的最优估计是后验分布的中位数。中位数对异常值不敏感，在厚尾分布下比均值更稳健。

- **0-1损失**：\(L(\hat{\theta}, \theta) = \mathbf{1}_{\hat{\theta} \neq \theta}\)。它的最优估计是后验分布的众数，即最大后验概率估计（MAP）。MAP 估计在计算上最简便（只需找峰值，无需积分），但也最容易过拟合。

下面通过一个二项分布的例子，在先验与后验之间建立起精确的解析关系——后验均值是先验均值与数据比例的凸组合，权重由先验强度与数据量共同决定。

**第二个模块：Monte Carlo 方法。**

在绝大多数实际问题中，后验分布 \(P(\theta|X)\) 的归一化常数 \(P(X) = \int P(X|\theta)P(\theta)d\theta\) 是一个高维积分，解析不可行。这是贝叶斯方法从理论走向应用的"最后一公里"。

解决方案是 **Monte Carlo 方法**——又称采样方法（Sample-based Method）或模拟方法（Simulation-based Method）。其核心思想是：既然无法解析计算后验的积分，那就从后验分布中抽取大量样本，用样本统计量替代解析积分。

具体来说：

- **L2损失（条件期望）**：用样本均值 \(\frac{1}{N}\sum_{i=1}^{N}\theta^{(i)}\) 逼近后验均值；
- **L1损失（中位数）**：将样本排序，取中间位置的值作为中位数的估计；
- **MAP估计（众数）**：先画出样本的直方图，再用光滑曲线拟合，最后寻找曲线的峰值点。

这种方法的美妙之处在于：**把微积分问题转化成了抽样问题。** 一旦掌握了从目标分布中抽取样本的能力，各种积分运算（求均值、求中位数、求分位数、求边缘分布）都可以通过样本统计量来逼近。

本讲在 Monte Carlo 方法中将重点介绍三类基础伪随机数生成器：**均匀分布**、**指数分布**和**高斯分布**的伪随机数生成算法。下面将详细推导它们背后的数学原理——从线性同余发生器到逆变换法，从 Box-Muller 变换到 Polar 方法。这些基础采样器是所有更复杂采样技术的"原子操作"，虽然在实际工程中很少直接靠它们从任意后验分布中采样，但理解它们的原理是通往更高级方法（如马尔可夫链蒙特卡洛 MCMC、粒子滤波等）的必经之路。

**本讲的定位**：如果说上一讲建立了贝叶斯的世界观（"贝叶斯之道"），那么本讲就是在搭建工具箱的第一层——决策工具和基础采样工具。在本讲奠定的基础上，后续课程中逐步引入从任意分布采样的高级方法。

---
## 2. 统计决策理论（Statistical Decision Theory）

统计决策理论，也称为**风险理论**（Risk Theory），是统计学中一个基础性的理论框架。它研究的核心问题是：**在不确定性环境下，如何做出最优的决策？**

更具体地说，给定观测数据 \( X \) 和未知参数 \( \theta \)，需要构造一个决策规则（即估计量）\( \hat{\theta}(X) \)，使得某个预定的损失函数 \( d(\theta, \hat{\theta}(X)) \) 尽可能小。统计决策理论提供了一套系统的语言——损失函数、风险函数、决策规则、最优性准则——来分析和比较不同的估计方法。

**频率学派与贝叶斯学派如何看待统计决策理论？**

两大学派在统计决策理论上的分歧，根源在于对参数 \( \theta \) 的理解不同。

在**频率学派**的视角下，\( \theta \) 是一个**未知但确定的常数**。损失函数 \( d(\theta, \hat{\theta}(X)) \) 的随机性完全来自于数据 \( X \)，因为 \( X \) 是随机的（每次观测都会波动），而 \( \theta \) 是固定的。因此，频率学派定义**风险函数**（Risk Function）为：

\[
R(\theta, \hat{\theta}) = \mathbb{E}_{X \mid \theta}\left[ d(\theta, \hat{\theta}(X)) \right] = \int d(\theta, \hat{\theta}(X)) P(X \mid \theta) \, dX
\]

注意这里期望只对 \( X \) 取，\( \theta \) 被视为给定的常数。风险函数 \( R(\theta, \hat{\theta}) \) 是 \( \theta \) 的函数——不同的 \( \theta \) 值对应不同的风险。频率学派通常无法在 \( \theta \) 上取平均（因为 \( \theta \) 没有分布），因此他们只能比较不同 \( \theta \) 下的风险，或者寻找在所有 \( \theta \) 上都表现良好的估计量（如 minimax 准则、容许性准则等）。

在**贝叶斯学派**的视角下，\( \theta \) 本身是一个**随机变量**，具有先验分布 \( P(\theta) \)。因此，损失函数的随机性来自两个来源：数据的随机性（\( X \)）和参数本身的随机性（\( \theta \)）。贝叶斯学派将损失函数对 \( X \) 和 \( \theta \) 的联合分布取期望，得到**贝叶斯风险**（Bayesian Risk）：

\[
\text{Bayes Risk} = \mathbb{E}_{X, \theta}\left[ d(\theta, \hat{\theta}(X)) \right] = \int \int d(\theta, \hat{\theta}(X)) P(X, \theta) \, dX \, d\theta
\]

贝叶斯风险是一个单一的数值（不是 \( \theta \) 的函数），因此贝叶斯学派可以直接在所有决策规则中寻找使贝叶斯风险最小的那一个。这正是贝叶斯统计决策理论的核心：**最小化后验期望损失**。

两大学派在统计决策理论上的立场差异，本质上是"是否允许在参数 \( \theta \) 上取平均"的差异。频率学派认为在常数上取平均没有意义，因此只能逐点分析风险；贝叶斯学派认为在参数上取平均正是对认知不确定性的量化，因此贝叶斯风险是自然的决策准则。这也解释了为什么贝叶斯方法在决策理论中有着更为统一的数学结构——给定一个损失函数，最优决策规则是唯一的，由后验分布直接给出。

下面从贝叶斯学派的角度来研究统计决策理论。

---

### 2.1 度量介绍

假设感兴趣的参数是 \( \theta \in \mathbb{R}^m \)，观测数据为 \( X \in \mathbb{R}^n \)。统计推断的目标是构造一个估计量 \( \hat{\theta}(X) \)，使得它尽可能接近真实的参数 \( \theta \)。

为了度量"接近"的程度，引入一个**距离函数**或**损失函数** \( d(\theta, \hat{\theta}(X)) \)。这个函数衡量当真实参数为 \( \theta \)、用 \( \hat{\theta}(X) \) 去估计它时所付出的代价。常见的距离函数包括：

- **平方距离**：\( d(\theta, \hat{\theta}) = \|\theta - \hat{\theta}\|_2^2 \)
- **绝对距离**：\( d(\theta, \hat{\theta}) = \|\theta - \hat{\theta}\|_1 \)
- **0-1 距离**：\( d(\theta, \hat{\theta}) = \mathbf{1}_{\theta \neq \hat{\theta}} \)

由于 \( X \) 和 \( \theta \) 都是随机的，\( d(\theta, \hat{\theta}(X)) \) 本身是一个随机变量。要研究一个随机变量，通常需要考察它的数字特征——这里关心它的期望值（即平均损失）。因此，目标是寻找一个估计量 \( \hat{\theta}(X) \)，使得在联合分布 \( (X, \theta) \) 下的期望损失最小：

\[
\boxed{
\theta^* = \arg\min_{\hat{\theta}(X)} \; \mathbb{E}_{X, \theta}\left[ d(\theta, \hat{\theta}(X)) \right]
}
\tag{15.1}
\]

其中 \( \mathbb{E}_{X, \theta} \) 表示在 \( X \) 和 \( \theta \) 的联合分布 \( P(X, \theta) \) 下取期望。这里的期望需要对所有可能的 \( X \) 和 \( \theta \) 进行积分（或求和）：

\[
\mathbb{E}_{X, \theta}\left[ d(\theta, \hat{\theta}(X)) \right]
= \int_{\mathbb{R}^m} \int_{\mathbb{R}^n} d(\theta, \hat{\theta}(X)) \, P(X, \theta) \, dX \, d\theta
\tag{15.2}
\]

这个期望可以写成两种等价的形式，它们分别对应了频率学派和贝叶斯学派的不同视角。

---

#### 2.1.1 频率学派视角：以似然表示

利用联合分布的分解 \( P(X, \theta) = P(X \mid \theta) P(\theta) \)，（15.2）式可以写成：

\[
\boxed{
\mathbb{E}_{X, \theta}\left[ d(\theta, \hat{\theta}(X)) \right]
= \int_{\mathbb{R}^m} \underbrace{\left( \int_{\mathbb{R}^n} d(\theta, \hat{\theta}(X)) \, P(X \mid \theta) \, dX \right)}_{\text{似然（给定 }\theta\text{ 下的期望损失）}} P(\theta) \, d\theta
}
\tag{15.3}
\]

在这个形式中，内层积分是在给定参数 \( \theta \) 的条件下对数据 \( X \) 取期望——这正是频率学派的风险函数 \( R(\theta, \hat{\theta}) \)。外层积分则是在先验分布 \( P(\theta) \) 下对这个风险取平均。这个形式直观地展示了频率学派的视角：固定 \( \theta \)，看数据波动带来的平均损失，然后在 \( \theta \) 上加权平均。

**但这个形式有一个问题：似然 \( P(X|\theta) \) 在实际中往往是无法知道的。**

为什么？因为 \( P(X|\theta) \) 描述的是"在给定参数 \( \theta \) 下，数据 \( X \) 是如何生成的"——它是一个完整的概率模型。在现实问题中，虽然可以假设数据的生成机制（比如"数据服从高斯分布"），但并不知道 \( \theta \) 的真实值，因此无法真正计算这个积分。如果 \( \theta \) 是已知的，那就不需要做统计推断——正是因为 \( \theta \) 未知，这个积分才无法直接计算。

因此，如果硬要用（15.3）式来优化，就会陷入循环困境：要知道 \( \theta \) 才能计算内层积分，但 \( \theta \) 正是要估计的未知量。

**工程上的替代方案：经验风险最小化。**

由于（15.3）中的内层积分无法计算，转而用数据本身的样本来近似它。这就是**经验风险最小化**（Empirical Risk Minimization, ERM）的基本思想：用样本均值替代总体期望。

\[
\arg\min_{\hat{\theta}(X)} \; \frac{1}{N} \sum_{i=1}^{N} d(\theta_i, \hat{\theta}(X_i))
\tag{15.4}
\]

注意这里的 \( \theta_i \) 和 \( X_i \) 是联合分布 \( P(X,\theta) \) 下的 \( N \) 个样本。当样本量 \( N \) 足够大时，根据大数定律，样本均值会收敛到总体期望（15.3）。经验风险最小化是机器学习中许多算法的理论基础——当真实的参数 \( \theta \) 未知时，无法直接计算期望损失，但可以用观测到的数据来近似它。

---

#### 2.1.2 贝叶斯学派视角：以后验表示

利用联合分布的另一种分解 \( P(X, \theta) = P(\theta \mid X) P(X) \)，（15.2）式可以写成：

\[
\boxed{
\mathbb{E}_{X, \theta}\left[ d(\theta, \hat{\theta}(X)) \right]
= \int_{\mathbb{R}^n} \underbrace{\left( \int_{\mathbb{R}^m} d(\theta, \hat{\theta}(X)) \, P(\theta \mid X) \, d\theta \right)}_{\text{后验期望损失}} P(X) \, dX
}
\tag{15.5}
\]

这个形式与（15.3）有本质的区别。内层积分是在给定数据 \( X \) 的条件下，对参数 \( \theta \) 取期望——这正是**后验期望损失**。外层积分则是在数据 \( X \) 的边缘分布 \( P(X) \) 下对这个期望取平均。

**这个形式展示了贝叶斯决策理论的核心结构。**

如果数据 \( X \) 已经观测到了（在工程中，数据是确定的、已知的），那么 \( P(X) \) 就是一个常数，外层积分只是将内层后验期望损失乘以一个常数因子。因此，最小化总期望损失等价于最小化**后验期望损失**：

\[
\boxed{
\theta^* = \arg\min_{\hat{\theta}(X)} \; \int d(\theta, \hat{\theta}(X)) \, P(\theta \mid X) \, d\theta
}
\tag{15.6}
\]

这也就是后验期望损失 \( E_{\theta|X} \) 的后验期望最小化。

**这个形式有一个非常重要的直观含义：在数据已经取得之后，数据反而变成了确定的，参数才是真正随机的。**

这正是贝叶斯视角的精髓。在没有看到数据之前，参数和数据的地位是对等的——两者都是随机的。但一旦数据被观测到，它就变成确定的了。此时，决策只依赖于后验分布 \( P(\theta \mid X) \)——这是对参数 \( \theta \) 的全部认知。数据的作用，是通过贝叶斯公式将先验更新为后验，然后在后验下做决策。

两种形式（15.3）和（15.5）的对比，清晰地揭示了频率学派和贝叶斯学派的思维差异：

| | 用似然表示（频率学派视角） | 用后验表示（贝叶斯学派视角） |
| :--- | :--- | :--- |
| **随机性来源** | 数据 \( X \) | 参数 \( \theta \) |
| **需要计算的积分** | \( \int d(\cdot) P(X\mid\theta) dX \) | \( \int d(\cdot) P(\theta\mid X) d\theta \) |
| **关键问题** | \( \theta \) 未知，无法计算积分 | 后验分布已知（至少形式上），可以计算 |
| **工程实现** | 经验风险最小化（用样本近似） | 后验期望最小化（需处理后验积分） |

这两种形式也说明了：贝叶斯决策理论虽然概念上优雅——给定后验即可做决策——但它将困难转移到了"计算后验积分"这个问题上。当 \( \theta \) 是高维向量时，\( \int d(\theta, \hat{\theta})P(\theta\mid X)d\theta \) 同样不可解析计算。这正是本讲后续要引入 Monte Carlo 方法的原因。

---

#### 2.1.3 小结

本节建立了统计决策理论的基本框架：

1. **损失函数** \( d(\theta, \hat{\theta}(X)) \) 量化了估计误差的代价；
2. **贝叶斯风险** 是在联合分布 \( P(X,\theta) \) 下的期望损失；
3. 贝叶斯风险有两种等价写法：用似然表示（频率学派视角）和用后验表示（贝叶斯学派视角）；
4. 在数据已经观测到的条件下，贝叶斯最优决策等价于最小化后验期望损失（15.6）。

在接下来的几节中，将针对具体的损失函数——平方损失、绝对损失和 0-1 损失——分别求解（15.6）式，得到对应的最优估计量：后验均值、后验中位数和后验众数。

### 2.2 常见的度量指标

在上一节中，建立了统计决策理论的一般框架：给定损失函数 \( d(\theta, \hat{\theta}(X)) \)，最优估计量是最小化后验期望损失的 \( \hat{\theta}(X) \)。不同的损失函数对应不同的最优估计量。本节将针对三种最常用的损失函数，分别推导它们的最优解。

---

#### 2.2.1 均方误差（MSE）

均方误差（Mean Squared Error）是最常用的回归损失函数，其定义为预测值与真实值之间差的平方：

\[
d(\theta, \hat{\theta}(X)) = (\theta - \hat{\theta}(X))^2
\tag{15.7}
\]

将 (15.7) 代入 (15.6)，得到后验期望损失：

\[
R(\theta, \hat{\theta}(X)) = \int_{\mathbb{R}^m} (\theta - \hat{\theta}(X))^2 \, P(\theta \mid X) \, d\theta
\tag{15.8}
\]

注意这里 \( \hat{\theta}(X) \) 是决策变量，\( P(\theta \mid X) \) 是后验分布。为了找到使 \( R \) 最小的 \( \hat{\theta}(X) \)，将 \( R \) 视为 \( \hat{\theta} \) 的函数，对其求梯度并令其为零。由于 \( R \) 是 \( \hat{\theta} \) 的二次函数（开口向上），梯度为零的点就是全局最小值点。

为了看得更清楚，将 (15.8) 展开：

\[
\begin{aligned}
R(\hat{\theta}) &= \int (\theta^2 - 2\theta\hat{\theta} + \hat{\theta}^2) \, P(\theta \mid X) \, d\theta \\
&= \int \theta^2 P(\theta \mid X) d\theta - 2\hat{\theta} \int \theta P(\theta \mid X) d\theta + \hat{\theta}^2 \int P(\theta \mid X) d\theta
\end{aligned}
\]

由于 \( \int P(\theta \mid X) d\theta = 1 \)，令 \( \mathbb{E}[\theta \mid X] = \int \theta P(\theta \mid X) d\theta \)，则：

\[
R(\hat{\theta}) = \mathbb{E}[\theta^2 \mid X] - 2\hat{\theta} \mathbb{E}[\theta \mid X] + \hat{\theta}^2
\]

对 \( \hat{\theta} \) 求导：

\[
\frac{\partial R}{\partial \hat{\theta}} = -2\mathbb{E}[\theta \mid X] + 2\hat{\theta}
\]

令导数为零：

\[
-2\mathbb{E}[\theta \mid X] + 2\hat{\theta} = 0
\]

解得：

\[
\boxed{
\theta^* = \mathbb{E}[\theta \mid X] = \int_{\mathbb{R}^m} \theta \, P(\theta \mid X) \, d\theta
}
\tag{15.9}
\]

**结论：在均方误差损失下，最优估计量是后验分布的均值（条件期望）。**

这个结论与第一卷中反复遇到的结果完全一致——在均方意义下的最优估计就是条件期望。只不过第一卷中是在数据和信号之间做条件期望 \( \mathbb{E}[S \mid X] \)，现在是在数据和参数之间做条件期望 \( \mathbb{E}[\theta \mid X] \)。两者的数学结构完全相同，只是变量的含义不同。

---

#### 2.2.2 平均绝对误差（MAE）

平均绝对误差（Mean Absolute Error）是另一种常用的回归损失函数，其定义为：

\[
d(\theta, \hat{\theta}(X)) = |\theta - \hat{\theta}(X)|
\tag{15.10}
\]

将 (15.10) 代入 (15.6)，得到后验期望损失：

\[
R(\theta, \hat{\theta}(X)) = \int_{\mathbb{R}^m} |\theta - \hat{\theta}(X)| \, P(\theta \mid X) \, d\theta
\tag{15.11}
\]

这个优化问题比 MSE 复杂，因为绝对值函数在 \( \theta = \hat{\theta} \) 处不可导，不能直接对整个函数求梯度。在前面的章节中，曾用次梯度法处理过类似的问题——只要 0 落在次梯度区间内，就算找到了最优解。但本节的推导将采用另一种更直接的技巧：**分段去掉绝对值符号**，将积分拆成两段分别处理。

---

#### 2.2.2a 预备知识：莱布尼茨积分求导法则

在推导 MAE 最优解的过程中，需要对形如 \( \int_{a(\theta)}^{b(\theta)} f(x, \theta) dx \) 的积分关于参数 \( \theta \) 求导。这类积分上下限中含有参数，不能简单地将求导号移入积分号内，而需要使用**莱布尼茨积分法则**（Leibniz Integral Rule）。

莱布尼茨法则给出了如下公式：

\[
\boxed{
\frac{\partial}{\partial \theta} \int_{a(\theta)}^{b(\theta)} f(x, \theta) \, dx
= \int_{a(\theta)}^{b(\theta)} \frac{\partial f}{\partial \theta}(x, \theta) \, dx
+ f(b(\theta), \theta) \cdot b'(\theta)
- f(a(\theta), \theta) \cdot a'(\theta)
}
\tag{15.12}
\]

这个公式的含义是：当积分上下限都是 \( \theta \) 的函数时，对 \( \theta \) 的导数由三部分组成——被积函数本身对 \( \theta \) 的偏导数在区间上的积分，加上在上限处的函数值乘以上限的导数，再减去在下限处的函数值乘以下限的导数。

下面逐项说明这个公式的由来。

**公式 (15.12) 的推导思路：**

考虑一个更一般的情形。设 \( g(\theta) = \int_{a(\theta)}^{b(\theta)} f(x, \theta) dx \)，需要计算 \( g'(\theta) \)。

首先将积分拆分为一个固定起点到 \( a(\theta) \) 的积分加上从 \( a(\theta) \) 到 \( b(\theta) \) 的积分，或者利用微积分基本定理。更直接的方法是对 (15.12) 做如下验证：当 \( b(\theta) \) 和 \( a(\theta) \) 都是常数时，公式退化为 \( \frac{\partial}{\partial \theta} \int_a^b f(x,\theta) dx = \int_a^b \frac{\partial f}{\partial \theta} dx \)，这是标准的"积分号下求导"。当上下限随 \( \theta \) 变化时，微积分基本定理表明，上限的变化贡献了 \( f(b(\theta),\theta) \cdot b'(\theta) \) 项，下限的变化贡献了 \( -f(a(\theta),\theta) \cdot a'(\theta) \) 项。

下面推导几种特殊形式的莱布尼茨法则，它们将在 MAE 的求导中直接使用。

---

**第一种形式：下限为常数 0，上限为 \( f(\theta) \)**

设 \( F(\theta) = \int_0^{f(\theta)} g(x, \theta) \, dx \)。

此时 \( a(\theta) = 0 \)，\( b(\theta) = f(\theta) \)，\( a'(\theta) = 0 \)，\( b'(\theta) = f'(\theta) \)。代入 (15.12)：

\[
\boxed{
\frac{\partial}{\partial \theta} \int_0^{f(\theta)} g(x, \theta) \, dx
= \int_0^{f(\theta)} \frac{\partial g}{\partial \theta}(x, \theta) \, dx
+ g(f(\theta), \theta) \cdot f'(\theta)
}
\tag{15.13}
\]

---

**第二种形式：上限为常数 0，下限为 \( h(\theta) \)**

设 \( F(\theta) = \int_{h(\theta)}^{0} g(x, \theta) \, dx \)。

此时 \( a(\theta) = h(\theta) \)，\( b(\theta) = 0 \)，\( a'(\theta) = h'(\theta) \)，\( b'(\theta) = 0 \)。代入 (15.12)：

\[
\boxed{
\frac{\partial}{\partial \theta} \int_{h(\theta)}^{0} g(x, \theta) \, dx
= \int_{h(\theta)}^{0} \frac{\partial g}{\partial \theta}(x, \theta) \, dx
- g(h(\theta), \theta) \cdot h'(\theta)
}
\tag{15.14}
\]

---

**第三种形式：一般情形，下限为 \( h(\theta) \)，上限为 \( f(\theta) \)**

设 \( F(\theta) = \int_{h(\theta)}^{f(\theta)} g(x, \theta) \, dx \)。

由 (15.12) 直接给出：

\[
\boxed{
\frac{\partial}{\partial \theta} \int_{h(\theta)}^{f(\theta)} g(x, \theta) \, dx
= \int_{h(\theta)}^{f(\theta)} \frac{\partial g}{\partial \theta}(x, \theta) \, dx
+ g(f(\theta), \theta) \cdot f'(\theta)
- g(h(\theta), \theta) \cdot h'(\theta)
}
\tag{15.15}
\]

---

下面回到 MAE 的优化问题。

为了去掉绝对值，将积分拆成两段：当 \( \theta \le \hat{\theta} \) 时，\( |\theta - \hat{\theta}| = \hat{\theta} - \theta \)；当 \( \theta \ge \hat{\theta} \) 时，\( |\theta - \hat{\theta}| = \theta - \hat{\theta} \)。

\[
R(\hat{\theta}) = \int_{-\infty}^{\hat{\theta}} (\hat{\theta} - \theta) \, P(\theta \mid X) \, d\theta
+ \int_{\hat{\theta}}^{\infty} (\theta - \hat{\theta}) \, P(\theta \mid X) \, d\theta
\tag{15.16}
\]

现在，将 \( R \) 视为 \( \hat{\theta} \) 的函数（注意：\( \hat{\theta} \) 是决策变量，同时也是积分上限/下限的边界点），对其求导并令其为零。

需要对 (15.16) 中的两项分别求导。

**第一项 \( R_1(\hat{\theta}) = \int_{-\infty}^{\hat{\theta}} (\hat{\theta} - \theta) P(\theta \mid X) \, d\theta \) 对 \( \hat{\theta} \) 的导数：**

这符合莱布尼茨法则的形式：下限为常数 \( -\infty \)，上限为 \( f(\hat{\theta}) = \hat{\theta} \)，被积函数 \( g(\theta, \hat{\theta}) = (\hat{\theta} - \theta) P(\theta \mid X) \)。

应用 (15.13)：

\[
\frac{\partial R_1}{\partial \hat{\theta}}
= \int_{-\infty}^{\hat{\theta}} \frac{\partial}{\partial \hat{\theta}} \left[ (\hat{\theta} - \theta) P(\theta \mid X) \right] d\theta
+ (\hat{\theta} - \hat{\theta}) P(\hat{\theta} \mid X) \cdot 1
\]

第二项中 \( (\hat{\theta} - \hat{\theta}) = 0 \)，所以第二项为零。第一项中被积函数对 \( \hat{\theta} \) 的偏导数为：

\[
\frac{\partial}{\partial \hat{\theta}} (\hat{\theta} - \theta) P(\theta \mid X) = P(\theta \mid X)
\]

（注意 \( P(\theta \mid X) \) 不依赖于 \( \hat{\theta} \)，因为它是后验分布，而 \( \hat{\theta} \) 是要选择的决策变量，不是后验分布的参数。）

因此：

\[
\frac{\partial R_1}{\partial \hat{\theta}} = \int_{-\infty}^{\hat{\theta}} P(\theta \mid X) \, d\theta
\tag{15.17}
\]

**第二项 \( R_2(\hat{\theta}) = \int_{\hat{\theta}}^{\infty} (\theta - \hat{\theta}) P(\theta \mid X) \, d\theta \) 对 \( \hat{\theta} \) 的导数：**

这符合莱布尼茨法则的形式：下限为 \( h(\hat{\theta}) = \hat{\theta} \)，上限为常数 \( \infty \)，被积函数 \( g(\theta, \hat{\theta}) = (\theta - \hat{\theta}) P(\theta \mid X) \)。

应用 (15.14)：

\[
\frac{\partial R_2}{\partial \hat{\theta}}
= \int_{\hat{\theta}}^{\infty} \frac{\partial}{\partial \hat{\theta}} \left[ (\theta - \hat{\theta}) P(\theta \mid X) \right] d\theta
- (\hat{\theta} - \hat{\theta}) P(\hat{\theta} \mid X) \cdot 1
\]

同样，第二项中 \( (\hat{\theta} - \hat{\theta}) = 0 \)，所以第二项为零。第一项中被积函数对 \( \hat{\theta} \) 的偏导数为：

\[
\frac{\partial}{\partial \hat{\theta}} (\theta - \hat{\theta}) P(\theta \mid X) = -P(\theta \mid X)
\]

因此：

\[
\frac{\partial R_2}{\partial \hat{\theta}} = - \int_{\hat{\theta}}^{\infty} P(\theta \mid X) \, d\theta
\tag{15.18}
\]

将 (15.17) 和 (15.18) 相加：

\[
\frac{\partial R}{\partial \hat{\theta}} = \frac{\partial R_1}{\partial \hat{\theta}} + \frac{\partial R_2}{\partial \hat{\theta}}
= \int_{-\infty}^{\hat{\theta}} P(\theta \mid X) \, d\theta - \int_{\hat{\theta}}^{\infty} P(\theta \mid X) \, d\theta
\]

令梯度为零：

\[
\int_{-\infty}^{\hat{\theta}} P(\theta \mid X) \, d\theta - \int_{\hat{\theta}}^{\infty} P(\theta \mid X) \, d\theta = 0
\]

\[
\boxed{
\int_{-\infty}^{\hat{\theta}} P(\theta \mid X) \, d\theta = \int_{\hat{\theta}}^{\infty} P(\theta \mid X) \, d\theta
}
\tag{15.19}
\]

这个等式说明：在最优解 \( \hat{\theta} \) 处，后验分布被分成了左右两部分，且两部分的概率质量相等。换句话说，\( \hat{\theta} \) 将后验分布的面积等分为两半。

这正是**中位数**的定义：**中位数是将概率密度函数分成左右两个相等面积的点**，即 \( F(\hat{\theta}) = 1/2 \)，其中 \( F \) 是后验分布的累积分布函数（CDF）。

因此：

\[
\boxed{
\theta^* = \text{中位数}(\theta \mid X)
}
\tag{15.20}
\]

**结论：在平均绝对误差损失下，最优估计量是后验分布的中位数。** 与均值相比，中位数对异常值不敏感，在厚尾分布下更稳健。

---

#### 2.2.3 0-1 损失

0-1 损失是一种用于分类问题的损失函数，其核心思想是：只要估计值足够接近真实值（在某个小邻域内），损失就为零；一旦超出这个邻域，损失就为一。它的定义如下：

\[
d(\theta, \hat{\theta}(X)) =
\begin{cases}
0, & |\hat{\theta}(X) - \theta| \le \delta \\
1, & |\hat{\theta}(X) - \theta| > \delta
\end{cases}
\tag{15.21}
\]

其中 \( \delta > 0 \) 是一个预先设定的小常数，表示"容忍范围"。这个损失函数的含义是：在容忍范围内，风险很小；超出容忍范围，风险很大。当 \( \delta \to 0 \) 时，0-1 损失退化为要求精确相等。

将 (15.21) 代入 (15.6)，得到后验期望损失：

\[
R(\hat{\theta}) = \int d(\theta, \hat{\theta}) \, P(\theta \mid X) \, d\theta
\]

由于损失函数在 \( |\theta - \hat{\theta}| > \delta \) 时取 1，在 \( |\theta - \hat{\theta}| \le \delta \) 时取 0，可以写成：

\[
\boxed{
R(\hat{\theta}) = 1 - \int_{\hat{\theta} - \delta}^{\hat{\theta} + \delta} P(\theta \mid X) \, d\theta
}
\tag{15.22}
\]

这个表达式的含义是：风险等于 1 减去后验分布在区间 \( [\hat{\theta} - \delta, \hat{\theta} + \delta] \) 上的概率质量。为了最小化风险，需要**最大化**这个区间内的后验概率质量。

当 \( \delta \) 充分小时，区间 \( [\hat{\theta} - \delta, \hat{\theta} + \delta] \) 非常窄。在这个小区间内，后验密度 \( P(\theta \mid X) \) 近似为常数 \( P(\hat{\theta} \mid X) \)，积分的近似为：

\[
\int_{\hat{\theta} - \delta}^{\hat{\theta} + \delta} P(\theta \mid X) \, d\theta \approx 2\delta \cdot P(\hat{\theta} \mid X)
\tag{15.23}
\]

因此 (15.22) 可以近似为：

\[
R(\hat{\theta}) \approx 1 - 2\delta \, P(\hat{\theta} \mid X)
\tag{15.24}
\]

由于 \( 1 \) 和 \( 2\delta \) 都是常数（与 \( \hat{\theta} \) 无关），最小化 \( R(\hat{\theta}) \) 等价于最大化 \( P(\hat{\theta} \mid X) \)：

\[
\arg\min_{\hat{\theta}} R(\hat{\theta}) \Longleftrightarrow \arg\max_{\hat{\theta}} P(\hat{\theta} \mid X)
\tag{15.25}
\]

\[
\boxed{
\theta^* = \arg\max_{\theta} P(\theta \mid X)
}
\tag{15.26}
\]

**结论：在 0-1 损失下，最优估计量是后验分布的众数，即最大后验概率估计（Maximum A Posteriori, MAP）。**

---

**什么是众数（Mode）？**

众数是概率分布中最常见的值——也就是概率密度函数取最大值时对应的点。在离散分布中，众数是出现频率最高的那个值；在连续分布中，众数是密度函数 \( P(\theta \mid X) \) 达到最大值的位置。

后验众数估计量 MAP 有一个重要的性质：**它依赖于后验密度函数的形状，而不依赖于具体的数值，因此计算相对简单**，不需要像计算后验均值那样需要积分。

当 \( \delta \to 0 \) 时，0-1 损失要求估计值与真实值精确相等。此时，最优估计量就是后验密度最大的点——因为只有在这个点附近（任意小的邻域内），后验概率质量才能被区间 \( [\hat{\theta} - \delta, \hat{\theta} + \delta] \) 捕获。换句话说，MAP 是"最大化后验概率"的估计量。

如果先验是均匀分布（平坦先验），则 MAP 估计退化为最大似然估计（MLE）。因此，MLE 可以看作是 MAP 在无信息先验下的特例。

---

#### 2.2.4 三种损失函数与最优估计量的对比

本节推导的三个结果汇总如下：

| 损失函数 | 名称 | 最优估计量 | 计算需求 |
| :--- | :--- | :--- | :--- |
| \( (\theta - \hat{\theta})^2 \) | 均方误差（MSE） | 后验均值 \( \mathbb{E}[\theta \mid X] \) | 需计算一阶矩积分 |
| \( \|\theta - \hat{\theta}\| \) | 平均绝对误差（MAE） | 后验中位数 | 需累积分布函数 |
| \( \mathbf{1}_{\|\theta - \hat{\theta}\| > \delta} \) | 0-1 损失 | 后验众数（MAP） | 只需找到密度峰值 |

这三个结果揭示了一个重要的模式：**不同的损失函数对应了后验分布的不同特征点——均值、中位数、众数。** 这三者在高斯分布（对称单峰）下重合，但在非对称或多峰分布下截然不同。因此，在实际问题中，选择哪个估计量不应该由方便性决定，而应该由愿意承担什么样的决策风险来决定。

---

### 2.3 综合运用的例子：二项分布与 Beta 共轭

前两节分别推导了三种损失函数下的最优估计量——后验均值、中位数和众数。这些结论给出了"从后验中取哪个点"的答案，但前提是已经知道了后验分布的具体形式。在实际问题中，后验分布的计算本身往往就是挑战。

本节通过一个完整的例子，展示贝叶斯推断的完整流程：从先验假设出发，写出似然函数，计算后验分布，最后求出后验均值。更重要的是，将在这个过程中揭示一个深刻的数学结构——**共轭先验**——它使得先验和后验保持相同的形式，让贝叶斯更新变得简洁而优美。

#### 2.3.1 问题设定与先验

假设有一个未知的成功概率 \( \theta \in [0, 1] \)，观测到了 \( n \) 次独立重复试验中的成功次数 \( X \)。这是一个典型的二项分布模型：

\[
X \mid \theta \sim \text{Binomial}(n, \theta)
\]

即：

\[
P(X = x \mid \theta) = \binom{n}{x} \theta^x (1-\theta)^{n-x}, \quad x = 0, 1, \cdots, n
\tag{15.27}
\]

目标是通过观测到的 \( X \) 来推断 \( \theta \)。

首先，假设先验 \( \theta \) 服从均匀分布：

\[
\theta \sim U(0, 1), \quad P(\theta) = 1, \quad 0 \le \theta \le 1
\tag{15.28}
\]

均匀分布表达的是"无信息先验"——在观测数据之前，所有可能的 \( \theta \) 值都是等可能的。这个假设看似合理，但实际上它并非真正的"无信息"，稍后会讨论这一点。

---

#### 2.3.2 联合分布

根据贝叶斯公式，后验分布正比于似然函数乘以先验。联合分布为：

\[
P(x, \theta) = P(x \mid \theta) P(\theta) = \binom{n}{x} \theta^x (1-\theta)^{n-x} \cdot 1
\tag{15.29}
\]

为了得到后验分布，需要计算归一化常数 \( P(x) \)：

\[
P(x) = \int_0^1 P(x, \theta) \, d\theta = \binom{n}{x} \int_0^1 \theta^x (1-\theta)^{n-x} \, d\theta
\tag{15.30}
\]

这个积分需要用到 Beta 分布的知识，先来回顾一下。

---

#### 2.3.3 Beta 分布回顾

Beta 分布是定义在 \( [0, 1] \) 区间上的连续概率分布，由两个正参数 \( a, b > 0 \) 控制其形状：

\[
\boxed{
f(\theta; a, b) = \frac{1}{B(a, b)} \theta^{a-1} (1-\theta)^{b-1}, \quad 0 \le \theta \le 1
}
\tag{15.31}
\]

其中 \( B(a, b) \) 是 Beta 函数，定义为：

\[
B(a, b) = \frac{\Gamma(a) \Gamma(b)}{\Gamma(a+b)}
\tag{15.32}
\]

Gamma 函数 \( \Gamma(x) = \int_0^\infty t^{x-1} \exp(-t) dt \) 满足一个重要的性质：当 \( n \in \mathbb{N} \) 时，\( \Gamma(n) = (n-1)! \)。因此，当 \( a, b \) 为正整数时：

\[
B(a, b) = \frac{(a-1)! (b-1)!}{(a+b-1)!}
\tag{15.33}
\]

于是 Beta 分布可以写成：

\[
f(\theta; a, b) = \frac{(a+b-1)!}{(a-1)! (b-1)!} \theta^{a-1} (1-\theta)^{b-1}
\tag{15.34}
\]

Beta 分布的一个重要性质是它的归一化积分：

\[
\int_0^1 \theta^{a-1} (1-\theta)^{b-1} \, d\theta = B(a, b) = \frac{(a-1)! (b-1)!}{(a+b-1)!}
\tag{15.35}
\]

这个公式将在后面的计算中反复使用。

---

#### 2.3.4 均匀先验下的后验分布

现在回到 (15.30)。利用 (15.35) 式，令 \( a = x+1 \)，\( b = n-x+1 \)：

\[
\int_0^1 \theta^x (1-\theta)^{n-x} \, d\theta = \frac{x! (n-x)!}{(n+1)!}
\tag{15.36}
\]

代入 (15.30)：

\[
P(x) = \binom{n}{x} \frac{x! (n-x)!}{(n+1)!}
= \frac{n!}{x! (n-x)!} \cdot \frac{x! (n-x)!}{(n+1)!}
= \frac{1}{n+1}
\tag{15.37}
\]

因此，在均匀先验下，观测到 \( X = x \) 的边缘概率恒为 \( 1/(n+1) \)，与 \( x \) 的取值无关。这是一个有趣的结果：在没有先验偏好的情况下，所有的成功次数 \( x = 0, 1, \cdots, n \) 都是等可能的。

后验分布为：

\[
P(\theta \mid x) = \frac{P(x, \theta)}{P(x)}
= \binom{n}{x} \theta^x (1-\theta)^{n-x} \cdot (n+1)
\tag{15.38}
\]

\[
\boxed{
P(\theta \mid x) = \frac{(n+1)!}{x! (n-x)!} \theta^x (1-\theta)^{n-x}
}
\tag{15.39}
\]

这是一个参数为 \( a = x+1 \)，\( b = n-x+1 \) 的 Beta 分布。验证一下：将 \( a = x+1 \)，\( b = n-x+1 \) 代入 (15.34)：

\[
f(\theta; x+1, n-x+1) = \frac{(n+1)!}{x! (n-x)!} \theta^x (1-\theta)^{n-x}
\]

与 (15.39) 完全一致。

---

#### 2.3.5 后验均值（MSE 最优估计）

在均方误差下，最优估计量是后验均值。计算：

\[
\mathbb{E}[\theta \mid x] = \int_0^1 \theta \, P(\theta \mid x) \, d\theta
= \int_0^1 \theta \cdot \frac{(n+1)!}{x! (n-x)!} \theta^x (1-\theta)^{n-x} \, d\theta
\]

\[
= \frac{(n+1)!}{x! (n-x)!} \int_0^1 \theta^{x+1} (1-\theta)^{n-x} \, d\theta
\]

再次利用 (15.35)，令 \( a = x+2 \)，\( b = n-x+1 \)：

\[
\int_0^1 \theta^{x+1} (1-\theta)^{n-x} \, d\theta
= \frac{(x+1)! (n-x)!}{(n+2)!}
\]

代入得：

\[
\mathbb{E}[\theta \mid x]
= \frac{(n+1)!}{x! (n-x)!} \cdot \frac{(x+1)! (n-x)!}{(n+2)!}
= \frac{(n+1)!}{(n+2)!} \cdot \frac{(x+1)!}{x!}
= \frac{x+1}{n+2}
\]

\[
\boxed{
\mathbb{E}[\theta \mid x] = \frac{x+1}{n+2}
}
\tag{15.40}
\]

**这个结果可以写成加权平均（凸组合）的形式：**

\[
\boxed{
\mathbb{E}[\theta \mid x] = \frac{n}{n+2} \cdot \underbrace{\frac{x}{n}}_{\text{数据比例 } \hat{p}} + \frac{2}{n+2} \cdot \underbrace{\frac{1}{2}}_{\text{先验均值 } \mathbb{E}[\theta]}
}
\tag{15.41}
\]

这个形式完美地展示了贝叶斯更新的本质：后验均值是"数据的信息"和"先验的信息"的加权平均。

- **数据比例** \( x/n \) 是样本中成功的频率，这是数据对 \( \theta \) 的直接估计；
- **先验均值** \( 1/2 \) 是均匀先验下 \( \theta \) 的期望；
- **权重**由 \( n \)（数据量）和 \( 2 \)（先验的"有效样本量"）决定。

当 \( n \) 很大时，数据权重 \( n/(n+2) \to 1 \)，后验均值几乎完全由数据决定；当 \( n \) 很小时，先验权重 \( 2/(n+2) \) 不可忽略，后验均值被先验拉向 \( 1/2 \)。

---

#### 2.3.6 无信息先验的困境

上面的推导在数学上是完美的——从均匀先验出发，得到了 Beta 后验，并求出了后验均值的闭式解。

**但这里有一个问题：均匀先验真的是"无信息"的吗？**

均匀先验 \( \theta \sim U(0,1) \) 似乎表达了"对 \( \theta \) 没有任何先验知识"的立场。但仔细想想，它其实包含了一个隐含的假设：**在 \( [0,1] \) 区间上，每个 \( \theta \) 值的先验权重是相等的。**

如果将参数做非线性变换，比如定义 \( \phi = \theta^2 \)，那么 \( \theta \) 上的均匀分布会诱导出 \( \phi \) 上的非均匀分布。这意味着"均匀先验"依赖于参数化的选择——它不是参数化不变的。换句话说，**均匀先验并非真正的"无信息"，它只是在一个特定的参数化下看起来是平坦的。**

更有甚者，在某些参数化下，均匀先验可能导致不合理的推断。例如，如果将成功概率 \( \theta \) 用对数几率 \( \log(\theta/(1-\theta)) \) 来表示，\( \theta \) 上的均匀先验会诱导出对数几率上的一个复杂先验，其形状强烈偏向两端。

因此，**一个好的先验应该满足两个条件**：
1. 它应该反映在观测数据之前的真实信念（如果有的话）；
2. 如果要表达"无信息"，它应该具有参数化不变性（如 Jeffreys 先验 \( P(\theta) \propto \sqrt{I(\theta)} \)）。

但真正想要的是：**先验和后验在形式上保持一致，这样贝叶斯更新就会变得极其简洁**。这就引出了下面的概念。

---

#### 2.3.7 共轭先验（Conjugate Prior）

**什么是共轭先验？**

设似然函数为 \( P(x \mid \theta) \)。如果先验分布 \( P(\theta) \) 和后验分布 \( P(\theta \mid x) \) 属于同一分布族，则称该先验为该似然的**共轭先验**。

换句话说，**共轭先验是"在贝叶斯更新下封闭"的分布族**——先验和后验具有相同的数学形式，只是参数被更新了。

共轭先验的好处是巨大的：
1. **计算简便**：后验分布可以直接通过更新先验参数得到，无需复杂的数值积分；
2. **解析可处理**：后验的均值、方差、分位数等都可以通过解析公式计算；
3. **直观的序贯更新**：每来一批新数据，只需要更新参数，不需要重新计算整个后验。

对于二项分布（Bernoulli 似然），Beta 分布是它的共轭先验：

\[
\underbrace{\text{Beta}(a, b)}_{\text{先验}} + \underbrace{\text{Binomial}(n, \theta)}_{\text{似然}} \longrightarrow \underbrace{\text{Beta}(a + x, b + n - x)}_{\text{后验}}
\tag{15.42}
\]

也就是说，如果先验是 Beta\( (a, b) \)，观测到 \( x \) 次成功和 \( n-x \) 次失败后，后验是 Beta\( (a+x, b+n-x) \)。

---

#### 2.3.8 Beta 先验下的后验推导

下面将先验从均匀分布推广到一般的 Beta 分布：

\[
P(\theta) = \frac{(a+b-1)!}{(a-1)! (b-1)!} \theta^{a-1} (1-\theta)^{b-1}, \quad 0 \le \theta \le 1
\tag{15.43}
\]

其中 \( a, b > 0 \) 是先验参数，控制了对 \( \theta \) 的先验信念：
- \( a \) 可以理解为"先验成功次数"；
- \( b \) 可以理解为"先验失败次数"；
- 先验均值 \( \mathbb{E}[\theta] = a/(a+b) \)；
- 先验的"有效样本量"为 \( a+b \)。

联合分布为：

\[
P(x, \theta) = P(x \mid \theta) P(\theta)
= \binom{n}{x} \theta^x (1-\theta)^{n-x} \cdot \frac{(a+b-1)!}{(a-1)! (b-1)!} \theta^{a-1} (1-\theta)^{b-1}
\]

合并 \( \theta \) 的幂次：

\[
P(x, \theta) = \binom{n}{x} \frac{(a+b-1)!}{(a-1)! (b-1)!}
\theta^{x+a-1} (1-\theta)^{n-x+b-1}
\tag{15.44}
\]

归一化常数 \( P(x) \)：

\[
P(x) = \binom{n}{x} \frac{(a+b-1)!}{(a-1)! (b-1)!}
\int_0^1 \theta^{x+a-1} (1-\theta)^{n-x+b-1} \, d\theta
\]

利用 (15.35)，令 \( a' = x+a \)，\( b' = n-x+b \)：

\[
\int_0^1 \theta^{x+a-1} (1-\theta)^{n-x+b-1} \, d\theta
= \frac{(x+a-1)! (n-x+b-1)!}{(n+a+b-1)!}
\tag{15.45}
\]

因此：

\[
P(x) = \binom{n}{x} \frac{(a+b-1)!}{(a-1)! (b-1)!}
\cdot \frac{(x+a-1)! (n-x+b-1)!}{(n+a+b-1)!}
\tag{15.46}
\]

后验分布为：

\[
P(\theta \mid x) = \frac{P(x, \theta)}{P(x)}
= \frac{(n+a+b-1)!}{(x+a-1)! (n-x+b-1)!}
\theta^{x+a-1} (1-\theta)^{n-x+b-1}
\tag{15.47}
\]

这正是一个 Beta 分布，参数为：

\[
\boxed{
\theta \mid x \sim \text{Beta}(a + x, \; b + n - x)
}
\tag{15.48}
\]

---

#### 2.3.9 后验均值

后验均值为：

\[
\mathbb{E}[\theta \mid x] = \int_0^1 \theta \, P(\theta \mid x) \, d\theta
\]

\[
= \frac{(n+a+b-1)!}{(x+a-1)! (n-x+b-1)!}
\int_0^1 \theta^{x+a} (1-\theta)^{n-x+b-1} \, d\theta
\]

利用 (15.35)，令 \( a'' = x+a+1 \)，\( b'' = n-x+b \)：

\[
\int_0^1 \theta^{x+a} (1-\theta)^{n-x+b-1} \, d\theta
= \frac{(x+a)! (n-x+b-1)!}{(n+a+b)!}
\]

代入得：

\[
\mathbb{E}[\theta \mid x]
= \frac{(n+a+b-1)!}{(x+a-1)! (n-x+b-1)!}
\cdot \frac{(x+a)! (n-x+b-1)!}{(n+a+b)!}
\]

\[
= \frac{(n+a+b-1)!}{(n+a+b)!} \cdot \frac{(x+a)!}{(x+a-1)!}
= \frac{x+a}{n+a+b}
\]

\[
\boxed{
\mathbb{E}[\theta \mid x] = \frac{x+a}{n+a+b}
}
\tag{15.49}
\]

---

#### 2.3.10 凸组合形式——奥妙所在

将 (15.49) 写成加权平均的形式：

\[
\boxed{
\mathbb{E}[\theta \mid x] = \frac{n}{n+a+b} \cdot \underbrace{\frac{x}{n}}_{\text{数据比例 } \hat{p}}
+ \frac{a+b}{n+a+b} \cdot \underbrace{\frac{a}{a+b}}_{\text{先验均值 } \mu_0}
}
\tag{15.50}
\]

这个形式揭示了贝叶斯更新的全部奥妙。

**结构一：权重由"有效样本量"决定。**

数据权重为 \( n/(n+a+b) \)，先验权重为 \( (a+b)/(n+a+b) \)。这里的 \( a+b \) 是 Beta 先验的"有效样本量"——它表达了对先验的信念强度。如果 \( a+b \) 很大（即先验很强），先验在后验中占的权重就大；如果 \( n \) 很大（即数据很多），数据占的权重就大。

**结构二：后验参数是先验参数与数据的直接叠加。**

\[
\text{后验参数} = \text{先验参数} + \text{数据统计量}
\]

\[
a_{\text{post}} = a_{\text{prior}} + x, \quad b_{\text{post}} = b_{\text{prior}} + (n-x)
\tag{15.51}
\]

这揭示了共轭先验的本质：**先验参数可以解释为"虚拟观测值"**。Beta 先验中的 \( a \) 和 \( b \) 可以理解为在看到真实数据之前，已经"虚拟地"观察到了 \( a \) 次成功和 \( b \) 次失败。新的数据来了之后，只需要在虚拟观测值上加上真实观测值即可。这种"数据叠加"的简洁性正是共轭先验的奥妙所在。

**结构三：无论先验参数 \( a, b \) 取什么值，后验均值总是先验均值与数据比例的凸组合。**

这是一个普遍性的结论——不依赖于具体的先验参数取值。只要先验是 Beta 分布（二项分布的共轭先验），后验均值就具有凸组合的形式。这完美地呼应了前面章节中"精度加权平均"的结论：在这个离散分布的例子里，精度就是"有效样本量"——数据的样本量越大，数据在加权平均中的权重越大；先验的有效样本量越大，先验在加权平均中的权重越大。

---

#### 2.3.11 两种先验的对比

为了更清楚地看到共轭先验的优势，将均匀先验和 Beta 先验的结果并排对比：

| | 均匀先验 \( U(0,1) \) | Beta 先验 \( \text{Beta}(a, b) \) |
| :--- | :--- | :--- |
| **先验分布** | \( \text{Beta}(1, 1) \) | \( \text{Beta}(a, b) \) |
| **先验均值** | \( 1/2 \) | \( a/(a+b) \) |
| **后验分布** | \( \text{Beta}(x+1, n-x+1) \) | \( \text{Beta}(a+x, b+n-x) \) |
| **后验均值** | \( \frac{x+1}{n+2} \) | \( \frac{x+a}{n+a+b} \) |
| **凸组合** | \( \frac{n}{n+2}\hat{p} + \frac{2}{n+2}\cdot\frac{1}{2} \) | \( \frac{n}{n+a+b}\hat{p} + \frac{a+b}{n+a+b}\cdot\frac{a}{a+b} \) |

均匀先验 \( U(0,1) \) 实际上是 Beta 先验在 \( a = b = 1 \) 时的特例。Beta 先验的参数 \( a \) 和 \( b \) 提供了两个自由度，可以独立地控制先验的**位置**（通过 \( a/(a+b) \)）和**强度**（通过 \( a+b \)）。而均匀先验固定了位置为 \( 1/2 \)，强度为 \( 2 \)，完全没有灵活性。

---

#### 2.3.12 数值例子

假设进行了 \( n = 20 \) 次试验，观察到 \( x = 8 \) 次成功，数据比例 \( \hat{p} = 0.4 \)。

使用 Beta 先验，取 \( a = 2 \)，\( b = 2 \)（先验均值为 \( 0.5 \)，强度为 \( 4 \)）：

\[
\mathbb{E}[\theta \mid x] = \frac{8 + 2}{20 + 4} = \frac{10}{24} \approx 0.417
\]

后验均值 0.417 介于数据比例 0.4 和先验均值 0.5 之间，略微偏向数据，因为数据量 \( n=20 \) 大于先验强度 \( a+b=4 \)。

如果增加先验强度，取 \( a = 5 \)，\( b = 5 \)（先验均值仍为 0.5，强度为 10）：

\[
\mathbb{E}[\theta \mid x] = \frac{8 + 5}{20 + 10} = \frac{13}{30} \approx 0.433
\]

更强的先验将后验均值拉得更靠近 0.5。

如果取 \( a = 1 \)，\( b = 3 \)（先验均值为 0.25，强度为 4）：

\[
\mathbb{E}[\theta \mid x] = \frac{8 + 1}{20 + 4} = \frac{9}{24} = 0.375
\]

先验均值 0.25 将后验均值拉向了更小的方向。

---

#### 2.3.13 小结

本节通过二项分布与 Beta 共轭先验的例子，揭示了贝叶斯推断的几个核心结构：

1. **共轭先验使得后验分布保持与先验相同的形式**，参数更新规则极其简洁：\( a \leftarrow a + x \)，\( b \leftarrow b + n - x \)。
2. **先验参数可以解释为"虚拟观测值"**，贝叶斯更新本质上是在先验的虚拟数据之上叠加真实观测数据。
3. **后验均值是先验均值与数据比例的凸组合**，权重由"有效样本量"决定——这正是"精度加权平均"在离散分布上的表现形式。
4. **共轭先验的选择提供了灵活性**：\( a/(a+b) \) 控制先验的位置（\( \theta \) 在哪里），\( a+b \) 控制先验的强度（对这个位置有多大把握）。

这个例子虽然是离散分布的特殊情况，但它揭示的结构——先验参数作为虚拟数据、后验参数等于先验参数加数据、后验均值是凸组合——实际上是一般贝叶斯推断的缩影。在后续章节中，当后验无法解析计算时，将用 Monte Carlo 方法来逼近这些量，但"先验与数据的加权平均"这一基本认知，将始终指导理解贝叶斯更新的行为。

---

## 3. 后验分布的计算：Monte Carlo 方法

在上一节中，通过二项分布与 Beta 共轭先验的例子，展示了后验分布可以解析求解的完美情况——先验和后验同属一个分布族，后验均值和方差都有闭式表达式。然而，这种"完美情况"在现实中极其罕见。

当面对更复杂的模型时——比如多峰分布、高维参数空间、非共轭先验、复杂的似然函数——后验分布 \( P(\theta \mid X) \) 的归一化常数 \( P(X) = \int P(X \mid \theta)P(\theta)d\theta \) 往往是一个无法解析计算的高维积分。没有归一化常数，就无法得到完整的后验密度函数，更不用说计算后验均值、中位数或众数了。

本节将介绍应对这一问题的核心方法：**Monte Carlo 方法**。它的基本思想是：既然解析计算行不通，就用随机采样来近似。

---

### 3.1 后验计算的 Monte Carlo 框架

在贝叶斯推断中，后验分布 \( P(\theta \mid X) \) 是对参数 \( \theta \) 的全部认知。一旦得到了后验，需要从中提取各种统计量来回答问题。

在均方误差（MSE）准则下，需要后验均值：

\[
\mathbb{E}[\theta \mid X] = \int_{-\infty}^{\infty} \theta \, P(\theta \mid X) \, d\theta
\tag{15.52}
\]

在平均绝对误差（MAE）准则下，需要后验中位数：

\[
\text{Median}(\theta \mid X) = \text{满足 } \int_{-\infty}^{m} P(\theta \mid X) d\theta = \frac{1}{2} \text{ 的 } m
\tag{15.53}
\]

在 0-1 损失下，需要后验众数：

\[
\text{Mode}(\theta \mid X) = \arg\max_{\theta} P(\theta \mid X)
\tag{15.54}
\]

这些量都涉及对后验分布的积分或优化。

**问题在于：这些积分通常是积不出来的。**

为什么会积不出来？原因有两个层面。

**第一，归一化常数不可计算。**

贝叶斯公式给出：

\[
P(\theta \mid X) = \frac{P(X \mid \theta) P(\theta)}{\int P(X \mid \theta) P(\theta) \, d\theta}
\]

分母 \( P(X) = \int P(X \mid \theta) P(\theta) d\theta \) 是一个积分。当 \( \theta \) 的维度 \( m \) 较大时，这是一个高维积分。在高维空间中，积分区域的体积随维度指数增长，数值积分方法（如网格积分）会遇到"维数灾难"——所需的网格点数随维度指数增加，计算量迅速变得不可行。

**第二，后验分布往往没有标准形式。**

即使不需要归一化常数（比如用 MCMC 方法直接采样），后验分布 \( P(\theta \mid X) \) 本身通常不是熟悉的标准分布（高斯、Beta、Gamma 等），无法直接用解析公式写出其均值和方差。它可能是多峰的、扭曲的、高维的，其统计量只能通过数值方法来估计。

因此，需要一种不依赖于解析积分的方法来近似这些统计量。这就是 **Monte Carlo 方法**（Monte Carlo Method）——也称为 **Sample Based 方法**（基于采样的方法）或 **Simulation Based 方法**（基于模拟的方法）。

---

**Monte Carlo 方法的基本思想：**

> 如果无法解析地计算一个分布 \( P(\theta) \) 的积分，那就从这个分布中抽取大量样本 \( \{\theta^{(1)}, \theta^{(2)}, \ldots, \theta^{(N)}\} \)，然后用样本统计量来近似真实的分布统计量。

这种方法的神奇之处在于：**把微积分问题（积分）转化成了抽样问题（统计）**。大数定律保证，当样本量 \( N \to \infty \) 时，样本统计量会收敛到真实的分布统计量。

Monte Carlo 方法在贝叶斯推断中的工作流程可以分为两步：

**第一步：从后验分布中生成随机样本。**

这是整个流程中最核心、最困难的一步。需要从目标分布 \( P(\theta \mid X) \) 中抽取独立同分布的样本 \( \{\theta^{(1)}, \theta^{(2)}, \ldots, \theta^{(N)}\} \)。

**第二步：处理样本，计算所需统计量。**

有了样本之后，各种统计量的计算就变成了简单的"数据处理"问题：

- **后验均值（MSE 最优估计）**：用样本均值近似
  \[
  \mathbb{E}[\theta \mid X] \approx \frac{1}{N} \sum_{i=1}^{N} \theta^{(i)}
  \tag{15.55}
  \]

- **后验中位数（MAE 最优估计）**：将样本排序，取中间位置的值
  \[
  \text{Median}(\theta \mid X) \approx \theta_{\left(\frac{N+1}{2}\right)}
  \tag{15.56}
  \]

- **后验众数（MAP 估计）**：先画出样本的直方图，然后用光滑曲线拟合直方图，找到曲线的峰值点
  \[
  \text{Mode}(\theta \mid X) \approx \arg\max_{\theta} \hat{f}(\theta)
  \tag{15.57}
  \]
  其中 \( \hat{f}(\theta) \) 是核密度估计（Kernel Density Estimation）得到的密度曲线。

Monte Carlo 方法的美妙之处在于：**一旦掌握了从任意分布中采样的方法，所有的贝叶斯推断问题——无论是求均值、中位数、众数，还是求分位数、可信区间——都可以通过简单的样本统计量来解决。**

---

**最重要的问题：如何从目标分布 \( F \) 中采样？**

在第一节（1.2）中简要提到过，对于一个任意分布 \( F \)，有一个标准的采样框架：

\[
\boxed{
X \sim F \quad \Longleftrightarrow \quad X = F^{-1}(U), \quad U \sim U(0, 1)
}
\tag{15.58a}
\]

这就是 **逆变换方法**（Inverse Transform Method）。它的核心思想是：如果 \( U \) 是 \( [0,1] \) 上的均匀随机变量，那么 \( F^{-1}(U) \) 的分布函数就是 \( F \)。

下面验证一下这个结论。设 \( Y = F^{-1}(U) \)，计算 \( Y \) 的累积分布函数：

\[
\begin{aligned}
F_Y(y) &= P(Y \le y) \\
&= P(F^{-1}(U) \le y) \\
&= P(U \le F(y)) \\
&= F(y)
\end{aligned}
\tag{15.58b}
\]

这里的推导用到了两个关键事实：第一，分布函数 \( F \) 是单调非递减的，所以 \( F^{-1}(U) \le y \iff U \le F(y) \)；第二，\( U \) 是 \( [0,1] \) 上的均匀分布，其累积分布函数为 \( P(U \le u) = u \)。

上面的推导验证了 \( F_Y(y) = F(y) \)。因此，\( Y \) 的累积分布函数与目标分布 \( F \) 完全相同，也就是说 \( Y \sim F \)。

逆变换方法的美妙之处在于：**只要能生成均匀分布的随机数，就能通过求分布函数 \( F \) 的反函数来生成任意分布的随机数。** 这个方法把"从任意分布采样"这个看似困难的问题，归结为了两个子问题：生成均匀随机数，以及求目标分布的反函数。

因此，整个采样的基础就归结为一个问题：**如何生成均匀分布 \( U(0,1) \) 的伪随机数？**

这就是下一节的主题。

---

### 3.2 生成均匀分布的伪随机数

均匀分布 \( U(0,1) \) 的随机数是所有随机数生成器的"原子操作"。一旦有了高质量的均匀随机数，就可以通过逆变换、Box-Muller 变换、接受-拒绝采样等方法生成其他分布的随机数。

在实际计算中，无法产生真正的随机数——计算机是确定性的机器。只能产生**伪随机数**（Pseudo-Random Numbers）：由确定性算法生成的一串数，它们看起来像是从均匀分布中独立抽取的，但实际上是确定性的。

#### 3.2.1 线性同余发生器（LCG）

最经典的均匀随机数生成方法是 **线性同余发生器**（Linear Congruential Generator, LCG），由 Lehmer 在 1951 年提出。其递推公式为：

\[
\boxed{
X_n = (a X_{n-1} + c) \bmod m
}
\tag{15.59}
\]

其中：
- \( X_0 \) 是**种子**（Seed），即初始值；
- \( a \) 是乘数（Multiplier）；
- \( c \) 是增量（Increment）；
- \( m \) 是模数（Modulus）。

生成的 \( X_n \) 是 \( \{0, 1, \ldots, m-1\} \) 中的整数。通常将其除以 \( m \) 来得到 \( (0,1) \) 区间上的均匀伪随机数：

\[
U_n = \frac{X_n}{m}
\tag{15.60}
\]

**LCG 示例：**

取 \( a = 7 \)，\( c = 3 \)，\( m = 10 \)，\( X_0 = 1 \)：

\[
\begin{aligned}
X_1 &= (7 \times 1 + 3) \bmod 10 = 10 \bmod 10 = 0 \\
X_2 &= (7 \times 0 + 3) \bmod 10 = 3 \\
X_3 &= (7 \times 3 + 3) \bmod 10 = 24 \bmod 10 = 4 \\
X_4 &= (7 \times 4 + 3) \bmod 10 = 31 \bmod 10 = 1
\end{aligned}
\]

序列为：1, 0, 3, 4, 1, 0, 3, 4, ... 周期为 4，非常短。

**LCG 的缺陷：**

1. **周期性问题**：LCG 生成的序列最终一定会进入循环。周期的长度最多为 \( m \)，但实际周期通常远小于 \( m \)。如果周期太短，生成的随机数序列会在重复中失去随机性。

2. **高位/低位相关性**：LCG 生成的随机数在低维空间中往往存在明显的结构。例如，如果将连续的三个随机数作为三维空间中的点 \( (U_n, U_{n+1}, U_{n+2}) \)，它们会落在少数几个平面上，而不是均匀地填充三维空间。这就是所谓的"格结构"（Lattice Structure）问题。

3. **对参数敏感**：LCG 的随机性质量高度依赖于参数 \( a, c, m \) 的选择。选择不当会导致极短的周期或严重的相关性。

4. **不满足严格的随机性测试**：许多 LCG 实现无法通过现代随机性测试套件（如 Die Hard 测试）。

#### 3.2.2 梅森旋转算法（Mersenne Twister）

为了克服 LCG 的缺陷，Makoto Matsumoto 和 Takuji Nishimura 在 1998 年提出了 **梅森旋转算法**（Mersenne Twister），这是目前最广泛使用的伪随机数生成器之一。

梅森旋转算法的核心思想是：维护一个大的状态向量，通过线性反馈移位寄存器（Linear Feedback Shift Register, LFSR）生成新的状态，然后从状态中输出随机数。它的名字来源于其周期——梅森素数 \( 2^{19937} - 1 \)。

你给出的递推公式是梅森旋转算法的一种简化表示：

\[
\begin{aligned}
X_n &= (A_1 + B_1 X_{n-1}) \bmod C_1 \\
Y_n &= (1 + L^{16})(1 + R^{16}) X_n \bmod C_2 \\
Z_n &= (1 + R^{16})(1 + L^{16}) Y_n \bmod C_3 \\
D_n &= (X_n + Y_n + Z_n) \bmod C_4
\end{aligned}
\tag{15.61}
\]

这里的核心思想是：**通过多次不同的线性变换和模运算，将不同来源的随机性混合在一起，消除单一 LCG 可能存在的相关性。**

- \( L^{16} \) 和 \( R^{16} \) 表示左移和右移 16 位的位运算操作；
- 四个不同的模数 \( C_1, C_2, C_3, C_4 \) 提供了不同的"随机源"；
- 最后通过模加运算将三者合并，使得输出的 \( D_n \) 不再依赖于单个生成器的内部结构。

这种方法被称为 **组合生成器**（Combined Generator）或 **混合生成器**——通过组合多个独立的随机数生成器，它们的缺陷（相关性、周期性、低维结构）会被相互抵消。

**MT19937 的关键特性：**

| 特性 | 数值 |
| :--- | :--- |
| 周期 | \( 2^{19937} - 1 \approx 10^{6001} \) |
| 状态大小 | 19937 比特（约 2.5 KB） |
| 输出维度 | 623 维均匀分布（在 623 维以内，点均匀填充超立方体） |
| 速度 | 非常快（用位运算实现） |

梅森旋转算法通过了几乎所有已知的随机性测试套件，包括 **Die Hard 测试**——这是由 George Marsaglia 设计的、由 16 个统计测试组成的随机性检验标准。

#### 3.2.3 现代方法

尽管梅森旋转算法仍然是许多科学计算软件（如 Python、R、MATLAB）的默认随机数生成器，但近年来的研究表明它也存在一些细微的缺陷（如对某些统计测试的失败率偏高）。目前更先进的方法包括：

- **PCG**（Permuted Congruential Generator）：在 LCG 的基础上增加了输出置换，具有更好的统计性质；
- **Xorshift** 家族：基于异或和移位操作，速度极快；
- **Cryptographic PRNGs**：如 ChaCha20，具有密码学级别的随机性（但速度较慢）。

在大多数工程应用中，梅森旋转算法已经足够好。重要的是理解一点：**无论使用哪种伪随机数生成器，得到的都只是对真正随机性的近似**。在后续的 MCMC 方法中，这种近似的影响可以通过增加采样长度和进行收敛诊断来控制和评估。

---

**小结：**

- Monte Carlo 方法将贝叶斯推断中的积分问题转化为采样问题；
- 一旦获得后验分布的样本，均值、中位数、众数都可以通过简单的样本统计量近似；
- 逆变换方法给出了从任意分布采样的通用框架：\( F^{-1}(U) \)；
- 采样的基础是均匀分布伪随机数；
- 线性同余发生器是最基础的方法，但存在周期性和相关性问题；
- 梅森旋转算法通过大的状态空间和位运算解决了这些问题，是目前最广泛使用的生成器；
- 本讲的目的是建立"采样替代积分"的基本认知，下一讲将介绍从任意分布中采样的更高级方法（MCMC）。

---

### 3.3 生成指数分布的伪随机数
原理
$
F(y) = 1-\exp(-\lambda y) \implies F^{-1} (y) = -\frac{1}{\lambda} \log(1-y) \\
\implies X_k = -\frac{1}{\lambda} \log(1 - u_k) \sim Exp(\lambda)
$

这个原理仔细推导一下，解释一下。
### 3.4 生成高斯分布的伪随机数

均匀分布的伪随机数生成之后，下一个关键问题是：如何生成高斯分布（正态分布）的伪随机数？高斯分布在信号处理中无处不在——噪声建模、先验分布、似然函数，几乎处处都有它的身影。

面临的问题是：高斯分布的累积分布函数 \( \Phi(x) = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{x} \exp(-t^2/2) dt \) 没有解析的闭式表达式，其反函数 \( \Phi^{-1} \) 也没有解析形式。因此，**逆变换法不能直接用于高斯分布**，需要另辟蹊径。

本节将介绍两种方法：一种是基于中心极限定理的近似方法（直观但代价大），另一种是 **Box-Muller 变换**（精确且高效）。下面从几何直观出发，用极坐标变换详细推导 Box-Muller 变换的每一步。

---

#### 3.4.1 方法一：基于中心极限定理的近似（不推荐）

中心极限定理指出：若 \( U_1, U_2, \ldots, U_n \overset{\text{i.i.d.}}{\sim} U(0,1) \)，则当 \( n \to \infty \) 时：

\[
\frac{\sum_{i=1}^{n} U_i - n/2}{\sqrt{n/12}} \xrightarrow{d} \mathcal{N}(0, 1)
\tag{15.62}
\]

也就是说，多个独立均匀随机变量的和，经过标准化后近似服从标准高斯分布。

**为什么不推荐？**

为了获得较好的近似，通常需要 \( n \ge 12 \)。这意味着每生成一个高斯随机数，就需要生成 12 个均匀随机数，计算效率很低。更重要的是，这种方法只能产生有界范围内的随机数（因为均匀分布的和是有界的），无法真正覆盖高斯分布的尾部，在需要极端值的应用中有本质缺陷。因此，不考虑这种方法。

---

#### 3.4.2 方法二：Box-Muller 变换——几何直观

现在介绍 Box-Muller 变换，它通过**极坐标变换**将两个独立均匀随机数精确地转化为两个独立的标准高斯随机数。

**几何直观：**

考虑二维平面上的一个随机点 \( (X_1, X_2) \)，其中 \( X_1 \) 和 \( X_2 \) 独立同分布，均服从标准高斯分布 \( \mathcal{N}(0, 1) \)。它们的联合概率密度为：

\[
p(x_1, x_2) = \frac{1}{2\pi} \exp\left(-\frac{x_1^2 + x_2^2}{2}\right)
\tag{15.63}
\]

这个密度函数有一个非常重要的性质：**它只依赖于 \( x_1^2 + x_2^2 \)**，即只依赖于点到原点的距离，而与方向无关。换句话说，这个二维高斯分布是**旋转对称的**——它的等高线是同心圆。

因此，如果用极坐标 \( (X_1, X_2) \to (R, \Theta) \) 来表示这个随机点：

\[
X_1 = R \cos \Theta, \quad X_2 = R \sin \Theta
\tag{15.64}
\]

那么：
- \( \Theta \) 应该在 \( (0, 2\pi) \) 上均匀分布（因为旋转对称性）；
- \( R \) 的分布只依赖于距离，与角度无关。

这个几何观察表明：**生成一个二维高斯分布，等价于生成一个随机方向（均匀角度）和一个随机半径（瑞利分布）。** 一旦有了 \( (R, \Theta) \)，就可以通过 (15.64) 得到两个独立的高斯随机变量。

下面的任务就变成了：证明 \( R \) 服从瑞利分布，并求出如何从均匀随机数生成瑞利分布。

---

#### 3.4.3 极坐标变换与变量替换

设 \( X_1, X_2 \overset{\text{i.i.d.}}{\sim} \mathcal{N}(0, 1) \)。它们的联合概率密度为：

\[
p_{X_1, X_2}(x_1, x_2) = \frac{1}{2\pi} \exp\left(-\frac{x_1^2 + x_2^2}{2}\right)
\tag{15.65}
\]

进行极坐标变换：

\[
x_1 = r \cos \theta, \quad x_2 = r \sin \theta, \quad r \ge 0, \quad 0 \le \theta < 2\pi
\tag{15.66}
\]

这个变换的**雅可比行列式**（Jacobian determinant）为：

\[
J = \left| \det \begin{pmatrix} 
\frac{\partial x_1}{\partial r} & \frac{\partial x_1}{\partial \theta} \\
\frac{\partial x_2}{\partial r} & \frac{\partial x_2}{\partial \theta}
\end{pmatrix} \right|
= \left| \det \begin{pmatrix} 
\cos \theta & -r \sin \theta \\
\sin \theta & r \cos \theta
\end{pmatrix} \right|
= |r \cos^2 \theta + r \sin^2 \theta| = r
\tag{15.67}
\]

多变量概率密度函数的变量替换公式为：

\[
p_{R, \Theta}(r, \theta) = p_{X_1, X_2}(r\cos\theta, r\sin\theta) \cdot |J|
\tag{15.68}
\]

代入 (15.65) 和 (15.67)：

\[
p_{R, \Theta}(r, \theta) = \frac{1}{2\pi} \exp\left(-\frac{r^2}{2}\right) \cdot r
\tag{15.69}
\]

其中 \( r \ge 0 \)，\( 0 \le \theta < 2\pi \)。

观察 (15.69)，\( p_{R, \Theta}(r, \theta) \) 可以分解为两个独立因子的乘积——一个只依赖于 \( r \)，另一个只依赖于 \( \theta \)：

\[
\boxed{
p_{R, \Theta}(r, \theta) = \underbrace{r \exp\left(-\frac{r^2}{2}\right)}_{\text{只依赖于 } r} \cdot \underbrace{\frac{1}{2\pi}}_{\text{只依赖于 } \theta}
}
\tag{15.70}
\]

因此，\( R \) 和 \( \Theta \) 是**相互独立的**，且：
- \( \Theta \sim U(0, 2\pi) \)（均匀分布）；
- \( R \) 的概率密度为 \( p_R(r) = r \exp(-r^2/2) \)，\( r \ge 0 \)，这正是瑞利分布。

---

#### 3.4.4 瑞利分布的逆变换法

现在有了 \( \Theta \) 的分布——均匀分布，可以直接用均匀随机数生成。但 \( R \) 的分布是瑞利分布，需要从瑞利分布中采样。

瑞利分布的概率密度函数为：

\[
\boxed{
p_R(r) = r \exp\left(-\frac{r^2}{2}\right), \quad r \ge 0
}
\tag{15.71}
\]

它的累积分布函数为：

\[
F_R(r) = \int_{0}^{r} t \exp\left(-\frac{t^2}{2}\right) dt
\tag{15.72}
\]

令 \( u = \frac{t^2}{2} \)，则 \( du = t \, dt \)。当 \( t = 0 \) 时 \( u = 0 \)，当 \( t = r \) 时 \( u = r^2/2 \)：

\[
F_R(r) = \int_{0}^{r^2/2} \exp(-u) \, du = 1 - \exp\left(-\frac{r^2}{2}\right)
\tag{15.73}
\]

因此：

\[
\boxed{
F_R(r) = 1 - \exp\left(-\frac{r^2}{2}\right), \quad r \ge 0
}
\tag{15.74}
\]

求反函数：设 \( u = F_R(r) \)，则 \( u = 1 - \exp(-r^2/2) \)，即 \( \exp(-r^2/2) = 1 - u \)，\( r^2/2 = -\log(1-u) \)，所以：

\[
\boxed{
F_R^{-1}(u) = \sqrt{-2 \log(1-u)}, \quad u \in (0, 1)
}
\tag{15.75}
\]

由于 \( U \sim U(0,1) \) 时，\( 1-U \sim U(0,1) \) 也服从均匀分布，可将 (15.75) 简化为：

\[
\boxed{
R = \sqrt{-2 \log U}, \quad U \sim U(0, 1)
}
\tag{15.76}
\]

---

#### 3.4.5 Box-Muller 变换的最终形式

下面把所有的部分组合在一起：

1. 生成两个独立的均匀随机数 \( U_1, U_2 \sim U(0, 1) \)；
2. 从均匀角度得到 \( \Theta = 2\pi U_2 \)；
3. 从瑞利分布得到 \( R = \sqrt{-2 \log U_1} \)；
4. 通过极坐标变换得到两个独立的标准高斯随机数：

\[
\boxed{
X_1 = \sqrt{-2 \log U_1} \cdot \cos(2\pi U_2), \quad
X_2 = \sqrt{-2 \log U_1} \cdot \sin(2\pi U_2)
}
\tag{15.77}
\]

这就是 **Box-Muller 变换**的最终形式。

---

#### 3.4.6 几何直观的完整图景

下面用几何语言来总结整个过程：

1. **二维高斯分布是旋转对称的**：它的密度函数 \( p(x_1, x_2) = \frac{1}{2\pi}\exp(-(x_1^2+x_2^2)/2) \) 只依赖于半径 \( r \)，不依赖于角度 \( \theta \)。

2. **旋转对称性的含义**：这意味着，如果在二维平面上随机放置一个点，其落点服从二维高斯分布，那么：
   - **方向是均匀随机的**——点出现在任何方向上的概率相同；
   - **半径是独立的**——点离原点有多远，与它朝哪个方向无关。

3. **半径的分布是瑞利分布**：\( p_R(r) = r\exp(-r^2/2) \)。这个分布表示"离原点越远的点，其概率密度先增后减"——这是二维高斯分布中"大多数点落在单位圆附近"的数学表达。

4. **逆变换法的应用**：瑞利分布的累积分布函数有闭式解，因此可以通过逆变换法高效采样。

5. **最终一步**：有了独立的半径和角度，通过 \( X_1 = R\cos\Theta, X_2 = R\sin\Theta \) 回到直角坐标系，就得到了两个独立的高斯随机数。

这个几何视角揭示了 Box-Muller 变换的本质：**它不是直接在高斯分布上做逆变换，而是利用二维高斯分布的旋转对称性，将一维的"高斯问题"转化为一维的"角度问题"（均匀分布）和一维的"半径问题"（瑞利分布）。** 后者有闭式解，因此整个问题变得可解。

---

**小结：**

- 高斯分布的 CDF 没有闭式形式，逆变换法不能直接使用；
- Box-Muller 变换利用二维高斯分布的旋转对称性，将问题转化为角度（均匀分布）和半径（瑞利分布）的采样；
- 瑞利分布的 CDF 有闭式形式 \( F_R(r) = 1 - \exp(-r^2/2) \)，其反函数为 \( \sqrt{-2\log(1-u)} \)；
- 最终得到 \( X_1 = \sqrt{-2\log U_1}\cos(2\pi U_2) \)，\( X_2 = \sqrt{-2\log U_1}\sin(2\pi U_2) \)；
- 该方法精确、高效、无需近似，是生成高斯伪随机数的标准方法。

---


## 4. 课后总结

### 4.1 核心逻辑链：从决策理论到 Monte Carlo 采样

本讲建立了一条从贝叶斯统计决策到伪随机数生成方法的完整链条：

1. **两种决策视角**：频率学派以**似然**为中心，最小化频率风险 $R(f) = \mathbb{E}_{X,Y}[\ell(f(X), Y)]$；贝叶斯学派以**后验**为中心，最小化后验期望损失 $\mathbb{E}_{\theta|X}[\ell(a, \theta)]$。两者形式上不同但最优解等价——贝叶斯决策规则 $a^* = \arg\min_a \mathbb{E}_{\theta|X}[\ell(a,\theta)]$ 是唯一同时最小化后验风险与频率风险的策略。

2. **Monte Carlo 求后验**：后验计算的核心问题是高维积分 $\mathbb{E}_{\theta|X}[g(\theta)] = \int g(\theta) p(\theta|X) d\theta$ 通常无闭式解。Monte Carlo 方法用样本均值近似积分：$\frac{1}{M}\sum_{i=1}^M g(\theta_i) \to \mathbb{E}_{\theta|X}[g(\theta)]$，收敛速度为 $O(1/\sqrt{M})$（与维度无关）。

3. **从均匀分布生成任意分布**：Monte Carlo 的前提是能从目标分布中采样。本讲从均匀分布 $U(0,1)$ 出发，通过 Box-Muller 变换精确生成高斯伪随机数。核心技巧是利用二维高斯分布的旋转对称性，将问题分解为角度（均匀）和半径（瑞利）两个独立子问题，后者可通过逆变换法求解。

### 4.2 频率学派 vs 贝叶斯学派对比

| 维度 | 频率学派 | 贝叶斯学派 |
| :--- | :--- | :--- |
| **核心量** | 似然 $P(X\|\theta)$ | 后验 $P(\theta\|X)$ |
| **参数观** | $\theta$ 是固定未知量 | $\theta$ 是随机变量 |
| **风险函数** | $R(f) = \mathbb{E}_{X,Y}[\ell]$ | $\mathbb{E}_{\theta\|X}[\ell(a,\theta)]$ |
| **优化方式** | 对 $p(X,Y)$ 求期望 | 对 $p(\theta\|X)$ 求期望 |
| **需要先验** | 不需要 | 需要 $P(\theta)$ |

### 4.3 重点公式

**贝叶斯决策规则：**
$$
\boxed{a^*(X) = \arg\min_a \int \ell(a, \theta) \, p(\theta | X) \, d\theta}
$$

**Monte Carlo 近似：**
$$
\boxed{\mathbb{E}_{\theta|X}[g(\theta)] \approx \frac{1}{M}\sum_{i=1}^M g(\theta_i), \quad \theta_i \sim p(\theta|X)}
$$

**Box-Muller 变换：**
$$
\boxed{X_1 = \sqrt{-2\ln U_1}\cos(2\pi U_2), \quad X_2 = \sqrt{-2\ln U_1}\sin(2\pi U_2)}
$$


## 5. 学习检查清单：自测核心知识点掌握情况

- [ ] 能区分频率学派和贝叶斯学派对"参数"的不同理解
- [ ] 能分别写出频率风险 $R(f) = \mathbb{E}_{X,Y}[\ell(f(X), Y)]$ 和后验风险 $\mathbb{E}_{\theta|X}[\ell(a,\theta)]$ 的定义
- [ ] 能解释为什么贝叶斯决策规则同时最小化两种风险
- [ ] 能推导后验分布的 Bayes 公式，并说明配分函数（归一化常数）的计算难点
- [ ] 能写出 Monte Carlo 估计的公式和收敛速度 $O(1/\sqrt{M})$
- [ ] 能解释 Monte Carlo 方法为何不依赖维数（与数值积分对比）
- [ ] 能阐述 Box-Muller 变换的几何直观：旋转对称性 → 角度均匀 + 半径瑞利
- [ ] 能推导瑞利分布 $p_R(r) = r\exp(-r^2/2)$ 的 CDF 及其反函数
- [ ] 能写出 Box-Muller 变换的完整推导过程（极坐标变换 → 雅可比行列式 → 逆变换法）
- [ ] 能对比中心极限定理近似和 Box-Muller 变换在精度和效率上的优劣


## 6. 思考题：拓展与挑战

1. **贝叶斯风险 vs 频率风险的数值关系**：在平方损失 $\ell(a, \theta) = (a - \theta)^2$ 下，贝叶斯最优决策是后验均值 $\hat{\theta}_B = \mathbb{E}[\theta|X]$。频率风险 $R(\hat{\theta}_B) = \mathbb{E}_{X,\theta}[(\hat{\theta}_B(X) - \theta)^2]$ 是否总等于贝叶斯风险？它们之间的一般关系是什么？

2. **Monte Carlo 的维度无关性边界**：Monte Carlo 的 $O(1/\sqrt{M})$ 收敛速度号称与维度 $d$ 无关，但方差常数 $\sigma_g^2 = \text{Var}[g(\theta)]$ 随 $d$ 增大而增加是常见的。在实际的 $d=100$ 维后验采样中，需要多少样本才能获得 $< 1\%$ 的相对误差？这个样本量在计算上是否实际可行？

3. **Box-Muller vs Marsaglia 极坐标法**：Box-Muller 需要计算 $\sin$ 和 $\cos$，在硬件上较慢。Marsaglia 极坐标法通过拒绝采样的方式避免了三角函数——生成单位圆内的均匀点然后缩放。试推导 Marsaglia 方法的公式，并对比两种方法：哪种更适合 GPU 实现？哪种的随机数质量更优？

4. **逆变换法何时可用**：Box-Muller 的精髓在于瑞利分布的 CDF 有解析逆。但许多重要分布（如 Beta、Gamma、Student-t）的 CDF 逆没有闭式。这时有哪些替代策略？试分类讨论：数值反演、查找表近似、以及更复杂的专用方法（如 Ahrens-Dieter 方法生成 Gamma）。

5. **从采样到统计决策的闭环**：Gibbs 采样、Metropolis-Hastings 等 MCMC 方法是 Monte Carlo 在高维后验计算中的标准工具。在这些方法中，样本是**相关的**（Markov 链），如何调整后验期望的估计公式？有效样本量（ESS）是什么？它如何影响 Monte Carlo 误差估计？


## 7. 实验设计：Box-Muller 实现与 Monte Carlo 收敛性验证

### 7.0 实验总览：从均匀随机数到统计推断的完整实践

本实验实现 Box-Muller 变换生成高斯随机数，验证 Monte Carlo 积分的收敛性，并在一个简单的贝叶斯决策场景中应用。

#### 实验矩阵

| 实验编号 | 实验名称 | 对比维度 | 预期结论 |
| :--- | :--- | :--- | :--- |
| 7.1 | Box-Muller 变换的实现与统计检验 | 正态性检验 | Box-Muller 生成的随机数通过 KS/QQ 检验 |
| 7.2 | Monte Carlo 积分收敛速度验证 | 样本数 $M$ | 误差按 $1/\sqrt{M}$ 衰减，与维度无关 |
| 7.3 | 贝叶斯决策：高斯均值估计的风险对比 | 决策策略 | 贝叶斯后验均值在 MSE 意义上最优 |

---

### 7.1 实验 1：Box-Muller 变换的实现与统计检验

**目的**：实现 Box-Muller 变换，用多种统计检验验证生成的随机数确实服从标准正态分布。

**步骤**：
1. 使用 `numpy.random.rand` 生成 $M = 10^4$ 个 $U(0,1)$ 随机数对 $(U_1, U_2)$
2. 通过 Box-Muller 变换 $X = \sqrt{-2\ln U_1}\cos(2\pi U_2)$ 生成 $M/2$ 个高斯样本
3. 绘制直方图并叠加标准正态 PDF 理论曲线
4. 执行 Kolmogorov-Smirnov 检验（`scipy.stats.kstest`），检验 $p$ 值
5. 绘制 QQ-plot 直观检查尾部偏差
6. 计算样本偏度（skewness）和峰度（kurtosis），与理论值 0 和 3 对比
7. 作为对比，用中心极限定理近似（12 个均匀数之和标准化）生成同样数量的样本，重复上述检验

**定量指标**：
- KS 检验的 $p$ 值（$p < 0.05$ 拒绝正态性假设）
- 偏度估计值（期望 0，标准差 $\approx \sqrt{6/M}$）
- 峰度估计值（期望 3，标准差 $\approx \sqrt{24/M}$）

**预期结果**：
- Box-Muller 生成的样本通过所有正态性检验（$p \gg 0.05$）
- 中心极限定理近似在尾部有明显偏差（QQ-plot 的尾部偏离直线）
- 验证 Box-Muller 是"精确"方法而中心极限定理是"近似"方法

---

### 7.2 实验 2：Monte Carlo 积分收敛速度验证

**目的**：数值验证 Monte Carlo 估计的 $O(1/\sqrt{M})$ 收敛速度，以及收敛速度与维度的无关性。

**步骤**：
1. 构造目标积分 $I_d = \int_{[0,1]^d} \sum_{i=1}^d x_i^2 \, dx_1\cdots dx_d = d/3$（该积分的解析值已知）
2. 分别取 $d = 1, 5, 20$，对每种维度用 Monte Carlo 估计：
   - 样本数 $M$ 从 $10^2$ 指数增长到 $10^6$
   - 用均匀分布采样 $x_i \sim U(0,1)$，计算 $\hat{I}_d = \frac{1}{M}\sum_{j=1}^M \sum_{i=1}^d x_{i,j}^2$
3. 对每个 $M$，重复 $T = 50$ 次独立估计，计算均方根误差 RMSE $= \sqrt{\frac{1}{T}\sum_t (\hat{I}^{(t)}_d - d/3)^2}$
4. 在 log-log 图上画出 RMSE vs $M$ 的三条曲线（$d = 1, 5, 20$），叠加理论线 $\propto 1/\sqrt{M}$

**定量指标**：
- log-log 图中各曲线的斜率（期望 $\approx -0.5$）
- 固定 $M = 10^4$ 时，$d=1$ 和 $d=20$ 的 RMSE 之比

**预期结果**：
- 三条曲线斜率为 $-1/2$（在 log-log 图上），验证 $O(1/\sqrt{M})$ 收敛
- 收敛速度几乎不随 $d$ 变化（三条曲线接近平行）
- 与梯形法则对比：梯形法则的收敛速度为 $O(M^{-2/d})$，在高维时远不如 Monte Carlo

---

### 7.3 实验 3：贝叶斯决策——高斯均值估计中的风险对比

**目的**：在一个简单的贝叶斯决策场景中，对比 MLE 估计和贝叶斯后验均值估计的 MSE 性能。

**步骤**：
1. 构造数据生成模型：$X_i \sim \mathcal{N}(\theta, \sigma^2)$，$i = 1,\dots,N$，其中 $\sigma^2 = 1$ 已知
2. 真实参数 $\theta_{\text{true}} = 2$，施加先验 $\theta \sim \mathcal{N}(0, \tau^2)$
3. 对固定样本量 $N = 5$，采样 $T = 1000$ 组数据
4. 对每组数据计算两种估计：
   - MLE：$\hat{\theta}_{\text{MLE}} = \bar{X} = \frac{1}{N}\sum X_i$
   - 贝叶斯后验均值：$\hat{\theta}_{\text{Bayes}} = \frac{\tau^2}{\tau^2 + \sigma^2/N}\bar{X}$（已知后验 $\mathcal{N}\left(\frac{\tau^2}{\tau^2+\sigma^2/N}\bar{X}, \frac{\tau^2\sigma^2/N}{\tau^2+\sigma^2/N}\right)$）
5. 计算两种估计的 MSE $= \frac{1}{T}\sum_t (\hat{\theta}^{(t)} - \theta_{\text{true}})^2$
6. 分别取 $\tau^2 = 1$（弱先验）和 $\tau^2 = 0.1$（强先验），重复实验

**定量指标**：
- MLE 的理论 MSE：$\text{MSE}_{\text{MLE}} = \sigma^2/N$
- 贝叶斯估计的理论 MSE：取决于先验方差 $\tau^2$ 和真实 $\theta_{\text{true}}$ 的偏差
- 实验 MSE vs 理论 MSE 的吻合度

**预期结果**：
- MLE 是无偏的，但方差较大（$1/N$）
- 贝叶斯估计有偏（向 0 收缩），但方差更小
- 当 $\tau^2$ 较小时，贝叶斯估计的 MSE 优于 MLE（偏差增量小于方差减少量）
- 当 $\theta_{\text{true}}$ 偏离先验均值 0 较远时，MMSE 估计的偏差代价大于方差收益



