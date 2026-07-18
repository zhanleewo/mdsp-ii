<div style="page-break-before: always; padding: 8% 8% 0 8%;">
 <h1 id="第一讲-MUSIC算法" style="text-align: center; margin-bottom: 2rem; border-bottom: none; display: block;">第一讲 MUSIC算法</h1> 
 <div style="display: flex; align-items: center; justify-content: center; gap: 12px; margin: 1.8rem auto;">
  <span style="flex:1; max-width:80px; height:1px; background: linear-gradient(to right, transparent, #888);"></span>
  <span style="display:inline-block; width:6px; height:6px; background:#38bdf8; border-radius:50%;"></span>
  <span style="flex:1; max-width:80px; height:1px; background: linear-gradient(to left, transparent, #888);"></span>
 </div>
</div>

<!-- # 第一讲 MUSIC算法 -->

> 多重信号分类算法与子空间测向
> 
> Multisensor Signal Classification Algorithm and Subspace Beamforming
> 
> 奠基人: Robert O. Schmidt
> 
> Multiple Emitter Location and Signal Parameter Estimation

进入 MUSIC 算法之前，先明确"阵列"是什么。

《现代数字信号处理 I》处理的对象是时间序列（单个传感器，采集一串时间样本）。《现代数字信号处理 II》的开篇——最优阵列处理——将视角从时间轴切换到空间轴（多个传感器在同一时刻采集空间样本）。

这是从一维信号处理到多维信号处理的扩展。下面先建立阵列和阵列处理的基本概念。

---

## 1. 阵列信号处理基础

### 1.1 阵列的定义与物理直觉

#### 1.1.1 物理直觉：多传感器协同与空间采样

想象一下，你闭上一只眼睛看世界，你只能看到二维画面，很难判断物体的远近。但当你睁开两只眼睛（双目视觉），大脑通过两只眼睛看到的微小差异（视差），就能瞬间算出物体的距离和深度。

阵列就是干这个用的——在空间中按特定几何位置排列多个传感器（天线/麦克风），利用它们之间的位置差来获取信号的空间信息，最直接的就是方向信息。

#### 1.1.2 典型应用案例：双耳定位、相控阵雷达、麦克风阵列

**例子 1：人类的双耳效应**
你闭着眼睛也能判断身后汽车喇叭声是从左边还是右边来的——因为两只耳朵收到的声音有时间差和强度差。大脑利用这个相位差算出方向，这就是一个天然的 2 元阵列。

**例子 2：相控阵雷达**
军事上的相控阵雷达是阵列的典型应用。天线面上排列着成千上万个发射/接收单元，雷达不转动，而是通过控制每个单元发射信号的相位，让波束在空间中瞬间指向任意方向——这就是电子扫描，比机械扫描快得多。

**例子 3：麦克风阵列**
智能音箱顶部有一圈麦克风。当你喊唤醒词时，它通过计算声音到达各麦克风的时间差来定位你的方向，然后增强该方向的接收灵敏度，抑制其他方向的噪声。

#### 1.1.3 数学定义：阵元、阵列流型与空间采样系统

在专业信号处理的语言中，阵列是这样定义的：

> **阵列（Array）** 是由多个相同的传感器（阵元）按照一定的空间几何布局（直线、圆环、平面或任意形状）组成的**空间采样系统**。它用于接收（或发射）在空间中传播的波（电磁波或声波）。

3.  **阵列流型（Array Manifold）：** 这是阵列处理最核心的"指纹"。阵列流型是一个函数，它描述了**对于来自某个特定方向 $\theta$ 的单一平面波，阵列上每个阵元接收到的信号的复数表达式（包含幅度和相位）**。该方向的信号到达不同阵元会走过不同的路程，产生不同的相位延迟，这种差异构成的向量，就是该方向的"标签"。只要标签唯一，我们就能逆向定位。

**一句话：阵列是"空间中的频谱分析仪"——傅里叶变换分析信号的频率成分，阵列分析信号的角度成分。**

---

### 1.2 阵列处理的两大核心任务：波束形成与到达角估计

把多个阵元摆在空间里，我们自然要问两个基本的问题。这两个问题构成了阵列处理的全部版图。

#### 1.2.1 直观理解：空间滤波（波束形成）与方向估计（测向）

阵列处理就干两件事：

1.  **把想要的方向的信号"捞"出来，把不想要的方向的干扰"滤"掉。**（就像在嘈杂的晚宴上，把听筒对准说话的人，屏蔽周围碰杯的声音）
2.  **找出信号到底是从哪个方向来的。**（就像飞机上的黑匣子搜寻器，必须精确找到水下信标的方向）

#### 1.2.2 典型场景：体育场 5G 基站天线阵列

**场景：** 大型体育场的 5G 基站天线阵列。

-   **波束形成（阵列的第一件事）：** 球馆里有上万名观众同时用手机看直播。基站上的大规模天线阵列（Massive MIMO）通过调节每个天线的权值，形成数十个极窄的"聚光灯"一样的波束，一个波束对准一群观众，定向发射信号。这么做不仅功率集中、传输距离远，还不会干扰其他波束里的用户。**这就是"空间滤波"——在物理空间上把人群分开了。**
-   **到达角估计（阵列的第二件事）：** 假设球馆里有一个来历不明的无线电干扰源干扰了通信。基站阵列接收信号后，利用算法算出该信号的来向（比如北偏东 23°），然后基站向指挥中心报告："干扰源在那个方向。"**这就是"测向"——估计空间中的角度参数。**

#### 1.2.3 理论框架与方法演进

在《最优阵列处理》的理论框架下，上述两件事分别对应以下两个核心分支：

##### 分支一：波束形成（Beamforming）—— 最优空域滤波

波束形成是**空间维度的滤波器设计**。时间滤波（维纳滤波，参考数字信号处理I）用的是抽头延迟线，空间滤波用的是**阵元空间加权**。

> **波束形成的本质：** 对各个阵元接收到的信号进行**复加权（幅度和相位调整）**并求和。通过设计不同的权向量 $\mathbf{w}$，使得阵列的方向图在期望信号方向增益最大，在干扰方向形成"零陷"（Nulling）。

这一分支的演进路径极为清晰，呼应了数字信号处理 I 中的线性最优估计：

-   **常规波束形成（CBF / Bartlett）：** 这是"时域匹配滤波"的空间版本。权向量就是导向矢量的共轭。它只做相位补偿，不处理干扰，分辨率受限于阵列孔径（类似时间域的 Rayleigh 限）。
-   **最优波束形成（MVDR / Capon）：** 这是"维纳滤波"的空间版本。在保证期望方向增益为 1 的约束下，最小化输出总功率（即抑制干扰和噪声）。它利用数据的协方差矩阵进行自适应滤波。
-   **稳健波束形成（Robust Beamforming）：** 当阵列流型存在误差或样本不足时，MVDR 会退化（类似数字信号处理I中提到的协方差矩阵病态问题）。于是引入"对角加载（Diagonal Loading）"——这和我们前面学过的**岭回归（正则化）**在数学上完全相同，都是在协方差矩阵上加一个 $\lambda \mathbf{I}$ 来换取稳健性。

##### 分支二：到达角估计（DOA Estimation）—— 空间谱估计

到达角估计是**空间频域的谱分析**。我们既然在"谱分析"一章花了大力气学周期图、AR 谱、Capon 谱，那么到了阵列里，面对的不是"时间频率 $\omega$"，而是"空间角度 $\theta$"。我们同样要把能量在角度域展开。

这一分支的演进路径同样沿着"非参数 -> 参数/超分辨"：

-   **延迟-相加法（常规 DOA）：** 对应周期图法。扫描整个角度空间，哪个角度输出的能量大，就认为信号从哪里来。分辨率受限于"瑞利限"（$\Delta\theta \propto \lambda / D$，其中 $D$ 是孔径）。
-   **Capon 谱估计（空间版本）：** 和数字信号处理I中完全一样，用在角度域的扫描，得到比常规法更尖的谱峰。
-   **子空间类方法（MUSIC / ESPRIT）：** 利用信号子空间与噪声子空间的正交性，对协方差矩阵进行特征分解，在角度域构造伪谱。它突破了瑞利限，实现了超分辨测向——即使两个信号的角度差远小于波束宽度，也能在谱上分开。

##### 分支三：阵列校准与参数估计（工程基础）

真实阵列存在阵元位置误差、通道不一致性（幅相误差）和互耦效应。如果算法基于理想的阵列流型而实际流型有误差，MUSIC 或 MVDR 也会失效。因此，阵列处理必须包含校准技术——利用已知方位的校正源来估计和补偿阵列流型误差。

---

### 1.3 知识迁移：数字信号处理 I（时域/频域）→ II（空域/角度域）

将两门课的对应关系整理如下：

| 维度 | 数字信号处理 I（时域/频域） | 数字信号处理 II（空域/角度域） |
| :--- | :--- | :--- |
| **信号载体** | 时间 $t$（采样间隔 $T_s$） | 空间位置 $d$（阵元间距 $d$） |
| **变换域** | 频率 $\omega$ | 角度 $\theta$（等价于波数 $k$） |
| **线性滤波器** | 维纳滤波（抽头延迟线） | **波束形成**（阵元空间加权） |
| **谱估计** | 周期图 / AR 谱（功率随频率变化） | **空间谱估计**（功率随角度变化） |
| **超分辨率** | AR 模型 / 最大熵 | **MUSIC / ESPRIT** |
| **核心数学** | 自相关 + Toeplitz 矩阵 | 协方差矩阵 + 特征分解（Hermitian） |
| **正则化** | 岭回归 / 对角加载 | 对角加载（稳健波束形成） |


![array-processing-basic-concepts](assets/01/03.png)
这张图展示了一个典型的**传感器阵列信号处理系统**的工作原理。它描述了如何从复杂的电磁环境中，利用多个传感器协同工作来接收并处理目标信号。

以下是图中各个部分的详细职责及其配合工作的流程：

##### ① 环境中的信号源（输入端）
图的上方展示了空间中存在的各种波场，它们是系统的输入对象：

- **Signal #1 & Signal #2 (信号源)**：
    - **职责**：这是系统想要接收或监测的目标信号（例如来自特定方向的无线电波、声波等）。
    - **特征**：它们以特定的方向传播向传感器阵列。
- **Noise Field (噪声场)**：
    - **职责**：代表背景干扰或无用的随机信号。
    - **特征**：图中的绿色箭头表示噪声可能来自四面八方，弥漫在整个空间中，会混杂在目标信号中一起被接收。
- **Interference (干扰)**：
    - **职责**：这是一个特定的强干扰源（紫色箭头），它不是随机的背景噪声，而是一个具体的、可能阻碍信号接收的外部因素。

##### ② Sensor Array (传感器阵列)
这是系统的"耳朵"或"眼睛"，位于图的中间部分。

- **组成**：由多个按特定几何形状排列的独立传感器（蓝色圆点）组成。
- **职责**：
    - **空间采样**：每个传感器独立地接收空间中的混合波形（信号+噪声+干扰）。
    - **捕捉差异**：由于波到达不同传感器的时间和相位存在微小差异（取决于波的方向），阵列能够捕捉到这些空间上的细微差别。这是后续进行方向估计和波束形成的基础。

##### ③ 传输链路
连接传感器和处理器的线条。

- **职责**：将每个传感器接收到的模拟电信号实时传输给后端的处理器。每一根线代表一个独立的通道。

##### ④ Array Signal Processor (阵列信号处理器)
这是系统的"大脑"，位于图的右下角。

- **职责**：对来自所有传感器的数据进行复杂的数学运算和算法处理。
- **核心功能**：
    - **波束形成**：通过调整各通道信号的权重（延时和相加），增强来自特定方向（如Signal #1或#2）的信号，同时抑制其他方向的噪声。
    - **干扰抑制**：识别并滤除那个紫色的"Interference"干扰源。
    - **信噪比提升**：从混杂的"Noise Field"中提取出清晰的有用信号。
    - **参数估计**：计算出信号到来的方向、频率等特征。

---

##### 总结：它们是如何配合工作的？

整个系统的工作流程可以概括为：**"多点接收 -> 数据传输 -> 联合处理"**。

1.  **混合接收**：空间中的**信号**、**噪声**和**干扰**同时传播到**传感器阵列**。每个传感器都会接收到一份包含所有成分的"混合数据"，但由于位置不同，每份数据的相位和幅度略有不同。
2.  **并行传输**：这些混合数据通过导线被同步传送到**阵列信号处理器**。
3.  **智能分离**：处理器利用算法（如波束成形技术），分析这些数据之间的时空相关性。它能够"算出"哪个方向是想要的信号，哪个方向是干扰。
4.  **最终输出**：处理器通过数学手段"指向"目标信号方向（增强Signal #1或#2），并在反方向形成"零陷"以抵消干扰，最终输出一个清晰、纯净的目标信号结果（图右侧的输出箭头）。

这种架构广泛应用于雷达探测、声呐系统、5G/6G通信基站以及射电天文望远镜等领域。

---

这里提到的"波束形成"是实现抗干扰的关键技术，需要我展开讲讲它的原理吗？

![array-processing-basic-concepts](assets/01/01.webp)
这张图非常直观地展示了**传感器阵列（Sensor Array）的核心能力——波束形成（Beamforming）**。

简单来说，这张图表达了：**通过控制多个传感器协同工作，可以将能量或接收灵敏度集中在一个特定的方向上，就像手电筒的光束一样，而不是像灯泡那样向四面八方发散。**

以下是图中各个部分的详细解读：

##### ① Sensor Array (传感器阵列 - 左下角的绿色点阵)
- **是什么**：这是一组排列整齐的传感器（比如天线、麦克风或声纳探头）。
- **作用**：它们是系统的物理基础。单个传感器的信号通常比较弱且没有方向性（像个普通灯泡），但当它们组合在一起并经过精确的信号处理（调整相位和幅度）后，就能产生"定向"效果。

##### ② Radiation Pattern (辐射方向图 - 中间的蓝紫色花瓣状图形)
这是整张图最核心的概念，它展示了阵列工作的**结果**。
- **主瓣 (Main Lobe - 那个最长的紫红色尖峰)**：
    - **含义**：这是能量最集中、信号最强的方向。
    - **比喻**：就像手电筒聚光后的光束，或者你把手拢在嘴边喊话时声音传得最远的那个方向。在这个方向上，系统发射的能量最大，或者接收信号的灵敏度最高。
- **旁瓣/副瓣 (Side Lobes - 周围那些较小的蓝色花瓣)**：
    - **含义**：这是能量泄漏的方向。虽然我们希望能量全部集中在主瓣，但物理上总会有一些能量分散到其他地方。这些小的波瓣就是干扰或能量浪费的区域，通常在工程设计中需要尽量抑制它们。

##### ③ Target (目标 - 右上角的黑色物体)
- **是什么**：这是系统想要探测、通信或打击的对象（例如飞机、潜艇、或者手机基站）。
- **关系**：注意看，**Radiation Pattern 的主瓣（紫红色尖峰）精准地指向了 Target**。这意味着系统已经"锁定"了目标，正在把绝大部分能量发射给它，或者正全神贯注地监听来自那个方向的微弱回波。

##### ④ 波纹线 (Wavefronts - 连接阵列和目标的弧线)
- **含义**：这代表正在传播的电磁波或声波。
- **动态过程**：这些波纹从传感器阵列发出，汇聚成一股能量流，径直冲向目标。

---

##### 总结：这张图到底想说什么？

这张图用可视化的方式解释了**相控阵雷达**或**定向通信技术**的基本原理：

1.  **能量聚焦**：它展示了如何将分散的能量"捏"成一股绳，指向特定目标，从而实现远距离探测或高效通信。
2.  **空间滤波**：只有在这个紫色尖峰指向的角度（Target所在的角度），信号才是最强的；其他角度的信号会被抑制。这就是为什么雷达能分辨出目标的具体方位。
3.  **电子扫描**：虽然图是静止的，但在实际应用中，不需要转动天线，只需要改变传感器阵列内部的电子信号相位，这个紫色的尖峰（主瓣）就可以瞬间指向天空中的任何一个目标。

![array-processing-basic-concepts](assets/01/06.png)
这张图揭示了**波束形成（Beamforming）技术背后的数学原理和硬件实现逻辑**。它展示了如何通过物理上的"延时"和数学上的"加权"，让多个天线接收到的信号在特定方向上叠加增强。

简单来说，它的核心思想是：**利用波到达不同天线的微小时间差，通过调整每个通道的权重，让来自目标方向的信号"步调一致"地相加，从而获得最强的输出。**

以下是图中各个关键部分的详细解读：

##### ① 几何结构与波前 (Wave Front & Geometry)
- **天线阵列**：图上方的一排倒三角形代表均匀排列的天线单元。
- **间距 $d$**：相邻两个天线之间的距离是固定的，标记为 $d$。
- **波前 (Wave Front)**：那条斜着的虚线代表平面波的波前。这意味着信号是从某个特定的角度 $\theta_0$ 射入的。
- **路径差**：这是最关键的几何关系。由于波是斜着过来的，它到达每个天线的时间是不一样的。
    - 到达第2个天线比第1个多走的距离是 $d \sin\theta_0$。
    - 到达第3个天线比第1个多走的距离是 $2d \sin\theta_0$。
    - 到达第 $N$ 个天线多走的距离是 $(N-1)d \sin\theta_0$。
    - **物理意义**：这段距离差导致了信号到达每个天线时存在**相位差**。如果不处理，直接把这些信号加起来，它们可能会因为相位不同而互相抵消，而不是增强。

##### ② 加权器 (Weights: $w^*$)
- **圆圈里的叉号 ($\times$)**：这代表乘法器。
- **$w_1^*, w_2^*, \dots, w_N^*$**：这些是复数权重系数。
    - **作用**：这是波束形成的核心算法部分。为了补偿上面提到的"路径差"带来的相位滞后，我们需要给每个天线的信号乘以一个特定的权重。
    - **原理**：这个权重通常包含一个相移项，用来抵消信号传播带来的相位差。这就好比为了让一队人齐步走，让走得慢的人（信号晚到的）先迈腿，或者给他们的步伐加一个修正值，最终让大家在同一点汇合时步调完全一致。

##### ③ 求和器 (Summer: $\Sigma$)
- **长方形框 $\Sigma$**：代表将所有经过加权处理后的信号进行相加。
- **建设性干涉**：当权重 $w$ 设置正确（即针对角度 $\theta_0$ 进行了匹配）时，来自 $\theta_0$ 方向的所有信号在求和时会相位相同，发生**建设性干涉**，信号幅度达到最大（如果是 $N$ 个天线，幅度理论上变为 $N$ 倍）。
- **破坏性干涉**：而对于其他方向的噪声或干扰信号，由于它们的相位差没有被补偿，加起来后会杂乱无章，甚至互相抵消，从而被抑制。

##### ④ 波束形成器输出 (Beamformer Output)
- 这是最终处理好的信号。在这个输出端，来自目标方向 $\theta_0$ 的信号被极大地增强了，而其他方向的干扰被滤除了。这就实现了我们在上一张图中看到的"像手电筒一样指向特定方向"的效果。

---

##### 总结三张图的逻辑联系
结合你之前发的两张图，我们可以梳理出一个完整的知识链条：

1.  **第一张图（宏观场景）**：展示了**为什么要用阵列**——为了在复杂的噪声和干扰环境中提取出微弱的目标信号。
2.  **第二张图（物理效果）**：展示了**阵列做到了什么**——形成了一个高增益的辐射方向图（主瓣），像探照灯一样精准地指向目标。
3.  **第三张图（微观原理，即本图）**：展示了**具体是怎么做到的**——通过计算波程差，利用加权求和的数学手段，在电路层面实现了信号的定向增强。

这三张图从应用场景、物理现象到数学原理，完整地解释了传感器阵列与波束形成技术。

---

#### 知识迁移小结：时域/频域工具到空域/角度域的完整映射

阵列处理的数学工具与数字信号处理 I 一一对应：

- 时域采样定理 → 空间采样（阵元间距 ≤ λ/2）
- 傅里叶变换 → 波数变换（角度-相位关系）
- 自适应滤波（LMS/RLS）→ 自适应波束形成
- CRLB → 角度估计的 CRLB

理解了这些对应关系，就可以用数字信号处理 I 中建立的数学直觉（投影、特征分解、正则化、谱峰搜索）来处理空间域的问题。
## 2. 均匀线性阵列（Uniform Linear Array, ULA）

### 2.1 基本概念

均匀线性阵列（ULA）是最简单、最经典、理论最成熟的阵列几何结构：\( N \) 个相同的阵元等间距 \( d \) 排列在一条直线上。该结构因其阵元位置构成等差数列，导向矢量具有范德蒙德结构，使得后续的波束形成、MUSIC、ESPRIT 等算法在 ULA 下均存在高效的闭式实现。

---

### 2.2 信号传播模型：远场平面波假设与窄带近似

#### 2.2.1 假设一：平面波（远场假设）

![plane wave](assets/01/08.jpg)

**（1）为什么能够假定为平面波——远场条件**

电磁波或声波源离阵列足够远时，波前从球面退化为平面。当信号源距离 \( R \) 远大于阵列孔径 \( D \) 时，到达不同阵元的波前曲率可以忽略。

工程上常用的远场判据为瑞利距离：

$$
R > \frac{2D^2}{\lambda}
\tag{1.1} $$

满足此条件时，实际球面波前与理想平面波前的相位误差小于 \( \pi/8 \)，平面波假设成立。

若信号源位于近场，波前呈球面，不同阵元接收到的信号不仅存在时延差，幅度也有差异（传播损耗 \( 1/R^2 \) 不同），且时延差与距离的平方有关。此时问题是距离-角度联合估计，属于近场源定位，本单元不涉及。

**（2）同时采样与空间差异**

阵列在时间采样上通常假定完全并行，即所有阵元在同一时刻完成采样。同时采样导致不同阵元的采样值差异完全来源于空间位置差异，而不混入时间先后因素。这就是"空间快拍"的含义——同一时刻不同空间点的数据。

设信号为 \( s(t) \)，以第一个阵元为参考点：

- 阵元 1 接收信号为 \( s(t) \)；
- 阵元 2 接收信号为 \( s(t+T) \)；
- 阵元 3 接收信号为 \( s(t+2T) \)；
- 阵元 \( N \) 接收信号为 \( s(t+(N-1)T) \)。

**（3）时延公式 \( T = d\cos\theta/c \) 的含义**

此处 \( \theta \) 定义为信号传播方向与阵列轴线（阵元排列方向）的夹角。相邻阵元间距为 \( d \)。平面波沿与轴线夹角 \( \theta \) 的方向传播时，到达相邻阵元所需额外传播路程为 \( d\cos\theta \)。

因此时延为：

$$
T = \frac{d\cos\theta}{c}
\tag{1.2} $$

其中：
- \( d \)：阵元间距（单位：米）；
- \( \cos\theta \)：将阵元间距投影到波传播方向上的因子；
- \( c \)：波速（电磁波为光速 \( 3\times 10^8 \text{ m/s} \)，声波为 \( 340 \text{ m/s} \)）。

**关于角度定义的说明：** 部分文献将 \( \theta \) 定义为与阵列法线的夹角，此时公式变为 \( T = d\sin\theta/c \)。两种定义本质相同，本书统一采用与轴线夹角定义，后续所有推导均基于此，不再另行说明。

**（4）核心推理链：从数据差异到角度估计**

这里出现了一个关键问题：\( \theta \) 是未知的。无论是通信还是雷达，这个角度正是想要搞清楚的目标。\( \theta \) 决定了时延 \( T \)，时延 \( T \) 导致了各阵元接收数据之间的差异。阵列处理把这个过程反过来：利用数据的差异算出时延，利用时延算出角度。这一推理链条可概括如下：

```mermaid
graph LR
    A[各阵元并行采样<br>获得数据向量] --> B[利用算法估计<br>阵元间时延差 T]
    B --> C["利用几何关系<br>$$ T = d\cos\theta/c $$"]
    C --> D["反解出到达角 $$ \theta $$"]
    style A fill:#f9f,stroke:#333,stroke-width:2px
    style D fill:#9cf,stroke:#333,stroke-width:2px
```


#### 2.2.2 假设二：窄带信号（Narrowband）

![narrowband](assets/01/09.jpg)

**（1）什么是窄带信号**

窄带信号是指信号带宽 \( B \) 远小于中心频率（载频）\( f_c \)：

$$
B \ll f_c
\tag{1.3} $$

**（2）高载频的原因——天线尺寸与带宽**

通信与雷达系统使用高载频（2.4 GHz、5 GHz、10 GHz、35 GHz 等），主要受天线尺寸和带宽需求的驱动。

天线要高效辐射电磁波，其物理尺寸必须与波长 \( \lambda = c/f \) 同量级（通常为 \( \lambda/2 \) 或 \( \lambda/4 \)）。若载频 100 MHz（\( \lambda = 3 \) m），天线需长约 1.5 m，无法集成于便携设备。载频升至 2.4 GHz（\( \lambda \approx 12.5 \) cm），天线仅需约 6 cm，可内置于手机中（如 IFA、PIFA 等微带天线）。

此外，高载频在满足窄带假设（相对带宽小）的同时仍能提供较大的绝对带宽。例如，相对带宽 1% 在 \( f_c = 100 \) MHz 时仅为 1 MHz，在 \( f_c = 100 \) GHz 时则为 1 GHz。注意：绝对带宽由数据速率与调制方式决定，与载频高低无直接因果关系。

**（3）窄带信号的数学表示**

窄带信号可写为：

$$
s(t) = a(t) \cdot \exp(j 2\pi f t)
\tag{1.4} $$

其中 \( a(t) \) 是慢变的复包络（包含幅度和相位调制信息），\( \exp(j 2\pi f t) \) 是快变的载波。\( a(t) \) 是基带信号，经过调制后成为带通信号 \( s(t) \)。

**（4）时延 \( T \) 下的窄带近似**

现在考察经过时延 \( T \) 后的信号：

$$
s(t+T) = a(t+T) \cdot \exp\big(j 2\pi f (t+T)\big)
\tag{1.5} $$

展开载波项：

$$
s(t+T) = a(t+T) \cdot \exp(j 2\pi f t) \cdot \exp(j 2\pi f T)
\tag{1.6} $$

由于窄带条件 \( B \ll f_c \) 且时延 \( T \) 很小（\( T \ll 1/B \)），包络 \( a(t) \) 在该时间尺度内变化可忽略：

$$
a(t+T) \approx a(t)
\tag{1.7} $$

于是得到窄带近似的核心公式：

$$
s(t+T) \approx a(t) \cdot \exp(j 2\pi f t) \cdot \exp(j 2\pi f T) = s(t) \cdot \exp(j 2\pi f T)
\tag{1.8} $$

窄带假设下，\( s(t+T) \) 不再是一般性的时移操作，而是退化为原信号乘以一个相位因子 \( \exp(j 2\pi f T) \)。换言之，时域时移在窄带条件下等价于复平面上的旋转。

**（5）为什么 \( a(t+T) \approx a(t) \) 是合理的——两个层面的解释**

**确定性信号层面：** 根据傅里叶变换的基本性质，频域窄（\( B \) 小）意味着时域宽。窄带信号的频谱集中在 \( f_c \) 附近很窄的区间 \( [f_c - B/2,\ f_c + B/2] \)，其时域包络的持续时间约为 \( 1/B \)，变化极为缓慢。当时间偏移 \( T \) 远小于 \( 1/B \) 时，包络几乎没有变化，即 \( a(t+T) \approx a(t) \)。这从时域直接证明了上述近似的合理性。

**随机信号层面（相关函数的平坦性）：** 对于平稳随机信号，自相关函数 \( R(\tau) = E[x(t)x^*(t-\tau)] \) 与功率谱密度 \( S(f) \) 构成傅里叶变换对：

$$
R(\tau) \xleftrightarrow{\mathscr{F}} S(f)
\tag{1.9} $$

窄带信号的功率谱 \( S(f) \) 仅在 \( [f_c - B/2,\ f_c + B/2] \) 内非零，且 \( B \) 很小。傅里叶变换的对偶性质表明：**频域越窄（\( B \) 小），时域相关函数 \( R(\tau) \) 的主瓣越宽（\( \approx 1/B \)）**。

- 若 \( B \) 很大（宽带信号），相关函数主瓣很窄，\( \tau \) 稍有偏离，\( R(\tau) \) 便迅速衰减至零，说明前后时刻几乎不相关。
- 若 \( B \) 很小（窄带信号），相关函数主瓣很宽，即使 \( \tau \) 有较大偏移，\( R(\tau) \) 仍接近 \( R(0) \)（信号功率），说明前后时刻高度相关。

时延 \( T \) 远小于 \( 1/B \)，位于相关函数主瓣内部，故 \( R(T) \approx R(0) \)。相关函数平坦意味着相关度高，相关度高意味着 \( x(t) \) 与 \( x(t+T) \) 几乎相等。这从统计意义上精确说明了窄带信号包络可近似相等的理论依据。

**（6）窄带假设的工程意义**

若不采用窄带假设，时延 \( T \) 包含在 \( s(t+T) \) 内部。估计 \( T \) 需要知道信号的具体波形，或进行互相关、反卷积等操作，且处理方法依赖于信号类型。

窄带假设将 \( T \) 从时域函数内部提取出来，变为一个独立的相移因子：

$$
s(t+T) = s(t) \cdot \exp(j\phi), \quad \phi = 2\pi f T\tag{1.10} $$

于是，估计时延 \( T \) 转化为估计相位 \( \phi \)。相位估计只需通过阵列各阵元的相位差来完成，不需要知道信号的具体波形。因此，只要满足窄带条件，算法与信号类型无关。

**窄带假设将非线性时延估计转化为线性相位估计，这是经典阵列处理算法的数学前提。**

---
### 2.3 方向矢量（Steering Vector）

方向矢量是阵列信号处理的核心概念，描述了方向信息在阵列数据中的编码方式。

假设阵列共有 \( M \) 个阵元。以第一个阵元为参考点，接收信号为 \( s(t) \)。根据前文的窄带假设，第 \( i \) 个阵元接收到的信号为 \( s(t + iT) \)，其中 \( T = d\cos\theta/c \) 是相邻阵元间的时延差。于是整个阵列在同一时刻接收到的信号向量为：

$$
\big( s(t),\ s(t+T),\ \cdots,\ s(t+(M-1)T) \big)
\tag{1.11} $$

利用窄带近似的核心结论 \( s(t+T) = s(t) \cdot \exp(j 2\pi f T) \)，可将上述向量改写为：

$$
\begin{aligned}
& \big( s(t),\ s(t)\exp(j 2\pi f T),\ \cdots,\ s(t)\exp(j 2\pi f (M-1)T) \big) \\
= & s(t) \cdot \underbrace{\big( 1,\ \exp(j 2\pi f T),\ \cdots,\ \exp(j 2\pi f (M-1)T) \big)}_{\text{方向矢量（Steering Vector）}}
\end{aligned}
\tag{1.12} $$

记方向矢量为：

$$
s(t) \big( \phi_0(T),\ \phi_1(T),\ \cdots,\ \phi_{M-1}(T) \big)
\tag{1.13} $$

其中：

$$
\phi_i(T) = \exp(j 2\pi f \cdot iT), \quad i = 0, 1, \cdots, M-1\tag{1.14} $$

这里 \( \phi_0(T) = 1 \) 作为参考相位。

上述表达式说明，在窄带假设下，信号 \( s(t) \) 与空间相位实现了完全解耦。

回顾窄带信号 \( s(t) = a(t)\exp(j 2\pi f t) \)，\( a(t) \) 是慢变复包络，承载信号的调制信息和时域结构。经过时延 \( T \) 后：

$$
s(t+T) = a(t+T)\exp(j 2\pi f (t+T))
\tag{1.15} $$

窄带条件 \( B \ll f \) 保证包络在时延 \( T \) 内近似不变（\( a(t+T) \approx a(t) \)），于是：

$$
s(t+T) \approx a(t)\exp(j 2\pi f t)\exp(j 2\pi f T) = s(t)\exp(j 2\pi f T)
\tag{1.16} $$

时延 \( T \) 的作用仅体现在复指数因子 \( \exp(j 2\pi f T) \) 中。信号内容 \( a(t) \) 在时延前后不变，被提取为所有阵元的公共因子 \( s(t) \)。

对整个阵列：各阵元的接收信号由一个公共因子 \( s(t) \) 乘以一个仅依赖于阵元位置 \( i \) 的复指数因子构成。\( s(t) \) 包含信号的全部时域特征，而阵元间的差异——即反推角度 \( \theta \) 的信息——全部集中在 \( \exp(j 2\pi f \cdot iT) \) 中。公共因子 \( s(t) \) 在后续处理中作为标量提取，不参与方向估计。

因此，方向矢量 \( \big( 1,\ \exp(j 2\pi f T),\ \cdots,\ \exp(j 2\pi f (M-1)T) \big) \) 完全由空间几何参数决定——载频 \( f \)、阵元间距 \( d \)、波速 \( c \)、来波方向 \( \theta \)——与信号的具体内容无关。只要窄带条件成立且载频相同，方向矢量与信号类型无关。

上述推导基于 ULA，相邻阵元间距相等，相位序列构成等差数列。若阵列为非均匀布局或位于二维/三维空间，方向矢量不再具有等差数列结构，但其本质不变——始终是各阵元相对于参考阵元的空间相位延迟向量。采用 \( \phi_i(T) \) 记号可以统一涵盖任意阵列几何：对于非均匀阵列，将 \( \phi_i(T) \) 替换为 \( \phi_i(\mathbf{p}_i, \mathbf{k}) = \exp(-j\mathbf{k}\cdot\mathbf{p}_i) \) 即可。


### 2.4 多信号源的阵列接收模型与方向矩阵

![multiple signals](assets/01/10.png)

现考虑空间中有 \( N \) 个信号源同时存在的一般情形。假定信号在空间中服从线性叠加原理——这是电磁波传播的基本物理定律——则阵列接收到的总场是各源信号的线性叠加。

设第 \( k \) 个信号源为 \( s_k(t) \)，其对应的相邻阵元时延为 \( T_k \)，则该信号在阵列上的方向矢量为：

$$
\big( 1,\ \exp(j 2\pi f T_k),\ \cdots,\ \exp(j 2\pi f (M-1) T_k) \big)
\tag{1.17} $$

其中每个 \( T_k \) 都包含了该信号源方向 \( \theta_k \) 的全部信息（\( T_k = d\cos\theta_k/c \)）。

于是，第 \( k \) 个信号源单独存在时，阵列接收向量为：

$$
s_k(t) \cdot \big( 1,\ \exp(j 2\pi f T_k),\ \cdots,\ \exp(j 2\pi f (M-1) T_k) \big)
\tag{1.18} $$

各信号源的内容 \( s_k(t) \) 互不相同——它们可能来自不同的发射机，携带不同的信息——但方向信息仍然只编码于各自的相位项 \( \exp(j 2\pi f \cdot iT_k) \) 中。2.3 节所揭示的"信号与相位解耦"性质，在多个信号源的场景下依然成立：每个信号源的内容只贡献了一个共同因子，而方向估计所需的信息全部包含在其对应的相位序列之中。

然而，上述以"每个信号源"为单位的组织形式虽然概念清晰，却不便于直接用于分析实际问题。原因在于，阵列的真实采样过程并非按信号源维度进行——实际数据采集是以阵元为基本单位的。每个阵元上采集到的瞬时数据，是所有信号源在该阵元位置处的场强叠加。

以第 1 个阵元（参考阵元）为例，其接收信号为所有信号源在同一时刻的瞬时值之和：

$$
X_1(t) = s_1(t) + s_2(t) + \cdots + s_N(t)
\tag{1.19} $$

以第 2 个阵元为例，其接收信号为每个信号源经过一个时延后的叠加。利用窄带近似将时延转换为相移：

$$
\begin{aligned}
X_2(t) &= s_1(t+T_1) + s_2(t+T_2) + \cdots + s_N(t+T_N) \\
&= s_1(t)\exp(j 2\pi f T_1) + s_2(t)\exp(j 2\pi f T_2) + \cdots + s_N(t)\exp(j 2\pi f T_N)
\end{aligned}
\tag{1.20} $$

依此类推，第 \( M \) 个阵元接收到的信号为：

$$
\begin{aligned}
X_M(t) =& s_1(t+(M-1)T_1) + s_2(t+(M-1)T_2) + \cdots + s_N(t+(M-1)T_N) \\
=& s_1(t)\exp(j 2\pi f (M-1)T_1) + s_2(t)\exp(j 2\pi f (M-1)T_2) \\
& + \cdots + s_N(t)\exp(j 2\pi f (M-1)T_N)
\end{aligned}
\tag{1.21} $$

将以上所有阵元的接收数据排列成一个列向量 \( \mathbf{X}(t) \)：

$$
\mathbf{X}(t) =
\begin{pmatrix}
X_1(t) \\
X_2(t) \\
\vdots \\
X_M(t)
\end{pmatrix}
\tag{1.22} $$

将各信号源的复包络排列成列向量 \( \mathbf{S}(t) \)：

$$
\mathbf{S}(t) =
\begin{pmatrix}
s_1(t) \\
s_2(t) \\
\vdots \\
s_N(t)
\end{pmatrix}
\tag{1.23} $$

将各信号源对应的方向矢量按列排列成矩阵，得到**方向矩阵（Steering Matrix）**：

$$
\small
\mathbf{A}(\boldsymbol{\theta}) = 
\begin{pmatrix}
1 & 1 & \cdots & 1 \\
\exp(j 2\pi f T_1) & \exp(j 2\pi f T_2) & \cdots & \exp(j 2\pi f T_N) \\
\vdots & \vdots & \ddots & \vdots \\
\exp(j 2\pi f (M-1) T_1) & \exp(j 2\pi f (M-1) T_2) & \cdots & \exp(j 2\pi f (M-1) T_N)
\end{pmatrix}
\tag{1.24} $$

方向矩阵的每一列都是一个方向矢量，对应一个信号源。列与列之间的差异来源于不同的时延 \( T_k \)，本质上即不同的来波方向 \( \theta_k \)。由于不同方向对应的方向矢量各不相同，方向矩阵 \( \mathbf{A}(\boldsymbol{\theta}) \) 完整地编码了全部 \( N \) 个信号源的空间信息。

实际接收信号中不可避免地包含噪声。设噪声向量为：

$$
\mathbf{N}(t) =
\begin{pmatrix}
n_1(t) \\
n_2(t) \\
\vdots \\
n_M(t)
\end{pmatrix}
\tag{1.25} $$

于是，阵列信号处理的完整数据模型可写为：

$$
\boxed{\mathbf{X}(t) = \mathbf{A}(\boldsymbol{\theta})\, \mathbf{S}(t) + \mathbf{N}(t)}\tag{1.26} $$


**至此，基本模型已经建立：**

$$
\mathbf{X} = \mathbf{A}(\boldsymbol{\theta}) \mathbf{S} + \underbrace{\mathbf{N}}_{\text{Noise}}
\tag{1.27} $$

该模型中，方向矩阵 \( \mathbf{A} \)（含角度参数 \( \theta_k \)）、信号波形 \( \mathbf{S} \)、噪声 \( \mathbf{N} \) 均未知，仅观测数据 \( \mathbf{X} \) 已知。这是一个欠定问题——方程数远少于未知数，直接求解不存在唯一解。

突破口在于信号与噪声的统计特性差异：信号部分在空间上具有确定结构——各信号源的导向矢量张成一个低维信号子空间，维度等于信号源个数 \( N \)；白噪声在各方向上均匀分布，无方向选择性。MUSIC 利用这一差异，对协方差矩阵进行特征分解，将观测空间拆分为信号子空间和噪声子空间，然后利用两者的正交性估计方向。

MUSIC 不需要信号波形 \( \mathbf{S} \) 和噪声 \( \mathbf{N} \) 的具体值，只需利用协方差矩阵的特征结构判断哪个方向的导向矢量与噪声子空间正交——该方向即为信号来向。


## 3. MUSIC 算法：多重信号分类与子空间测向

多信号分辨算法


### 3.1 算法动机：从 PCA 到子空间直觉

**多维分布存在显著的分布方向**

将每次快拍 \( \mathbf{x}(t) \in \mathbb{C}^M \) 视为 \( M \) 维复空间中的一个点。只有白噪声时，这些点在各个方向上的散布均匀——协方差矩阵为 \( \sigma^2 \mathbf{I} \)，各方向特征值相等。

当空间中有 \( K \) 个信号源时，在无噪声理想情况下，所有快拍 \( \mathbf{x}(t) = \sum_{k=1}^K \mathbf{a}(\theta_k) s_k(t) \) 均落在 \( K \) 个导向矢量张成的子空间内，数据在高维空间中被压缩到一个 \( K \) 维子空间（\( K < M \)）。加入噪声后，数据被均匀地扩展到整个 \( M \) 维空间——噪声在所有方向上贡献相同方差 \( \sigma^2 \)，信号仅在 \( K \) 个特定方向上贡献额外方差。

因此，协方差矩阵 \( \mathbf{R}_{XX} \) 的特征值呈现可辨识的结构：\( K \) 个大特征值（信号占优方向），\( M-K \) 个相等的小特征值（纯噪声方向）。PCA 正是通过特征分解来识别这些方差较大的方向。

**PCA 即协方差矩阵的特征分解**

将 \( \mathbf{R}_{XX} \) 进行特征分解：

$$
\mathbf{R}_{XX} = \sum_{i=1}^{M} \lambda_i \mathbf{e}_i \mathbf{e}_i^H, \quad \lambda_1 \geq \lambda_2 \geq \cdots \geq \lambda_M
\tag{1.28} $$

其中 \( \lambda_i \) 是特征值，\( \mathbf{e}_i \) 是对应的特征向量。PCA 告诉我们：数据在 \( \mathbf{e}_1 \) 方向上的方差最大，在 \( \mathbf{e}_2 \) 方向上其次，以此类推。前 \( K \) 个特征向量张成信号子空间，后 \( M-K \) 个特征向量张成噪声子空间。

PCA 本身不直接给出角度——它只能指出方差大的方向，但未将这些方向与物理角度 \( \theta \) 关联。MUSIC 的核心思路是：信号子空间恰好由导向矢量 \( \{\mathbf{a}(\theta_1), \cdots, \mathbf{a}(\theta_K)\} \) 张成。因此，真实角度 \( \theta_k \) 对应的导向矢量 \( \mathbf{a}(\theta_k) \) 必然在信号子空间中，与噪声子空间正交。

这一正交性即 MUSIC 的判据：扫描所有可能的 \( \theta \)，计算 \( \mathbf{a}(\theta) \) 到噪声子空间的距离，距离为零时对应的 \( \theta \) 即为信号来向。

### 3.2 噪声与信号的统计假设

为了从 $\mathbf{x}(t)$ 中提取角度信息，对噪声和信号做以下标准假设：

**噪声假设：**
1. **零均值**：$\mathbb{E}[\mathbf{n}(t)] = \mathbf{0}$。
2. **空间白化**：$\mathbb{E}[\mathbf{n}(t) \mathbf{n}^H(t)] = \sigma^2 \mathbf{I}_N$。各阵元噪声功率相同，互不相关。
3. **与信号独立**：$\mathbb{E}[\mathbf{s}(t) \mathbf{n}^H(t)] = \mathbf{0}_{K \times N}$。若噪声与信号相关，说明存在阵元互耦等非理想因素。

**信号假设：**
1. **零均值**：$\mathbb{E}[\mathbf{s}(t)] = \mathbf{0}$。
2. **满秩协方差**：$\mathbf{R}_{SS} \triangleq \mathbb{E}[\mathbf{s}(t) \mathbf{s}^H(t)]$ 正定。若各信号互不相关，$\mathbf{R}_{SS} = \operatorname{diag}(P_1, P_2, \dots, P_K)$。

### 3.3 协方差矩阵：从欠定到可解

当前问题的结构：

> **已知**：阵列快拍 $\mathbf{x}(t)$（$N$ 维向量）。
> **未知**：方向矩阵 $\mathbf{A}(\theta)$（含 $K$ 个角度）、信号波形 $\mathbf{s}(t)$、噪声 $\mathbf{n}(t)$。
> **目标**：估计 $\theta_1, \dots, \theta_K$。

单快拍层面，已知量远少于未知量，无法求解。

但采集**多个快拍**（$t = 1, 2, \dots, L$）后，可以利用二阶统计量——协方差矩阵 $\mathbf{R}_{XX}$：

$$
\boxed{\mathbf{R}_{XX} \triangleq \mathbb{E}[\mathbf{x}(t) \mathbf{x}^H(t)] = \mathbf{A}(\theta) \, \mathbf{R}_{SS} \, \mathbf{A}^H(\theta) + \sigma^2 \mathbf{I}_N}\tag{1.29} $$

推导（代入 → 展开四项 → 信号与噪声独立消去交叉项 → 合并）：

$$
\begin{aligned}
\mathbf{R}_{XX} &= \mathbb{E}\big[(\mathbf{A}\mathbf{s} + \mathbf{n})(\mathbf{A}\mathbf{s} + \mathbf{n})^H\big] \\
&= \mathbb{E}[\mathbf{A}\mathbf{s}\mathbf{s}^H\mathbf{A}^H] + \underbrace{\mathbb{E}[\mathbf{A}\mathbf{s}\mathbf{n}^H]}_{= \mathbf{0}} + \underbrace{\mathbb{E}[\mathbf{n}\mathbf{s}^H\mathbf{A}^H]}_{= \mathbf{0}} + \mathbb{E}[\mathbf{n}\mathbf{n}^H] \\
&= \mathbf{A} \, \underbrace{\mathbb{E}[\mathbf{s}\mathbf{s}^H]}_{\mathbf{R}_{SS}} \, \mathbf{A}^H + \underbrace{\mathbb{E}[\mathbf{n}\mathbf{n}^H]}_{\sigma^2 \mathbf{I}_N}
\end{aligned}
\tag{1.30} $$

> $\mathbf{A} \mathbf{R}_{SS} \mathbf{A}^H$ 秩为 $K$（$\mathbf{A}$ 列满秩，$\mathbf{R}_{SS}$ 正定）。$\sigma^2 \mathbf{I}_N$ 满秩。$\mathbf{R}_{XX}$ 的特征分解将 $N$ 维空间分为 $K$ 维**信号子空间**和 $(N-K)$ 维**噪声子空间**，两者**相互正交**。

问题转化为：

> **如何从 $\mathbf{R}_{XX}$ 的特征结构中反推出 $\mathbf{A}(\theta)$ 中的角度 $\theta_1, \dots, \theta_K$。**

MUSIC 的做法是：利用信号子空间与噪声子空间的正交性——当且仅当 $\theta$ 等于真实信源方向时，$\mathbf{a}(\theta)$ 与噪声子空间正交。在角度域构造伪谱：

$$
P_{\text{MUSIC}}(\theta) = \frac{1}{\|\mathbf{a}^H(\theta) \mathbf{U}_N\|^2}\tag{1.31} $$

$\mathbf{U}_N$ 是噪声子空间的基。当扫描到真实方向时，分母趋于零，出现谱峰。MUSIC 突破了常规波束形成的瑞利限，实现了超分辨测向。下面展开完整推导。

### 3.4 MUSIC 算法的推导与实现

现在我们从阵列模型和特征分解两个角度，将协方差矩阵的结构完全展开，逐步推导出 MUSIC 算法的核心表达式。

阵列接收数据的协方差矩阵为：

$$
\mathbf{R}_X = \mathbb{E}\left[\mathbf{x}(t)\mathbf{x}^H(t)\right]
\tag{1.32} $$

将数据模型 \( \mathbf{x}(t) = \mathbf{A}(\boldsymbol{\theta})\mathbf{s}(t) + \mathbf{n}(t) \) 代入：

$$
\begin{aligned}
\mathbf{R}_X &= \mathbb{E}\left[(\mathbf{A}\mathbf{s} + \mathbf{n})(\mathbf{A}\mathbf{s} + \mathbf{n})^H\right] \\
&= \mathbb{E}\left[\mathbf{A}\mathbf{s}\mathbf{s}^H\mathbf{A}^H\right] + \mathbb{E}\left[\mathbf{n}\mathbf{n}^H\right] + \mathbb{E}\left[\mathbf{A}\mathbf{s}\mathbf{n}^H\right] + \mathbb{E}\left[\mathbf{n}\mathbf{s}^H\mathbf{A}^H\right]
\end{aligned}
\tag{1.33} $$

利用噪声与信号独立的假设（\( \mathbb{E}[\mathbf{s}\mathbf{n}^H] = \mathbf{0} \)，\( \mathbb{E}[\mathbf{n}\mathbf{s}^H] = \mathbf{0} \)），以及白噪声假设 \( \mathbb{E}[\mathbf{n}\mathbf{n}^H] = \sigma^2 \mathbf{I}_M \)，得到：

$$
\mathbf{R}_X = \mathbf{A} \mathbf{R}_S \mathbf{A}^H + \sigma^2 \mathbf{I}_M
\tag{1.34} $$

其中 \( \mathbf{R}_S = \mathbb{E}[\mathbf{s}(t)\mathbf{s}^H(t)] \) 是信号的协方差矩阵。

---

#### 3.4.1 特征分解：从数据域到子空间域

对 \( \mathbf{R}_X \) 进行特征分解，设其特征值为 \( \lambda_1 \geq \lambda_2 \geq \cdots \geq \lambda_M \)，对应的特征向量为 \( \mathbf{u}_1, \mathbf{u}_2, \cdots, \mathbf{u}_M \)，则：

$$
\mathbf{R}_X = \sum_{i=1}^{M} \lambda_i \mathbf{u}_i \mathbf{u}_i^H = \mathbf{U} \boldsymbol{\Lambda} \mathbf{U}^H
\tag{1.35} $$

其中：

$$
\mathbf{U} = [\mathbf{u}_1, \mathbf{u}_2, \cdots, \mathbf{u}_M] \in \mathbb{C}^{M \times M}, \quad \mathbf{U}\mathbf{U}^H = \mathbf{U}^H\mathbf{U} = \mathbf{I}_M
\tag{1.36} $$

$$
\boldsymbol{\Lambda} = \text{diag}(\lambda_1, \lambda_2, \cdots, \lambda_M), \quad \lambda_1 \geq \lambda_2 \geq \cdots \geq \lambda_M \geq 0
\tag{1.37} $$

由于协方差矩阵是 Hermitian 矩阵（\( \mathbf{R}_X^H = \mathbf{R}_X \)），其特征值全部为实数且非负，特征向量构成一组标准正交基。

---

#### 3.4.2 阵列模型与特征分解的统一视角

现在将阵列模型的表达式与特征分解的表达式并置：

$$
\mathbf{R}_X = \mathbf{U} \boldsymbol{\Lambda} \mathbf{U}^H = \mathbf{A} \mathbf{R}_S \mathbf{A}^H + \sigma^2 \mathbf{I}_M
\tag{1.38} $$

将 \( \sigma^2 \mathbf{I}_M \) 移到等号左边：

$$
\mathbf{U} \boldsymbol{\Lambda} \mathbf{U}^H - \sigma^2 \mathbf{I}_M = \mathbf{A} \mathbf{R}_S \mathbf{A}^H
\tag{1.39} $$

由于 \( \mathbf{U}\mathbf{U}^H = \mathbf{I}_M \)，可以将 \( \sigma^2 \mathbf{I}_M \) 写成 \( \sigma^2 \mathbf{U}\mathbf{U}^H \)：

$$
\mathbf{U} (\boldsymbol{\Lambda} - \sigma^2 \mathbf{I}_M) \mathbf{U}^H = \mathbf{A} \mathbf{R}_S \mathbf{A}^H
\tag{1.40} $$

记 \( \boldsymbol{\Lambda}_S = \boldsymbol{\Lambda} - \sigma^2 \mathbf{I}_M \)，则：

$$
\mathbf{U} \boldsymbol{\Lambda}_S \mathbf{U}^H = \mathbf{A} \mathbf{R}_S \mathbf{A}^H
\tag{1.41} $$

其中：

$$
\boldsymbol{\Lambda}_S = \text{diag}(\lambda_1 - \sigma^2, \lambda_2 - \sigma^2, \cdots, \lambda_M - \sigma^2)
\tag{1.42} $$

---

#### 3.4.3 信号子空间与噪声子空间的分离

现在考察特征值 \( \lambda_i \) 与噪声方差 \( \sigma^2 \) 的关系。

对于噪声特征向量 \( \mathbf{u}_i \)（\( i > K \)），其特征值等于噪声功率：

$$
\lambda_i = \sigma^2, \quad i = K+1, K+2, \cdots, M
\tag{1.43} $$

这意味着：

$$
\lambda_i - \sigma^2 = 0, \quad i = K+1, K+2, \cdots, M
\tag{1.44} $$

对于信号特征向量 \( \mathbf{u}_i \)（\( i \leq K \)），其特征值大于噪声功率：

$$
\lambda_i > \sigma^2, \quad i = 1, 2, \cdots, K
\tag{1.45} $$

因此，\( \boldsymbol{\Lambda}_S \) 的结构为：

$$
\boldsymbol{\Lambda}_S = 
\begin{pmatrix}
\lambda_1 - \sigma^2 & & & & \\
& \lambda_2 - \sigma^2 & & & \\
& & \ddots & & \\
& & & \lambda_K - \sigma^2 & \\
& & & & 0 \\
& & & & & \ddots \\
& & & & & & 0
\end{pmatrix}
=
\begin{pmatrix}
\boldsymbol{\Lambda}_s & \mathbf{0} \\
\mathbf{0} & \mathbf{0}
\end{pmatrix}
\tag{1.46} $$

其中 \( \boldsymbol{\Lambda}_s = \text{diag}(\lambda_1 - \sigma^2, \lambda_2 - \sigma^2, \cdots, \lambda_K - \sigma^2) \) 是 \( K \times K \) 的正定对角矩阵。

相应地，将特征向量矩阵 \( \mathbf{U} \) 分块为信号子空间和噪声子空间：

$$
\mathbf{U} = [\mathbf{U}_s \quad \mathbf{U}_N]
\tag{1.47} $$

其中：

$$
\mathbf{U}_s = [\mathbf{u}_1, \mathbf{u}_2, \cdots, \mathbf{u}_K] \in \mathbb{C}^{M \times K}, \quad \text{张成信号子空间}
\tag{1.48} $$

$$
\mathbf{U}_N = [\mathbf{u}_{K+1}, \mathbf{u}_{K+2}, \cdots, \mathbf{u}_M] \in \mathbb{C}^{M \times (M-K)}, \quad \text{张成噪声子空间}
\tag{1.49} $$

于是：

$$
\mathbf{U} \boldsymbol{\Lambda}_S \mathbf{U}^H = 
[\mathbf{U}_s \quad \mathbf{U}_N]
\begin{pmatrix}
\boldsymbol{\Lambda}_s & \mathbf{0} \\
\mathbf{0} & \mathbf{0}
\end{pmatrix}
\begin{pmatrix}
\mathbf{U}_s^H \\
\mathbf{U}_N^H
\end{pmatrix}
= \mathbf{U}_s \boldsymbol{\Lambda}_s \mathbf{U}_s^H
\tag{1.50} $$

因此，结合 \( \mathbf{U} \boldsymbol{\Lambda}_S \mathbf{U}^H = \mathbf{A} \mathbf{R}_S \mathbf{A}^H \)，得到：

$$
\mathbf{U}_s \boldsymbol{\Lambda}_s \mathbf{U}_s^H = \mathbf{A} \mathbf{R}_S \mathbf{A}^H
\tag{1.51} $$

---

#### 3.4.4 信号子空间与噪声子空间的正交性

##### 预备知识：span、ker 与正交补空间

在推导信号与噪声子空间的正交关系之前，先统一几个在后续推导中反复出现的线性代数基本概念。

###### 1. 张成空间（span）

**定义：** 给定一组向量 \( \{\mathbf{v}_1, \mathbf{v}_2, \cdots, \mathbf{v}_r\} \subset \mathbb{C}^M \)，它们的**张成空间**（span）是指这组向量的所有线性组合构成的集合：

$$
\text{span}\{\mathbf{v}_1, \cdots, \mathbf{v}_r\} = \left\{ \sum_{i=1}^r c_i \mathbf{v}_i \;\Big|\; c_i \in \mathbb{C} \right\}
$$

当向量组拼成矩阵 \( \mathbf{V} = [\mathbf{v}_1, \cdots, \mathbf{v}_r] \in \mathbb{C}^{M \times r} \) 时，span 等价于该矩阵的**列空间**（column space / range）：

$$
\text{span}\{\mathbf{v}_1, \cdots, \mathbf{v}_r\} = \text{range}(\mathbf{V}) = \{\mathbf{V}\mathbf{c} : \mathbf{c} \in \mathbb{C}^r\}
$$

直观理解：span 就是这组向量"够得着"的所有方向——把它们按任意系数组合，得到的全体向量就构成一个子空间。

###### 2. 核空间（kernel / nullspace）

**定义：** 设 \( \mathbf{B} \in \mathbb{C}^{M \times M} \) 是一个方阵，它的**核空间**（kernel，也称零空间 nullspace）定义为所有被 \( \mathbf{B} \) 映射到零向量的输入构成的集合：

$$
\ker(\mathbf{B}) = \{\mathbf{x} \in \mathbb{C}^M : \mathbf{B}\mathbf{x} = \mathbf{0}\}
$$

**秩-零化度定理（Rank–Nullity Theorem）：**

$$
\dim(\text{range}(\mathbf{B})) + \dim(\ker(\mathbf{B})) = M
$$

直观理解：一个矩阵作用于全空间，有些方向会被"保留"（range），有些方向会被"消灭"（kernel）。两者的维度之和恰好等于全空间的维度。

**kernel 与 range 的对偶关系（一般形式）：** 对任意矩阵 \( \mathbf{B} \in \mathbb{C}^{M \times N} \)，有：

$$
\text{range}(\mathbf{B})^\perp = \ker(\mathbf{B}^H)
\tag{*}
$$

**证明：** \( \mathbf{y} \in \ker(\mathbf{B}^H) \iff \mathbf{B}^H \mathbf{y} = \mathbf{0} \iff \forall \mathbf{x}, \langle \mathbf{B}\mathbf{x}, \mathbf{y} \rangle = \mathbf{x}^H \mathbf{B}^H \mathbf{y} = 0 \iff \mathbf{y} \perp \text{range}(\mathbf{B}) \iff \mathbf{y} \in \text{range}(\mathbf{B})^\perp \)。

两端取正交补即得对偶形式：\( \text{range}(\mathbf{B}) = (\ker(\mathbf{B}^H))^\perp \)。后面"第二步"中将用此关系，从 kernel 相等推出 range 相等。

###### 3. 正交补空间（orthogonal complement）

**定义：** 设 \( V \) 是 \( \mathbb{C}^M \) 的线性子空间，则 \( V \) 的**正交补空间**（记为 \( V^\perp \)）定义为：

$$
V^\perp = \{\mathbf{x} \in \mathbb{C}^M : \langle \mathbf{x}, \mathbf{v} \rangle = 0,\ \forall \mathbf{v} \in V\}
\tag{1.51a} $$

即 \( V^\perp \) 是 \( \mathbb{C}^M \) 中所有与 \( V \) 内每一个向量都正交的向量组成的集合。

**正交补的基本性质：**

1. **直和分解：** \( V \oplus V^\perp = \mathbb{C}^M \)，即全空间可以唯一地分解为 \( V \) 与 \( V^\perp \) 的直和。
2. **维数关系：** \( \dim(V) + \dim(V^\perp) = M \)。
3. **互为补：** \( (V^\perp)^\perp = V \)，正交补是相互的——若 \( W = V^\perp \)，则 \( V = W^\perp \)。
4. **正交投影：** \( \mathbb{C}^M \) 中的任一向量可以唯一地写为 \( \mathbf{x} = \mathbf{v} + \mathbf{w} \)，其中 \( \mathbf{v} \in V \)，\( \mathbf{w} \in V^\perp \)，且 \( \langle \mathbf{v}, \mathbf{w} \rangle = 0 \)。

通俗地说，正交补就是把整个空间"切成"互相垂直的两半，每个子空间中的向量都和另一个子空间中的所有向量正交。

##### 为什么正交补恰好等于核空间？——矩阵的 range–kernel 分解

将上述抽象概念落实到矩阵运算上，有一个极其重要的定理：

> **定理（Hermitian 矩阵的 range–kernel 正交分解）**
> 
> 设 \( \mathbf{B} \in \mathbb{C}^{M \times M} \) 是 Hermitian 矩阵（\( \mathbf{B}^H = \mathbf{B} \)），则：
> $$ \mathbb{C}^M = \text{range}(\mathbf{B}) \oplus \ker(\mathbf{B}) $$
> 且该直和是**正交直和**：\( \text{range}(\mathbf{B}) \perp \ker(\mathbf{B}) \)。

**证明思路：**

1. **正交性：** 任取 \( \mathbf{y} = \mathbf{B}\mathbf{x} \in \text{range}(\mathbf{B}) \) 和 \( \mathbf{z} \in \ker(\mathbf{B}) \)，则 \( \langle \mathbf{y}, \mathbf{z} \rangle = \mathbf{z}^H \mathbf{B} \mathbf{x} = (\mathbf{B}^H \mathbf{z})^H \mathbf{x} = (\mathbf{B} \mathbf{z})^H \mathbf{x} = \mathbf{0}^H \mathbf{x} = 0 \)。关键一步利用了 \( \mathbf{B}^H = \mathbf{B} \) 和 \( \mathbf{B}\mathbf{z} = \mathbf{0} \)。

2. **直和：** 任取 \( \mathbf{w} \in \mathbb{C}^M \)，将其分解为 \( \mathbf{w} = \underbrace{\mathbf{B} \mathbf{B}^+ \mathbf{w}}_{\in\ \text{range}(\mathbf{B})} + \underbrace{(\mathbf{I} - \mathbf{B} \mathbf{B}^+) \mathbf{w}}_{\in\ \ker(\mathbf{B})} \)，其中 \( \mathbf{B}^+ \) 是 \( \mathbf{B} \) 的 Moore–Penrose 伪逆。这证明了 \( \mathbb{C}^M = \text{range}(\mathbf{B}) + \ker(\mathbf{B}) \)。

3. 结合维数公式 \( \dim(\text{range}(\mathbf{B})) + \dim(\ker(\mathbf{B})) = M \)（秩-零化度定理），知和为直和，且 \( \ker(\mathbf{B}) = (\text{range}(\mathbf{B}))^\perp \)，即**核空间恰好就是值域的正交补**。

**在 MUSIC 中的应用：** 这个定理直接告诉我们，一旦确定了某个投影矩阵的 range，它的 kernel 自然就是 range 的正交补。这对理解信号/噪声子空间的关系至关重要——参见下文的"第二步"。

---

##### 第一步：信号子空间与噪声子空间的正交性

回顾特征分解的结果（式 1.35）：协方差矩阵 \( \mathbf{R}_X \) 的 \( M \) 个特征向量构成酉矩阵

$$
\mathbf{U} = [\mathbf{u}_1, \mathbf{u}_2, \cdots, \mathbf{u}_M] \in \mathbb{C}^{M \times M}
$$

酉矩阵的定义是 \( \mathbf{U}^H \mathbf{U} = \mathbf{I}_M \)，展开即：

$$
\mathbf{u}_i^H \mathbf{u}_j = \delta_{ij} =
\begin{cases}
1, & i = j \\
0, & i \neq j
\end{cases}
\tag{1.51b} $$

这意味着：\( \mathbf{U} \) 的所有列向量构成 \( \mathbb{C}^M \) 的一组**标准正交基**。

在此基础上，按特征值大小将 \( \mathbf{U} \) 分为两块（式 1.47–1.49）：

$$
\mathbf{U} = [\mathbf{U}_s \mid \mathbf{U}_N]
$$

- \( \mathbf{U}_s = [\mathbf{u}_1, \cdots, \mathbf{u}_K] \) 对应大特征值，列向量张成**信号子空间** \( \text{span}\{\mathbf{U}_s\} \)。
- \( \mathbf{U}_N = [\mathbf{u}_{K+1}, \cdots, \mathbf{u}_M] \) 对应小特征值，列向量张成**噪声子空间** \( \text{span}\{\mathbf{U}_N\} \)。

由于 \( \mathbf{U}_s \) 的每一列与 \( \mathbf{U}_N \) 的每一列来自同一个标准正交基，它们彼此正交，因此：

$$
\mathbf{U}_s^H \mathbf{U}_N = \mathbf{0}_{K \times (M-K)}, \quad \mathbf{U}_N^H \mathbf{U}_s = \mathbf{0}_{(M-K) \times K}
\tag{1.52} $$

这意味着 \( \text{span}\{\mathbf{U}_s\} \) 和 \( \text{span}\{\mathbf{U}_N\} \) 是互相正交的子空间。

又因为：

- \( \dim(\text{span}\{\mathbf{U}_s\}) = K \)
- \( \dim(\text{span}\{\mathbf{U}_N\}) = M - K \)
- \( K + (M - K) = M \)

对比正交补的定义可知，**噪声子空间正是信号子空间的正交补**，反之亦然：

$$
\text{span}\{\mathbf{U}_N\} = \bigl(\text{span}\{\mathbf{U}_s\}\bigr)^\perp, \qquad
\text{span}\{\mathbf{U}_s\} = \bigl(\text{span}\{\mathbf{U}_N\}\bigr)^\perp
\tag{1.53} $$

于是整个空间 \( \mathbb{C}^M \) 有一个正交直和分解：

$$
\boxed{\mathbb{C}^M = \text{span}\{\mathbf{U}_s\} \oplus \text{span}\{\mathbf{U}_N\}}
\tag{1.54}
$$

###### 与 range–kernel 分解的对应关系

这个直和分解并非巧合，它是前述 range–kernel 定理的直接体现。考虑**噪声子空间投影矩阵**：

$$
\mathbf{P}_N = \mathbf{U}_N \mathbf{U}_N^H
$$

\( \mathbf{P}_N \) 是 Hermitian 矩阵（\( \mathbf{P}_N^H = \mathbf{P}_N \)）且幂等（\( \mathbf{P}_N^2 = \mathbf{P}_N \)），它是一个正交投影算子。对其应用 range–kernel 分解定理：

- \( \text{range}(\mathbf{P}_N) = \text{span}\{\mathbf{U}_N\} \) —— 噪声子空间（\(M-K\) 维）
- \( \ker(\mathbf{P}_N) = \text{span}\{\mathbf{U}_s\} \) —— 信号子空间（\(K\) 维）

定理断言 \( \text{range}(\mathbf{P}_N) \perp \ker(\mathbf{P}_N) \) 且二者直和为 \( \mathbb{C}^M \)，即：

$$
\mathbb{C}^M = \text{range}(\mathbf{P}_N) \oplus \ker(\mathbf{P}_N) = \text{span}\{\mathbf{U}_N\} \oplus \text{span}\{\mathbf{U}_s\}
\tag{1.54a}
$$

等价地，对**信号子空间投影矩阵** \( \mathbf{P}_s = \mathbf{U}_s \mathbf{U}_s^H \) 也有：

$$
\text{range}(\mathbf{P}_s) = \text{span}\{\mathbf{U}_s\}, \quad \ker(\mathbf{P}_s) = \text{span}\{\mathbf{U}_N\}
$$

**关键理解：** 信号子空间和噪声子空间并非两个"独立找来"的子空间，而是**同一个 Hermitian 投影算子的 range 和 kernel**——一个是矩阵"能到达"的方向，另一个是矩阵"归零"的方向。二者由秩-零化度定理天然配对，维度互补且正交。

---

##### 第二步：信号子空间与方向矩阵张成同一子空间

上一节式 (1.51) 给出了关键等式：

$$
\mathbf{U}_s \boldsymbol{\Lambda}_s \mathbf{U}_s^H = \mathbf{A} \mathbf{R}_S \mathbf{A}^H
\tag{1.51}
$$

目标是证明：\( \text{span}\{\mathbf{U}_s\} = \text{span}\{\mathbf{A}\} \)。

###### 证明策略：走 kernel 路径

直接比较两边的 range 容易含混。更清晰的路径是：**先证明两边的核空间相等，再利用 range 与 kernel 的对偶关系推出 span 相等。**

需要用到的对偶关系（预备知识中已建立，此处重申）：

$$
\text{range}(\mathbf{B})^\perp = \ker(\mathbf{B}^H), \qquad
\text{range}(\mathbf{B}) = \bigl(\ker(\mathbf{B}^H)\bigr)^\perp
\tag{1.55a}
$$

因此，只要证明 \( \ker(\mathbf{U}_s^H) = \ker(\mathbf{A}^H) \)，两侧取正交补即得 \( \text{span}\{\mathbf{U}_s\} = \text{span}\{\mathbf{A}\} \)。

###### 已知条件回顾

- \( \mathbf{U}_s \in \mathbb{C}^{M \times K} \)，列满秩，\( \mathbf{U}_s^H \mathbf{U}_s = \mathbf{I}_K \)；
- \( \boldsymbol{\Lambda}_s = \text{diag}(\lambda_1 - \sigma^2, \cdots, \lambda_K - \sigma^2) \succ \mathbf{0} \)（正定对角矩阵，可逆）；
- \( \mathbf{A} = [\mathbf{a}(\theta_1), \cdots, \mathbf{a}(\theta_K)] \in \mathbb{C}^{M \times K} \)，列满秩（不同方向的导向矢量线性无关）；
- \( \mathbf{R}_S \in \mathbb{C}^{K \times K} \) 是信号协方差矩阵，正定（\( \mathbf{R}_S \succ \mathbf{0} \)），可逆。

###### 方向一：\( \ker(\mathbf{U}_s^H) \subseteq \ker(\mathbf{A}^H) \)

设 \( \mathbf{y} \in \ker(\mathbf{U}_s^H) \)，即 \( \mathbf{U}_s^H \mathbf{y} = \mathbf{0} \)。

将式 (1.51) 两边右乘 \( \mathbf{y} \)：

$$
\mathbf{U}_s \boldsymbol{\Lambda}_s (\mathbf{U}_s^H \mathbf{y}) = \mathbf{A} \mathbf{R}_S (\mathbf{A}^H \mathbf{y})
$$

左边 \( \mathbf{U}_s \boldsymbol{\Lambda}_s \cdot \mathbf{0} = \mathbf{0} \)，故：

$$
\mathbf{A} \mathbf{R}_S (\mathbf{A}^H \mathbf{y}) = \mathbf{0}
\tag{1.55b}
$$

左乘 \( (\mathbf{A}^H \mathbf{y})^H \)：

$$
(\mathbf{A}^H \mathbf{y})^H \mathbf{A} \mathbf{R}_S (\mathbf{A}^H \mathbf{y}) = \mathbf{0}
$$

令 \( \mathbf{z} = \mathbf{A}^H \mathbf{y} \in \mathbb{C}^K \)，则上式即 \( \mathbf{z}^H \mathbf{R}_S \mathbf{z} = 0 \)。

由于 \( \mathbf{R}_S \succ \mathbf{0} \)（正定），二次型 \( \mathbf{z}^H \mathbf{R}_S \mathbf{z} = 0 \) 当且仅当 \( \mathbf{z} = \mathbf{0} \)。因此：

$$
\mathbf{A}^H \mathbf{y} = \mathbf{0} \quad \Longrightarrow \quad \mathbf{y} \in \ker(\mathbf{A}^H)
$$

即 \( \ker(\mathbf{U}_s^H) \subseteq \ker(\mathbf{A}^H) \)。

###### 方向二：\( \ker(\mathbf{A}^H) \subseteq \ker(\mathbf{U}_s^H) \)

设 \( \mathbf{y} \in \ker(\mathbf{A}^H) \)，即 \( \mathbf{A}^H \mathbf{y} = \mathbf{0} \)。

将式 (1.51) 两边右乘 \( \mathbf{y} \)：

$$
\mathbf{U}_s \boldsymbol{\Lambda}_s (\mathbf{U}_s^H \mathbf{y}) = \underbrace{\mathbf{A} \mathbf{R}_S (\mathbf{A}^H \mathbf{y})}_{\mathbf{A} \mathbf{R}_S \cdot \mathbf{0} = \mathbf{0}} = \mathbf{0}
$$

即：

$$
\mathbf{U}_s \boldsymbol{\Lambda}_s (\mathbf{U}_s^H \mathbf{y}) = \mathbf{0}
\tag{1.55c}
$$

左乘 \( \mathbf{U}_s^H \)：

$$
\mathbf{U}_s^H \mathbf{U}_s \boldsymbol{\Lambda}_s (\mathbf{U}_s^H \mathbf{y}) = \mathbf{0}
$$

利用 \( \mathbf{U}_s^H \mathbf{U}_s = \mathbf{I}_K \)：

$$
\boldsymbol{\Lambda}_s (\mathbf{U}_s^H \mathbf{y}) = \mathbf{0}
$$

由于 \( \boldsymbol{\Lambda}_s \succ \mathbf{0} \)（可逆），必有：

$$
\mathbf{U}_s^H \mathbf{y} = \mathbf{0} \quad \Longrightarrow \quad \mathbf{y} \in \ker(\mathbf{U}_s^H)
$$

即 \( \ker(\mathbf{A}^H) \subseteq \ker(\mathbf{U}_s^H) \)。

###### 结论

两个方向结合，得到：

$$
\boxed{\ker(\mathbf{U}_s^H) = \ker(\mathbf{A}^H)}
\tag{1.55d}
$$

利用 kernel–range 对偶关系 (1.55a)，两侧取正交补：

$$
\bigl(\ker(\mathbf{U}_s^H)\bigr)^\perp = \bigl(\ker(\mathbf{A}^H)\bigr)^\perp
\quad \Longrightarrow \quad
\text{range}(\mathbf{U}_s) = \text{range}(\mathbf{A})
$$

即：

$$
\boxed{\text{span}\{\mathbf{U}_s\} = \text{span}\{\mathbf{A}\}}
\tag{1.55}
$$

**通俗理解：** 式 (1.51) 这个矩阵等式，从 kernel 角度看，意味着"被 \( \mathbf{U}_s^H \) 消灭的向量"与"被 \( \mathbf{A}^H \) 消灭的向量"完全一致。由于 range 与 kernel 互为对方的正交补，range 自然也就相同——**信号特征向量张成的子空间，就是方向矩阵列向量张成的子空间。** 这正是它被命名为"信号子空间"的根本原因。

由于 \( \mathbf{A} \) 的每一列就是某个信号源方向的导向矢量 \( \mathbf{a}(\theta_k) \)，(1.55) 的直接推论是：

$$
\mathbf{a}(\theta_k) \in \text{span}\{\mathbf{U}_s\}, \quad k = 1, 2, \cdots, K
\tag{1.56}
$$

---

#### 第三步：导向矢量与噪声子空间的正交性——MUSIC 的核心

现在将前两步的结论组合在一起：

1. **由 (1.53)/(1.54)：** 噪声子空间 \( \text{span}\{\mathbf{U}_N\} \) 是信号子空间 \( \text{span}\{\mathbf{U}_s\} \) 的正交补。
2. **由 (1.56)：** 每一个真实的导向矢量 \( \mathbf{a}(\theta_k) \) 都落在信号子空间 \( \text{span}\{\mathbf{U}_s\} \) 中。

根据正交补的定义——信号子空间中任意向量都与噪声子空间中任意向量正交，立即得出：

$$
\boxed{\mathbf{U}_N^H \mathbf{a}(\theta_k) = \mathbf{0}_{(M-K) \times 1}, \quad k = 1, 2, \cdots, K}
\tag{1.57}
$$

逐行展开就是：

$$
\mathbf{u}_{K+1}^H \mathbf{a}(\theta_k) = 0,\quad \mathbf{u}_{K+2}^H \mathbf{a}(\theta_k) = 0,\quad \cdots,\quad \mathbf{u}_M^H \mathbf{a}(\theta_k) = 0
\tag{1.58}
$$

这 \( M-K \) 个标量等式共同表明：每一个真实信号源对应的导向矢量 \( \mathbf{a}(\theta_k) \) **与整个噪声子空间的每一维都正交**。

---

#### 小结：逻辑链回顾

整个推导的逻辑链条如下。注意第二条线的结构：**先证 kernel 相等，再利用 kernel–range 对偶关系推出 span 相等**——这恰好对应预备知识中建立的 $\text{range}(\mathbf{B}) = (\ker(\mathbf{B}^H))^\perp$。

$$
\begin{CD}
\mathbf{U} \text{ 是酉矩阵} @>>> \mathbf{U}_s \perp \mathbf{U}_N \text{ 且 } \dim(\mathbf{U}_s) + \dim(\mathbf{U}_N) = M \\
@VVV \\
\text{span}\{\mathbf{U}_N\} = \bigl(\text{span}\{\mathbf{U}_s\}\bigr)^\perp \quad \text{(噪声子空间=信号子空间的正交补)}
\end{CD}
$$

$$
\begin{CD}
\mathbf{U}_s \boldsymbol{\Lambda}_s \mathbf{U}_s^H = \mathbf{A} \mathbf{R}_S \mathbf{A}^H @>>> \ker(\mathbf{U}_s^H) = \ker(\mathbf{A}^H) \\
@VVV \\
\text{span}\{\mathbf{U}_s\} = \text{span}\{\mathbf{A}\} \\
@VVV \\
\mathbf{a}(\theta_k) \in \text{span}\{\mathbf{U}_s\}
\end{CD}
$$

两条线交汇：

$$
\begin{CD}
\mathbf{a}(\theta_k) \in \text{span}\{\mathbf{U}_s\} \quad , \quad \text{span}\{\mathbf{U}_N\} \perp \text{span}\{\mathbf{U}_s\} \\
@VVV \\
\mathbf{U}_N^H \mathbf{a}(\theta_k) = \mathbf{0}
\end{CD}
$$

这一结论构成了 MUSIC 算法谱函数设计的理论基础。下一节将以此构建用于估计信号到达角的谱峰搜索函数。

---

#### 3.4.5 MUSIC 谱函数的导出

对于任意一个试探方向 \( \theta \)，其对应的导向矢量为 \( \mathbf{a}(\theta) \)。如果 \( \theta \) 恰好等于某个真实信号源的方向 \( \theta_k \)，则 \( \mathbf{a}(\theta) \) 与噪声子空间正交，即：

$$
\mathbf{U}_N^H \mathbf{a}(\theta_k) = \mathbf{0} \quad \Longrightarrow \quad \| \mathbf{U}_N^H \mathbf{a}(\theta_k) \|^2 = 0
\tag{1.59} $$

如果 \( \theta \) 不是真实信号源的方向，则 \( \mathbf{a}(\theta) \) 在噪声子空间上存在非零投影，即：

$$
\| \mathbf{U}_N^H \mathbf{a}(\theta) \|^2 > 0
\tag{1.60} $$

因此，定义 MUSIC 空间谱为：

$$
\boxed{P_{\text{MUSIC}}(\theta) = \frac{1}{\| \mathbf{U}_N^H \mathbf{a}(\theta) \|^2}}
\tag{1.61} $$

或者等价地：

$$
P_{\text{MUSIC}}(\theta) = \frac{1}{\mathbf{a}^H(\theta) \mathbf{U}_N \mathbf{U}_N^H \mathbf{a}(\theta)}
\tag{1.62} $$

当 \( \theta \) 等于真实信号方向时，分母为零，谱函数出现尖峰（理论上趋于无穷大）；当 \( \theta \) 偏离真实信号方向时，分母非零，谱函数取有限值。

在实际实现中，由于有限样本估计存在误差，分母不会严格为零，而是在真实信号方向处取得极小值，谱函数在该处出现尖锐峰值。**这些峰值的位置就是信号源到达角 \( \theta_k \) 的估计值**：

$$
\hat{\theta}_k = \arg\max_{\theta} P_{\text{MUSIC}}(\theta), \quad k = 1, 2, \cdots, K
\tag{1.63} 
$$

**MUSIC的谱反应的是方向矢量与噪声子空间的正交程度**
**因此叫伪峰、伪谱**


---

#### 3.4.6 算法流程总结

**输入：** 阵列快拍数据 \( \{\mathbf{x}(t)\}_{t=1}^{N_s} \)，信号源个数 \( K \)

**步骤 1：** 估计协方差矩阵

$$
\hat{\mathbf{R}}_X = \frac{1}{N_s} \sum_{t=1}^{N_s} \mathbf{x}(t) \mathbf{x}^H(t)
\tag{1.64} $$

**步骤 2：** 对 \( \hat{\mathbf{R}}_X \) 进行特征分解

$$
\hat{\mathbf{R}}_X = \sum_{i=1}^{M} \lambda_i \mathbf{u}_i \mathbf{u}_i^H, \quad \lambda_1 \geq \lambda_2 \geq \cdots \geq \lambda_M
\tag{1.65} $$

**步骤 3：** 确定噪声子空间

取后 \( M-K \) 个特征向量构成噪声子空间：

$$
\hat{\mathbf{U}}_N = [\mathbf{u}_{K+1}, \mathbf{u}_{K+2}, \cdots, \mathbf{u}_M]
\tag{1.66} $$

**步骤 4：** 扫描角度，计算 MUSIC 谱

对每个试探角度 \( \theta \in [-90^\circ, 90^\circ] \)，计算：

$$
P_{\text{MUSIC}}(\theta) = \frac{1}{\mathbf{a}^H(\theta) \hat{\mathbf{U}}_N \hat{\mathbf{U}}_N^H \mathbf{a}(\theta)}
\tag{1.67} $$

**步骤 5：** 峰值检测

找出 \( P_{\text{MUSIC}}(\theta) \) 的前 \( K \) 个峰值，其对应的角度即为信号源的到达角估计值 \( \hat{\theta}_1, \hat{\theta}_2, \cdots, \hat{\theta}_K \)。

**输出：** \( \hat{\theta}_1, \hat{\theta}_2, \cdots, \hat{\theta}_K \)


#### 3.4.7 MUSIC 与 Capon 的数学联系与性能对比

MUSIC 和 Capon（MVDR）是阵列处理中两种主要的高分辨谱估计方法。它们的谱函数形式类似——均为分母趋零时产生尖峰——这源于两者之间的数学联系。

##### Capon 谱的回顾

Capon 谱（空间版本的 MVDR 波束形成器）定义为：

$$
\boxed{P_{\text{Capon}}(\theta) = \frac{1}{\mathbf{a}^H(\theta) \mathbf{R}_X^{-1} \mathbf{a}(\theta)}}
\tag{1.68}
$$

它来源于一个约束优化问题：在保证期望方向增益为 1 的条件下，最小化阵列输出总功率。与 MUSIC 并置比较：

$$
P_{\text{MUSIC}}(\theta) = \frac{1}{\mathbf{a}^H(\theta) \, \mathbf{U}_N \mathbf{U}_N^H \, \mathbf{a}(\theta)},
\qquad
P_{\text{Capon}}(\theta) = \frac{1}{\mathbf{a}^H(\theta) \, \mathbf{R}_X^{-1} \, \mathbf{a}(\theta)}
$$

MUSIC 分母中用的是**噪声子空间投影矩阵** \( \mathbf{U}_N \mathbf{U}_N^H \)，Capon 用的是**协方差矩阵的逆** \( \mathbf{R}_X^{-1} \)。两者看似不同，但存在精确的代数关系。

##### 数学桥梁：\( \mathbf{R}_X^{-1} \) 的谱分解

将 \( \mathbf{R}_X \) 的特征分解代入逆矩阵：

$$
\mathbf{R}_X = \mathbf{U}_s \boldsymbol{\Lambda}_s \mathbf{U}_s^H + \sigma^2 \mathbf{U}_N \mathbf{U}_N^H
\tag{1.69}
$$

两边取逆（利用 \( \mathbf{U} \) 的酉性质）：

$$
\boxed{\mathbf{R}_X^{-1} = \mathbf{U}_s \boldsymbol{\Lambda}_s^{-1} \mathbf{U}_s^H + \frac{1}{\sigma^2} \mathbf{U}_N \mathbf{U}_N^H}
\tag{1.70}
$$

这个展开说明：**\( \mathbf{R}_X^{-1} \) 包含噪声子空间投影项 \( \mathbf{U}_N \mathbf{U}_N^H \)**，与信号子空间项 \( \mathbf{U}_s \boldsymbol{\Lambda}_s^{-1} \mathbf{U}_s^H \) 加权叠加。

##### 高 SNR 下的趋同

当信噪比很高时（\( \sigma^2 \ll \lambda_i,\ i=1,\dots,K \)），\( \boldsymbol{\Lambda}_s^{-1} \) 的每个对角元素 \( 1/(\lambda_i - \sigma^2) \) 都是有界的小量，而 \( 1/\sigma^2 \) 趋近于无穷大。因此：

$$
\mathbf{R}_X^{-1} \;\xrightarrow{\text{高 SNR}}\; \frac{1}{\sigma^2} \mathbf{U}_N \mathbf{U}_N^H
\tag{1.71}
$$

代入 Capon 谱：

$$
\boxed{P_{\text{Capon}}(\theta) \approx \sigma^2 \cdot P_{\text{MUSIC}}(\theta) \quad (\text{高 SNR})}
\tag{1.72}
$$

**高 SNR 极限下，Capon 谱和 MUSIC 谱仅差一个常数因子 \( \sigma^2 \)**，谱峰位置完全一致。

##### 几何解释

MUSIC 和 Capon 解决同一个问题的路径不同：

| | MUSIC | Capon |
|:---|:---|:---|
| **策略** | 显式分离信号/噪声子空间 | 直接对协方差矩阵求逆 |
| **消去信号的方式** | 特征分解后丢弃 \( \mathbf{U}_s \)，保留 \( \mathbf{U}_N \) | 逆矩阵自动赋予噪声方向更大权重 |
| **核心运算** | \( \mathbf{U}_N \mathbf{U}_N^H \)（正交投影） | \( \mathbf{R}_X^{-1} \)（白化+逆滤波） |
| **需要先验知识** | 信号源个数 \( K \) | 不需要 \( K \)，但需要 \( \mathbf{R}_X \) 满秩 |
| **理论分辨率** | 超分辨（突破瑞利限） | 超分辨，但略逊于 MUSIC |

两种方法的共同点：**信号方向的导向矢量 \( \mathbf{a}(\theta_k) \) 在信号子空间中，与噪声子空间正交。MUSIC 直接检验这种正交性，Capon 则通过 \( \mathbf{R}_X^{-1} \) 隐式实现——噪声方向被放大（除以 \( \sigma^2 \)），信号方向被抑制（除以大特征值）。**

##### 有限 SNR 下的分歧

在有限信噪比下，二者的差异开始显现。将 \( \mathbf{R}_X^{-1} \) 代入 Capon 谱分母：

$$
\mathbf{a}^H(\theta) \mathbf{R}_X^{-1} \mathbf{a}(\theta) =
\underbrace{\mathbf{a}^H(\theta) \mathbf{U}_s \boldsymbol{\Lambda}_s^{-1} \mathbf{U}_s^H \mathbf{a}(\theta)}_{\text{信号子空间贡献，}>0}
+ \frac{1}{\sigma^2}
\underbrace{\mathbf{a}^H(\theta) \mathbf{U}_N \mathbf{U}_N^H \mathbf{a}(\theta)}_{\text{MUSIC 的分母}}
$$

当 \( \theta = \theta_k \) 时，右边第二项为零（MUSIC 的判定条件），但**第一项并不为零**——真实导向矢量在信号子空间中仍有非零投影。这意味着 Capon 谱峰的高度受信号子空间残余项的制约，而 MUSIC 不存在这个"自污染"问题。因此：

- **MUSIC 的谱峰理论上可以趋于无穷**（分母严格为零），分辨率更高；
- **Capon 的谱峰高度有限**，当信号源角度接近时，谱峰合并得更早。

这解释了为什么在工程实践中，若已知信号源个数 \( K \)，MUSIC 的分辨能力通常优于 Capon；但 Capon 不需要已知 \( K \)，在信号源个数不确定的场景下更为稳健。

##### 统一视角：协方差矩阵的"逆"是关键

MUSIC 和 Capon 的关系可以总结为一条简洁的公式链：

$$
\mathbf{R}_X^{-1}
\;=\;
\underbrace{\mathbf{U}_s \boldsymbol{\Lambda}_s^{-1} \mathbf{U}_s^H}_{\text{信号方向权重小}}
\;+\;
\underbrace{\frac{1}{\sigma^2} \mathbf{U}_N \mathbf{U}_N^H}_{\text{噪声方向权重大}}
\;\xrightarrow{\sigma^2 \to 0}\;
\frac{1}{\sigma^2} \mathbf{U}_N \mathbf{U}_N^H
$$

> **Capon 是"隐式的 MUSIC"**——协方差矩阵求逆这个操作，本身就内在地重加权了信号子空间和噪声子空间，只不过它把二者混在了一起。MUSIC 的贡献在于把这种混合"拆开"，通过特征分解将信号和噪声彻底分离，从而获得了更纯粹的正交性检验和更高的分辨率。

### 3.8 工程实践中的退化现象与对策

前面几节推导了 MUSIC 和 Capon 的理想模型——假定协方差矩阵精确已知、噪声是理想白噪声、信号源互不相干。但在实际系统中，这些假设往往不成立，导致谱估计出现各种退化现象。本节梳理最常见的几类问题及其成因。

#### 3.8.1 伪峰问题：谱峰数目多于真实信号源数

理想情况下，MUSIC 谱在 \(K\) 个真实信号方向上各产生一个尖峰。但在实践中，经常出现**多于 \(K\) 个峰**的情况：

**成因 1：协方差矩阵估计误差**

实际中 \( \mathbf{R}_X \) 由有限快拍数 \(N\) 估计得到：

$$
\hat{\mathbf{R}}_X = \frac{1}{N} \sum_{t=1}^{N} \mathbf{x}(t) \mathbf{x}^H(t)
$$

当 \(N\) 不够大时，\( \hat{\mathbf{R}}_X \) 与真实 \( \mathbf{R}_X \) 存在偏差。体现在特征结构上：噪声特征值不再严格等于 \( \sigma^2 \)，而是散布在 \( \sigma^2 \) 附近。这导致 \( \mathbf{U}_N \mathbf{U}_N^H \mathbf{a}(\theta) \) 在非信号方向上也不严格为零——原本平坦的噪声基底出现起伏，某些角度上的起伏恰好形成"伪峰"。

**成因 2：信号源个数 \(K\) 估计错误**

MUSIC 需要预先知道 \(K\) 来划分子空间。若 \(K\) 估计偏大，会将部分噪声特征向量错误地划入噪声子空间（实际已是信号子空间），使得噪声子空间的维度变小，导致谱中产生多余的零点/极点，表现为伪峰。

**为什么 \(K\) 估计偏大 → 噪声子空间维度变小 → 出现更多极点（伪峰）？**

关键在于 MUSIC 谱的定义：

\[
P_{MUSIC}(\theta) = \frac{1}{\mathbf{a}^H(\theta) \mathbf{U}_N \mathbf{U}_N^H \mathbf{a}(\theta)}
\]

谱的**极点**（尖峰）出现在分母为零的地方，即 \(\mathbf{U}_N^H \mathbf{a}(\theta) = \mathbf{0}\)，也就是 \(\mathbf{a}(\theta)\) 完全落在噪声子空间里。

现在分两种情况：

**1. \(K\) 估计正确（\(K_{est} = K_{true}\)）**

噪声子空间维度 = \(M - K_{true}\)，它包含了所有与信号子空间正交的方向。真实信号方向满足正交条件，谱出现 \(K\) 个尖峰。其他方向不满足，谱值平坦。

**2. \(K\) 估计偏大（\(K_{est} > K_{true}\)）**

此时我们把前 \(K_{est}\) 个特征向量划为"信号子空间"，但其中第 \(K_{true}+1\) 到第 \(K_{est}\) 个实际上是**噪声特征向量**（特征值 ≈ \(\sigma^2\)）。

这意味着：
- "信号子空间"被扩大了，混入了一些本该属于噪声子空间的向量
- 相应地，"噪声子空间"\(\mathbf{U}_N\) 的维度从 \(M - K_{true}\) **缩小**为 \(M - K_{est}\)

**噪声子空间维度变小意味着什么？**

\(\mathbf{U}_N\) 的列数变少了，它能"张成"的子空间变小了。原来需要落在一个高维子空间里才能满足 \(\mathbf{U}_N^H \mathbf{a}(\theta) = \mathbf{0}\)，现在只需落在一个**低维**子空间里就够了。

换句话说：**正交条件的"门槛"降低了**——有更多的 \(\mathbf{a}(\theta)\) 方向能够满足分母为零的条件，即产生更多的极点。

**直观类比：**

想象噪声子空间是一个"过滤网"，只有与这个网完全正交的方向才能产生尖峰。网越大（维度越高），能通过过滤的方向越少（只有真实的 \(K\) 个信号方向）。网变小了（维度降低），过滤变松了，更多的方向"漏过去"产生伪峰。


**成因 3：阵列校准误差**

导向矢量模型假设阵元位置、增益、相位精确已知。实际中，阵元位置偏差、通道不一致性等会使真实的导向矢量偏离模型 \( \mathbf{a}(\theta) \)。MUSIC 扫描使用的是理想 \( \mathbf{a}(\theta) \)，它与真实的噪声子空间不再严格正交，谱峰的位置和数量都会出现偏差。

**成因 4：相干信号（多径）**

这是最严重也最常见的问题。当存在多径传播时，若干到达波是同一信号的相干副本：

- 直达波：\( \mathbf{a}(\theta_1) s(t) \)
- 反射波：\( \alpha \cdot \mathbf{a}(\theta_2) s(t) \)，\( \alpha \) 为复反射系数

此时信号协方差矩阵 \( \mathbf{R}_S \) **秩亏**（rank < \(K\)），信号子空间维度小于实际信号源数。特征分解后，部分相干信号的导向矢量会"泄漏"到噪声子空间中，导致 \( \mathbf{U}_N^H \mathbf{a}(\theta_k) \neq \mathbf{0} \)——真实方向的谱峰消失或减弱，同时在错误位置出现伪峰。

**对策：**
- **高信噪比** (MUSIC的缺陷)
- 空间平滑（Spatial Smoothing）：将阵列划分为重叠子阵，对各子阵协方差矩阵求平均，恢复 \( \mathbf{R}_S \) 的秩
- 增加快拍数，改善 \( \hat{\mathbf{R}}_X \) 的估计质量
- 使用稳健的 \(K\) 估计方法（MDL、AIC 等信息论准则）
- 阵列校准补偿

#### 3.8.2 谱峰抑制：真实信号峰衰减或消失

**成因 1：低 SNR**

回顾 3.4.5 节的几何解释：当噪声功率 \( \sigma^2 \) 增大，协方差矩阵的信号特征值 \( \lambda_i \) 和噪声特征值 \( \sigma^2 \) 之间的"差距"缩小，子空间界限变得模糊。在极限情况下（\( \sigma^2 \to \infty \)），所有特征值趋于相等，信号子空间和噪声子空间无法区分，MUSIC 谱退化为平坦曲线——所有峰消失。

**成因 2：角度间隔过小**

当两个信号源的角度 \( \theta_1, \theta_2 \) 非常接近时，对应的导向矢量 \( \mathbf{a}(\theta_1), \mathbf{a}(\theta_2) \) 几乎共线。尽管理论上 MUSIC 具有超分辨能力，但在有限快拍和有限 SNR 下，两个峰可能合并为一个，或被噪声基底淹没。

**成因 3：快拍数不足**

有限快拍导致样本协方差矩阵的噪声特征值散布不均（特征值扩散）。这等效于噪声子空间的估计质量下降，\( \mathbf{U}_N^H \mathbf{a}(\theta) \) 不再在真实方向上精确归零，谱峰高度降低。极端情况下（\( N < M \)），样本协方差矩阵不满秩，子空间估计彻底失效。

**综合退化规律：**

| 因素 | 对 MUSIC 谱的影响 | 对 Capon 谱的影响 |
|:---|:---|:---|
| 低 SNR | 峰降低、展宽 | 峰降低（Capon 谱峰本身的动态范围就小于 MUSIC） |
| 快拍数少 | 伪峰增多、真实峰降低 | 旁瓣抬高、分辨率下降 |
| 信号角度接近 | 峰合并 | 峰合并（比 MUSIC 更早合并） |
| 相干信号 | 峰消失、位置偏移 | 受影响但程度略轻（Capon 对相干信号有一定鲁棒性） |
| \(K\) 估计错误 | 伪峰或漏峰 | 不受影响（Capon 不需要 \(K\)） |

#### 3.8.3 Capon 的独特退化模式：相关信号对消

Capon 谱有一个独特的退化模式：当信号之间存在相关性时，Capon 波束形成器会在信号方向**相互对消**——不是产生伪峰，而是**真实峰被抑制**。这与 MUSIC 的退化机制不同：

- **MUSIC（相干信号）**：部分信号导向矢量"掉出"信号子空间，落入噪声子空间 → 正交条件被破坏 → 峰消失或偏移
- **Capon（相关信号）**：波束形成器的最优权重在约束条件下 \( \mathbf{w}^H \mathbf{a}(\theta_1)=1 \) 的同时，会在 \( \theta_2 \) 方向自动形成一个零点以最小化输出功率 → \( \theta_2 \) 处的谱峰被对消

因此，在面对相干/相关信号时，Capon 和 MUSIC **都会退化**，但退化方式不同——理解这一区别有助于在实际系统中联合使用两种方法进行交叉验证。

#### 3.8.4 工程实践建议汇总

1. **快拍数**：尽量保证 \( N \gg 2M \)（至少 \( N > M \)），以获得稳定的协方差矩阵估计
2. **\(K\) 的确定**：在 MUSIC 之前先用 MDL/AIC 准则估计信号源个数
3. **相干信号**：优先使用空间平滑预处理，或改用基于空间平滑的 MUSIC（SS-MUSIC）
4. **阵列校准**：定期进行阵列校准，补偿通道幅相不一致
5. **交叉验证**：结合 MUSIC 和 Capon 的估计结果——若二者一致，置信度高；若不一致，提示存在相干信号或模型失配

## 4. 课后总结

### 4.1 核心逻辑链：从阵列数据到子空间角度估计

本讲从时间域信号处理过渡到空间域阵列处理。核心方法如下：

1. **阵列**：阵列是空间中的频谱分析仪——将多个传感器在空间中的位置差转化为相位差，从而将角度估计问题转化为相位估计问题。

2. **窄带假设**：窄带假设 \( B \ll f_c \) 将时延 \( T \) 从信号内部提取为独立的相位因子 \( \exp(j 2\pi f T) \)。信号内容与方向信息解耦——方向矢量 \( \mathbf{a}(\theta) \) 仅取决于几何参数（载频、阵元间距、来波方向），与信号波形无关。**窄带假设是经典阵列处理理论的数学前提。**

3. **MUSIC**：阵列数据位于 \( M \) 维空间，但信号仅占据 \( K \) 维子空间。通过对协方差矩阵进行特征分解，将空间分为信号子空间与噪声子空间，利用二者正交性，通过扫描角度检验导向矢量与噪声子空间的正交程度，实现超分辨测向。

4. **子空间方法**：**用低维结构的几何特性（正交性）来解决高维数据的参数估计问题。** 不直接求解欠定的未知量，而是利用协方差矩阵的特征结构间接推断角度。

### 4.2 三种测向方法对比：CBF / Capon / MUSIC

| 维度 | 常规波束形成（CBF） | Capon（MVDR） | MUSIC |
| :--- | :--- | :--- | :--- |
| **核心运算** | 相位补偿 + 求和 | 约束优化（权重自适应） | 特征分解 + 正交性检验 |
| **谱的表达式** | \( \mathbf{a}^H \hat{\mathbf{R}}_X \mathbf{a} \) | \( 1 / (\mathbf{a}^H \hat{\mathbf{R}}_X^{-1} \mathbf{a}) \) | \( 1 / \|\mathbf{U}_N^H \mathbf{a}\|^2 \) |
| **分辨率** | 受瑞利限约束 | 超分辨（优于 CBF） | **超分辨（最优）** |
| **需要先验知识** | 无 | 无 | **需要信号源个数 \( K \)** |
| **理论谱峰高度** | 有限 | 有限 | **理论上无穷大** |
| **对相干信号** | 不受影响 | 信号对消 | **秩亏 → 失效**（需空间平滑） |
| **计算复杂度** | 低 | 中（需矩阵求逆） | 中高（需特征分解） |

### 4.3 重点概念总结

#### 4.3.1 阵列数据模型

$$
\boxed{\mathbf{X}(t) = \mathbf{A}(\boldsymbol{\theta}) \mathbf{S}(t) + \mathbf{N}(t)}
$$

- \( \mathbf{A}(\boldsymbol{\theta}) \in \mathbb{C}^{M \times K} \)：方向矩阵，每列为一个导向矢量 \( \mathbf{a}(\theta_k) \)
- \( \mathbf{S}(t) \in \mathbb{C}^{K \times 1} \)：信号复包络
- \( \mathbf{N}(t) \in \mathbb{C}^{M \times 1} \)：空间白噪声，\( \mathbb{E}[\mathbf{N}\mathbf{N}^H] = \sigma^2 \mathbf{I}_M \)

#### 4.3.2 协方差矩阵的特征结构

$$
\mathbf{R}_X = \mathbf{A} \mathbf{R}_S \mathbf{A}^H + \sigma^2 \mathbf{I}_M
$$

- 特征值：\( \lambda_1 \geq \cdots \geq \lambda_K > \lambda_{K+1} = \cdots = \lambda_M = \sigma^2 \)
- 信号子空间：前 \( K \) 个特征向量 \( \mathbf{U}_s = [\mathbf{u}_1, \cdots, \mathbf{u}_K] \)
- 噪声子空间：后 \( M-K \) 个特征向量 \( \mathbf{U}_N = [\mathbf{u}_{K+1}, \cdots, \mathbf{u}_M] \)
- 正交性：\( \mathbf{U}_s \perp \mathbf{U}_N \)，\( \text{span}\{\mathbf{U}_s\} \oplus \text{span}\{\mathbf{U}_N\} = \mathbb{C}^M \)

#### 4.3.3 MUSIC 的核心定理

$$
\text{span}\{\mathbf{U}_s\} = \text{span}\{\mathbf{A}\}
\quad\Longrightarrow\quad
\boxed{\mathbf{U}_N^H \mathbf{a}(\theta_k) = \mathbf{0},\quad k=1,\dots,K}
$$

**证明路线**：先证 \( \ker(\mathbf{U}_s^H) = \ker(\mathbf{A}^H) \)，再利用 kernel–range 对偶关系推出 span 相等。

#### 4.3.4 MUSIC 伪谱

$$
P_{\text{MUSIC}}(\theta) = \frac{1}{\mathbf{a}^H(\theta) \mathbf{U}_N \mathbf{U}_N^H \mathbf{a}(\theta)}
$$

- 当 \( \theta = \theta_k \) 时，分母趋零 → 谱峰（理论上无穷大）
- 当 \( \theta \neq \theta_k \) 时，分母非零 → 有限值

#### 4.3.5 MUSIC 与 Capon 的关系

$$
\mathbf{R}_X^{-1} = \mathbf{U}_s \boldsymbol{\Lambda}_s^{-1} \mathbf{U}_s^H + \frac{1}{\sigma^2} \mathbf{U}_N \mathbf{U}_N^H
\quad\xrightarrow{\sigma^2 \to 0}\quad
\frac{1}{\sigma^2} \mathbf{U}_N \mathbf{U}_N^H
$$

- **高 SNR 下**：Capon 谱与 MUSIC 谱仅差常数因子，二者等价
- **有限 SNR 下**：Capon 存在信号子空间"自污染"项，谱峰高度有限；MUSIC 谱峰理论上可趋于无穷
- **关系**：Capon 可视为"隐式的 MUSIC"——协方差矩阵求逆本身就在重加权子空间，MUSIC 则将这一混合显式分离开

#### 4.3.6 实际应用中的退化与对策

| 问题 | 成因 | 对 MUSIC 的影响 | 对策 |
| :--- | :--- | :--- | :--- |
| **伪峰增多** | \( K \) 估计偏大（噪声子空间维度变小，正交门槛降低） | 非信号方向产生尖峰 | MDL/AIC 准则定阶 |
| | 有限快拍（噪声特征值散布） | 噪声基底起伏形成伪峰 | 增加快拍数 \( N \gg 2M \) |
| **真实峰消失** | 相干信号（\( \mathbf{R}_S \) 秩亏） | 信号导向矢量泄漏到噪声子空间 | 空间平滑预处理 |
| | 低 SNR | 信号/噪声特征值界限模糊 | 提高 SNR 或增加阵元数 |
| **谱峰合并** | 信号角度过于接近 | 两峰合为一个宽峰 | 增加阵元数（扩大孔径） |

---

## 5. 学习检查清单：自测核心知识点掌握情况

- [ ] 能用一句话解释"阵列是什么"，以及阵列处理与时间域信号处理的核心对应关系
- [ ] 能描述均匀线性阵列（ULA）的远场平面波假设和窄带信号假设，以及它们各自为什么必要
- [ ] 能写出窄带近似 \( s(t+T) \approx s(t) \cdot \exp(j 2\pi f T) \)，并解释为什么这是阵列处理的数学前提
- [ ] 能写出多信号源的阵列数据模型 \( \mathbf{X}(t) = \mathbf{A}\mathbf{S}(t) + \mathbf{N}(t) \)，说明各矩阵的维度含义
- [ ] 能从协方差矩阵 \( \mathbf{R}_X = \mathbf{A}\mathbf{R}_S\mathbf{A}^H + \sigma^2\mathbf{I}_M \) 出发，说明其特征值的分布结构及其几何意义
- [ ] 能完整推导 \( \text{span}\{\mathbf{U}_s\} = \text{span}\{\mathbf{A}\} \)（两个方向：kernel 包含 + range-kernel 对偶）
- [ ] 能由此推出 \( \mathbf{U}_N^H \mathbf{a}(\theta_k) = \mathbf{0} \)，并写出 MUSIC 伪谱表达式
- [ ] 能解释为什么 MUSIC 的谱叫"伪谱"——它反映的是正交程度，而非功率
- [ ] 能推导 MUSIC 与 Capon 在高 SNR 下的等价关系（\( \mathbf{R}_X^{-1} \) 的谱分解）
- [ ] 能分析 \( K \) 估计偏大时为何伪峰增多（噪声子空间维度变小 → 正交门槛降低）
- [ ] 能分析相干信号为何使 MUSIC 失效（\( \mathbf{R}_S \) 秩亏 → 信号子空间维度不足）
- [ ] 能写出 MUSIC 算法的完整流程（5 步：估计协方差 → 特征分解 → 确定噪声子空间 → 扫描 → 峰值检测）

## 6. 思考题：拓展与挑战

1. **窄带假设的边界**：如果信号不满足窄带假设（如超宽带 UWB 信号），时延 \( T \) 不再能近似为相位因子。此时 MUSIC 算法还能用吗？如果不能，哪些核心推导步骤会失效？有没有可能推广 MUSIC 到宽带场景？

2. **子空间分离的"灰色地带"**：理想白噪声下，噪声特征值严格相等（\( \lambda_{K+1} = \cdots = \lambda_M = \sigma^2 \)），信号与噪声的界限清晰可辨。但在实际数据中，有限快拍导致噪声特征值不再严格相等。此时你如何判断 \( K \)？仅靠"看特征值的拐点"是否可靠？MDL 和 AIC 准则各自偏好什么场景？

3. **\( K \) 估计错误的非对称性**：\( K \) 估计偏大将部分噪声特征向量划入信号子空间，导致噪声子空间维度变小，出现伪峰。那 \( K \) 估计偏小呢？此时噪声子空间维度变大，对 MUSIC 谱会产生什么影响？为什么实践中 \( K \) 估计偏小比偏大更"灾难性"？

4. **导向矢量与特征向量的关系**：MUSIC 的核心定理是 \( \text{span}\{\mathbf{U}_s\} = \text{span}\{\mathbf{A}\} \)，但信号特征向量 \( \mathbf{u}_1, \cdots, \mathbf{u}_K \) 并不等于导向矢量 \( \mathbf{a}(\theta_1), \cdots, \mathbf{a}(\theta_K) \)——仅仅是张成同一子空间。这意味着什么？能否直接从 \( \mathbf{U}_s \) 读出角度？为什么 MUSIC 还需要谱峰搜索？

5. **MUSIC 与 Capon 的选择**：在信号源个数 \( K \) 不确定的场景下，Capon 比 MUSIC 更稳健（不需要 \( K \)）。但既然 Capon 可以通过 \( \mathbf{R}_X^{-1} \) 隐式地实现类似效果，为什么还需要 MUSIC？请从分辨率极限、对正交性检验的纯粹度、以及"信号子空间自污染"三个角度分析。

6. **相干信号与空间平滑**：两个相干信号对应同一个信号源的不同多径副本。此时 \( \mathbf{R}_S \) 的秩降为 1（而非 2）。为什么空间平滑能恢复秩？这个操作牺牲了什么？如果将阵列划分为子阵，有效孔径变小，这对分辨率有何影响？

## 7. 实验设计：基于 ESP32 天线阵列的 WiFi 信号 DOA 估计

### 7.0 实验总览：利用 WiFi CSI 在 ESP32 阵列上实现 MUSIC 测向

本实验使用多块 ESP32 开发板搭建真实射频天线阵列，通过提取 WiFi 信号的 **CSI（Channel State Information，信道状态信息）**，在物理世界中验证 MUSIC 算法的来波方向估计能力。

#### 核心思路

WiFi 信号（802.11n @ 2.4 GHz）在室内传播时，接收端每根天线收到的信号之间存在由路径差导致的**相位差**——这正是 MUSIC 所需的空间相位信息。ESP32 提供 CSI 提取接口，可直接读取出每个子载波上复数形式的信道响应 \( H_k = |H_k| e^{j\phi_k} \)，天然构成 MUSIC 算法的输入数据。

```
┌───────────────────────────────────────────────────────────┐
│               WiFi 发射源 (ESP32 Beacon)                    │
│                   θ = 待测角度                              │
│                      │                                     │
│                 d sin θ                                    │
│                   ↙                                       │
│   ════╪════════╪════════╪════════╪══  阵列法线 (θ=0°)     │
│     Ant0    Ant1     Ant2     Ant3       ← 接收阵列        │
│     ┌─┐     ┌─┐     ┌─┐     ┌─┐                          │
│     │ESP│   │ESP│   │ESP│   │ESP│   4 个接收节点           │
│     │32 │   │32 │   │32 │   │32 │                          │
│     └─┬─┘   └─┬─┘   └─┬─┘   └─┬─┘                          │
│       │USB    │USB    │USB    │USB   CSI 数据汇聚到 PC      │
│       └───────┴───┬───┴───────┘                           │
│                   ▼                                        │
│           ┌──────────────┐                                 │
│           │   PC Python  │   MUSIC 算法离线处理             │
│           │ 协方差矩阵   │                                  │
│           │ 特征分解     │  →  空间谱 → θ_est              │
│           └──────────────┘                                 │
└───────────────────────────────────────────────────────────┘
```

#### 为什么 MUSIC + WiFi CSI 适合做这个实验

1. **窄带假设成立**：802.11n 单个子载波带宽 312.5 kHz，载频 2.412 GHz，带宽比 \( B/f_c \approx 1.3 \times 10^{-4} \)，窄带假设成立，时延可精确简化为相位旋转。
2. **天然的多快拍**：一次 WiFi 数据帧在 20 MHz 带宽下提供 64 个子载波，每个子载波的 CSI 是复数，64 个子载波 = 64 个快拍，足够估计协方差矩阵。
3. **相位信息精度高**：CSI 提供的是基带等效信道响应，相位信息直接编码了不同天线之间的路径差，无需额外同步机制。

#### 物理参数计算

| 参数 | 符号 | 值 | 计算依据 |
| :--- | :--- | :--- | :--- |
| WiFi 信道 1 中心频率 | \( f_c \) | 2.412 GHz | 802.11 标准 |
| 波长 | \( \lambda = c/f_c \) | **12.44 cm** | \( c = 3 \times 10^8 \) m/s |
| 阵列半波长间距 | \( d = \lambda/2 \) | **6.22 cm** | 避免栅瓣 |
| 子载波间隔 | \( \Delta f \) | 312.5 kHz | 20 MHz / 64 |
| 窄带比 | \( B/f_c \) | \( \approx 1.3 \times 10^{-4} \) | 远小于 1 |

---

### 7.1 物料清单（BOM）

| 器材 | 规格/型号 | 数量 | 用途 | 参考单价 |
| :--- | :--- | :-: | :--- | :--- |
| **ESP32-S3-DevKitC** | 双核 240 MHz，带 IPEX 天线座 | 5 | 1 个发射 + 4 个接收 | ¥25—40 |
| **2.4 GHz 外置棒状天线** | SMA 公头，增益 3 dBi | 4 | 接收阵列的阵元 | ¥3—8 |
| **IPEX-SMA 转接线** | IPEX 1 代转 SMA 母座，15 cm | 4 | 将 ESP32 板载 IPEX 引出为 SMA | ¥2—5 |
| **铝型材导轨** | 2020 欧标，长度 30 cm | 1 | 阵列支架 | ¥8 |
| **T 型螺母 + 螺丝** | M4，适配 2020 型材 | 8 套 | 固定 ESP32 板 | ¥5 |
| **亚克力/3D 打印托盘** | ESP32 板卡固定座 | 4 块 | 方便调节方向和间距 | 自备 |
| **USB Hub (4 口)** | USB 2.0，带外部供电 | 1 | 4 块 ESP32 同时连接 PC | ¥15 |
| **USB Type-C 数据线** | 带数据传输，50 cm | 5 | 供电 + 串口通信 | ¥5/根 |
| **量角器 + 激光测距仪** | — | 1 套 | 标定信号源角度 | 自备 |
| **移动电源 / 锂电池** | 5V 输出 | 1 | 给发射端 ESP32 独立供电（避免 USB 线缆影响） | 自备 |

> 总预算约 **¥120—200**（不含 PC），远低于一套矢量网络分析仪。

#### 天线选型要点

- **必须使用外置天线**：ESP32 板载 PCB 天线位置固定在板上，无法调节阵元间距。必须通过 IPEX 座引出到外置 SMA 天线。
- **天线必须一致**：4 根天线应为同一型号、同一批次，以保证幅相一致性。否则阵元间的固有幅相失配会引入系统误差，需要额外的校准步骤。
- **极化方向对齐**：所有天线垂直安装（垂直极化），与 WiFi 发射天线的极化方向一致。

---

### 7.2 硬件搭建与系统集成

#### 7.2.1 天线阵列的机械安装

1. 在 2020 铝型材导轨上，按间距 \( d = 6.2 \) cm 标记 4 个安装位
2. 每个安装位固定一块 ESP32-S3 板，使 IPEX 插座朝向一致
3. 通过 IPEX-SMA 转接线将外置棒状天线垂直安装在导轨上方
4. **关键**：确保 4 根天线的相位中心在同一水平线上，天线杆体互相平行

```
        ↕ 天线 (垂直极化)
        │    │    │    │
        ■    ■    ■    ■    ← SMA 天线座
        │    │    │    │
   ═════╪════╪════╪════╪══  2020 铝型材导轨
    ●   ●   ●   ●   ●       ← ESP32 板卡 (用螺母固定于型材槽中)
   Ant0 Ant1 Ant2 Ant3
   ├─ d ─┤
   ├───── (M-1)d = 3d = 18.7 cm 孔径 ────┤
```

- 阵列法线方向为导轨的垂直平分线方向（即天线平面前方）
- 将整个阵列置于桌面边缘或三脚架上，确保前方 2 米内无障碍物
- 阵列高度约 1—1.5 米（模拟典型 AP 或基站高度）

#### 7.2.2 发射源的设置

**方案 A（单源实验）**：1 块 ESP32 作为 WiFi Beacon 发射源

- 烧录 Arduino 固件，配置为 `WiFi.softAP("DOA_Test", NULL, 1)` （信道 1，2.412 GHz）
- 设置 Beacon 间隔为 100 ms（`esp_wifi_set_config` 的 `beacon_interval` 参数），确保每秒 10 次 CSI 快照机会
- 用移动电源供电，放置在高约 1.2 m 的三脚架或支架上
- 通过激光测距仪和量角器标定其相对于阵列中心的距离 \( R \) 和角度 \( \theta \)

**方案 B（双源实验）**：2 块 ESP32 各自创建一个 SSID 不同的 AP，放置在两个已知角度位置

- 两个 AP 设置相同信道（信道 1），各自以 100 ms 间隔发送 Beacon
- 两个 AP 之间的距离应满足：阵列处看向两者的夹角 ≥ 瑞利限

**距离要求**：发射源距离阵列 \( R \geq \frac{2D^2}{\lambda} \)（远场条件），其中 \( D = (M-1)d = 18.7 \) cm 为阵列孔径。计算得 \( R \geq \frac{2 \times 0.187^2}{0.1244} \approx 0.56 \) m。实际操作中取 \( R \geq 1.5 \) m 即可保证远场平面波近似成立。

#### 7.2.3 USB 连接与供电

```
┌─────────┐     ┌──────────┐
│ 发射ESP32│     │USB Hub   │ ← 外部 5V 供电
│ (充电宝) │     │          │
└─────────┘     ├─ Port1 ── ESP32-RX0 (天线 0)
                ├─ Port2 ── ESP32-RX1 (天线 1)
                ├─ Port3 ── ESP32-RX2 (天线 2)
                ├─ Port4 ── ESP32-RX3 (天线 3)
                └─ Upstream ── PC USB
```

- 4 块接收端 ESP32 全部通过 USB Hub 连接到同一台 PC
- PC 端使用 `pyserial` 同时打开 4 个 COM 口，分别读取各节点的 CSI 数据
- **发射端 ESP32 必须用独立移动电源供电**——不能和接收端共用 USB Hub，否则发射信号会通过 USB 地线耦合到接收端，产生虚假的强相关

---

### 7.3 CSI 数据采集与 MUSIC 快拍构造

#### 7.3.1 CSI 是什么

ESP32 的 WiFi 驱动在每收到一个 WiFi 帧时，可以报告该帧的 CSI——即物理层对每个子载波的信道估计：

$$H(k) = |H(k)| \cdot e^{j\phi(k)}, \quad k = 0, 1, \ldots, N_{\text{sub}}-1$$

- 802.11n 20 MHz 模式下，\( N_{\text{sub}} = 64 \) 个子载波（含导频和空子载波，有效约 52 个）
- \( |H(k)| \)：幅度响应（对应 RSSI 在每个子载波上的分解）
- \( \phi(k) \)：相位响应（**包含不同天线之间的路径差信息**）

#### 7.3.2 MUSIC 的输入构造

对于第 \( k \) 个子载波，4 块接收 ESP32 各自提取的 CSI 值构成一个 \( 4 \times 1 \) 的复向量——这正是 MUSIC 所需的一个"快拍"：

$$\mathbf{x}_k = \begin{bmatrix} H_0(k) \\ H_1(k) \\ H_2(k) \\ H_3(k) \end{bmatrix} \in \mathbb{C}^{4 \times 1}$$

其中 \( H_m(k) \) 为第 \( m \) 根天线在子载波 \( k \) 上的 CSI 值。

对 64 个子载波（去掉导频和空子载波后约 52 个有效数据子载波），我们获得约 **52 个快拍**，足以估计 \( 4 \times 4 \) 的协方差矩阵。

#### 7.3.3 时间同步要求

4 块接收 ESP32 必须在**同一 WiFi 帧**上提取 CSI。由于 Beacon 帧是广播帧，4 块 ESP32 均在微秒级时间窗口内收到同一帧——这个级别的同步对 MUSIC 来说足够。

**实际做法**：
- 所有接收 ESP32 侦听同一 WiFi 信道（信道 1），开启混杂模式（`esp_wifi_set_promiscuous(true)`）
- 每收到一个帧，提取其源 MAC 地址和序列号
- 用 **MAC + 序列号** 作为帧的唯一标识，在 PC 端对齐 4 路数据
- 丢弃未对齐的帧（某路丢帧的情况）

---

### 7.4 实验 1：单源来波方向估计（基础验证）

**目的**：验证 ESP32 天线阵列 + MUSIC 能否正确估计一个 WiFi 信号源的来波方向。

**步骤**：
1. 将单个发射 ESP32（AP 模式）放置在阵列前方约 2 m 处，角度 \( \theta = 0^\circ \)（正对阵列法线）
2. 4 块接收 ESP32 各采集 200 个 Beacon 帧的 CSI 数据（约 20 秒），通过 USB 发送到 PC
3. PC 端用 52 个有效子载波的 CSI 构造 \( \mathbf{x}_k \) 向量，计算 \( 4 \times 4 \) 样本协方差矩阵 \( \hat{\mathbf{R}}_X \)
4. 特征分解，取最小特征值对应的特征向量构成噪声子空间 \( \mathbf{U}_N \)（此时 \( K=1 \)）
5. 扫描 \( \theta \in [-90^\circ, 90^\circ] \)，步长 \( 0.5^\circ \)，计算 MUSIC 伪谱
6. 改变发射源角度为 \( \theta = -45^\circ, -30^\circ, 0^\circ, +30^\circ, +45^\circ \)，重复实验
7. 绘制各角度下的 MUSIC 空间谱叠加图，标注谱峰位置与真实角度的偏差

**定量指标**：
- 测角误差：\( \Delta\theta = |\hat{\theta}_{\text{MUSIC}} - \theta_{\text{true}}| \)
- 谱峰宽度：半功率波束宽度（3 dB 宽度），反映角度分辨率
- RMS 误差：每个角度采集 10 组数据后的均方根误差

**预期结果**：
- 谱峰位置与真实角度基本吻合（误差 ≤ 5°）
- 正前方（\( 0^\circ \)）估计最准确，大角度（\( \pm 60^\circ \) 以上）精度下降（导向矢量非线性加剧 + 天线方向图衰减）
- 在某些角度可能出现镜像峰（前后模糊，ULA 本身无法区分 ±θ，需借助天线方向性或阵列朝向排除）

**思考**：如果测量误差显著大于预期，可能的原因有哪些？列出至少三个硬件层面的误差源（天线一致性、多径反射、近场效应、通道间幅相失配等），并讨论哪些可以通过校准消除、哪些是本质限制。

---

### 7.5 实验 2：双源超分辨能力（核心实验）

**目的**：验证 MUSIC 在角度差小于瑞利限时，是否仍能分辨两个 WiFi 信号源——这是 MUSIC 超分辨能力的最直接证明。

**背景**：
- 瑞利限（瑞利分辨率判据）：两个等功率信号源可分辨的最小角度差为
  $$\Delta\theta_{\text{Rayleigh}} \approx \frac{0.891 \lambda}{M d \cos\theta} \cdot \frac{180^\circ}{\pi}$$
  对于 \( M=4, d=\lambda/2, \theta=0^\circ \)：\( \Delta\theta_{\text{Rayleigh}} \approx 25.5^\circ \)
- 如果 MUSIC 能分辨 \( 10^\circ—15^\circ \) 的角度差（远小于 25.5°），即验证了超分辨能力

**步骤**：
1. 准备两个发射 ESP32，分别创建 SSID 为 `SRC_A` 和 `SRC_B` 的 AP，同信道（信道 1）
2. 使用量角器和激光测距仪标定两个发射源位置：
   - **组 1（超分辨测试）**：\( \theta_A = -10^\circ, \theta_B = +10^\circ \)（角度差 20°，略小于瑞利限）
   - **组 2（极限测试）**：\( \theta_A = -7^\circ, \theta_B = +7^\circ \)（角度差 14°，远小于瑞利限）
   - **组 3（大角度间隔）**：\( \theta_A = -30^\circ, \theta_B = +30^\circ \)（角度差 60°，轻松分辨的参考组）
3. 两个发射源均距离阵列约 2 m
4. 接收阵列采集 CSI 数据，按源 MAC 地址分离两个源的 CSI
5. 构造 MUSIC 协方差矩阵时，将两个源的 CSI 向量**叠加**（模拟两个源同时存在）：
   $$\mathbf{x}_k = \mathbf{x}_{k}^{(A)} + \mathbf{x}_{k}^{(B)}$$
6. 设定 \( K=2 \)，取前 2 大特征值对应的特征向量构成信号子空间
7. 扫描 MUSIC 伪谱，检测谱峰数量和位置
8. 同组数据也用 CBF（Bartlett 波束形成）生成空间谱，与 MUSIC 谱并列绘制对比

**定量指标**：
- 能否分辨？（两峰之间存在 ≥ 3 dB 的凹陷）
- 角度估计误差：每个峰的 \( \Delta\theta \)
- 谱峰分辨率对比：MUSIC vs CBF 的分辨临界角

**预期结果**：
- 组 1（20° 角差）：MUSIC 大概率成功分辨两个峰；CBF 谱可能只呈现一个宽峰或肩峰
- 组 2（14° 角差）：MUSIC 可能仍能分辨（取决于 SNR 和多径环境），CBF 必然失败
- 组 3（60° 角差）：两者均能分辨，但 MUSIC 的峰更尖锐
- 室内多径可能导致出现额外的伪峰——这正是真实环境的复杂性

**思考**：
1. 实际环境中，两个 WiFi 源的发射功率可能不同（不对称功率）。当 \( P_A \gg P_B \) 时，弱信号源的 MUSIC 谱峰会如何变化？低功率源被高功率源"淹没"的机理是什么？
2. 实验中用 MAC 地址分离源来获取"干净"的 CSI——如果两个源恰好同时发射（Beacon 碰撞），此时无法按 MAC 分离，你如何处理这种相干场景？

---

### 7.6 实验 3：室内多径环境下的 MUSIC 表现

**目的**：观察真实室内多径反射对 MUSIC 算法的影响——这是纯仿真无法获得的体验。

**步骤**：
1. 使用单个发射 ESP32，放置于 \( \theta_{\text{true}} = 20^\circ \)，距离 2 m
2. 在三种环境下各采集一组 CSI 数据：
   - **环境 A（吸波/空旷）**：阵列前方为开阔空间，周围无大型金属物体
   - **环境 B（单反射面）**：在阵列另一侧约 1 m 处置一块大型金属板（如白板、铝箔覆盖的硬纸板），制造一个可预测的反射路径
   - **环境 C（典型室内）**：正常实验室/教室环境，有墙壁、家具、金属柜等
3. 对每组数据运行 MUSIC（设定 \( K=2 \) 或 \( K=3 \) 以尝试捕获反射路径）
4. 对比三种环境下的 MUSIC 空间谱

**预期结果**：
- 环境 A：单峰，位置接近 \( 20^\circ \)
- 环境 B：除 \( 20^\circ \) 的直达峰外，在镜像方向（金属板对面的对称角度）出现一个较弱但可辨识的峰——这是相干信号场景（直达 + 反射发自同一源），MUSIC 可能因秩亏而失效。**此时可进一步验证空间平滑的效果。**
- 环境 C：可能出现多个峰，其中一些是真实的反射路径，一些是 \( K \) 估计不当导致的伪峰

**思考**：
1. 在环境 B 中，如果反射峰和直达峰同时被 MUSIC 检测到，它们的相位关系是什么？反射路径的额外传播距离是否会导致两者在子载波维度上呈现相位旋转差异？
2. WiFi 信号在室内传播时，墙壁反射的衰减约为 3—6 dB（取决于材料）。如果反射路径的功率低于直达路径 10 dB 以上，MUSIC 还能检测到它吗？这个检测极限由什么因素决定？

---

### 7.7 实验 4：相位测角 vs 幅度测角（MUSIC vs RSSI）

**目的**：对比基于相位信息的 MUSIC 测角与基于幅度信息的 RSSI 测角，直观感受为什么"相位才是阵列处理的灵魂"。

**步骤**：
1. 使用单发射源，在 \( \theta = -60^\circ, -40^\circ, -20^\circ, 0^\circ, 20^\circ, 40^\circ, 60^\circ \) 的 7 个角度上分别采集数据
2. 从 CSI 数据中提取两种信息：
   - **幅度**：\( |H_m(k)| \)，取所有子载波的平均值作为每根天线的 RSSI 代理
   - **相位**：\( \angle H_m(k) \)，用于 MUSIC
3. **RSSI 方法**：用 4 根天线的 RSSI 比值查表或最小二乘拟合来估计角度（这是很多简单 IoT 定位系统的做法）
4. **MUSIC 方法**：利用完整的 CSI 相位信息
5. 绘制两种方法的 **估计角度 vs 真实角度** 散点图，计算各自的 RMSE

**预期结果**：
- RSSI 方法：大角度时误差急剧增大（天线方向图非线性），小角度灵敏度极低（4 根天线间距仅 6 cm，RSSI 差异太小），整体 RMSE 在 \( 10^\circ—30^\circ \) 量级
- MUSIC 方法：全角度范围内 RMSE 在 \( 2^\circ—8^\circ \) 量级
- 直观说明：**阵列测角的核心是相位信息而非幅度**，MUSIC 正是因为利用了相位信息才有超分辨能力

**思考**：为什么间距 6 cm 的两根天线收到的信号幅度几乎相同，但相位却可以有显著差异？用 \( d \sin\theta / \lambda \) 计算 \( \theta = 30^\circ \) 时相邻天线的相位差，体会"相位比幅度对微小空间位移更敏感"这一物理事实。

---

### 7.8 实验报告要求

每个实验需提交：
1. **硬件搭建照片**：阵列整体照片、天线间距测量特写、发射源标定场景照
2. **MUSIC 空间谱图**：每个实验条件下清晰标注的角度—谱曲线（dB 坐标）
3. **定量结果表**：峰位置、偏差、RMS 误差、检测率
4. **与 CBF/理论极限的对比**：在谱图上叠加 CBF 曲线和瑞利限标记线
5. **异常分析**：对实验中出现伪峰、峰消失、偏移等现象的物理分析和改进建议

---

---

## 附录：子空间方法的几何直觉与理论延伸

### A.1 高维数据中的低维结构：子空间的几何直觉

MUSIC 的核心是一个几何命题：**阵列接收数据位于高维空间（维度 = 阵元数 \( M \)），但信号仅占据一个低维子空间（维度 = 信号源数 \( K \)）。**

考虑无噪声情形下的阵列快拍：

$$
\mathbf{x}(t) = \mathbf{A}(\boldsymbol{\theta}) \mathbf{s}(t) = \sum_{k=1}^{K} \mathbf{a}(\theta_k) s_k(t)
$$

关键点在于：**无论 \( s_k(t) \) 如何随时间变化，\( \mathbf{x}(t) \) 始终是 \( K \) 个固定向量 \( \{\mathbf{a}(\theta_1), \cdots, \mathbf{a}(\theta_K)\} \) 的线性组合。** 这 \( K \) 个导向矢量张成信号子空间，所有快拍数据都落在这个子空间内，而非散布在整个 \( M \) 维空间。

从几何角度看：数据在高维空间中被压缩到一个低维平面上，该平面的方向由导向矢量 \( \mathbf{a}(\theta_k) \) 决定，而导向矢量由物理角度 \( \theta_k \) 决定。因此，**找到子空间的方向，就等于找到了信号源的角度。** 角度估计问题可转化为子空间识别问题。

### A.2 噪声干扰下的信号/噪声子空间分裂

当噪声存在时，数据向量变为：

$$
\mathbf{x}(t) = \underbrace{\sum_{k=1}^{K} \mathbf{a}(\theta_k) s_k(t)}_{\text{信号部分，落在 }K\text{ 维子空间内}} + \underbrace{\mathbf{n}(t)}_{\text{噪声，充满全空间}}
$$

白噪声在所有方向上具有相同方差，各方向互不相关。因此，噪声会将数据点从信号子空间推出，散布到整个 \( M \) 维空间，但这种散布是各向同性的——噪声不偏向任何特定方向。

于是，整个 \( M \) 维空间可以被分解为两个正交的子空间：

- **信号子空间 \( \mathcal{S} \)**：由 \( K \) 个导向矢量张成，维度为 \( K \)。这是信号贡献的主要方向。
- **噪声子空间 \( \mathcal{N} \)**：信号子空间的正交补空间，维度为 \( M - K \)。这是纯噪声占主导的方向。

这两个子空间互为正交补空间：\( \mathcal{S} \perp \mathcal{N} \)，\( \mathcal{S} \oplus \mathcal{N} = \mathbb{C}^M \)。

这一结构是 MUSIC 算法成立的全部几何基础。

### A.3 子空间辨识的工具：特征分解

问题在于：我们只有离散的快拍数据 \( \{\mathbf{x}(t)\}_{t=1}^{N_s} \)，并不知道哪些方向属于信号子空间、哪些属于噪声子空间。

答案是协方差矩阵的特征分解。

当噪声为白噪声时，协方差矩阵具有如下结构：

$$
\mathbf{R}_X = \mathbf{A} \mathbf{R}_S \mathbf{A}^H + \sigma^2 \mathbf{I}_M
$$

这一结构使得特征分解具有特殊的性质。设 \( \{\lambda_i\}_{i=1}^{M} \) 为 \( \mathbf{R}_X \) 的特征值，按降序排列：

$$
\lambda_1 \geq \lambda_2 \geq \cdots \geq \lambda_K > \lambda_{K+1} = \lambda_{K+2} = \cdots = \lambda_M = \sigma^2
$$

前 \( K \) 个特征值（大于 \( \sigma^2 \)）对应的特征向量张成信号子空间；后 \( M-K \) 个特征值（等于 \( \sigma^2 \)）对应的特征向量张成噪声子空间。

**特征分解将数据的总方差按方向重要性排序。** 信号子空间对应的方向方差大（信号贡献额外能量），噪声子空间对应的方向方差小且相等（仅有噪声的均匀贡献）。

### A.4 子空间方法的优势：角度估计无需已知信号波形

子空间方法绕开了传统参数估计中信号波形未知的问题。

传统的最大似然估计需要同时估计角度和信号波形——信号波形 \( \mathbf{s}(t) \) 的维度随快拍数线性增长，这是一个高维非线性优化问题，计算量大且容易陷入局部最优。

子空间方法不同。特征分解将信号子空间和噪声子空间分离后，不需要信号波形 \( \mathbf{s}(t) \) 的具体值，只需利用正交性条件：

$$
\mathbf{a}(\theta_k) \perp \mathcal{N} \quad \Longleftrightarrow \quad \mathbf{U}_N^H \mathbf{a}(\theta_k) = \mathbf{0}
$$

**信号波形完全从问题中消失了**——无论 \( s_k(t) \) 是正弦波、脉冲还是随机信号，只要满足窄带条件，导向矢量与噪声子空间的正交性始终成立。

这就是 MUSIC 实现超分辨的原因：它将角度估计从高维非线性问题转化为一维的扫描搜索问题。

### A.5 非线性视角

子空间识别本质上是非线性问题。尽管使用了线性代数工具，但涉及的步骤是非线性的：

1. 计算协方差矩阵（二次型）；
2. 求解特征值问题；
3. 根据特征值大小划分子空间（阈值判定）。

此外，导向矢量 \( \mathbf{a}(\theta) \) 中 \( \theta \) 出现在指数函数的自变量中（\( \exp(j 2\pi d \cos\theta / \lambda) \)），从数据到角度的映射是强非线性的。

MUSIC 的非线性体现在三个层面：

**层面一：模型定阶。** 实际中不知道信号源个数 \( K \)，需要通过 AIC/MDL 等准则判断，这是非线性决策问题。

**层面二：谱函数的非线性。** \( P_{\text{MUSIC}}(\theta) = 1 / \|\mathbf{U}_N^H \mathbf{a}(\theta)\|^2 \) 是 \( \theta \) 的高度非线性函数，分母为零处产生尖峰。这种尖锐的谱结构带来了超分辨能力。

**层面三：投影算子的数据依赖性。** \( \mathbf{U}_N \mathbf{U}_N^H \) 从数据中估计得到，而非预先给定，整个处理流程是一个自适应非线性系统。

**子空间方法用线性代数工具解决了非线性的几何结构提取问题。**

### A.6 子空间方法的延伸：ESPRIT、Root-MUSIC 与统一框架

子空间思想是多种阵列信号处理算法的共同基础：

- **ESPRIT（Estimation of Signal Parameters via Rotational Invariance Techniques）**：利用阵列的平移不变性，通过两个子阵之间的旋转关系直接求解角度，避免了 MUSIC 的谱峰搜索，将问题进一步简化为特征值分解。

- **WSF（Weighted Subspace Fitting）**：将信号子空间的估计与导向矢量张成的子空间进行加权拟合，在统计最优意义下估计角度。

- **Root-MUSIC**：利用 ULA 导向矢量的范德蒙德结构，将谱峰搜索转化为多项式求根问题，以更低的计算量获得相同精度。

- **子空间 DOA 估计的统一框架**：所有子空间方法共享同一个基本前提——信号子空间与噪声子空间的正交性。不同算法的区别在于如何利用这一正交性来提取角度信息。

子空间方法的核心可概括为：**用低维结构的几何特性（正交性）来解决高维数据的参数估计问题。** 它从数据中提取的不是具体数值，而是一个几何结构（子空间），然后从结构中读出物理参数。这种"结构优先、参数其次"的思路，使其在低信噪比、少快拍、信号相关等困难条件下仍能有效工作。

<div style="page-break-before: always;"></div>