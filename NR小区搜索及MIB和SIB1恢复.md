## NR小区搜索及MIB和SIB1恢复

* *本示例演示了如何使用5G Toolbox™对实时gNodeB信号进行同步、解调和解码*
* *该示例解码MIB和第一个系统信息块（SIB1）*
* *解码MIB和SIB1需要一个全面的接收机，能够解调和解码大部分下行链路信道和信号。*

---

### 简介

在用户设备（UE）与网络通信之前，它必须执行小区搜索和选择程序，并获得初始系统信息。该过程的第一步是获取帧同步、找出小区标识以及解码MIB和SIB1。本示例展示了如何使用5G工具箱执行这些步骤。

您可以使用本示例捕获的I/Q采样波形，也可以使用nrWaveformGenerator生成包含同步信号（SS）突发和SIB1的本地波形。对于本地生成的波形，示例执行以下步骤：

* **波形生成：**
  * 使用5G工具箱中的下行链路波形发生器配置和生成同步信号突发，其中包含**MIB、CORESET 0、PDCCH和携带SIB1的PDSCH。**
  * 发射机可以提高一个SS块的信噪比，但不进行波束成形。有关SSB波束成形的更多信息，请参阅NR SSB Beam Sweeping。

* **AWGN：**对波形施加加性白高斯噪声（AWGN）
* **接收器：** 对接收波形进行各种同步和解调处理，以确定系统帧号、小区标识和 SSB，并解码 MIB。这些都是对 PDCCH 中的下行链路控制信息（DCI）进行盲解码所需的信息。接收器利用 DCI 配置 PDSCH 解调器，解码 DL-SCH 并最终恢复 SIB1。

**这些图显示了接收器内部的处理步骤：**

---

![img](https://picture-cloud-storage-pyp.oss-cn-beijing.aliyuncs.com/img/202307201514196.png)

![img](https://picture-cloud-storage-pyp.oss-cn-beijing.aliyuncs.com/img/202307191129578.png)

---

### 接收器配置

主要代码：

```matlab
  % 生成包含 SS burst 和 SIB1 的波形
    % 配置单元格标识
    config = struct();
```

---

**struct数组：**

*结构体数组*是使用名为*字段*的数据容器将相关数据组合在一起的数据类型。每个字段都可以包含任意类型的数据。可以使用 structName.fieldName格式的圆点表示法来访问字段中的数据。

---

```matlab
    config.NCellID = 102;
    
    % 配置 SS 突发
    config.BlockPattern = 'Case B';         % FR1: 'Case A','Case B','Case C'. FR2: 'Case D','Case E'
    config.TransmittedBlocks = ones(1,8);   % 传输的 SS 数据块的位图
    config.SubcarrierSpacingCommon = 15;    % SIB1 子载波间隔（千赫）（FR1 为 15 或 30，FR2 为 60 或 120
    config.EnableSIB1 = 1;                  % 设置为 0 时禁用 SIB1

    % 设置 NR 频段所需的最小信道带宽，以便在 FR1 中配置 CORESET 0（见 TS 38.101-1 表 5.3.5-1）。
    % 在 FR1 中配置 CORESET 0（见 TS 38.101-1 表 5.3.5-1）
    config.MinChannelBW = 5; % 5, 10, 40 MHz

    % 通过提高功率（从而提高信噪比）引入波束成形增益。
    % 提高一个 SSB 及相关 SIB1 PDCCH 和 PDSCH 的功率（从而提高 SNR）
    boost = 6; % SNR boost in dB
    config.Power = zeros(size(config.TransmittedBlocks));
    config.Power(1) = boost; % boost the first SSB

    % 配置并生成包含 SS 脉冲串和 SIB1 的波形
    wavegenConfig = hSIB1WaveformConfiguration(config);
    [txWaveform,waveInfo] = nrWaveformGenerator(wavegenConfig);
    txOfdmInfo = waveInfo.ResourceGrids(1).Info;
```

---

**hSIB1WaveformConfiguration：**

* WAVEGENCONFIG = hSIB1WaveformConfiguration(CONFIG)会创建一个nrDLCarrierConfig 配置对象

* 用于生成携带主信息块的 SS 突发和携带第一系统的**控制和数据通道**

* **控制和数据通道**携带第一个系统信息块 (SIB1)的控制和数据通道

---

**nrDLCarrierConfig：**

**解释：**nrDLCarrierConfig 对象用于设置单分量载波 5G 下行链路波形的参数。调用 nrWaveformGenerator 函数时，使用该对象配置 5G 下行链路波形生成。

**该对象定义了下行链路波形的这些方面：**

* Frequency range**（频率范围）** 
* Channel bandwidth**（信道带宽）**
* Cell identity**（小区标识）**

* Waveform duration**（波形持续时间）**

* Subcarrier spacing (SCS) carriers**（子载波间隔【SCS】载波）**

* Bandwidth parts (BWPs)**（带宽部件【BWP】）**

* Synchronization signal (SS) burst**（同步信号（SS）突发）**

* Control resource sets (CORESETs)**（控制资源集）（CORESET）**

* Search spaces**（搜索空间）**

* Physical downlink control channel (PDCCH) and PDCCH demodulation reference signal (DM-RS) **（物理下行链路控制信道 (PDCCH) 和 PDCCH 解调参考信号 (DM-RS)）**

* Physical downlink shared channel (PDSCH), PDSCH DM-RS, and PDSCH phase-tracking reference signal (PT-RS) **（物理下行链路共享信道 (PDSCH)、PDSCH DM-RS 和 PDSCH 相位跟踪参考信号 (PT-RS)）**

* Channel state information reference signal (CSI-RS) **（信道状态信息参考信号 (CSI-RS)）**

---

**nrWaveformGenerator：**

**解释：**[wave,info] = nrWaveformGenerator(cfg) 为指定的配置 cfg 生成 5G NR 波形波。输入 cfg 指定了单个或多个子载波间隔（SCS）载波和带宽部分（BWP）的下行或上行配置参数。

* 如果 cfg 是 nrDLCarrierConfig 对象，配置还会指定同步信号 (SS) 突发、控制资源集 (CORESET)、搜索空间、物理下行链路控制信道 (PDCCH) 和相关解调参考信号 (DM-RS)、物理下行链路共享信道 (PDSCH) 和相关 DM-RS 和相位跟踪参考信号 (PT-RS)，以及信道状态信息参考信号 (CSI-RS)。
* 该函数还会返回一个结构信息，其中包含有关资源网格和波形资源的信息。

函数中txWaveform为波形信息，可以在matlab中绘制出函数图像，由于txWaveform是复数信息，故绘制出的仅为时间关于幅度的函数。

```matlab
x = 1 : 1 : size(txWaveform)
plot(x,txWaveform)
xlabel('time')
ylabel('txWaveform')
```

**不过绘制的时候系统自动忽略了一些参数：**![image-20230720142409175](https://picture-cloud-storage-pyp.oss-cn-beijing.aliyuncs.com/img/202307201424203.png)

**结果图如下：**

![image-20230720142433705](https://picture-cloud-storage-pyp.oss-cn-beijing.aliyuncs.com/img/202307201424731.png)

**waveInfo：**就是返回的结构信息

![image-20230720142611699](https://picture-cloud-storage-pyp.oss-cn-beijing.aliyuncs.com/img/202307201426721.png)

代码中的**waveInfo.ResourceGrids(1).Info**信息也可以在控制台查看：

![image-20230720142714216](https://picture-cloud-storage-pyp.oss-cn-beijing.aliyuncs.com/img/202307201427244.png)

---

```matlab
% 在波形中添加白高斯噪声。
% 请注意，信噪比仅适用于增强的 SSB / SIB1
rng('default'); % 重置随机数生成器
SNRdB = 20; % AWGN 的信噪比
rxWaveform = awgn(txWaveform,SNRdB-boost,-10*log10(double(txOfdmInfo.Nfft)));
```

---

**awgn：**

* Y = awgn(X,snr,signalpower) 接受以 dBW 为单位的输入信号功率值。

* 要在添加噪声前测量 X 的功率，请将 signalpower 指定为 "实测"。

* 如果输入信号功率因衰减而随时间变化，且信道相干时间大于输入持续时间，那么在循环中重复调用 awgn 函数时，"测量 "选项不会生成所需的平均信噪比。

---

```matlab
    % 配置接收器
    % 采样率
    sampleRate = txOfdmInfo.SampleRate;

    % 符号相位补偿频率（赫兹）。
    % 函数nrWaveformGenerator 不对生成的波形进行符号相位补偿。
    % 生成波形。
    fPhaseComp = 0; %载波中心频率（赫兹）

    % 最小通道带宽（兆赫）
    minChannelBW = config.MinChannelBW;

    % 在接收器上配置必要的突发参数
    refBurst.BlockPattern = config.BlockPattern;
    % n = numel(A) 返回数组 A 中的元素数目 n 等同于 prod(size(A))。
    refBurst.L_max = numel(config.TransmittedBlocks); 

% 从配置的突发和接收机参数中获取 OFDM 信息
nrbSSB = 20;
scsSSB = hSSBurstSubcarrierSpacing(refBurst.BlockPattern); % 这是一个获取子载波间隔的函数
rxOfdmInfo = nrOFDMInfo(nrbSSB,scsSSB,'SampleRate',sampleRate);
```

---

**nrOFDMInfo：**

* info = nrOFDMInfo(载波) 为指定的载波配置参数提供与正交频分复用（OFDM）调制相关的尺寸信息。

* info = nrOFDMInfo(nrb,scs) 为指定的资源块数量 nrb 和子载波间隔 scs 提供 OFDM 信息。
* info = nrOFDMInfo(,名称,值) 除前面任何语法中的输入参数外，使用一个或多个名称-值对参数指定选项。

---

```matlab
% 显示接收波形的频谱图
figure;
nfft = rxOfdmInfo.Nfft;
spectrogram(rxWaveform(:,1),ones(nfft,1),0,nfft,'centered',sampleRate,'yaxis','MinThreshold',-130);
title('Spectrogram of the Received Waveform')
```

---

整段代码运行结果图：

![image-20230720145839483](https://picture-cloud-storage-pyp.oss-cn-beijing.aliyuncs.com/img/202307201458525.png)

---

### PSS 搜索和频率偏移校正

接收器按照以下步骤进行 PSS 搜索和粗频率偏移估算：

* 用候选频率偏移对接收到的波形进行移频。候选偏移间隔为半个子载波。使用 searchBW 控制频率偏移搜索带宽。

* 将频率偏移后的接收波形与三个可能的 PSS 序列（NID2）逐一**相关**，并提取最强的相关峰值。**参考 PSS 序列的频率居中**。因此，最强相关峰值提供了相对于载波中心频率的**粗频率偏移量**。
* 该峰值还表明在接收波形中检测到了三个 PSS（NID2）中的哪一个，以及最佳信道条件的时间瞬间。

* 通过将 SSB 中每个 OFDM 符号的**循环前缀**与 OFDM 符号的相应有用部分相关联，估算出低于半个子载波的频率偏移。这种相关性的相位与波形中的频率偏移成正比。

```matlab
disp(' -- Frequency correction and timing estimation --')

% 以千赫为单位指定频率偏移搜索带宽
searchBW = 6*scsSSB;
% SS 脉冲串波形的频率偏移校正
[rxWaveform,freqOffset,NID2] = hSSBurstFrequencyCorrect(rxWaveform,refBurst.BlockPattern,sampleRate,searchBW);
disp([' Frequency offset: ' num2str(freqOffset,'%.0f') ' Hz'])
```

**结果图：**

![image-20230720151253912](https://picture-cloud-storage-pyp.oss-cn-beijing.aliyuncs.com/img/202307201512942.png)

![image-20230720153430549](https://picture-cloud-storage-pyp.oss-cn-beijing.aliyuncs.com/img/202307201534569.png)

---

### 时间同步和 OFDM 解调

* 接收机利用频率搜索过程中检测到的参考 PSS 序列，估算出最强 SS 块的定时偏移。

* 经过频率偏移校正后，接收器可以认为参考 PSS 和接收波形的中心频率是一致的。

* 最后，接收器 OFDM 对同步波形进行解调，提取 SS 块。

```matlab
% 利用检测到的 PSS 创建用于时序估计的参考网格。
% PSS放在参考网格的第二个 OFDM 符号中，以避免第一个 OFDM 符号的特殊 CP 长度。
refGrid = zeros([nrbSSB*12 2]);
refGrid(nrPSSIndices,2) = nrPSS(NID2); % 正确 CP 长度的第二个 OFDM 符号
```

---

**nrPSS:**sym = nrPSS(ncellid) 返回物理层小区标识号 ncellid 的主同步信号（PSS）符号。该函数实现了 TS 38.211 第 7.4.2.2 节 [1]。

nrPSSIndices这个变量是PSS对应子载波的序号数，从下图中可以看出这一点。

![image-20230720164729852](https://picture-cloud-storage-pyp.oss-cn-beijing.aliyuncs.com/img/202307201647888.png)

这张图的子载波数是从0开始标记序号的，在matlab中是从1开始标记序号的，因此是从57 - 183。为了更好的观看这里将列向量转换成行向量进行观察。

![image-20230720164904289](https://picture-cloud-storage-pyp.oss-cn-beijing.aliyuncs.com/img/202307201649330.png)

最终变量refGrid中也仅有第二列的这些序号位置有值，如下图所示：

![image-20230720165005400](https://picture-cloud-storage-pyp.oss-cn-beijing.aliyuncs.com/img/202307201650425.png)

---

```matlab
% 时序估计。这是在检测到 SSB 之前的 OFDM 符号的时序偏移。
% 由于参考网格的内容而导致的检测到的 SSB 的时间偏移。
nSlot = 0;
timingOffset = nrTimingEstimate(rxWaveform,nrbSSB,scsSSB,nSlot,refGrid,'SampleRate',sampleRate);
```

---

[offset,mag] = nrTimingEstimate(carrier,waveform,refGrid)：

* 通过将输入波形与参考波形进行交叉相关来执行实用的定时估计。
* 函数通过使用正交频分复用（OFDM）调制参考资源网格 refGrid 来获取参考波形。
* 函数返回输入波形中每个接收天线的估计时序偏移（offset）和估计脉冲响应幅度（mag）。

[offset,mag] = nrTimingEstimate(waveform,nrb,scs,initialNSlot,refGrid)： 

* 通过调制参考资源网格 refGrid 并使用 OFDM 调制，以子载波间隔 scs 和初始时隙数 initialNSlot 跨 nrb 资源块，获得参考波形。

---

```matlab
% 同步、OFDM 解调和提取最强 SS 块
rxGrid = nrOFDMDemodulate(rxWaveform(1+timingOffset:end,:),nrbSSB,scsSSB,nSlot,'SampleRate',sampleRate);
rxGrid = rxGrid(:,2:5,:);
```

---

grid = nrOFDMDemodulate（carrier,waveform）通过解调波形（OFDM 调制波形）来恢复载波配置参数 carrier 的载波资源阵列。

grid = nrOFDMDemodulate(waveform,nrb,scs,initialNSlot) 为 nrb、指定的资源块数、子载波间隔 scs 和初始时隙数 initialNSlot 解调波形。

**rxGrid最开始为一个240 × 557的向量：**

![image-20230720175034350](https://picture-cloud-storage-pyp.oss-cn-beijing.aliyuncs.com/img/202307201750373.png)

**代码中仅获取了第2到5列：**

![image-20230720175134810](https://picture-cloud-storage-pyp.oss-cn-beijing.aliyuncs.com/img/202307201751828.png)

---

```matlab
% 以样本为单位显示定时偏移。
% 由于符号长度是以 FFT 样本测量的，因此要根据接收机的采样率。
srRatio = sampleRate/(scsSSB*1e3*rxOfdmInfo.Nfft);
firstSymbolLength = rxOfdmInfo.SymbolLengths(1)*srRatio;
str = sprintf(' Time offset to synchronization block: %%.0f samples (%%.%.0ff ms) \n',floor(log10(sampleRate))-3);
fprintf(str,timingOffset+firstSymbolLength,(timingOffset+firstSymbolLength)/sampleRate*1e3);
```

最后的目的是为了测量定时偏移：

* sampleRate是采样率，通过上述公式可以计算出srRatio，srRatio是一个OFDM符号在时域上是采样时间的几倍，换句话来讲就是一个OFDM符号需要的采样次数。
* firstSymbolLength是通过OFDM基本信息以及srRatio参数计算得出的第一个OFDM符号的点数。

![image-20230720175718479](https://picture-cloud-storage-pyp.oss-cn-beijing.aliyuncs.com/img/202307201757503.png)

* 总共的点数就是偏移量与第一个OFDM符号点数相加，这里是1644 + 556 = 2200.

![image-20230720175847581](https://picture-cloud-storage-pyp.oss-cn-beijing.aliyuncs.com/img/202307201758602.png)

* 采样一次需要1/sampleRate的时间，故总共产生 (timingOffset+firstSymbolLength)/sampleRate 时间的偏移，最后再乘以1e3换算成毫秒单位得到我们最终的结果

![image-20230720175922220](https://picture-cloud-storage-pyp.oss-cn-beijing.aliyuncs.com/img/202307201759240.png)

---

### SSS 搜索

* 接收器从接收到的网格中提取与 SSS 相关的资源元素，并将其与本地生成的每个可能的 SSS 序列相关联。

* 最强 PSS 序列和 SSS 序列的索引组合给出物理层小区标识，这是 PBCH DM-RS 和 PBCH 处理所必需的。
