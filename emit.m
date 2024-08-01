burst.BlockPattern = 'Case B';%选择模式B传输
burst.SSBPeriodicity = 20;%SSB实际发送周期可以为5ms、10ms或20ms，但不能大于20ms
burst.NFrame = 4;%5G NR的系统帧号（无线帧号）{0-1023}，默认为0。
burst.SSBTransmitted = [1 1 1 1 1 1 1 1];%SSB突发集中SSB最大的个数，上文提到的   ，此处SSB突发集中的个数为8。
burst.NCellID = 102;%小区ID，设置为102

gnb.SubcarrierSpacing = 30;%子载波间隔30KHz，我上传的代码这里写成15KHz，修改一下。
gnb.NRB = 24;%在传输带宽为10MHz时，可以传输24个RB，一个RB为12个子载波，其余带宽作为保护间隔。我上传的代码这里写成52，修改一下。
gnb.CyclicPrefix = 'Normal';%设置OFDM循环前缀的长度，15KHz为标准长度

burst.SampleRate = 15.36;%采样率 15.36M，这里其实如果各种琢磨，会比较难理解。根据上面的参数举例，10M/30KHz=333.33，所以采用512点FFT，采样率512×30KHz=15.36M。涉及到OFDM的相关知识，在下面放一个链接，讲解原理。
K = 24*12;%子载波数目=NRB * 12
burst.FrequencyPointA = -K/2 * gnb.SubcarrierSpacing * 1e3;%5G NR新定义的频域参考点。

burst.DMRSTypeAPosition = 2;%标准中定义的dmrs-TypeA-Position，DM-RS第一个符号的位置（下行同步机制中有讲解DM-RS在一个SSB中的位置）
burst.PDCCHConfigSIB1 = 17;%参考38.213 - section13 - UE procedure for monitoring Type0-PDCCH CSS sets
burst.CellBarred = 0;%广播MIB的两种方式
burst.IntraFreqReselection = 0;%内部频率选择问题
burst.SubcarrierSpacingCommon = 30;%子载波间隔30kHz
burst.DisplayBurst = true;%显示SSB

ntxants = 8;
nrxants = 2;

velocity = 30.0;
fc = 4e9;
c = physconst('lightspeed');
fd = (velocity*1000/3600)/c*fc;%最大多普勒频移
channel = nrTDLChannel;%matlab输入 help nrTDLChannel 查看使用方法
channel.Seed = 24;
channel.DelayProfile = 'TDL-C';
channel.DelaySpread = 300e-9;
channel.MaximumDopplerShift = fd;
channel.MIMOCorrelation = 'Medium';
channel.Polarization = 'Cross-Polar';
channel.NumTransmitAntennas = ntxants;
channel.NumReceiveAntennas = nrxants;

SNRdB = 10;

disp(burst);

[~,burstInfo] = hSSBurstInfo(burst);
disp(burstInfo);

burstGrid = hSSBurst(burst)

W = fft(eye(ntxants)) / sqrt(ntxants);

beamformedGrid = zeros([size(burstGrid) ntxants]);
blockSubcarriers = burstInfo.OccupiedSubcarriers;
for ssb = 1:length(burstInfo.SSBIndex)
    blockSymbols = burstInfo.OccupiedSymbols(ssb,:);
    block = burstGrid(blockSubcarriers,blockSymbols);
    Wssb = W(mod(ssb-1,ntxants)+1,:);
    beamformedBlock = reshape(block(:) * Wssb,[size(block) ntxants]);
    beamformedGrid(blockSubcarriers,blockSymbols,:) = beamformedBlock;
end

beamformedGrid = beamformedGrid(:,1:max(burstInfo.OccupiedSymbols(:))+1,:);
ofdmConfig.SubcarrierSpacing = burstInfo.SubcarrierSpacing;%子载波间隔
ofdmConfig.NRB = burstInfo.NRB;%RB个数
ofdmConfig.CyclicPrefix = burstInfo.CyclicPrefix;%循环前缀的长度
ofdmConfig.Windowing = 0;
[burstWaveform,ofdmInfo] = hOFDMModulate(ofdmConfig,beamformedGrid);%进行OFDM调制

rxWaveform = channel(burstWaveform);

rng('default');%产生伪随机数
% rxWaveform = awgn(rxWaveform,SNRdB,-10*log10(double(ofdmInfo.Nfft)));