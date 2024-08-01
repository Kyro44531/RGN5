pssIndices = nrPSSIndices;%PSS信号的位置    第一个OFDM符号，所在的子载波范围为[57,183]
pssGrid = zeros([240 4]);%用于存放时域4个OFDM符号，频域上4×240个子载波
refGrid = zeros([ofdmInfo.NSubcarriers ofdmInfo.SymbolsPerSlot]);%每个时隙固定14个OFDM符号
k = burstInfo.OccupiedSubcarriers;
l = burstInfo.OccupiedSymbols(1,:);
figure;
hold on;
peak_value = zeros(1,3);
peak_index = zeros(1,3);
for NID2 = [0 1 2]%NID2三种形式
    pssRef = nrPSS(NID2);%help nrPSS 查看功能
    pssGrid(pssIndices) = pssRef;%生成的PSS信号放入pssGrid
    refGrid(k,l) = pssGrid;%放入SSB中PSS对应的位置
    refWaveform = hOFDMModulate(ofdmConfig,refGrid);%做OFDM调制
    refWaveform = refWaveform(refWaveform~=0);   
    corr = zeros(size(rxWaveform));
    for r = 1:size(rxWaveform,2)%做自相关运算
        antcorr = xcorr(rxWaveform(:,r),refWaveform);
        corr(:,r) = antcorr(size(rxWaveform,1):end);
    end
    corr = sum(abs(corr),2);%相加之后得到自相关运算的结果（这里可以了解一下matlab做自相关的原理）
    [peak_value(NID2+1),peak_index(NID2+1)] = max(corr)%找出三个本地的PSS信号分别与接收的信号做自相关的峰值
    plot(corr);
end
% Plot PSS correlations
axis([1 length(rxWaveform(:,1)) 0 max(peak_value)*1.1]);
title('PSS Correlations (time domain)');
ylabel('Magnitude');
xlabel('Sample Index');
% Determine NID2 by finding the strongest correlation
NID2 = find(peak_value==max(peak_value)) - 1;%找出NID2
%至此已经找出NID2。
% Determine timing offset（找出最大值处的peak_index）
offset = peak_index(NID2+1) - 1;
% Plot selected NID2
plot(offset+1,peak_value(NID2+1),'kx','LineWidth',2,'MarkerSize',8);
lgd = legend;
lgd.Interpreter = 'latex';
legends = "$N_{ID}^{(2)}$ = " + num2cell(0:2);
legend([legends "$N_{ID}^{(2)}$ = " + num2str(NID2)]);
% Extract strongest burst
offset = offset - ofdmInfo.SymbolLengths(1);%SymbolLengths=info.CyclicPrefixLengths + info.Nfft
rxGrid = hOFDMDemodulate(ofdmConfig,rxWaveform(1+offset:end,:));
rxGrid = rxGrid(burstInfo.OccupiedSubcarriers,2:5,:);%解调出PSS信号。

% Extract the received SSS symbols from the SS/PBCH block
sssIndices = nrSSSIndices;%SSS信号在SSB中的索引
sssRx = nrExtractResources(sssIndices,rxGrid);%找出SSS信号  127×2
% Correlate received SSS symbols with each possible SSS sequence（检测SSS信号过程）
sssEst = zeros(1,336);
for NID1 = 0:335
    ncellid = (3*NID1) + NID2;
    sssRef = nrSSS(ncellid);
    sssEst(NID1+1) = sum(abs(mean(sssRx .* conj(sssRef),1)).^2);%在实际通信系统中也是这样检测的，336种NID2
end
% Plot SSS correlations
figure;
stem(0:335,sssEst,'o');
title('SSS Correlations (frequency domain)');
xlabel('$N_{ID}^{(1)}$','Interpreter','latex');
ylabel('Magnitude');
axis([-1 336 0 max(sssEst)*1.1]);
% Determine NID1 by finding the strongest correlation
NID1 = find(sssEst==max(sssEst)) - 1;
% Plot selected NID1
hold on;
plot(NID1,max(sssEst),'kx','LineWidth',2,'MarkerSize',8);
lgd = legend;
lgd.Interpreter = 'latex';
legend(["correlations" "$N_{ID}^{(1)}$ = " + num2str(NID1)]);
% Form overall cell identity from NID1 and NID2
ncellid = (3*NID1) + NID2;

% Extract the received PBCH DM-RS symbols from the SS/PBCH block
dmrsIndices = nrPBCHDMRSIndices(ncellid);%根据得到的小区ID来获取dmrsIndices 
%大家在单步调试查看dmrsIndices时，可以看到初始位置243（对应标准中的0+v, 4+v, ...）
[dmrsRx,dmrsRxIndices] = nrExtractResources(dmrsIndices,rxGrid);%提取出DM-RS信号和索引

% Correlate received DM-RS symbols with each possible DM-RS sequence
dmrsEst = zeros(1,8);
for ibar_SSB = 0:7
    dmrsRef = nrPBCHDMRS(ncellid,ibar_SSB);%nrPBCHDMRS生成标准中定义的所有DM-RS信号
    dmrsEst(ibar_SSB+1) = sum(abs(mean(dmrsRx .* conj(dmrsRef),1)).^2);%做相关检测的过程
    
end

% Plot PBCH DM-RS correlations
figure;
stem(0:7,dmrsEst,'o');
title('PBCH DM-RS Correlations (frequency domain)');
xlabel('$\overline{i}_{SSB}$','Interpreter','latex');
xticks(0:7);
ylabel('Magnitude');
axis([-1 8 0 max(dmrsEst)*1.1]);

% Record ibar_SSB for the strongest correlation
ibar_SSB = find(dmrsEst==max(dmrsEst)) - 1;%找出相关值最大的横坐标（索引）

% Plot selected ibar_SSB
hold on;
plot(ibar_SSB,max(dmrsEst),'kx','LineWidth',2,'MarkerSize',8);
lgd = legend;
lgd.Interpreter = 'latex';
legend(["correlations" "$\overline{i}_{SSB}$ = " + num2str(ibar_SSB)]);

%% Channel Estimation using PBCH DM-RS
% Now that the PBCH DM-RS sequence is known, a channel estimate for the
% SS/PBCH block can be created by estimating the channel in each PBCH DM-RS
% resource element location and interpolating across the SS/PBCH block. An
% estimate of the additive noise on the PBCH DM-RS is also performed.

% Channel estimation, using linear interpolation of the PBCH DM-RS
dmrsRef = nrPBCHDMRS(ncellid,ibar_SSB);
dmrsSubs = double(nrPBCHDMRSIndices(ncellid,'IndexStyle','subscript'));
hest = zeros([240 4 nrxants 1]);
[l_hest,k_hest] = meshgrid(1:size(hest,2),1:size(hest,1));
dmrsEsts = dmrsRx .* conj(dmrsRef);
for r = 1:nrxants
    f = scatteredInterpolant(dmrsSubs(:,2),dmrsSubs(:,1),dmrsEsts(:,r));
    hest(:,:,r) = f(l_hest,k_hest);
end

% Noise estimation, based on the difference between the PBCH DM-RS in
% symbols 2 and 4 (1-based) of the SS/PBCH block. This technique assumes
% that the channel itself does not change between the two symbols
dmrsEstsSym2 = dmrsEsts(dmrsSubs(:,2)==2,:);
dmrsEstsSym4 = dmrsEsts(dmrsSubs(:,2)==4,:);
dmrsEstsSym2and4 = [dmrsEstsSym2(:) dmrsEstsSym4(:)];
dmrsNoise = mean(dmrsEstsSym2and4,2) - dmrsEstsSym2and4;
nest = var(dmrsNoise(:)) * 2;

% Extract the received PBCH symbols from the SS/PBCH block
[pbchIndices,pbchIndicesInfo] = nrPBCHIndices(ncellid);%根据得到的小区ID来获取dmrsIndices pbchIndices&pbchIndicesInfo
pbchRx = nrExtractResources(pbchIndices,rxGrid);%根据pbchIndices，在接收的时频资源（一堆复数信号）中提取出PBCH信号

% Plot received PBCH constellation before equalization
figure;
plot(pbchRx,'o');
title('Received PBCH constellation');
m = max(abs([real(pbchRx(:)); imag(pbchRx(:))])) * 1.1;%取实部和虚部来画出星座图
axis([-m m -m m]);

% Configure 'v' for PBCH scrambling according to TS 38.211 Section 7.3.3.1
% 'v' is also the 2 LSBs of the SS/PBCH block index for L=4, or the 3 LSBs
% for L=8 or 64
%这部分时加扰PBCH加扰操作，解调过程中需要加扰相关信息，把标准中定义放在代码后。
if (burstInfo.L==4)
    v = mod(ibar_SSB,4);
else
    v = ibar_SSB;
end
ssbIndex = v;

% PBCH equalization and CSI calculation
pbchHest = nrExtractResources(pbchIndices,hest);
[pbchEq,csi] = nrEqualizeMMSE(pbchRx,pbchHest,nest);%这两步和上面相同
Qm = pbchIndicesInfo.G / pbchIndicesInfo.Gd;
csi = repmat(csi.',Qm,1);
csi = reshape(csi,[],1);

% Plot received PBCH constellation after equalization
%显示均衡之后的星座图
figure;
plot(pbchEq,'o');
title('Equalized PBCH constellation');
m = max(abs([real(pbchEq(:)); imag(pbchEq(:))])) * 1.1;
axis([-m m -m m]);

% PBCH demodulation
pbchBits = nrPBCHDecode(pbchEq,ncellid,v,nest);

% Calculate RMS PBCH EVM
pbchRef = nrPBCH(pbchBits<0,ncellid,v);
evm = comm.EVM;
evm_rms = evm(pbchEq,pbchRef);

% Display calculated EVM
disp(['RMS PBCH EVM = ' num2str(evm_rms,'%0.3f') '%']);

%% BCH Decoding
% The BCH bit estimates are weighted with Channel State Information (CSI)
% from the MMSE equalizer then BCH decoding is performed, consisting of
% rate recovery, polar decoding, CRC decoding, descrambling and separating
% the 24 BCH transport block bits from the 8 additional timing-related
% payload bits.

% Apply CSI
pbchBits = pbchBits .* csi;

% Perform BCH decoding including rate recovery, polar decoding and CRC
% decoding. PBCH descrambling and separation of the BCH transport block
% bits 'trblk' from 8 additional payload bits A...A+7 is also performed:
%   A ... A+3: 4 LSBs of System Frame Number
%         A+4: half frame number
% A+5 ... A+7: for L=64, 3 MSBs of the SS/PBCH block index
%              for L=4 or 8, A+5 is the MSB of the subcarrier offset k_SSB
polarListLength = 8;
[~,err,trblk,sfn4lsb,nHalfFrame,msbidxoffset] = ...
    nrBCHDecode(pbchBits,polarListLength,burstInfo.L,ncellid);%可以成功解码BCH

% Use 'msbidxoffset' value to set bits of 'k_SSB' or 'ssbIndex', depending
% on the number of SS/PBCH blocks in the burst
if (burstInfo.L==64)
    ssbIndex = ssbIndex + (bi2de(msbidxoffset.','left-msb') * 8);
    k_SSB = 0;
else
    k_SSB = msbidxoffset * 16;
end

% Display the BCH CRC
disp(['BCH CRC = ' num2str(err)]);

% Displaying the SSB index
disp(['SSB index = ' num2str(ssbIndex)]);

%% MIB Parsing
% Finally the 24 decoded BCH transport block bits are parsed into a
% structure which represents the MIB message fields. This includes
% reconstituting the 10-bit System Frame Number (SFN) |NFrame| from the 6
% MSBs in the MIB and the 4 LSBs in the PBCH payload bits. It also includes
% incorporating the MSB of the subcarrier offset |k_SSB| from the PBCH
% payload bits in the case of L=4 or 8 SS/PBCH blocks per burst.

% Create set of subcarrier spacings signaled by the 7th bit of the decoded
% MIB, the set is different for FR1 (L=4 or 8) and FR2 (L=64)
if (burstInfo.L==64)
    commonSCSs = [60 120];
else
    commonSCSs = [15 30];
end

% Create a structure of MIB fields from the decoded MIB bits. The BCH
% transport block 'trblk' is the RRC message BCCH-BCH-Message, consisting
% of a leading 0 bit then 23 bits corresponding to the MIB
mib.NFrame = bi2de([trblk(2:7); sfn4lsb] .','left-msb');
mib.SubcarrierSpacingCommon = commonSCSs(trblk(8) + 1);
mib.k_SSB = k_SSB + bi2de(trblk(9:12).','left-msb');
mib.DMRSTypeAPosition = 2 + trblk(13);
mib.PDCCHConfigSIB1 = bi2de(trblk(14:21).','left-msb');
mib.CellBarred = trblk(22);
mib.IntraFreqReselection = trblk(23);

% Display the MIB structure
disp(mib);