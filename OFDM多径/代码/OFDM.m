clear all;
close all;
clc;
%---------------- 参数设置部分----------------%
Nsp=52;             %系统子载波数（不包括直流载波）
Nfft=64;            % FFT 长度
Ncp=16;             % 循环前缀长度
Ns=Nfft+Ncp;        % 1个完整OFDM符号长度
noc=53;             % 包含直流载波的总的子载波数
Nd=6;               % 每帧包含的OFDM符号数(不包括训练符号)
M1=4;               % QPSK调制
sr=250000;          % OFDM符号速率
SNR=20;         	% 信噪比
ts=1/sr/Ns;                    % OFDM符号抽样时间间隔
t=0:ts:(Ns*(Nd+1)-1)*ts;       % 抽样时刻
fd=100;                        % 最大多普勒频移

%----------------三径信道的参数----------------%
h=rayleigh(fd,t);                   % 生成单径Rayleigh衰落信道
h1=sqrt(1/2)*h;                     % 第一径的功率是总功率的1/2
h2=sqrt(1/3)*h;                     % 第二径的功率是总功率的1/3
h3=sqrt(1/6)*h;                     % 第三径的功率是总功率的1/6
h2=[zeros(1,4) h2(1:end-4)];        % 第二径延时4个点
h3=[zeros(1,8) h2(1:end-8)];        % 第三径延时8个点


%-----------------产生训练序列-----------------%
%-训练符号频域数据,采用802.11a中的长训练符号数据-%
Preamble=[1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 ...
    1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1];
Preamble1=zeros(1,Nfft);                            % 1X64全零矩阵
Preamble1(2:27)=Preamble(27:end);                   % 交织
Preamble1(39:end)=Preamble(1:26);
preamble1=ifft(Preamble1);                          % 训练符号时域数据
preamble1=[preamble1(Nfft-Ncp+1:end) preamble1];    % 加入16个点的循环前缀

 %-------------------发射机部分-------------------%
 msg1=randsrc(Nsp,Nd,[0:M1-1]);         % QPSK信息数据（52X6的矩阵）
 
 data1=qammod(msg1,M1)/sqrt(2);         % QPSK调制并归一化
 
 data2=zeros(Nfft,Nd);                  % 根据FFT要求，对数据重排（64X6的矩阵）
 
 data2(2:27,:)=data1(27:end,:);
 data2(39:end,:)=data1(1:26,:);
 
 data2=ifft(data2);                          % IFFT变换
 
 data2=[data2(Nfft-Ncp+1:end,:);data2];      % 加入循环前缀
 
 spow1=norm(data2,'fro').^2/(Nsp*Nd);        % 计算符号能量
 
 %下面进行的是加入导频的工作，加入导频后，每一帧含有7个符号，每个符号中由于
 %加入了循环前缀，因此含有80个点
 data3=zeros(Ns,(Nd+1));                % 加入训练符号（80X7的矩阵）

data3(:,1)=preamble1.'; % 在每一帧的开头加入导频序列
data3(:,2:(Nd+1))=data2(:,1:Nd);
% 在导频序列之后加入6个符号构成一帧
 
 data3=reshape(data3,1,Ns*(Nd+1));      % 并串变换（1X560的矩阵）
 
 data31=zeros(1,length(data3));
 data32=zeros(1,length(data3));
 data31(5:end)=data3(1:end-4);               % 第二径接收到的符号
 data32(9:end)=data3(1:end-8);               % 第三径接收到的符号
 
 sigma1=sqrt(1/2*spow1/log2(M1)*10.^(-SNR/10)); % 根据SNR计算噪声标准差
 
 dd1=data3(1:Ns*(Nd+1)); %取出第一径每一帧的数据
 dd2=data31(1:Ns*(Nd+1));%取出第二径每一帧的数据
 dd3=data32(1:Ns*(Nd+1));%取出第三径每一帧的数据
 
 hh1=h1(1:Ns*(Nd+1));    % 当前帧的3径信道参数
 hh2=h2(1:Ns*(Nd+1));
 hh3=h3(1:Ns*(Nd+1));
 
 % 信号通过3径衰落信道，并加入高斯白噪声
 r1=hh1.*dd1+hh2.*dd2+hh3.*dd3+sigma1*(randn(1,length(dd1))+1i*randn(1,length(dd1)));
 
 %-------------------接收机部分-------------------%
 r1=reshape(r1,Ns,Nd+1);             % 串并变换(80X7的矩阵)
 
 r1=r1(Ncp+1:end,:);                 % 移除循环前缀
 
 R1=fft(r1);                         % fft运算
 
 R1=[R1(39:end,:);R1(2:27,:)];       % 数据重排（解交织）
 
 HH1=(Preamble.')./R1(:,1);          % LS信道估计
 
 HH1=HH1*ones(1,Nd);                 % 按一帧的规格拓展信道估计结果
 
 x1=R1(:,2:end).*HH1;                % 简单的信道均衡
 
 x1=qamdemod(x1.*sqrt(2),M1);        % 数据解调
 
 [neb1,ber1]=biterr(x1,msg1,log2(M1)); % 误比特率
 
%调制信号 
figure(1);
subplot(2,1,1);
plot(t,real(data3(1,:)));
xlabel('time');
ylabel('Ampulitude');
title('Re_ signal');
subplot(2,1,2);
plot(t,imag(data3(1,:)));
xlabel('time');
ylabel('Ampulitude');
title('Im_ signal');
%频谱图
freq_signal=fft(data3(1,:));
sig_signal=fftshift(freq_signal);
n=length(freq_signal);
fs=4*10^6;
f=(0:n-1)*fs/n;
psd_signal=10*log10(abs(sig_signal).^2/max(abs(sig_signal).^2));
figure(2);
plot(f,psd_signal);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('OFDM signal spectrum');
%发送、接收信号对比
figure(3)
subplot(2,1,1);
stem(msg1);
title('Transmitted Messages');
subplot(2,1,2);
stem(x1);
title('Received Messages');
xlabel(['BER=',num2str(ber1)]) 
