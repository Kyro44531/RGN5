clear all;
close all;
clc;
%---------------- �������ò���----------------%
Nsp=52;             %ϵͳ���ز�����������ֱ���ز���
Nfft=64;            % FFT ����
Ncp=16;             % ѭ��ǰ׺����
Ns=Nfft+Ncp;        % 1������OFDM���ų���
noc=53;             % ����ֱ���ز����ܵ����ز���
Nd=6;               % ÿ֡������OFDM������(������ѵ������)
M1=4;               % QPSK����
sr=250000;          % OFDM��������
SNR=20;         	% �����
ts=1/sr/Ns;                    % OFDM���ų���ʱ����
t=0:ts:(Ns*(Nd+1)-1)*ts;       % ����ʱ��
fd=100;                        % ��������Ƶ��

%----------------�����ŵ��Ĳ���----------------%
h=rayleigh(fd,t);                   % ���ɵ���Rayleigh˥���ŵ�
h1=sqrt(1/2)*h;                     % ��һ���Ĺ������ܹ��ʵ�1/2
h2=sqrt(1/3)*h;                     % �ڶ����Ĺ������ܹ��ʵ�1/3
h3=sqrt(1/6)*h;                     % �������Ĺ������ܹ��ʵ�1/6
h2=[zeros(1,4) h2(1:end-4)];        % �ڶ�����ʱ4����
h3=[zeros(1,8) h2(1:end-8)];        % ��������ʱ8����


%-----------------����ѵ������-----------------%
%-ѵ������Ƶ������,����802.11a�еĳ�ѵ����������-%
Preamble=[1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 ...
    1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1];
Preamble1=zeros(1,Nfft);                            % 1X64ȫ�����
Preamble1(2:27)=Preamble(27:end);                   % ��֯
Preamble1(39:end)=Preamble(1:26);
preamble1=ifft(Preamble1);                          % ѵ������ʱ������
preamble1=[preamble1(Nfft-Ncp+1:end) preamble1];    % ����16�����ѭ��ǰ׺

 %-------------------���������-------------------%
 msg1=randsrc(Nsp,Nd,[0:M1-1]);         % QPSK��Ϣ���ݣ�52X6�ľ���
 
 data1=qammod(msg1,M1)/sqrt(2);         % QPSK���Ʋ���һ��
 
 data2=zeros(Nfft,Nd);                  % ����FFTҪ�󣬶��������ţ�64X6�ľ���
 
 data2(2:27,:)=data1(27:end,:);
 data2(39:end,:)=data1(1:26,:);
 
 data2=ifft(data2);                          % IFFT�任
 
 data2=[data2(Nfft-Ncp+1:end,:);data2];      % ����ѭ��ǰ׺
 
 spow1=norm(data2,'fro').^2/(Nsp*Nd);        % �����������
 
 %������е��Ǽ��뵼Ƶ�Ĺ��������뵼Ƶ��ÿһ֡����7�����ţ�ÿ������������
 %������ѭ��ǰ׺����˺���80����
 data3=zeros(Ns,(Nd+1));                % ����ѵ�����ţ�80X7�ľ���

data3(:,1)=preamble1.'; % ��ÿһ֡�Ŀ�ͷ���뵼Ƶ����
data3(:,2:(Nd+1))=data2(:,1:Nd);
% �ڵ�Ƶ����֮�����6�����Ź���һ֡
 
 data3=reshape(data3,1,Ns*(Nd+1));      % �����任��1X560�ľ���
 
 data31=zeros(1,length(data3));
 data32=zeros(1,length(data3));
 data31(5:end)=data3(1:end-4);               % �ڶ������յ��ķ���
 data32(9:end)=data3(1:end-8);               % ���������յ��ķ���
 
 sigma1=sqrt(1/2*spow1/log2(M1)*10.^(-SNR/10)); % ����SNR����������׼��
 
 dd1=data3(1:Ns*(Nd+1)); %ȡ����һ��ÿһ֡������
 dd2=data31(1:Ns*(Nd+1));%ȡ���ڶ���ÿһ֡������
 dd3=data32(1:Ns*(Nd+1));%ȡ��������ÿһ֡������
 
 hh1=h1(1:Ns*(Nd+1));    % ��ǰ֡��3���ŵ�����
 hh2=h2(1:Ns*(Nd+1));
 hh3=h3(1:Ns*(Nd+1));
 
 % �ź�ͨ��3��˥���ŵ����������˹������
 r1=hh1.*dd1+hh2.*dd2+hh3.*dd3+sigma1*(randn(1,length(dd1))+1i*randn(1,length(dd1)));
 
 %-------------------���ջ�����-------------------%
 r1=reshape(r1,Ns,Nd+1);             % �����任(80X7�ľ���)
 
 r1=r1(Ncp+1:end,:);                 % �Ƴ�ѭ��ǰ׺
 
 R1=fft(r1);                         % fft����
 
 R1=[R1(39:end,:);R1(2:27,:)];       % �������ţ��⽻֯��
 
 HH1=(Preamble.')./R1(:,1);          % LS�ŵ�����
 
 HH1=HH1*ones(1,Nd);                 % ��һ֡�Ĺ����չ�ŵ����ƽ��
 
 x1=R1(:,2:end).*HH1;                % �򵥵��ŵ�����
 
 x1=qamdemod(x1.*sqrt(2),M1);        % ���ݽ��
 
 [neb1,ber1]=biterr(x1,msg1,log2(M1)); % �������
 
%�����ź� 
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
%Ƶ��ͼ
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
%���͡������źŶԱ�
figure(3)
subplot(2,1,1);
stem(msg1);
title('Transmitted Messages');
subplot(2,1,2);
stem(x1);
title('Received Messages');
xlabel(['BER=',num2str(ber1)]) 
