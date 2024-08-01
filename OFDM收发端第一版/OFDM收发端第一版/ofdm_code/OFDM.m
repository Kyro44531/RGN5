clear 
%% 参数设置
N_sc = 52;                            % 系统子载波数（不包括直流载波）、number of subcarrierA
N_fft = 64;                            % FFT 长度
N_cp =16;                             % 循环前缀长度、Cyclic prefix
N_symbo = N_fft + N_cp;    % 1个完整OFDM符号长度
N_c = 53;                              % 包含直流载波的总的子载波数、number of carriers
M = 4;                                   % 4PSK调制
SNR = 0:1:25;                        % 仿真信噪比
N_frm = 10;                           % 每种信噪比下的仿真帧数、frame
Nd = 6;                                  % 每帧包含的OFDM符号数
P_f_inter = 6;                         % 导频间隔
data_station = [];                   % 导频位置
L = 7;                                     % 卷积码约束长度
tblen = 6*L;                           % Viterbi译码器回溯深度
stage = 3;                              % m序列的阶数
ptap1 = [1 3];                        % m序列的寄存器连接方式
regi1 = [1 1 1];                      % m序列的寄存器初始值


%% 基带数据数据产生
P_data = randi([0 1],1,N_sc*Nd*N_frm);
% 基带数据总数 = 子载波数  * 每一帧下OFDM符号数 * 帧数

%% 信道编码（卷积码、或交织器）
% 卷积码：前向纠错非线性码
% 交织：使突发错误最大限度的分散化
trellis = poly2trellis(7,[133 171]);       % (2,1,7)卷积编码
% 个人认为并不是(2,1,7)卷积码，根据mathwork文文档实际上只有6个移位寄存器
% 下面是生成卷积编码的代码，卷积码可以采用维特比译码，具有一定的纠错能力，因此使用在本程序中较好。
code_data = convenc(P_data,trellis);
% 后续还可以采用交织代码，因为交织可以分散突发错误
%
%                             待后续补充...
% ----------------------------------------------------------

%% qpsk调制
data_temp1 = reshape(code_data,log2(M),[])';         % 以每组2比特进行分组，M=4
data_temp2 = bi2de(data_temp1);                            % 二进制转换为十进制数字，一共有0,1,2,3四种可能（因为M = 4）
modu_data = pskmod(data_temp2,M,pi/M);            % 调用matlab函数，实现QPSK调制
% 本质上就是将十进制数字代表的数据转换为复数，映射为二维坐标平面上的若干个对称点。这里最后一个参数pi/m = 2pi/M /2用来旋转，
% 使得星座点在图中呈对称状态。

% 绘图，观察星座点（其实没什么必要）
scatterplot(modu_data),grid;                  % 星座图(也可以取实部用plot函数)

%% 扩频
%————————————————————————————————————————————————————————%
% 扩频通信信号所占有的频带宽度远大于所传信息必需的最小带宽
% 根据香农定理，扩频通信就是用宽带传输技术来换取信噪比上的好处，这就是扩频通信的基本思想和理论依据。
% 扩频就是将一系列正交的码字与基带调制信号内积
% 扩频后数字频率变成了原来的m倍。码片数量 = 2（符号数）* m（扩频系数）
%————————————————————————————————————————————————————————%
code = mseq(stage,ptap1,regi1,N_sc);                                                % 生成扩频码，也是生成m序列
% m序列均为±1，因此这里需要将0和1转换为+1和-1
code = code * 2 - 1;                                                                             % 将1、0变换为1、-1
% 为每个子载波上生成对应数量的数据量
% 一共有N_sc行，因为一共有N_sc个子载波
modu_data = reshape(modu_data,N_sc,length(modu_data) / N_sc);
spread_data = spread(modu_data,code);                                            % 扩频
spread_data=reshape(spread_data,[],1);                                              % 这里进行并串变换
% 上述进行并串变换的原因是因为之后要进行IFFT / FFT操作
% 这些操作是基于2^n个数据进行的
% 因此要将子载波数（算上导频）扩展到2^n，这里选取的是64
% 所以要在串行数据后面进行填0操作，之后在进行串并转换

%% 插入导频 
P_f = 3+3*1i;                                                 % Pilot frequency
% 通过导频间隔确定导频位置
P_f_station = 1:P_f_inter:N_fft;                    % 导频位置
pilot_num = length(P_f_station);                   % 导频数量

% 默认第一个位置为导频，下面的步骤是将数据的索引值存储在data_station中
for img=1:N_fft                                             % 数据位置
    if mod(img,P_f_inter)~=1                          % mod(a,b)就是求的是a除以b的余数
        data_station=[data_station,img];
    end
end
data_row=length(data_station);                    % 行数为所需要使用到的数据数
% 列数上取整，这里为了让可传输容量大于等于数据量，否则会出问题
data_col=ceil(length(spread_data)/data_row);
pilot_seq=ones(pilot_num,data_col) * P_f;    % 将导频放入矩阵
data=zeros(N_fft,data_col);                           % 预设整个矩阵
data(P_f_station(1:end),:) = pilot_seq;           % 对pilot_seq按行取

% 为了进行快速傅里叶/逆变换，需要进行补零操作
if data_row*data_col >= length(spread_data)
    data2=[spread_data;zeros(data_row*data_col-length(spread_data),1)];% 将数据矩阵补齐，补0是虚载频~
end

%% 串并转换
data_seq=reshape(data2,data_row,data_col);
data(data_station(1:end),:)=data_seq;                       % 将导频与数据合并

%% IFFT
ifft_data=ifft(data); 

%% 插入保护间隔、循环前缀
Tx_cd=[ifft_data(N_fft-N_cp+1:end,:);ifft_data];         % 把ifft的末尾N_cp个数补充到最前面

%% 并串转换
Tx_data=reshape(Tx_cd,[],1);                                       % 由于传输需要

%% 信道（通过多经瑞利信道、或信号经过AWGN信道）
Ber=zeros(1,length(SNR));                                          % 译码后
Ber2=zeros(1,length(SNR));                                        % 译码前

for jj = 1:5:length(SNR)
    rx_channel = awgn(Tx_data,SNR(jj),'measured');    % 添加高斯白噪声
    
%% 串并转换
    Rx_data1 = reshape(rx_channel,N_fft+N_cp,[]);
    
%% 去掉保护间隔、循环前缀
    Rx_data2=Rx_data1(N_cp+1:end,:);

%% FFT
    fft_data=fft(Rx_data2);
    
%% 信道估计与插值（均衡）
    data3=fft_data(1:N_fft,:); 
    Rx_pilot=data3(P_f_station(1:end),:);                      % 接收到的导频
    h=Rx_pilot./pilot_seq;                                             % pilot_seq都是3 + 3i
    % 分段线性插值：插值点处函数值由连接其最邻近的两侧点的线性函数预测。对超出已知点集的插值点用指定插值方法计算函数值
    H=interp1( P_f_station(1:end)',h,data_station(1:end)','linear','extrap');

%% 信道校正
    data_aftereq=data3(data_station(1:end),:) ./ H;
%% 并串转换
% 首先将并行数据串行处理，即放入一个列向量中
    data_aftereq=reshape(data_aftereq,[],1);
% 进行去0操作，因为后面部分数据是填零处理
    data_aftereq=data_aftereq(1:length(spread_data));
% 将数据恢复成N_sc行的矩阵
    data_aftereq=reshape(data_aftereq,N_sc,length(data_aftereq)/N_sc);
    
%% 解扩
    demspread_data = despread(data_aftereq,code);       % 数据解扩
    figure
    hold on;
    xlabel('Real');
    ylabel('Imag');
    title(sprintf("信噪比为%d时的接收点星座图", jj));
    scatter(real(demspread_data),imag(demspread_data));
    
%% QPSK解调
% 解调
    demodulation_data=pskdemod(demspread_data,M,pi/M);    
% 转换成流向量
    De_data1 = reshape(demodulation_data,[],1);
% 十进制数据转换成二进制数据
    De_data2 = de2bi(De_data1);
% 转换成行向量
    De_Bit = reshape(De_data2',1,[]);

%% （解交织）
%% 信道译码（维特比译码）
    trellis = poly2trellis(7,[133 171]);
    rx_c_de = vitdec(De_Bit,trellis,tblen,'trunc','hard');   % 硬判决

%% 计算误码率
    [~,Ber2(jj)] = biterr(De_Bit(1:length(code_data)),code_data);% 译码前的误码率
    [err, Ber(jj)] = biterr(rx_c_de(1:length(P_data)),P_data);% 译码后的误码率
    
    figure
    subplot(2,1,1);
    x=0:1:30;
    stem(x,P_data(1:31));
    ylabel('amplitude');
    title('发送数据');
    legend('4PSK调制、卷积译码、有扩频');
    
    subplot(2,1,2);
    x=0:1:30;
    stem(x,rx_c_de(1:31));
    ylabel('amplitude');
    title(sprintf("信噪比为%d时的接收数据", jj));
    legend('4PSK调制、卷积译码、有扩频');
end
 figure
 semilogy(SNR,Ber2,'b-s');
 hold on;
 semilogy(SNR,Ber,'r-o');
 hold on;
 legend('4PSK调制、卷积码译码前（有扩频）','4PSK调制、卷积码译码后（有扩频）');
 hold on;
 xlabel('SNR');
 ylabel('BER');
 title('AWGN信道下误比特率曲线');
