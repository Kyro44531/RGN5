function [mout] = mseq(n, taps, inidata, num)
%% 说明部分
% n            : m序列的阶数n
% taps       : 反馈寄存器的连接位置（也可以理解为抽头位置）
% inidata   : 寄存器的初始值序列（初始值为多少都可以）
% num       : 输出的m序列的个数（和子载波数相关）
% mout      : 输出的m序列，如果num>1,则每一行为一个m序列
%% 实现部分
% 规范输出序列的大小，一共有num行，即一共需要输出num组数据
% 每一行共有2^n - 1列，即m序列的长度。我们需要获取完整的m序列才可以正常实现系统功能
mout = zeros(num,2^n-1);

% 除反馈端（图中最左端） 一共fpos个抽头位置
fpos = zeros(n,1);
% 在规定位置（m序列生成电路中有抽头的位置）规定为1，即这个地方参与循环移位的反馈。
fpos(taps) = 1;

% 进行循环位移操作生成m序列
% 由于对于相同长度，相同电路的m序列只有一个，故进行一次循环就可以
for ii=1:2^n-1
    mout(1,ii) = inidata(n);                               % 寄存器的输出值
    temp = mod(inidata*fpos,2);                     % 计算反馈数据，inidata * fpos就是反馈值，这里使用行向量乘以列向量，故相当于直接将各抽头值相加！
    inidata(2:n) = inidata(1:n-1);                      % 每一次输出一位结果后进行||寄存器移位一次||
    inidata(1) = temp;                                      % 更新第1个寄存器的值，为计算出的反馈数据  
end

% 如果要输出多个m序列，生成其他m序列，这里使用到了第一行的m序列进行循环移位的操作
% 下面的shift为循环位移函数
if num > 1                                          
    for ii=2:num
        mout(ii,:) = shift(mout(ii-1,:),1);
    end
end
