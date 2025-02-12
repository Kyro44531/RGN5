function [out] = spread(data, code)
%% 说明部分
%   data    : 输入数据序列
%   code   : 扩频码序列
%   out      : 扩频后的输出数据序列
%% 实现部分
switch nargin
case { 0 , 1 }                                  %如果输入参数个数不对，提示错误
    error('缺少输入参数');
end
%%
% 分别获取数据和扩频码的行数与列数
[hn,vn] = size(data);
[hc,vc] = size(code);

% 确保子载波数对应确定量的扩频码
if hn > hc                                      %如果扩频码数小于输入的待扩频的数据序列，提示错误
    error('缺少扩频码序列');
end

% 在正常情况下，hn和hc的数量应该是相等的
% 扩频码在一整个周期中才能起到作用，因此我们需要实现整周期分别相乘
% 因此这里一共使用了vn * vc 列数据
out = zeros(hn,vn*vc);

% 输出数据
% code(ii,:).'*data(ii,:)为矩阵，这个矩阵一共有vc行，即m序列的个数
% 这个矩阵每一列的数据为对应的m序列位的值与第（这列）个数据值相乘。
% 举个例子，假如m序列位1 -1 1 1 -1，数据为2，则这个矩阵为[2 -2 2 2 -2]'（这里的m序列可能不太恰当）
% 对于每一行（每一个hn），将上述矩阵转换为串行数据（这里也可以理解为并串转换）
% 最终还是输出矩阵的形式，每一个子载波的数据是经过扩频后的！
for ii=1:hn
    out(ii,:) = reshape(code(ii,:).'*data(ii,:),1,vn*vc);
end
