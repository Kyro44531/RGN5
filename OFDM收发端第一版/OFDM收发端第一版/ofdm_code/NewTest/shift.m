function [outregi] = shift(inregi,shiftr)
%% 说明部分
% inrege     : 输入序列
% shiftr       : 循环右移的位数
% outregi    : 输出序列
%% 实现部分
% 首先获取输入序列的长度，这对后面的循环移位非常必要
v  = length(inregi);
% 创建一个输出序列矩阵，因为是循环移位，因此长度和输入序列是相同的
outregi = inregi;
%  rem函数是一个取余数的函数，如果循环移位长度大于序列长度，这样就相当于移位了（rem(shifter,v)次）
shiftr = rem(shiftr,v);

% 分别判断右移和作左移的情况。
% 右移的最右边一位移动到左边第一位，左移也是如此
if shiftr > 0
    outregi(:,1:shiftr) = inregi(:,v-shiftr+1:v);     %循环移位
    outregi(:,1+shiftr:v) = inregi(:,1:v-shiftr);
elseif shiftr < 0
    outregi(:,1:v+shiftr) = inregi(:,1-shiftr:v);
    outregi(:,v+shiftr+1:v) = inregi(:,1:-shiftr);
end
