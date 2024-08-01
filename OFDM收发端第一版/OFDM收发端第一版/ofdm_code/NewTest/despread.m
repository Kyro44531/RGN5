function out = despread(data, code)
%% 说明部分
%   data   : 输入数据序列
%   code   : 解扩使用的扩频码序列
%   out     : 解扩后的输出数据序列
%% 实现部分
switch nargin                           %如果输入参数个数不对，提示错误
case { 0 , 1 }
    error('缺少输入参数');
end
%%
% 分别获取数据和扩频码的行数与列数
[hn,vn] = size(data);
[hc,vc] = size(code);                  

% 在正常情况下，hn和hc的数量应该是相等的
% 需要恢复成原始数据的数目，故需要vn / vc
out    = zeros(hc,vn/vc);                  

for ii=1:hc
    % 将data数据的其中某一行拿出
    % 并将这一行转换为vc 行 vn/vc列的矩阵
    xx=reshape(data(ii,:),vc,vn/vc);
    % 解扩操作
    out(ii,:)= code(ii,:)*xx/vc;
end

