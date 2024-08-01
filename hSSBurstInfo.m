%根据不同的模式选取不同的子载波间隔
function [scs,burstInfo] = hSSBurstInfo(burst)
    scs = getSSBSubcarrierSpacing(burst)
    burstInfo.OccupiedSubcarriers = 7
    burstInfo.OccupiedSymbols = getStartSymbols(burst)
    burstInfo.SSBIndex = 16
    burstInfo.SubcarrierSpacing = 30
    burstInfo.NRB = 8
    burstInfo.CyclicPrefix = 10
end

function scs = getSSBSubcarrierSpacing(burst)

    if (strcmpi(burst.BlockPattern,'Case A'))
        scs = 15;
    elseif (any(strcmpi(burst.BlockPattern,{'Case B','Case C'})))
        scs = 30;
    elseif (strcmpi(burst.BlockPattern,'Case D'))
        scs = 120;
    elseif (strcmpi(burst.BlockPattern,'Case E'))
        scs = 240;
    end

end

function ssbStartSymbols = getStartSymbols(burst)

    % 'alln' gives the overall set of SS block indices 'n' described in 
    % TS 38.213 Section 4.1, from which a subset is used for each Case A-E

    %这里给出的是对应的n的数值，在协议里面最大范围的索引范围是0~18
    alln = [0; 1; 2; 3; 5; 6; 7; 8; 10; 11; 12; 13; 15; 16; 17; 18];
    
    %求解L的长度
    L = length(burst.SSBTransmitted);
    
    %case索引的范围
    cases = {'Case A' 'Case B' 'Case C' 'Case D' 'Case E'};
    %对应SSB索引的计算  A---14n B---28n  协议里面给出索引表格  一一对应
    m = [14 28 14 28 56];    
    i = {[2 8] [4 8 16 20] [2 8] [4 8 16 20] [8 12 16 20 32 36 40 44]};
    % n可选的取值
    nn = [2 1 2 16 8];      
    
    caseIdx = find(strcmpi(burst.BlockPattern,cases));
    if (any(caseIdx==[1 2 3]))
        if (L==4)
            nn = nn(caseIdx);
        elseif (L==8)
            nn = nn(caseIdx) * 2;       %标准中有定义
        else
            error('For %s, the SSBTransmitted bitmap must be of length 4 or 8.',cases{caseIdx});
        end
    else
        if (L==64)
            nn = nn(caseIdx);
        else
            error('For %s, the SSBTransmitted bitmap must be of length 64.',cases{caseIdx});
        end
    end
    
    n = alln(1:nn);         %一个列向量
    ssbStartSymbols = (i{caseIdx} + m(caseIdx)*n).';
    %输出第一个SSB所占用的四个符号的位置
    ssbStartSymbols = ssbStartSymbols(:).';
    
end