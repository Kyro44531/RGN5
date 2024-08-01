%根据不同的模式选取不同的子载波间隔
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
