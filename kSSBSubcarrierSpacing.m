function scsKSSB = kSSBSubcarrierSpacing(scsCommon)
% Subcarrier spacing of k_SSB, as defined in TS 38.211 Section 7.4.3.1

    if scsCommon > 30  % FR2
        scsKSSB = scsCommon;
    else
        scsKSSB = 15;
    end

end