%hMCS Target code rate and modulation 
%   [TCR,MOD] = hMCS(IMCS) returns target code rate TCR and modulation MOD
%   corresponding to the modulation and coding scheme index IMCS according
%   to TS 38.214 Table 5.1.3.1-1.

%   Copyright 2020 The MathWorks, Inc.

function [tcr,modulation] = hMCS(imcs)
   
    % TS 38.214 Table 5.1.3.1-1: MCS index table 1 for PDSCH (qam64)
    table1 = [...
        0	2	120	0.2344
        1	2	157	0.3066
        2	2	193	0.3770
        3	2	251	0.4902
        4	2	308	0.6016
        5	2	379	0.7402
        6	2	449	0.8770
        7	2	526	1.0273
        8	2	602	1.1758
        9	2	679	1.3262
        10	4	340	1.3281
        11	4	378	1.4766
        12	4	434	1.6953
        13	4	490	1.9141
        14	4	553	2.1602
        15	4	616	2.4063
        16	4	658	2.5703
        17	6	438	2.5664
        18	6	466	2.7305
        19	6	517	3.0293
        20	6	567	3.3223
        21	6	616	3.6094
        22	6	666	3.9023
        23	6	719	4.2129
        24	6	772	4.5234
        25	6	822	4.8164
        26	6	873	5.1152
        27	6	910	5.3320
        28	6	948	5.5547];

    Qm = [2 4 6 8];
    modulation = {'QPSK' '16QAM' '64QAM' '256QAM'};
    
    row = table1(table1(:,1)==imcs,:);
    
    tcr = row(:,3) / 1024;
    modulation = modulation{Qm==row(:,2)};
    
end
