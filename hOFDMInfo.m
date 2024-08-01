% Generate the CP and frame numerology for the 1ms subframe
function info = nrNumerology(info,deltaf,slotduration)
end

% lteOFDMInfo functionality (LTE baseline numerology with any NRB)
function info = lteNumerology(enb)

    % Valid the parameters
    if ~any(strcmpi(enb.CyclicPrefix,{'Extended', 'Normal'}))
        error('The cyclic prefix type must be one Normal or Extended');
    end
  
    % Get FFT size (allowing NRB > 110)
    nFFT = power(2,ceil(log2(enb.NRB*12/0.85)));
    nFFT = max(128,nFFT);
    
    info.SamplingRate = nFFT * 15e3;
    info.Nfft = nFFT;

    ecp = strcmpi(enb.CyclicPrefix,'Extended');
     
    % The number of window samples is chosen in accordance with the maximum
    % values implied by TS 36.101, Tables F.5.3-1, and F.5.4-1.
    if (isfield(enb,'Windowing'))
        w = enb.Windowing;
    else
        w = 0;
        if (ecp)
            switch (nFFT)
                case 128
                    w = 4;
                case 256
                    w = 6;
                case 512
                    w = 4;
                case 1024
                    w = 6;
                case 2048
                    w = 8;
            end
        else
            switch (nFFT)
                case 128
                    w = 4;
                case 256
                    w = 6;
                case 512
                    w = 4;
                case 1024
                    w = 6;
                case 2048
                    w = 8;
            end
        end
        if w==0
            % Additional rule for other FFT sizes
            w = max(0,8-2*(11-(log2(nFFT))));  
        end
    end    
    info.Windowing = w;
      
    % CP lengths for 2048 point FFT per LTE slot (6 or 7 symbols)
    if ecp
        cpLengths = [512 512 512 512 512 512];
    else
        cpLengths = [160 144 144 144 144 144 144];
    end
    
    % Scale according to the FFT size and repeat the slots for a LTE subframe
    info.CyclicPrefixLengths = repmat(cpLengths*double(nFFT)/2048,1,2);
    
end