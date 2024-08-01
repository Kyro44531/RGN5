function [timingOffset,nLeadingFrames] = hTimingOffsetToFirstFrame(offset,burst,ssbIdx,nHalfFrame,sampleRate)

    % As the symbol lengths are measured in FFT samples, scale the symbol
    % lengths to account for the receiver sample rate. Non-integer delays
    % are approximated at the end of the process.
    scs = hSSBurstSubcarrierSpacing(burst.BlockPattern);
    ofdmInfo = nrOFDMInfo(1,scs,'SampleRate',sampleRate); % smallest FFT size for SCS-SR
    srRatio = sampleRate/(scs*1e3*ofdmInfo.Nfft);
    symbolLengths = ofdmInfo.SymbolLengths*srRatio;

    % Adjust timing offset to the start of the SS block. This step removes
    % the extra offset introduced in the reference grid during PSS search,
    % which contained the PSS in the second OFDM symbol.
    offset = offset + symbolLengths(1);

    % Timing offset is adjusted so that the received grid starts at the
    % frame head i.e. adjust the timing offset for the difference between
    % the first symbol of the strongest SSB, and the start of the frame
    burstStartSymbols = hSSBurstStartSymbols(burst.BlockPattern,burst.L_max); % Start symbols in SSB numerology
    ssbFirstSym = burstStartSymbols(ssbIdx+1); % 0-based

    % Adjust for whole subframes
    symbolsPerSubframe = length(symbolLengths);
    subframeOffset = floor(ssbFirstSym/symbolsPerSubframe);
    samplesPerSubframe = sum(symbolLengths);
    timingOffset = offset - (subframeOffset*samplesPerSubframe);

    % Adjust for remaining OFDM symbols
    symbolOffset = mod(ssbFirstSym,symbolsPerSubframe);
    timingOffset = timingOffset - sum(symbolLengths(1:symbolOffset));

    % The first OFDM symbol of the SSB is defined with respect to the
    % half-frame where it is transmitted. Adjust for half-frame offset
    timingOffset = timingOffset - nHalfFrame*5*samplesPerSubframe;

    % Adjust offset to the first frame in the waveform that is scheduled
    % for SSB transmission.
    repetitions = ceil(timingOffset/(20*samplesPerSubframe));
    timingOffset = round(timingOffset - repetitions*20*samplesPerSubframe);

    % Calculate the number of leading frames before the detected one
    nLeadingFrames = 2*repetitions;

end