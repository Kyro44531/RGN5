function plotResourceGrid(rxGrid,refBurst,systemInfo,nLeadingFrames,ssbIndex,nHalfFrame)

    % Extract SSB and common SCS from reference SS burst and initial system
    % information
    scsSSB = hSSBurstSubcarrierSpacing(refBurst.BlockPattern);
    scsCommon = systemInfo.SubcarrierSpacingCommon;
    scsRatio = scsSSB/scsCommon;

    % Number of subcarriers, symbols and frames.
    [K,L] = size(rxGrid);
    symbolsPerSubframe = 14*scsCommon/15;
    numFrames = ceil(L/(10*symbolsPerSubframe));

    % Define colors and auxiliary plotting function
    basePlotProps = {'LineStyle','-','LineWidth',1};
    occasionColor = 0.7*[1 1 1];
    detectionColor = [200,0,0]/255;

    frameBoundaryColor = 0.1*[1 1 1];
    boundingBox = @(y,x,h,w,varargin)rectangle('Position',[x+0.5 y-0.5 w h],basePlotProps{:},varargin{:});

    % Create figure and display resource grid
    figure;
    imagesc(abs(rxGrid(:,:,1))); axis xy; hold on;

    % Add vertical frame lines
    x = repmat((0:numFrames-1)*10*symbolsPerSubframe,3,1);
    x(3,:) = NaN;
    y = repmat([0;K;NaN],1,numFrames);
    plot(x(:),y(:),'Color',frameBoundaryColor);

    % Determine frequency origin of the SSB in common numerology
    ssbCenter = K/2;
    halfSSB = 10*12*scsRatio;
    scsKSSB = kSSBSubcarrierSpacing(scsCommon);
    kSSBFreqOff = systemInfo.k_SSB*scsKSSB/scsCommon;
    ssbFreqOrig = ssbCenter - halfSSB + kSSBFreqOff + 1;

    % Determine time origin of the SSB in common numerology
    ssbStartSymbols = hSSBurstStartSymbols(refBurst.BlockPattern,refBurst.L_max);
    ssbStartSymbols = ssbStartSymbols + 5*symbolsPerSubframe*nHalfFrame;
    ssbHeadSymbol = ssbStartSymbols(ssbIndex+1)/scsRatio;
    ssbTailSymbol = floor((ssbStartSymbols(ssbIndex+1)+4)/scsRatio)-1;

    % Draw bounding boxes around all SS/PBCH block occasions
    w = ssbTailSymbol - ssbHeadSymbol + 1;
    for i = 1:ceil(numFrames/2)
        s = ssbHeadSymbol + (i-1)*2*10*symbolsPerSubframe + 5*symbolsPerSubframe*nHalfFrame;
        if s <= (L - w)
            boundingBox(ssbFreqOrig,s,240*scsRatio,w,'EdgeColor',occasionColor);
        end
    end

    % Draw bounding box for detected SS/PBCH block
    s = ssbHeadSymbol + nLeadingFrames*10*symbolsPerSubframe + 5*symbolsPerSubframe*nHalfFrame;
    boundingBox(ssbFreqOrig,s,240*scsRatio,w,basePlotProps{:},'EdgeColor',detectionColor)

    % Add text next to detected SS/PBCH block
    str = sprintf('SSB#%d',ssbIndex);
    text(s+w+1,ssbFreqOrig+24,0,str,'FontSize',10,'Color','w')

    % Create legend. Since rectangles don't show up in legend, create a
    % placeholder for bounding boxes.
    plot(NaN,NaN,basePlotProps{:},'Color',occasionColor);
    plot(NaN,NaN,basePlotProps{:},'Color',detectionColor);
    legend('Frame boundary','Occasion','Detected')
    xlabel('OFDM symbol'); ylabel('Subcarrier');

    % Add title including frame numbers
    firstNFrame = systemInfo.NFrame - nLeadingFrames;
    nframes = mod(firstNFrame + (0:numFrames-1),1024);
    sfns = sprintf('(%d...%d)',nframes(1),nframes(end));
    title(['Received Resource Grid. System Frame Number: ' sfns]);

end