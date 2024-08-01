function plotResourceGridSIB1(slotGrid,carrier,pdcch,pdsch,tcr,K0)

    % Display the OFDM grid of the slot containing decoded PDCCH
    figure;
    imagesc(abs(slotGrid(:,:,1))); axis xy
    xlabel('OFDM symbol');
    ylabel('subcarrier');
    title('Slot Containing Decoded PDCCH');

    aggregationLevelIndex = log2(pdcch.AggregationLevel)+1;
    candidate = pdcch.AllocatedCandidate;

    % Define auxiliary plotting function
    color = [200,0,0]/255;
    boundingBox = @(y,x,h,w,varargin)rectangle('Position',[x+0.5 y-0.5 w h],'EdgeColor',color,varargin{:});

    % Highlight PDCCH in resource grid
    carrier.NSlot = carrier.NSlot - K0; % Substract slot offset K0 for subscripts calculations
    subsPdcch = nrPDCCHSpace(carrier,pdcch,'IndexStyle','Subs');
    subsPdcch = double(subsPdcch{aggregationLevelIndex}(:,:,candidate));
    x = min(subsPdcch(:,2))-1; X = max(subsPdcch(:,2))-x;
    y = min(subsPdcch(:,1)); Y = max(subsPdcch(:,1))-y+1;
    boundingBox(y,x,Y,X);
    str = sprintf('PDCCH \nAggregation Level: %d\nCandidate: %d',2.^(aggregationLevelIndex-1),candidate-1);
    text(x+X+1,y+Y/2,0,str,'FontSize',10,'Color','w')

    % Highlight PDSCH and PDSCH DM-RS in resource grid
    carrier.NSlot = carrier.NSlot + K0; % Add back slot offset K0 for subscripts calculations
    subsPdschSym = double(nrPDSCHIndices(carrier,pdsch,'IndexStyle','subscript'));
    subsPdschDmrs = double(nrPDSCHDMRSIndices(carrier,pdsch,'IndexStyle','subscript'));
    subsPdsch = [subsPdschSym;subsPdschDmrs];
    x = min(subsPdsch(:,2))-1; X = max(subsPdsch(:,2))-x;
    y = min(subsPdsch(:,1)); Y = max(subsPdsch(:,1))-y+1;
    boundingBox(y,x,Y,X);
    str = sprintf('PDSCH (SIB1) \nModulation: %s\nCode rate: %.2f',pdsch.Modulation,tcr);
    text(x+4,y+Y+60,0, str,'FontSize',10,'Color','w')

end