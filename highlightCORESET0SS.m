function highlightCORESET0SS(csetSubcarriers,monSlots,detSlot,pdcch,dciCRC)

    ssFirstSym = pdcch.SearchSpace.StartSymbolWithinSlot;
    csetDuration = pdcch.CORESET.Duration;

    % Define colors and plotting function
    basePlotProps = {'LineStyle','-','LineWidth',1};
    occasionColor = 0.7*[1 1 1];
    detectionColor = [200,0,0]/255;
    boundingBox = @(y,x,h,w,varargin)rectangle('Position',[x+0.5 y-0.5 w h],basePlotProps{:},varargin{:});

    % Highlight all CORESET 0/SS occasions related to the detected SSB
    k0 = csetSubcarriers(1);
    K = length(csetSubcarriers);
    ssSym = 14*monSlots + ssFirstSym ;
    for i = 1:length(ssSym)
        boundingBox(k0,ssSym(i),K,csetDuration,'EdgeColor',occasionColor);
    end

    if dciCRC == 0
        % Highlight decoded PDCCH
        ssSym = 14*detSlot + ssFirstSym;
        boundingBox(k0,ssSym,K,csetDuration,'EdgeColor',detectionColor);

        % Add text next to decoded PDCCH
        text(ssSym+csetDuration+1,k0+24,0,'PDCCH','FontSize',10,'Color','w')
    end

end