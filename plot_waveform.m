fl=dir('D:\WF\**\waveform.mat');
for f=fl'
    fstr=load(fullfile(f.folder,f.name));
    w=fstr.waveform;
    d=ceil(sqrt(size(w,1)));
    fh=figure('Color','w','Position',[100,100,800,800]);
    for sidx=1:size(w,1)
        subplot(d,d,sidx)
        onetrace=mean(squeeze(w{sidx,3}.data(1,:,:))');
        plot(onetrace);
        set(gca(),'XTick',[],'YTick',[])
        xlim([0,92])
        
        
        % Extracellular Spike Waveform Dissociates Four Functionally Distinct Cell Classes in Primate Cortex
        % Current Biology 2019 Earl K. Miller,Markus Siegel
        % excluded waveforms that satisfied any of three criteria for atypical shape:
        % (1) the amplitude of the main trough was smaller than the subsequent positive
        % peak (n = 41), (2) the trace was noisy, defined as > = 6 local maxima of magnitude > = 0.01 (n = 38), (3) there was one or more local
        % maxima in the period between the main trough and the subsequent peak (n = 35).
        
        % aligned to baseline during extraction
        
        [trough,idx_min]=min(onetrace);
        [late_peak,idx_late_peak]=max(onetrace(idx_min:end));
        
        if trough < late_peak
            
        end
        
        
        
    end
    pause
    close(fh);
end