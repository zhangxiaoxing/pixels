% motifs PSTH
% [14,18,22,33,34,68,100,102,114]


%% show case
for sessid=[68,100,102,114]
    currentTrls=[];
    for maxid=1:10
        [ttrial,~,meta]=wave.chain_loop_SC_spk(false,sessid,maxid,'skip_plot',true,'by_cover_rate',true);
        if ismember(meta(1),currentTrls) || isempty(ttrial)
            continue
        end
        currentTrls=[currentTrls,meta(1)];
        covered=zeros(1,1);
        tsmin=max(min(cell2mat(ttrial(:,4))));
        umotif=unique([ttrial{:,2}].');
        % sort


        fh=figure();
        fh.Position(1:2)=[100,100];
        tiledlayout(2,1)
        % motif raster
        nexttile()
        hold on
        % for each motif
        for mmid=1:size(ttrial,1)
            %extract ts
            ts=(ttrial{mmid,4}(:,2)-tsmin)./30;
            %scatter %optional alternative plot
            [~,umotifid]=ismember(ttrial{mmid,2}{1},umotif);
            scatter(ts,umotifid,'k','|')
            onset=floor(ts(1))+1;
            offset=ceil(ts(end));
            if offset>numel(covered)
                covered(end+1:offset)=0;
            end
            covered(onset:offset)=covered(onset:offset)+1;
        end

        xlim([0,meta(3)*1000]);
        xlabel('Time (ms)')
        ylim([0.5,size(umotif,1)+0.5])
        ylabel('Motifs #')
        edges = find(diff([0,covered~=0,0]==1));
        onset = edges(1:2:end-1);  % Start indices
        run_length = edges(2:2:end)-onset;  % Consecutive ones counts
        ymax=max(ylim());
        for ii=1:numel(onset)
            plot([onset(ii),onset(ii)+run_length(ii)],[ymax,ymax],'r-','LineWidth',4);
        end
        % PSTH
        nexttile()
        hold on
        covered(end+1:meta(3)*1000)=0;
        plot(250:500:meta(3)*1000-250,mean(reshape(covered,500,[])).*10,'r-');
        ylabel('Motif frequency (Hz)')
        xlabel('Time (ms)')
        title("S"+sessid+"M"+maxid+"T"+meta(1));
        appendfig('fn','motif_psth.pdf','close',true)
    end
end

%% statistics
% moved to chain_loop_stats->plot_motif_freq()

