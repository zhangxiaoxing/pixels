%% spike proportion
function motif_spike_total_spike_proportion()
load(fullfile('binary','delay_vs_iti_pivot_spike.mat'),'pivot_dict')
ucid=fieldnames(pivot_dict.delay);
ucid_chain=fieldnames(pivot_dict.delay_chain);
ucid_loop=fieldnames(pivot_dict.delay_loop);
usess=str2double(unique(regexp([ucid;ucid_chain;ucid_loop],'(?<=ses)\d{1,3}(?=s)','match','once')));
ratios=nan(numel(ucid),2);
ratios_chain=nan(numel(ucid_chain),2);
ratios_loop=nan(numel(ucid_loop),2);
for sessid=reshape(usess,1,[])
    [~,~,trials,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true,'jagged',true);
    sesssel=startsWith(ucid,"ses"+sessid+"s");

    for selidx=reshape(find(sesssel),1,[])
        curruid=regexp(ucid(selidx),'(?<=U)\d*$','match','once');
        currsample=str2double(regexp(ucid(selidx),'(?<=\ds)\d(?=d)','match','once'));
        currdelay=str2double(regexp(ucid(selidx),'(?<=s\dd)\d(?=U)','match','once'));
        motifspk=pivot_dict.delay.(ucid{selidx}).numEntries;
        if currsample==1
            samp=4;
        else
            samp=8;
        end
        trialsel=find(trials(:,5)==samp & trials(:,8)==currdelay & all(trials(:,9:10),2));
        idsel=strcmp(FT_SPIKE.label,curruid);
        allspk=nnz(ismember(FT_SPIKE.trial{idsel},trialsel) & FT_SPIKE.time{idsel}>1 & FT_SPIKE.time{idsel}<=currdelay+1);
        ratios(selidx,:)=[motifspk,allspk];

        chainidx=strcmp(ucid_chain,ucid{selidx});
        if any(chainidx)
            chainspk=pivot_dict.delay_chain.(ucid{selidx}).numEntries;
            ratios_chain(chainidx,:)=[chainspk,allspk];
        end

        loopidx=strcmp(ucid_loop,ucid{selidx});
        if any(loopidx)
            loopspk=pivot_dict.delay_loop.(ucid{selidx}).numEntries;
            ratios_loop(loopidx,:)=[loopspk,allspk];
        end
    end
end

if false
    spkratio.all_motif_spikes=table(replace(ucid,'ses','session'),ratios(:,1),ratios(:,2),'VariableNames',{'NeuronID','Motif_spikes','Delay_spikes'});
    spkratio.chain_spikes=table(replace(ucid_chain,'ses','session'),ratios_chain(:,1),ratios_chain(:,2),'VariableNames',{'NeuronID','Chain_spikes','Delay_spikes'});
    spkratio.loop_spikes=table(replace(ucid_loop,'ses','session'),ratios_loop(:,1),ratios_loop(:,2),'VariableNames',{'NeuronID','Loop_spikes','Delay_spikes'});
    
    fid=fopen(fullfile('binary','upload','F2O_motif_spike_over_WM_delay_spike.json'),'w');
    fprintf(fid,jsonencode(spkratio));
    fclose(fid)

end
fh=figure();
tiledlayout(1,3)
for dd={ratios,ratios_chain,ratios_loop}
    nexttile();
    currratios=dd{1}(:,1)./dd{1}(:,2);
    bh=boxplot(currratios,'Colors','k','Whisker',inf);
    qtrs=prctile(currratios,[25,50,75]);
    ylim([0,0.25])
    % set(gca(),'YScale','log')
    title(sprintf('max%.2f',max(currratios)))
end
sgtitle("L->R:chn+loop,chain,loop")

save(fh,fullfile("binary","motif_spike_total_spike_proportion.fig"))

