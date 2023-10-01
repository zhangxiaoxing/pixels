%% spike proportion
function motif_spike_total_spike_proportion(opt)
arguments
    opt.denovo (1,1) logical = false
end
if opt.denovo
    shuf_ratios=cell(100,1);
    shuf_chain=cell(100,1);
    shuf_loop=cell(100,1);
    blame=vcs.blame();
    for shufidx=0:100
        if rem(shufidx,10)==0
            disp("shuf#"+shufidx)
        end
        if shufidx==0
            pivot_dict=wave.replay.delay_vs_iti_motif_per_spike(skip_save=true,shuf=false,delay_only=true);
        else
            pivot_dict=wave.replay.delay_vs_iti_motif_per_spike(skip_save=true,shuf=true,shufidx=shufidx,delay_only=true);
        end
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
                if shufidx==0
                    ratios(selidx,:)=[motifspk,allspk];
                else
                    shuf_ratios{shufidx}(selidx,:)=[motifspk,allspk];
                end

                chainidx=strcmp(ucid_chain,ucid{selidx});
                if any(chainidx)
                    chainspk=pivot_dict.delay_chain.(ucid{selidx}).numEntries;
                    if shufidx==0
                        ratios_chain(chainidx,:)=[chainspk,allspk];
                    else
                        shuf_chain{shufidx}(chainidx,:)=[chainspk,allspk];
                    end
                end

                loopidx=strcmp(ucid_loop,ucid{selidx});
                if any(loopidx)
                    loopspk=pivot_dict.delay_loop.(ucid{selidx}).numEntries;
                    if shufidx==0
                        ratios_loop(loopidx,:)=[loopspk,allspk];
                    else
                        shuf_loop{shufidx}(loopidx,:)=[loopspk,allspk];
                    end
                end
            end
        end
        save(fullfile('binary','motif_spike_total_spike_proportion.mat'),'shuf_ratios','shuf_loop','shuf_chain','ratios','ratios_chain','ratios_loop','blame')
    end
else
    load(fullfile('binary','motif_spike_total_spike_proportion.mat'),'shuf_ratios','shuf_loop','shuf_chain','ratios','ratios_chain','ratios_loop')

end

if false
    % spkratio.all_motif_spikes=table(replace(ucid,'ses','session'),ratios(:,1),ratios(:,2),'VariableNames',{'NeuronID','Motif_spikes','Delay_spikes'});
    % spkratio.chain_spikes=table(replace(ucid_chain,'ses','session'),ratios_chain(:,1),ratios_chain(:,2),'VariableNames',{'NeuronID','Chain_spikes','Delay_spikes'});
    % spkratio.loop_spikes=table(replace(ucid_loop,'ses','session'),ratios_loop(:,1),ratios_loop(:,2),'VariableNames',{'NeuronID','Loop_spikes','Delay_spikes'});
    spkratio.observed_motif_spike_and_memory_delay_spike=ratios;
    spkratio.shuffled_motif_spike_and_memory_delay_spike=shuf_ratios;

    spkratio.observed_chain_spike_and_memory_delay_spike=ratios_chain;
    spkratio.shuffled_chain_spike_and_memory_delay_spike=shuf_chain;

    spkratio.observed_loop_spike_and_memory_delay_spike=ratios_loop;
    spkratio.shuffled_loop_spike_and_memory_delay_spike=shuf_loop;

    fid=fopen(fullfile('binary','upload','F2O_SF6A_motif_spike_over_WM_delay_spike.json'),'w');
    fprintf(fid,jsonencode(spkratio));
    fclose(fid)
end
obdata={ratios,ratios_chain,ratios_loop};
shufdata={shuf_ratios,shuf_chain,shuf_loop};
fh=figure();
tiledlayout(1,3)
for didx=1:3
    nexttile();
    dd=obdata{didx};
    ss=cell2mat(shufdata{didx});
    currratios=[dd(:,1)./dd(:,2),ones(size(dd,1),1);...
        ss(:,1)./ss(:,2),2.*ones(size(ss,1),1)];

    bh=boxplot(currratios(:,1),currratios(:,2),'Colors','k','Whisker',inf);
    qtrs=prctile(currratios,[25,50,75]);
    ylim([0,0.25])
    % set(gca(),'YScale','log')
    title(sprintf('max%.2f,%.2f',...
        max(currratios(currratios(:,2)==1,1)),...
        max(currratios(currratios(:,2)==2,1))))
end
sgtitle("L->R:chn+loop,chain,loop")

savefig(fh,fullfile("binary","motif_spike_total_spike_proportion.fig"))

shufmm=cell2mat(cellfun(@(x) nanmean(x(:,1)./x(:,2)), shuf_ratios, 'UniformOutput', false));
smm=mean(shufmm);
sstd=std(shufmm);
zscore=(mean(ratios(:,1)./ratios(:,2))-smm)./sstd;
disp(zscore)