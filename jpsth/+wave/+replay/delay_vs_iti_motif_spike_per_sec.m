function composite_spk_per_sec=delay_vs_iti_motif_spike_per_sec(chain_replay,ring_replay,trials_dict,opt)
arguments
    chain_replay = []
    ring_replay = []
    trials_dict = []
    opt.skip_save (1,1) logical = true
    opt.per_unit_motif (1,1) logical = false
end

if isempty(ring_replay)
    load(fullfile('binary','motif_replay.mat'),'ring_replay');
end

if isempty(chain_replay)
    load(fullfile('binary','motif_replay.mat'),'chain_replay');
end

if isempty(trials_dict)
    load(fullfile('binary','trials_dict.mat'),'trials_dict');
end

sps=30000;
% per session
composite_spk_per_sec=[];
for sess=reshape(unique([ring_replay.session;chain_replay.session]),1,[])
    disp(sess)
    session_tick=wave.replay.sessid2length(sess);
    trials=cell2mat(trials_dict(sess));
    out_task_sec=(trials(1,1)./sps-60)+((session_tick-trials(end,2))./sps-60-1);


    for onewave=["s1d3","s1d6","s2d3","s2d6"]
        chain_sel=chain_replay.session==sess & chain_replay.wave==onewave;
        ring_sel=ring_replay.session==sess & ring_replay.wave==onewave;

        chain_su=[chain_replay.meta{chain_sel,2}];
        ring_su=[ring_replay.meta{ring_sel,2}];

        if numel([chain_su,ring_su])==numel(unique([chain_su,ring_su]))
            continue
        end

        [delay_spikes,iti_spikes,out_task_spikes,npdelay_spikes,npiti_spikes]=deal([]);
        switch onewave
            case "s1d3"
                samp=4;delay=3;
            case "s1d6"
                samp=4;delay=6;
            case "s2d3"
                samp=8;delay=3;
            case "s2d6"
                samp=8;delay=6;
        end
        % chain ------------------------------------------------
        for chainii=reshape(find(chain_sel),1,[])
            % per preferred trial
            trl_align=chain_replay.trl_align{chainii};
            pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
            if opt.per_unit_motif
                delay_spikes=[delay_spikes;1,chainii,numel(chain_replay.meta{chainii,2}).*nnz(pref_delay)];
            else
                delay_spikes=unique([delay_spikes;reshape(chain_replay.meta{chainii,2}+chain_replay.ts{chainii}(pref_delay,:).*100000,[],1)]);

                nonpref_delay=all(trl_align(:,5:6)==1,2) & trl_align(:,3)~=samp & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1) & trl_align(:,4)==delay;
                npdelay_spikes=unique([npdelay_spikes;reshape(chain_replay.meta{chainii,2}+chain_replay.ts{chainii}(nonpref_delay,:).*100000,[],1)]);

                pref_succeed_iti=all(trl_align(:,5:7)==1,2)... % WT, pref
                    & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test/ 1s sample /1s test, 3s response?
                    & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));
                iti_spikes=unique([iti_spikes;reshape(chain_replay.meta{chainii,2}+100000.*chain_replay.ts{chainii}(pref_succeed_iti,:),[],1)]);

                nonpref_iti=all(trl_align(:,5:6)==1,2) & trl_align(:,3)~=samp & trl_align(:,4)==delay...
                    & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test/ 1s sample /1s test, 3s response?
                    & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));
                npiti_spikes=unique([npiti_spikes;reshape(chain_replay.meta{chainii,2}+100000.*chain_replay.ts{chainii}(nonpref_iti,:),[],1)]);

                %  corresponding network in pre task, post task
                pre_post_motif=(trl_align(:,8)==1 & trl_align(:,9)>60) ...
                    | (trl_align(:,8)<0 & trl_align(:,2)>(60+5+trl_align(:,4)));
                out_task_spikes=unique([out_task_spikes;reshape(chain_replay.meta{chainii,2}+100000.*chain_replay.ts{chainii}(pre_post_motif,:),[],1)]);
            end
        end

        % loops -------------------------------------------
        for cii=reshape(find(ring_sel),1,[])
            % per preferred trial
            trl_align=ring_replay.trl_align{cii};
            pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
            if opt.per_unit_motif
                delay_spikes=[delay_spikes;2,cii,numel(ring_replay.meta{cii,2}).*nnz(pref_delay)];
            else
                delay_spikes=unique([delay_spikes;cell2mat(ring_replay.ts_seq{cii}(pref_delay))+100000.*cell2mat(ring_replay.ts{cii}(pref_delay))]);

                nonpref_delay=all(trl_align(:,5:6)==1,2) & trl_align(:,3)~=samp & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1)  & trl_align(:,4)==delay;
                npdelay_spikes=unique([npdelay_spikes;cell2mat(ring_replay.ts_seq{cii}(nonpref_delay))+100000.*cell2mat(ring_replay.ts{cii}(nonpref_delay))]);

                pref_succeed_iti=all(trl_align(:,5:7)==1,2)... % WT, pref
                    & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test/ 1s sample /1s test, 3s response?
                    & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));
                iti_spikes=unique([iti_spikes;cell2mat(ring_replay.ts_seq{cii}(pref_succeed_iti))+100000.*cell2mat(ring_replay.ts{cii}(pref_succeed_iti))]);


                nonpref_iti=all(trl_align(:,5:6)==1,2) & trl_align(:,3)~=samp & trl_align(:,4)==delay...
                    & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test/ 1s sample /1s test, 3s response?
                    & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));
                npiti_spikes=unique([iti_spikes;cell2mat(ring_replay.ts_seq{cii}(nonpref_iti))+100000.*cell2mat(ring_replay.ts{cii}(nonpref_iti))]);

                % corresponding network in pre task, post task
                pre_post_motif=(trl_align(:,8)==1 & trl_align(:,9)>60) ...
                    | (trl_align(:,8)<0 & trl_align(:,2)>(60+5+trl_align(:,4)));

                out_task_spikes=unique([out_task_spikes;cell2mat(ring_replay.ts_seq{cii}(pre_post_motif))+100000.*cell2mat(ring_replay.ts{cii}(pre_post_motif))]);
            end
        end

        trl_sel=find(trials(:,5)==samp & trials(:,8)==delay & all(trials(:,9:10)>0,2));
        np_trl_sel=find(trials(:,5)==setdiff([4,8],samp) & trials(:,8)==delay & all(trials(:,9:10)>0,2));

        pref_delay_sec=sum(diff(trials(trl_sel,1:2),1,2)./sps-1);
        np_delay_sec=sum(diff(trials(np_trl_sel,1:2),1,2)./sps-1);

        trials(end+1,:)=trials(end,2)+14*sps;
        pref_succeed_iti_sec=sum((trials(trl_sel+1,1)-trials(trl_sel,2))./sps-4); % 1s test + 3s response
        npiti_sec=sum((trials(np_trl_sel+1,1)-trials(np_trl_sel,2))./sps-4); % 1s test + 3s response
        trials(end,:)=[];
        if opt.per_unit_motif
            composite_spk_per_sec=[composite_spk_per_sec;...
                cell2table({sess,onewave,[delay_spikes,repmat(pref_delay_sec,size(delay_spikes,1),1)]},...
                'VariableNames',{'sess','wave','delay'})];
        else
            composite_spk_per_sec=[composite_spk_per_sec;...
                cell2table({sess,onewave,[numel(delay_spikes),pref_delay_sec],[numel(iti_spikes),pref_succeed_iti_sec],[numel(out_task_spikes),out_task_sec],...
                [numel(npdelay_spikes),np_delay_sec],[numel(npiti_spikes),npiti_sec]},...
                'VariableNames',{'sess','wave','delay','iti','out_task','npdelay','npiti'})];
        end
    end
end

if ~opt.skip_save
    blame=vcs.blame();
    if opt.per_unit_motif
        save(fullfile("binary","chains_loops_freq.mat"),"composite_spk_per_sec","blame");
    else
        save(fullfile("binary","delay_iti_motif_spike_per_sec.mat"),"composite_spk_per_sec","blame");
    end
end

%% plot
if opt.per_unit_motif
    dmat=cell2mat(composite_spk_per_sec.delay);
    qtrs=prctile(dmat(:,3)./dmat(:,4),[25,50,75]);
    fh=figure();
    tiledlayout(1,3)
    nexttile()
    chainfreq=dmat(dmat(:,1)==1,3)./dmat(dmat(:,1)==1,4);
    boxplot(chainfreq,'Whisker',inf,'Colors','k','Widths',1)
    set(gca,'YScale','log','YLim',[1e-2,1e2],'XTick',[])
    title('Chains')
    ylabel('Motif spike rate (Hz, cross-motifs)')
    nexttile()
    ringfreq=dmat(dmat(:,1)==2,3)./dmat(dmat(:,1)==2,4);
    boxplot(ringfreq,'Whisker',inf,'Colors','k','Widths',1)
    set(gca,'YScale','log','YLim',[1e-2,1e2],'XTick',[])
    ylabel('Motif spike rate (Hz, cross-motifs)')
    title('Loops')
    nexttile()
    mfreq=dmat(:,3)./dmat(:,4);
    boxplot(mfreq,'Whisker',inf,'Colors','k','Widths',1)
    set(gca,'YScale','log','YLim',[1e-2,1e2],'XTick',[])
    ylabel('Motif spike rate (Hz, cross-motifs)')
    title('All motifs')
    savefig(fh,fullfile("binary","chains_loops_freq.fig"))

else
    dps=composite_spk_per_sec.delay;
    dprop=(dps(:,1)./dps(:,2));

    npdps=composite_spk_per_sec.npdelay;
    npdprop=(npdps(:,1)./npdps(:,2));

    ips=composite_spk_per_sec.iti;
    iprop=(ips(:,1)./ips(:,2));

    npips=composite_spk_per_sec.npiti;
    npiprop=(npips(:,1)./npips(:,2));

    sps=composite_spk_per_sec.out_task;
    sprop=(sps(:,1)./sps(:,2)); % Due to 4 waves

    mm=[mean(dprop),mean(npdprop),mean(iprop),mean(npiprop),mean(sprop)];
    sem=[std(dprop),std(npdprop),std(iprop),std(npiprop),std(sprop)]./sqrt([numel(dprop),numel(npdprop),numel(iprop),numel(npiprop),numel(sprop)]);


    fh=figure('Position',[100,100,400,300]);
    hold on
    bh=bar(mm,'FaceColor','none');
    errorbar(bh.XEndPoints,bh.YEndPoints,sem,'k.')
    set(gca,'XTick',1:5,'XTickLabelRotation',90,'XTickLabel',{'Delay','NP.Delay','ITI','NP.ITI','out_task'});
    ylabel('Motif spike frequency (Hz)');

    [~,piti]=ttest(dprop,iprop);
    [~,psur]=ttest(dprop,sprop);

    title(sprintf('%.4f,',piti,psur));
    savefig(fh,fullfile("binary","delay_iti_motif_spike_per_sec.fig"))
end
end

