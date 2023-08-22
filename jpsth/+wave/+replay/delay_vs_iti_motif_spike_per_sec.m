function composite_spk_per_sec=delay_vs_iti_motif_spike_per_sec(chain_replay,ring_replay_tbl,trials_dict)
sps=30000;
% per session
composite_spk_per_sec=[];
for sess=reshape(unique([ring_replay_tbl.session;chain_replay.session]),1,[])
    disp(sess)
    session_tick=wave.replay.sessid2length(sess);
    trials=cell2mat(trials_dict(sess));
    surround_sec=(trials(1,1)./sps-60)+((session_tick-trials(end,2))./sps-60-1);
    

    for onewave=["s1d3","s1d6","s2d3","s2d6"]
        chain_sel=chain_replay.session==sess & chain_replay.wave==onewave;
        ring_sel=ring_replay_tbl.session==sess & ring_replay_tbl.wave==onewave;
        if nnz(chain_sel)+nnz(ring_sel)<2 % TODO: should be 1 or 2?
            continue
        end
        
        [delay_spikes,iti_spikes,surround_spikes]=deal([]);

        % chain ------------------------------------------------
        for chainii=reshape(find(chain_sel),1,[])
            % per preferred trial
            trl_align=chain_replay.trl_align{chainii};
            pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
            delay_spikes=unique([delay_spikes;reshape(chain_replay.ts{chainii}(pref_delay,:),[],1)]);

            pref_succeed_iti=all(trl_align(:,5:7)==1,2)... % WT, pref
                & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test/ 1s sample /1s test, 3s response?
                & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));
            
            iti_spikes=unique([iti_spikes;reshape(chain_replay.ts{chainii}(pref_succeed_iti,:),[],1)]);

                %  corresponding network in pre task, post task
            pre_post_motif=(trl_align(:,8)==1 & trl_align(:,9)>60) ...
                | (trl_align(:,8)<0 & trl_align(:,2)>(60+5+trl_align(:,4)));
            surround_spikes=unique([surround_spikes;reshape(chain_replay.ts{chainii}(pre_post_motif,:),[],1)]);
        end

        % loops -------------------------------------------
        for cii=reshape(find(ring_sel),1,[])
            % per preferred trial
            trl_align=ring_replay_tbl.trl_align{cii};
            pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
            delay_spikes=unique([delay_spikes;cell2mat(ring_replay_tbl.ts{cii}(pref_delay))]);

            pref_succeed_iti=all(trl_align(:,5:7)==1,2)... % WT, pref
                & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test/ 1s sample /1s test, 3s response?
                & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));
            iti_spikes=unique([iti_spikes;cell2mat(ring_replay_tbl.ts{cii}(pref_succeed_iti))]);

                % corresponding network in pre task, post task
            pre_post_motif=(trl_align(:,8)==1 & trl_align(:,9)>60) ...
                | (trl_align(:,8)<0 & trl_align(:,2)>(60+5+trl_align(:,4)));

            surround_spikes=unique([surround_spikes;cell2mat(ring_replay_tbl.ts{cii}(pre_post_motif))]);
        end
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
        trl_sel=find(trials(:,5)==samp & trials(:,8)==delay & all(trials(:,9:10)>0,2));
        pref_delay_sec=sum(diff(trials(trl_sel,1:2),1,2)./sps-1);
        
        trials(end+1,:)=trials(end,2)+14*sps;
        pref_succeed_iti_sec=sum((trials(trl_sel+1,1)-trials(trl_sel,2))./sps-4); % 1s test + 3s response
        trials(end,:)=[];

        composite_spk_per_sec=[composite_spk_per_sec;...
            cell2table({sess,onewave,[numel(delay_spikes),pref_delay_sec],[numel(iti_spikes),pref_succeed_iti_sec],[numel(surround_spikes),surround_sec]},...
            'VariableNames',{'sess','wave','delay','iti','surround'})];
    end
end

%% plot

dps=composite_spk_per_sec.delay;
dprop=(dps(:,1)./dps(:,2));

ips=composite_spk_per_sec.iti;
iprop=(ips(:,1)./ips(:,2));

sps=composite_spk_per_sec.surround;
sprop=(sps(:,1)./sps(:,2)); % Due to 4 waves

mm=[mean(dprop),mean(iprop),mean(sprop)];
sem=[std(dprop),std(iprop),std(sprop)]./sqrt([numel(dprop),numel(iprop),numel(sprop)]);


fh=figure('Position',[100,100,400,300]);
hold on
bh=bar(mm,'FaceColor','none');
errorbar(bh.XEndPoints,bh.YEndPoints,sem,'k.')
set(gca,'XTick',1:3,'XTickLabelRotation',90,'XTickLabel',{'Delay','ITI','Surround'});
ylabel('Motif spike frequency (Hz)');

[~,piti]=ttest(dprop,iprop);
[~,psur]=ttest(dprop,sprop);

title(sprintf('%.4f,',piti,psur));







keyboard();
end

