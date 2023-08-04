% motif freq ramping before, after task

motif_replay=chain_replay;

before_task_sums=nan(0,15);
before_task_sums_long=nan(0,15);
after_task_sums=nan(0,15);
after_task_sums_long=nan(0,15);

for dd=reshape(fieldnames(motif_replay),1,[])
    for ww=reshape(fieldnames(motif_replay.(dd{1})),1,[])
        for cc=reshape(fieldnames(motif_replay.(dd{1}).(ww{1})),1,[])
            onechain=motif_replay.(dd{1}).(ww{1}).(cc{1});

            % before task
            before_sel=onechain.trl_align(:,1)<0 & onechain.trl_align(:,9)<=30;
            if any(before_sel)
                motif_ts=onechain.trl_align(before_sel,9);
                onehist=histcounts(motif_ts,0:2:30,'Normalization','probability');
                before_task_sums=[before_task_sums;onehist];
            end

            % before-long
            if any(onechain.trl_align(:,9)>600)
                before_sel=onechain.trl_align(:,1)<0 & onechain.trl_align(:,9)<=600;
                if any(before_sel)
                    motif_ts=onechain.trl_align(before_sel,9);
                    onehist=histcounts(motif_ts,0:20:600,'Normalization','probability');
                    before_task_sums_long=[before_task_sums_long;onehist];
                end
            end

            % after task
            after_onset=60+2+onechain.trl_align(:,4);
            after_sel=onechain.trl_align(:,8)<0 & onechain.trl_align(:,2)>after_onset & onechain.trl_align(:,2)<=after_onset+30;
            if any(after_sel)
                motif_ts=onechain.trl_align(after_sel,2);
                onehist=histcounts(motif_ts,after_onset:2:after_onset+30,'Normalization','probability');
                after_task_sums=[after_task_sums;onehist];
            end

            if any(onechain.trl_align(:,2)>after_onset+600)
                after_sel=onechain.trl_align(:,8)<0 & onechain.trl_align(:,2)>after_onset & onechain.trl_align(:,2)<=after_onset+600;
                if any(after_sel)
                    motif_ts=onechain.trl_align(after_sel,2);
                    onehist=histcounts(motif_ts,after_onset:20:after_onset+600,'Normalization','probability');
                    after_task_sums_long=[after_task_sums_long;onehist];
                end
            end

        end
    end
end

figure()
tiledlayout(1,2)
nexttile
hold on;

beforelh=plot(-590:20:-10,fliplr(mean(before_task_sums_long)),'r-');
afterlh=plot(10:20:590,mean(after_task_sums_long),'k-');
ylabel('Normalized motif spike rate')
xlabel('Time from task (s)')
nexttile
hold on
beforeh=plot(-29:2:-1,fliplr(mean(before_task_sums)),'r-');
afterh=plot(1:2:29,mean(after_task_sums),'k-');
ylabel('Normalized motif spike rate')
xlabel('Time from task (s)')