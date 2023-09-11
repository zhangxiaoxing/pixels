function sums=uniq_spike_freq(chain_replay,ring_replay,trials_dict)
arguments
    chain_replay = []
    ring_replay = []
    trials_dict = []
end
if isempty(chain_replay) || isempty (ring_replay)
    load(fullfile('binary','motif_replay.mat'),'ring_replay','chain_replay');
end
if isempty(trials_dict)
    load(fullfile('binary','trials_dict.mat'),'trials_dict');
end

sps=30000;

sums=cell2struct({[];[];[];[]},{'delay','iti','before','after'});


usess=union(chain_replay.session,ring_replay.session);
for sessid=reshape(usess,1,[])
    trials=cell2mat(trials_dict(sessid));
    session_tick=wave.replay.sessid2length(sessid);
    
    for wid=["s1d3","s2d3","s1d6","s2d6"]
        delay_motif_spk=cell(0,1);
        delay_motif_id=[];
        delay_motif_seq=cell(0,1);

        iti_motif_spk=cell(0,1);
        iti_motif_id=[];
        iti_motif_seq=cell(0,1);

        before_motif_spk=cell(0,1);
        before_motif_id=[];
        before_motif_seq=cell(0,1);

        after_motif_spk=cell(0,1);
        after_motif_id=[];
        after_motif_seq=cell(0,1);


        disp(wid+"of sess"+sessid)
        samp=str2double(regexp(wid,"(?<=s)\d(?=d)",'match','once')).*4;
        delay=str2double(regexp(wid,"(?<=d)\d$",'match','once'));

        ringsel=find(ring_replay.session==sessid & ring_replay.wave==wid);
        chainsel=find(chain_replay.session==sessid & chain_replay.wave==wid);

        trl_sel=find(trials(:,5)==samp & trials(:,8)==delay & all(trials(:,9:10)>0,2));
        pref_delay_sec=sum(diff(trials(trl_sel,1:2),1,2)./sps-1);
        trials(end+1,:)=trials(end,2)+14*sps;
        pref_succeed_iti_sec=sum((trials(trl_sel+1,1)-trials(trl_sel,2))./sps-4); % 1s test + 3s response
        trials(end,:)=[];
        before_sec=trials(1,1)./sps-60;
        after_sec=(session_tick-trials(end,2))./sps-60-1;

        for ridx=reshape(ringsel,1,[])
            trl_align=ring_replay.trl_align{ridx};

            pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
            pref_succeed_iti=all(trl_align(:,5:7)==1,2)... % WT, pref
                & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test/ 1s sample /1s test, 3s response?
                & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));
            before_task=(trl_align(:,8)==1 & trl_align(:,9)>60);
            after_task=trl_align(:,8)<0 & trl_align(:,2)>(60+5+trl_align(:,4));

            currts=ring_replay.ts{ridx}(pref_delay);
            if ~isempty(currts)
                delay_motif_spk=[delay_motif_spk;{currts}];
                delay_motif_id=[delay_motif_id;"r"+ridx];
                delay_motif_seq=[delay_motif_seq;{ring_replay.ts_seq{ridx}(pref_delay)}];
            end

            currts=ring_replay.ts{ridx}(pref_succeed_iti);
            if ~isempty(currts)
                iti_motif_spk=[iti_motif_spk;{currts}];
                iti_motif_id=[iti_motif_id;"r"+ridx];
                iti_motif_seq=[iti_motif_seq;{ring_replay.ts_seq{ridx}(pref_succeed_iti)}];
            end
            

            currts=ring_replay.ts{ridx}(before_task);
            if ~isempty(currts)
                before_motif_spk=[before_motif_spk;{currts}];
                before_motif_id=[before_motif_id;"r"+ridx];
                before_motif_seq=[before_motif_seq;{ring_replay.ts_seq{ridx}(before_task)}];
            end

            currts=ring_replay.ts{ridx}(after_task);
            if ~isempty(currts)
                after_motif_spk=[after_motif_spk;{currts}];
                after_motif_id=[after_motif_id;"r"+ridx];
                after_motif_seq=[after_motif_seq;{ring_replay.ts_seq{ridx}(after_task)}];
            end
        end

        for cidx=reshape(chainsel,1,[])
            trl_align=chain_replay.trl_align{cidx};

            pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
            pref_succeed_iti=all(trl_align(:,5:7)==1,2)... % WT, pref
                & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test/ 1s sample /1s test, 3s response?
                & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));
            before_task=(trl_align(:,8)==1 & trl_align(:,9)>60);
            after_task=trl_align(:,8)<0 & trl_align(:,2)>(60+5+trl_align(:,4));

            currts=chain_replay.ts{cidx}(pref_delay,:);
            if ~isempty(currts)
                delay_motif_spk=[delay_motif_spk;{mat2cell(currts.',size(currts,2),ones(size(currts,1),1)).'}];
                delay_motif_id=[delay_motif_id;"c"+cidx];
                delay_motif_seq=[delay_motif_seq;{repmat({chain_replay.meta{cidx,2}.'},nnz(pref_delay),1)}];
            end
            
            currts=chain_replay.ts{cidx}(pref_succeed_iti,:);
            if ~isempty(currts)
                iti_motif_spk=[iti_motif_spk;{mat2cell(currts.',size(currts,2),ones(size(currts,1),1)).'}];
                iti_motif_id=[iti_motif_id;"c"+cidx];
                iti_motif_seq=[iti_motif_seq;{repmat({chain_replay.meta{cidx,2}.'},nnz(pref_succeed_iti),1)}];
            end


            currts=chain_replay.ts{cidx}(before_task,:);
            if ~isempty(currts)
                before_motif_spk=[before_motif_spk;{mat2cell(currts.',size(currts,2),ones(size(currts,1),1)).'}];
                before_motif_id=[before_motif_id;"c"+cidx];
                before_motif_seq=[before_motif_seq;{repmat({chain_replay.meta{cidx,2}.'},nnz(before_task),1)}];
            end

            currts=chain_replay.ts{cidx}(after_task,:);
            if ~isempty(currts)
                after_motif_spk=[after_motif_spk;{mat2cell(currts.',size(currts,2),ones(size(currts,1),1)).'}];
                after_motif_id=[after_motif_id;"c"+cidx];
                after_motif_seq=[after_motif_seq;{repmat({chain_replay.meta{cidx,2}.'},nnz(after_task),1)}];
            end
        end

        %% plot
        delay_allspk=cell2mat(cellfun(@(x) cell2mat(x),delay_motif_spk,'UniformOutput',false));
        delay_allid=cell2mat(cellfun(@(x) cell2mat(x),delay_motif_seq,'UniformOutput',false));
        delay_uspkts=unique([delay_allspk,delay_allid],'rows');
        sums.delay=[sums.delay;size(delay_uspkts,1),pref_delay_sec];
        iti_allspk=cell2mat(cellfun(@(x) cell2mat(x),iti_motif_spk,'UniformOutput',false));
        iti_allid=cell2mat(cellfun(@(x) cell2mat(x),iti_motif_seq,'UniformOutput',false));
        iti_uspkts=unique([iti_allspk,iti_allid],'rows');
        sums.iti=[sums.iti;size(iti_uspkts,1),pref_succeed_iti_sec];

        before_allspk=cell2mat(cellfun(@(x) cell2mat(x),before_motif_spk,'UniformOutput',false));
        before_allid=cell2mat(cellfun(@(x) cell2mat(x),before_motif_seq,'UniformOutput',false));
        before_uspkts=unique([before_allspk,before_allid],'rows');
        sums.before=[sums.before;size(before_uspkts,1),before_sec];

        after_allspk=cell2mat(cellfun(@(x) cell2mat(x),after_motif_spk,'UniformOutput',false));
        after_allid=cell2mat(cellfun(@(x) cell2mat(x),after_motif_seq,'UniformOutput',false));
        after_uspkts=unique([after_allspk,after_allid],'rows');
        sums.after=[sums.after;size(after_uspkts,1),after_sec];
    end
end

dall=sums.delay(:,1)./sums.delay(:,2);
dhat=median(dall);
dci=bootci(1000,@(x) median(x), dall);
iall=sums.iti(:,1)./sums.iti(:,2);
ihat=median(iall);
ici=bootci(1000,@(x) median(x), iall);
ball=sums.before(:,1)./sums.before(:,2);
bhat=median(ball);
bci=bootci(1000,@(x) median(x), ball);
aall=sums.after(:,1)./sums.after(:,2);
ahat=median(aall);
aci=bootci(1000,@(x) median(x), aall);

pi=signrank(dall,iall);
pb=signrank(dall,ball);
pa=signrank(dall,aall);
mm=[dhat,ihat,bhat,ahat];
ci=[dci,ici,bci,aci];
figure()
hold on
bar(mm,'FaceColor','k')
errorbar(1:4,mm,ci(1,:)-mm,ci(2,:)-mm,'k.')
title(sprintf('p%.4f',pi,pb,pa));
ylabel('Unique motif spike frequencey (Hz)')
xlim([0.25,4.75])
set(gca,'XTick',1:4,'XTickLabel',{'Delay','ITI','Before','After'})


end





