% motifs PSTH

sessid=100;
maxid=4;
[mmax,ttrial]=wave.chain_loop_SC_spk(false,sessid,maxid,'skip_plot',true);
% [~,~,trials,~,~,~]=ephys.getSPKID_TS(sessid,'keep_trial',false);
% 
% if ismember(3,waveid)
%     trial_sel=find(trials(:,5)==8 & trials(:,8)==3 & all(trials(:,9:10)>0,2));
% elseif ismember(1,waveid)
%     trial_sel=find(trials(:,5)==4 & trials(:,8)==3 & all(trials(:,9:10)>0,2));
% end
ring_spikes=cell(0,numel(sess_ring));
for t=reshape(trial_sel,1,[])
    onetrial=cell(1,0);
    onecid=cell(1,0);
    for onefn=reshape(sess_ring(outperm),1,[])
        ts_id=pstats.congru.(onefn{1}).ts_id();
        tsel=ts_id(:,5)==t & ts_id(:,6)~=0; % practically the new ts segment in mmax etc
        time_trial={reshape(ts_id(tsel,4),1,[])};
        cid_trial={reshape(ts_id(tsel,2),1,[])};
        onetrial(1,end+1)=time_trial;
        onecid(1,end+1)=cid_trial;
    end
    ring_spikes=[ring_spikes;onetrial];
end
%% hebbian composite loop showcase
totalspk=arrayfun(@(x) numel([ring_spikes{x,:}]),1:size(ring_spikes,1));
[totalspk,sidx]=sort(totalspk,'descend');

tt=sidx(1);
covered=zeros(1,3000);
%         join_spk=[];
figure()
hold on
for rridx=1:size(ring_spikes,2)
    if isempty(ring_spikes{tt,rridx})
        continue
    end
    ts=ring_spikes{tt,rridx}-1;
    scatter(ts,rridx,'k','|')
    for ii=2:numel(ts)
        if ts(ii)-ts(ii-1)<0.01005
            onset=floor(ts(ii-1)*1000);
            offset=ceil(ts(ii)*1000);
            %                     disp([onset,offset]);
            covered(onset:offset)=1;
        end
    end
    %             join_spk=[join_spk,[repmat(rridx,1,numel(ring_spikes{tt,rridx}));ring_cid{tt,rridx};ring_spikes{tt,rridx}-1]];
end

xlim([0,3])
xlabel('Time (s)')
ylim([0.5,size(ring_spikes,2)+0.5])
%% currently working on this
ylabel('Loops #')
edges = find(diff([0,covered,0]==1));
onset = edges(1:2:end-1);  % Start indices
run_length = edges(2:2:end)-onset;  % Consecutive ones counts
ymax=max(ylim());
for ii=1:numel(onset)
    plot([onset(ii),onset(ii)+run_length(ii)]./1000,[ymax,ymax],'r-','LineWidth',2);
end
