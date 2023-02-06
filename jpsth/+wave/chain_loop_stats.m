% interate sessions
% load single spike chain data
% load burst spike chain data
function single_multi_chains()
ssc=load('chain_tag.mat','out');
perchaindur=struct();
[perchaindur.d6.size,perchaindur.d6.dur,perchaindur.d3.size,perchaindur.d3.dur,perchaindur.d6.int,perchaindur.d3.int]=deal([]);
for dur=reshape(fieldnames(out),1,[])
    perchaindur=struct();
    [perchaindur.size,perchaindur.dur]=deal([]);
    for wv=reshape(fieldnames(out.(dur{1})),1,[])
        for lp=reshape(fieldnames(out.(dur{1}).(wv{1})),1,[])
            perchaindur.size=[perchaindur.size;size(out.(dur{1}).(wv{1}).(lp{1}).ts,2)];
            perchaindur.dur=[perchaindur.dur;{diff(out.(dur{1}).(wv{1}).(lp{1}).ts(:,[1,end]),1,2)./30}];
        end
    end
    statss.("d"+dur)=perchaindur;
end
singlehist=histcounts([cell2mat(statss.dd3.dur);cell2mat(statss.dd6.dur)],0:10:100,'Normalization','pdf');


%% multi spike
msc=load('chain_sust_tag_600.mat','out');
perchaindur=struct();
[perchaindur.d6.size,perchaindur.d6.dur,perchaindur.d3.size,perchaindur.d3.dur,perchaindur.d6.int,perchaindur.d3.int]=deal([]);
for dur=reshape(fieldnames(out),1,[])
%     [perchaindur.size,perchaindur.dur]=deal([]);
    for wv=reshape(fieldnames(out.(dur{1})),1,[])
        for lp=reshape(fieldnames(out.(dur{1}).(wv{1})),1,[])
%             keyboard();
            perchaindur.(dur{1}).size=[perchaindur.(dur{1}).size,cellfun(@(x) size(x,1),out.(dur{1}).(wv{1}).(lp{1}).ts)];
            perchaindur.(dur{1}).dur=[perchaindur.(dur{1}).dur,cellfun(@(x) diff(x([1,end],3),1,1),out.(dur{1}).(wv{1}).(lp{1}).ts)./30];
            perchaindur.(dur{1}).int=[perchaindur.(dur{1}).int,cell2mat(cellfun(@(x) diff(x(:,3),1,1).',out.(dur{1}).(wv{1}).(lp{1}).ts,'UniformOutput',false))./30];
        end
    end
    statsm.("d"+dur)=perchaindur;
end

multihist=histcounts([statsm.dd3.d3.dur,statsm.dd6.d6.dur],[0:10:100,200:100:600],'Normalization','pdf');
end




%% load single spike loop data
function single_spike_composite()
stats=struct();
waveids={[1 5],[2 5],[3 6],[4 6],[1 7],[2 8],[3 7],[4 8]};
for wid=1:4
    waveid=waveids{wid};
    for sessid=reshape(unique(sess),1,[])
        allsessfn=fn(sess==sessid);
        sess_ring=cell(0);
        for onefn=reshape(allsessfn,1,[])
            if all(ismember(pstats.congru.(onefn{1}).rstats{1,4},waveid),'all')
                sess_ring=[sess_ring;onefn];
            end
        end
        if size(sess_ring,1)<2
            continue
        end
        trials=pstats.congru.(sess_ring{1}).trials;
        if ismember(3,waveid)
            trial_sel=find(trials(:,5)==8 & trials(:,8)==3 & all(trials(:,9:10)>0,2));
        elseif ismember(1,waveid)
            trial_sel=find(trials(:,5)==4 & trials(:,8)==3 & all(trials(:,9:10)>0,2));
        elseif ismember(2,waveid)
            trial_sel=find(trials(:,5)==4 & trials(:,8)==6 & all(trials(:,9:10)>0,2));
        elseif ismember(4,waveid)
            trial_sel=find(trials(:,5)==8 & trials(:,8)==6 & all(trials(:,9:10)>0,2));
        end
        ring_spikes=cell(0,numel(sess_ring));
%         ring_cid=cell(0,numel(sess_ring));
        for t=reshape(trial_sel,1,[])
            onetrial=cell(1,0);
%             onecid=cell(1,0);
            for onefn=reshape(sess_ring,1,[])
                ts_id=pstats.congru.(onefn{1}).ts_id();
                tsel=ts_id(:,5)==t & ts_id(:,6)~=0;
                time_trial={reshape(ts_id(tsel,4),1,[])};
%                 cid_trial={reshape(ts_id(tsel,2),1,[])};
                onetrial(1,end+1)=time_trial;
%                 onecid(1,end+1)=cid_trial;
            end
            ring_spikes=[ring_spikes;onetrial];
%             ring_cid=[ring_cid;onecid];
        end
        sfn=sprintf('w%ds%d',wid,sessid);
        stats.(sfn)=cell(0);
        for tt=1:size(ring_spikes,1)
            covered=zeros(1,3000);
            for rridx=1:size(ring_spikes,2)
                if isempty(ring_spikes{tt,rridx})
                    continue
                end
                ts=ring_spikes{tt,rridx}-1;
                for ii=2:numel(ts)
                    if ts(ii)-ts(ii-1)<0.01
                        onset=ceil(ts(ii-1)*1000);
                        offset=ceil(ts(ii)*1000);
                        covered(onset:offset)=1;
                    end
                end
            end

            edges = find(diff([0,covered,0]==1));
            onset = edges(1:2:end-1);  % Start indices
            run_length = edges(2:2:end)-onset;  % Consecutive ones counts
            stats.(sfn)=[stats.(sfn),{reshape(run_length,1,[])}];
        end

        %    max([stats{:}])
    end
end
end


%% load burst spike loop data
function burst_spike_composite()
dbfile=fullfile("bzdata","rings_wave_burst_600.db");
conn=sqlite(dbfile,"readonly");
keys=table2array(conn.fetch("SELECT name FROM sqlite_master WHERE type='table'"));
usess=str2double(unique(regexp(keys,'(?<=s)\d{1,3}(?=r[3-5])','match','once')));
% statistics
stats=struct();
for sessid=[4 9 14 18]%reshape(unique(usess),1,[])
    [~,~,trials,~,~,~]=ephys.getSPKID_TS(sessid,'keep_trial',false,'skip_spike',true);
    sessfn=contains(keys,"s"+num2str(sessid)+"r") & endsWith(keys,"_ts");
    for wid=1:4
        switch wid
            case 1
                wsel=startsWith(keys,"d3s1d3") | startsWith(keys,"d3olf_s1") | startsWith(keys,"d3dur_d3");
                trial_sel=find(trials(:,5)==4 & trials(:,8)==3 & all(trials(:,9:10)>0,2));
            case 2
                wsel=startsWith(keys,"d6s1d6") | startsWith(keys,"d6olf_s1") | startsWith(keys,"d6dur_d6");
                trial_sel=find(trials(:,5)==4 & trials(:,8)==6 & all(trials(:,9:10)>0,2));
            case 3
                wsel=startsWith(keys,"d3s2d3") | startsWith(keys,"d3olf_s2") | startsWith(keys,"d3dur_d3");
                trial_sel=find(trials(:,5)==8 & trials(:,8)==3 & all(trials(:,9:10)>0,2));
            case 4
                wsel=startsWith(keys,"d6s2d6") | startsWith(keys,"d6olf_s2") | startsWith(keys,"d6dur_d6");
                trial_sel=find(trials(:,5)==8 & trials(:,8)==6 & all(trials(:,9:10)>0,2));
        end
        wave_key=keys(sessfn & wsel);
        if size(wave_key,1)<2
            continue
        end
        ts_id=struct();
        burst_ts=struct();
        for onefn=reshape(wave_key,1,[])
            ts_id.(onefn)=table2array(sqlread(conn,replace(onefn,"_ts","_tsid")));
            burst_ts.(onefn)=table2array(sqlread(conn,onefn));
        end

        ring_spikes=cell(0,numel(wave_key));
        %         ring_cid=cell(0,numel(sess_ring));
        for t=reshape(trial_sel,1,[])
            onetrial=cell(1,0);
            %             onecid=cell(1,0);
            for onefn=reshape(wave_key,1,[])
                disp({t,onefn})
                tsel=false(size(ts_id.(onefn),1),1);
                for suidx=1:ts_id.(onefn)(end,3)
                    tsel(ts_id.(onefn)(:,5)==t & ts_id.(onefn)(:,3)==suidx & ismember(ts_id.(onefn)(:,1),burst_ts.(onefn)(burst_ts.(onefn)(:,2)==suidx,4)))=true;
                end
                time_trial={sortrows(ts_id.(onefn)(tsel,[3,4]),2)};
                %                 cid_trial={reshape(ts_id(tsel,2),1,[])};
                onetrial(1,end+1)=time_trial;
                %                 onecid(1,end+1)=cid_trial;
            end
            ring_spikes=[ring_spikes;onetrial]; % n Trial x trial
            %             ring_cid=[ring_cid;onecid];
        end
        sfn=sprintf('w%ds%d',wid,sessid);
        stats.(sfn)=cell(0);
        for tt=1:size(ring_spikes,1)
            if wid==1 || wid==3
                covered=zeros(1,3000);
            else
                covered=zeros(1,6000);
            end
            for rridx=1:size(ring_spikes,2)
                if isempty(ring_spikes{tt,rridx})
                    continue
                end
                ts=ring_spikes{tt,rridx}-1;
                for ii=2:size(ts,1)
                    if (ts(ii,1)~=ts(ii-1,1) && ts(ii,2)-ts(ii-1,2)<0.01)...
                            || (ts(ii,1)==ts(ii-1,1) && ts(ii,2)-ts(ii-1,2)<0.02)
                        onset=ceil(ts(ii-1,2)*1000);
                        offset=ceil(ts(ii,2)*1000);
                        covered(onset:offset)=1;
                    end
                end
            end

            edges = find(diff([0,covered,0]==1));
            onset = edges(1:2:end-1);  % Start indices
            run_length = edges(2:2:end)-onset;  % Consecutive ones counts
            stats.(sfn)=[stats.(sfn),{reshape(run_length,1,[])}];
        end
    end
end
conn.close();
end


% union all neuron spikes
% tag chained-loops spikes
% time constant stats
