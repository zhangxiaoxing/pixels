% interate sessions
% load single spike chain data
% load burst spike chain data
% TODO: probe memory, remove unused sessions if necessary

%% single spike chain
sschain=load('chain_tag.mat','out');
keys=[struct2cell(structfun(@(x) fieldnames(x), sschain.out.d6, 'UniformOutput', false));...
    struct2cell(structfun(@(x) fieldnames(x), sschain.out.d3, 'UniformOutput', false))];
keys=vertcat(keys{:});
ssc_sess=unique(str2double(regexp(keys,'(?<=s)\d{1,3}(?=c)','match','once')));
%% multi spike
bschain=load('chain_sust_tag_600.mat','out');
keys=[struct2cell(structfun(@(x) fieldnames(x), bschain.out.d6, 'UniformOutput', false));...
    struct2cell(structfun(@(x) fieldnames(x), bschain.out.d3, 'UniformOutput', false))];
keys=vertcat(keys{:});
bsc_sess=unique(str2double(regexp(keys,'(?<=s)\d{1,3}(?=c)','match','once')));
%% single spike loop
load(fullfile('bzdata','rings_spike_trial_tagged.mat'),'pstats');
pstats=rmfield(pstats,"nonmem");
ssl_sess=unique(str2double(regexp(fieldnames(pstats.congru),'(?<=s)\d{1,3}(?=r)','match','once')));
%% burst spike loop, keys only
dbfile=fullfile("bzdata","rings_wave_burst_600.db");
conn=sqlite(dbfile,"readonly");
bslkeys=table2array(conn.fetch("SELECT name FROM sqlite_master WHERE type='table'"));
close(conn);
bsl_sess=unique(str2double(regexp(bslkeys,'(?<=s)\d{1,3}(?=r)','match','once')));

usess=intersect(intersect(intersect(ssc_sess,bsc_sess),ssl_sess),bsl_sess);

% single spk chn:1, burst spk chn:2, single spk loop:4, burst spk loop:8
for sessid=18
    ucid=[];

    [~,~,trials,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true);
    FT_SPIKE.lc_tag=cell(size(FT_SPIKE.timestamp));
    for cidx=1:numel(FT_SPIKE.timestamp)
        FT_SPIKE.lc_tag{cidx}=zeros(size(FT_SPIKE.timestamp{cidx}));
    end

    %% single spike chain
    for dur=["d6","d3"]
        for wid=reshape(fieldnames(sschain.out.(dur)),1,[])
            for cc=reshape(fieldnames(sschain.out.(dur).(wid{1})),1,[])
                if ~startsWith(cc{1},['s',num2str(sessid),'c'])
                    continue
                end
                onechain=sschain.out.(dur).(wid{1}).(cc{1});
                for cidx=1:size(onechain.ts,2)
                    cid=onechain.meta{1}(cidx);
                    ucid=[ucid,cid];
                    cidsel=find(strcmp(FT_SPIKE.label,num2str(cid)));
                    totag=onechain.ts(:,cidx);
                    [ism,totagidx]=ismember(totag,FT_SPIKE.timestamp{cidsel});
                    FT_SPIKE.lc_tag{cidsel}(totagidx)=FT_SPIKE.lc_tag{cidsel}(totagidx)+1;
                end
            end
        end
    end
    clear sschain
%% multi spike chain
    for dur=["d6","d3"]
        for wid=reshape(fieldnames(bschain.out.(dur)),1,[])
            for cc=reshape(fieldnames(bschain.out.(dur).(wid{1})),1,[])
                if ~startsWith(cc{1},['s',num2str(sessid),'c'])
                    continue
                end
                onechain=bschain.out.(dur).(wid{1}).(cc{1});
                for cidx=1:size(onechain.meta{1},2)
                    cid=onechain.meta{1}(cidx);
                    ucid=[ucid,cid];
                    cidsel=find(strcmp(FT_SPIKE.label,num2str(cid)));
                    for bscid=1:numel(onechain.ts)
                        totag=onechain.ts{bscid}(onechain.ts{bscid}(:,1)==cidx,3);
                        [ism,totagidx]=ismember(totag,FT_SPIKE.timestamp{cidsel});
                        FT_SPIKE.lc_tag{cidsel}(totagidx)=FT_SPIKE.lc_tag{cidsel}(totagidx)+2;
                    end
                end
            end
        end
    end
    clear bschain

    %% single spike loop
    
    for cc=reshape(fieldnames(pstats.congru),1,[])
        if ~startsWith(cc{1},['s',num2str(sessid),'r'])
            continue
        end
        onechain=pstats.congru.(cc{1});
        for cidx=1:size(onechain.rstats{3},2)
            cid=onechain.rstats{3}(cidx);
            ucid=[ucid,cid];
            cidsel=find(strcmp(FT_SPIKE.label,num2str(cid)));
            totag=onechain.ts_id(onechain.ts_id(:,6)>0 & onechain.ts_id(:,2)==cid,1);
            [ism,totagidx]=ismember(totag,FT_SPIKE.timestamp{cidsel});
            FT_SPIKE.lc_tag{cidsel}(totagidx)=FT_SPIKE.lc_tag{cidsel}(totagidx)+4;
        end
    end
    clear pstats
    %% burst spike loop, keys only
    dbfile=fullfile("bzdata","rings_wave_burst_600.db");
    conn=sqlite(dbfile,"readonly");
    for cc=reshape(bslkeys,1,[])
        if ~contains(cc,['s',num2str(sessid),'r']) || ~endsWith(cc,'_ts')
            continue
        end
       
        onechain=table2array(conn.sqlread(cc));
        chainmeta=table2array(conn.sqlread(replace(cc,'_ts','_meta')));
        for cidx=1:numel(chainmeta)
            cid=chainmeta(cidx);
            ucid=[ucid,cid];
            cidsel=find(strcmp(FT_SPIKE.label,num2str(cid)));
            totag=onechain(onechain(:,2)==cidx,4);
            [ism,totagidx]=ismember(totag,FT_SPIKE.timestamp{cidsel});
            FT_SPIKE.lc_tag{cidsel}(totagidx)=FT_SPIKE.lc_tag{cidsel}(totagidx)+8;
        end
    end
    conn.close();
end
function quickdirty(FT_SPIKE)
total=0;
tagged=0;
ts_tag=[];
for ii=1:numel(FT_SPIKE.lc_tag)
    if any(FT_SPIKE.lc_tag{ii}>0,'all')
        total=total+numel(FT_SPIKE.lc_tag{ii});
        tagged=tagged+nnz(FT_SPIKE.lc_tag{ii}>0);
        ts_tag=[ts_tag;repmat(ii,numel(FT_SPIKE.timestamp{ii}),1),FT_SPIKE.timestamp{ii}.',FT_SPIKE.time{ii}.',FT_SPIKE.trial{ii}.',FT_SPIKE.lc_tag{ii}.'];
    end
end
disp([total,tagged,tagged./total])
end








%% =============No function below this line================================

function chains()

%% single spike chain
sschain=load('chain_tag.mat','out');
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

singlechainhist=histcounts([cell2mat(statss.dd3.dur);cell2mat(statss.dd6.dur)],0:10:100,'Normalization','pdf');


%% multi spike chain
bschain=load('chain_sust_tag_600.mat','out');
%%
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


multichainhist=histcounts([statsm.dd3.d3.dur,statsm.dd6.d6.dur],[0:10:100,200:100:600],'Normalization','pdf')
end

%% load single spike loop data
function single_spike_composite()
%%
load(fullfile('bzdata','rings_spike_trial_tagged.mat'),'pstats')
pstats=rmfield(pstats,"nonmem");

%%
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




function db_test()
dbfile=fullfile("bzdata","rings_wave_burst_600.db");
conn=sqlite(dbfile,"readonly");
keys=table2array(conn.fetch("SELECT name FROM sqlite_master WHERE type='table'"));
keys(1)
conn.fetch("SELECT COUNT(*) FROM "+keys(1));
end