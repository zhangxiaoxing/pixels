% TODO: confirm inter-loop-transition with spkTS-cid-map
% TODO: join selective and non-selective loops

% We found that larger proportions of spikes were associated with 
% such activity loops of same-memory neurons than that of non-memory
% neurons (fig. Sx). For all the recorded neurons, we observed xx% of
% spikes belong to some activity loops, which was much higher than
% the shuffled level (). 

function pstats=rings_time_constant(sums_all,opt)
arguments
    sums_all
    opt.load_file (1,1) logical = false
end
% function cong_dir=rings_su_tcom_order(sums_all) % modified from
% persistent sums_all
% if isempty(sums_all)
%     load(fullfile('bzdata','sums_ring_stats_all.mat'));
% end
if ~opt.load_file
    pstats=struct();
    pstats.congru=struct();
    pstats.nonmem=struct();

    su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
    wrs_mux_meta=ephys.get_wrs_mux_meta();

    rstats=cell(0,11);
    for rsize=3:5
        one_rsize=sums_all{rsize-2};
        curr_sess=-1;
        for ridx=1:size(one_rsize,1)
            if curr_sess~=one_rsize{ridx,1}
                curr_sess=one_rsize{ridx,1};
                sesscid=su_meta.allcid(su_meta.sess==curr_sess);
                sesswaveid=wrs_mux_meta.wave_id(su_meta.sess==curr_sess);
                sess_wave_map=containers.Map(num2cell(sesscid),num2cell(sesswaveid));
            end
            curr_waveid=cell2mat(sess_wave_map.values(num2cell(one_rsize{ridx,3})));
            %         if (~all(ismember(curr_part,{'CH','BS'}),'all'))
            %             continue
            %         end
            [rwid,seltype]=bz.rings.ring_wave_type(curr_waveid);

            if ~strcmp(rwid,'congru') && ~strcmp(rwid,'nonmem')
                continue
            end
            rstats=[rstats;one_rsize(ridx,:),curr_waveid,seltype,rsize,rwid];
        end
    end
%     keyboard()
    rstats=rstats(cell2mat(rstats(:,6))>0.1 & cellfun(@(x) numel(unique(x)),rstats(:,3))==cell2mat(rstats(:,9)),:);
    usess=unique(cell2mat(rstats(:,1)));

    for sessid=usess.'
        [spkID,spkTS,trials,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true);
        rid=find(cell2mat(rstats(:,1))==sessid);
        for ri=reshape(rid,1,[])
            ts_id=[];
            cids=rstats{ri,3};
            disp({sessid,ri});

            per_cid_spk_cnt=zeros(size(cids));
            for in_ring_pos=1:numel(cids) % TODO, 1:rsize
                one_ring_sel=spkID==cids(in_ring_pos);
                per_cid_spk_cnt(in_ring_pos)=nnz(one_ring_sel);
                rawts=spkTS(one_ring_sel);

                ft_sel=strcmp(FT_SPIKE.label,num2str(cids(in_ring_pos)));
                ft_ts=FT_SPIKE.timestamp{ft_sel};
                ft_trl_time=FT_SPIKE.time{ft_sel};
                ft_trl=FT_SPIKE.trial{ft_sel};

                [tt,tspos]=ismember(ft_ts,rawts);
                ext_time=repmat(-realmax,numel(rawts),1);
                ext_time(tspos)=ft_trl_time;

                ext_trl=repmat(-realmax,numel(rawts),1);
                ext_trl(tspos)=ft_trl;

                ts_id=cat(1,ts_id,[rawts,... % 1
                    repmat(cids(in_ring_pos),per_cid_spk_cnt(in_ring_pos),1),...  % 2
                    ones(per_cid_spk_cnt(in_ring_pos),1)*in_ring_pos,...  % 3
                    ext_time,...  % 4
                    ext_trl]); % 5
            end
            ts_id=sortrows(ts_id,1);
            ts_id=[ts_id,full(rstats{ri,5}.tags)]; % join TS, ring tag % 6
            
            uidtag=sprintf('s%dr%d',sessid,ri);
            if strcmp(rstats{ri,8},'none')
                d3=find(ismember(trials(:,5),[4 8]) & trials(:,8)==3 & all(trials(:,9:10)~=0,2));
                d6=find(ismember(trials(:,5),[4 8]) & trials(:,8)==6 & all(trials(:,9:10)~=0,2));
                tssel=(ismember(ts_id(:,5),d3) & ts_id(:,4)>=1 & ts_id(:,4)<4) ...
                    | (ismember(ts_id(:,5),d6) & ts_id(:,4)>=1 & ts_id(:,4)<7);
                rsums=ts_id(tssel,:);
                pstats.nonmem.(uidtag).ts_id=rsums;
                pstats.nonmem.(uidtag).rstats=rstats(ri,[1:3,7:10]);
                pstats.nonmem.(uidtag).trials=trials;
            else
                [pref3,pref6]=bz.rings.preferred_trials_rings(rstats{ri,7},trials);
                tssel=(ismember(ts_id(:,5),pref3) & ts_id(:,4)>=1 & ts_id(:,4)<4)...
                    | (ismember(ts_id(:,5),pref6) & ts_id(:,4)>=1 & ts_id(:,4)<7);
                rsums=ts_id(tssel,:);
                pstats.congru.(uidtag).ts_id=rsums;
                pstats.congru.(uidtag).rstats=rstats(ri,[1:3,7:10]);
                pstats.congru.(uidtag).trials=trials;
            end
        end
    end
%     keyboard();
    blame=vcs.blame();
    save(fullfile('bzdata','rings_spike_trial_tagged.mat'),'pstats','blame','-v7.3')
else
    load(fullfile('bzdata','rings_spike_trial_tagged.mat'),'pstats')
end
plotStats(pstats)
end


function plotStats(pstats)
fn=fieldnames(pstats.congru);
sess=str2double(string(regexp(fn,'(?<=s)\d+(?=r.*)','match')));
waveid=[3 6];

figure();pie(categorical(sess));
for sessid=33%[100,18,33]
    sess_ring=fn(sess==sessid);
    su_per_region=cell(0);
    for onefn=reshape(sess_ring,1,[])
        if any(ismember(waveid,pstats.congru.(onefn{1}).rstats{4}),'all')
            su_per_region(end+1)=pstats.congru.(onefn{1}).rstats(3);
        end
    end
    usus=unique([su_per_region{:}]);
%     sumap=nan(numel(sess_ring),numel(usus));
    sumap=[];
    for onefn=reshape(sess_ring,1,[])
         sumap=[sumap;ismember(usus,pstats.congru.(onefn{1}).rstats{3})];
    end
    Z=linkage(sumap);
    ZV=linkage(sumap.');
    figure();
    tiledlayout(1,3);
    nexttile()
    [H,T,outperm]=dendrogram(Z,0,'Orientation','left','ColorThreshold','default');
    nexttile();
    [HV,TV,outpermV]=dendrogram(ZV,0,'Orientation','left','ColorThreshold','default');
    nexttile()
    sumapV=sumap(:,outpermV);
    imagesc(~sumapV(outperm,:))
    colormap('gray')
    set(gca,'YDir','normal')
    
    trials=pstats.congru.(sess_ring{1}).trials;
    if ismember(3,waveid)
        trial_sel=find(trials(:,5)==8 & trials(:,8)==3 & all(trials(:,9:10)>0,2));
    elseif ismember(1,waveid)
        trial_sel=find(trials(:,5)==4 & trials(:,8)==3 & all(trials(:,9:10)>0,2));
    end
    ring_spikes=cell(0,numel(sess_ring));
    ring_cid=cell(0,numel(sess_ring));
    for t=reshape(trial_sel,1,[])
        onetrial=cell(1,0);
        onecid=cell(1,0);
        for onefn=reshape(sess_ring(outperm),1,[])
            ts_id=pstats.congru.(onefn{1}).ts_id();
            tsel=ts_id(:,5)==t & ts_id(:,6)~=0;
            time_trial={reshape(ts_id(tsel,4),1,[])};
            cid_trial={reshape(ts_id(tsel,2),1,[])};
            onetrial(1,end+1)=time_trial;
            onecid(1,end+1)=cid_trial;
        end
        ring_spikes=[ring_spikes;onetrial];
        ring_cid=[ring_cid;onecid];
    end
    %% showcase
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
            if ts(ii)-ts(ii-1)<0.01
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
   ylabel('Loops #')
   edges = find(diff([0,covered,0]==1));
   onset = edges(1:2:end-1);  % Start indices
   run_length = edges(2:2:end)-onset;  % Consecutive ones counts
   ymax=max(ylim());
   for ii=1:numel(onset)
       plot([onset(ii),onset(ii)+run_length(ii)]./1000,[ymax,ymax],'r-','LineWidth',2);
   end
end

%% statistics
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
blame=vcs.blame();
save(fullfile('bzdata','hebbian_ring.mat'),'stats','blame')
C=struct2cell(stats).';
expd=[C{:}];
expdd=[expd{:}];
figure();
histogram(expdd,'Normalization','probability');
set(gca(),'YScale','log','XScale','log')
ylim([1e-3,1])
end



