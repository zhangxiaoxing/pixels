% We found that larger proportions of spikes were associated with
% such activity loops of same-memory neurons than that of non-memory
% neurons (fig. Sx). For all the recorded neurons, we observed xx% of
% spikes belong to some activity loops, which was much higher than
% the shuffled level ().
classdef rings_time_constant <handle
    methods (Static)
        function [pstats,tag_compress]=stats(sums_all,sel_meta,opt)
            arguments
                sums_all  % w/ loop tags
                sel_meta
                opt.load_file (1,1) logical = false
                opt.skip_save (1,1) logical = false
                opt.odor_only (1,1) logical = false
                opt.compress (1,1) logical = false
            end
            % function cong_dir=rings_su_tcom_order(sums_all) % modified from
            % persistent sums_all
            % if isempty(sums_all)
            %     load(fullfile('bzdata','sums_ring_stats_all.mat'));
            % end
            tag_compress=struct();
            if ~opt.load_file
                pstats=struct();
                pstats.congru=struct();
                pstats.nonmem=struct();

                su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);

                rstats=cell(0,11);
                for rsize=3:5
                    one_rsize=sums_all{rsize-2};
                    curr_sess=-1;
                    for ridx=1:size(one_rsize,1)
                        if curr_sess~=one_rsize{ridx,1}
                            curr_sess=one_rsize{ridx,1};
                            sesscid=su_meta.allcid(su_meta.sess==curr_sess);
                            sesswaveid=sel_meta.wave_id(su_meta.sess==curr_sess);
                            sess_wave_map=containers.Map(num2cell(sesscid),num2cell(sesswaveid));
                        end
                        curr_waveid=cell2mat(sess_wave_map.values(num2cell(one_rsize{ridx,3})));

                        [rwid,seltype]=bz.rings.ring_wave_type(curr_waveid);

                        if (~strcmp(rwid,'congru') && ~strcmp(rwid,'nonmem'))...
                                ||(opt.odor_only && strcmp(seltype,'dur'))
                            continue
                        end

                        rstats=[rstats;one_rsize(ridx,:),curr_waveid,seltype,rsize,rwid];
                        %  ////////////////^^^^^^^^^^\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                        % {sessidx,ring_id,cids,per_cid_spk_cnt,ring_stats,coact_count./(ts_id(end,1)./30000)}
                    end
                end
                
                rstats=rstats(cell2mat(rstats(:,6))>0.1 & cellfun(@(x) numel(unique(x)),rstats(:,3))==cell2mat(rstats(:,9)),:);
                % ring spike rate > 0.1 Hz & no cross link inside loops
                usess=unique(cell2mat(rstats(:,1)));

                for sessid=usess.'
                    [spkID,spkTS,trials,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true);
                    rid=find(cell2mat(rstats(:,1))==sessid);
                    for ri=reshape(rid,1,[])
                        if opt.compress && (any(ismember(rstats{ri,7},[7 8]),'all') ||(any(rstats{ri,7}==0,'all') && ~all(rstats{ri,7}==0,'all')))
                            % no difference between 3s or 6s
                            continue
                        end

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
                            if opt.compress
                                % TODO: complete output codes
                                % TODO: copied from olf code, WIP
                                onets=arrayfun(@(x) ts_id(ts_id(:,6)==x,1),setdiff(unique(ts_id(:,6)),0),'UniformOutput',false);
                                onemeta=rstats(ri,[1 3 7]);
                                pref_samp='none';
                                % no dur pref
                                tag_compress.none.(pref_samp).(sprintf('s%dr%d',sessid,ri))=...
                                    cell2struct({onets;onemeta;trials},{'ts','meta','trials'});

                                %
                            else
                                d3=find(ismember(trials(:,5),[4 8]) & trials(:,8)==3 & all(trials(:,9:10)~=0,2));
                                d6=find(ismember(trials(:,5),[4 8]) & trials(:,8)==6 & all(trials(:,9:10)~=0,2));
                                % in-delay selection
                                % TODO: optional switch
                                tssel=(ismember(ts_id(:,5),d3) & ts_id(:,4)>=1 & ts_id(:,4)<4) ...
                                    | (ismember(ts_id(:,5),d6) & ts_id(:,4)>=1 & ts_id(:,4)<7);
                                rsums=ts_id(tssel,:);

                                for ii=reshape(setdiff(unique(rsums(:,6)),0),1,[])
                                    if nnz(rsums(:,6)==ii)<=numel(cids)
                                        rsums(rsums(:,6)==ii,6)=0;
                                    end
                                end
                                %<<<<<<<<<<<<<<<<<<
                                pstats.nonmem.(uidtag).ts_id=...
                                    table(uint32(rsums(:,1)),uint16(rsums(:,2)),uint8(rsums(:,3)),rsums(:,4),uint8(rsums(:,5)),uint16(rsums(:,6)),...
                                    'VariableNames',["TS","CID","REL_POS","Time","Trial","Loop_tag"]);

                                pstats.nonmem.(uidtag).rstats=rstats(ri,[1:3,7:10]);
                                pstats.nonmem.(uidtag).trials=trials;
                            end
                        else
                            % in-delay selection
                            if opt.compress
                                % TODO: complete output codes
                                onets=arrayfun(@(x) ts_id(ts_id(:,6)==x,1),setdiff(unique(ts_id(:,6)),0),'UniformOutput',false);
                                onemeta=rstats(ri,[1 3 7]);
                                % 3s 6s or both
                                if all(ismember(rstats{ri,7},[1 2 5]),'all')
                                    pref_samp='olf_s1';
                                elseif all(ismember(rstats{ri,7},[3 4 6]),'all')
                                    pref_samp='olf_s2';
                                else % TODO: verify and remove if unnecessary
                                    warning("Incongruent ring under congruent context")
                                    keyboard()
                                end
                                if all(ismember(rstats{ri,7},[1 3 5 6]),'all') % 3s
                                    tag_compress.d3.(pref_samp).(sprintf('s%dr%d',sessid,ri))=...
                                        cell2struct({onets;onemeta;trials},{'ts','meta','trials'});
                                end
                                if all(ismember(rstats{ri,7},[2 4 5 6]),'all') % 6s
                                    tag_compress.d6.(pref_samp).(sprintf('s%dr%d',sessid,ri))=...
                                        cell2struct({onets;onemeta;trials},{'ts','meta','trials'});
                                end
                            else
                                [pref3,pref6]=bz.rings.preferred_trials_rings(rstats{ri,7},trials);
                                tssel=(ismember(ts_id(:,5),pref3) & ts_id(:,4)>=1 & ts_id(:,4)<4)...
                                    | (ismember(ts_id(:,5),pref6) & ts_id(:,4)>=1 & ts_id(:,4)<7);
                                rsums=ts_id(tssel,:);
                                % remove trial-time cut-off loops
                                for ii=reshape(setdiff(unique(rsums(:,6)),0),1,[])
                                    if nnz(rsums(:,6)==ii)<=numel(cids)
                                        rsums(rsums(:,6)==ii,6)=0;
                                    end
                                end

                                pstats.congru.(uidtag).ts_id=...
                                    table(uint32(rsums(:,1)),uint16(rsums(:,2)),uint8(rsums(:,3)),rsums(:,4),uint8(rsums(:,5)),uint16(rsums(:,6)),...
                                    'VariableNames',["TS","CID","REL_POS","Time","Trial","Loop_tag"]);
                            end
                            %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                            pstats.congru.(uidtag).rstats=rstats(ri,[1:3,7:10]);
                            pstats.congru.(uidtag).trials=trials;
                        end
                    end
                end
                %     keyboard();
                if ~opt.skip_save
                    blame=vcs.blame();
                    save(fullfile('bzdata','rings_spike_trial_tagged.mat'),'pstats','blame','-v7.3')
                end
            else
                load(fullfile('bzdata','rings_spike_trial_tagged.mat'),'pstats','blame')
            end
%             plotStats(pstats)
        end


        function plotStats(pstats)
            fn=fieldnames(pstats.congru); %rings
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
            end

            %% statistics
            stats=struct();
            waveids={[1 5],[2 5],[3 6],[4 6],[1 7],[2 8],[3 7],[4 8]};
            for wid=1:8
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
                    for t=reshape(trial_sel,1,[])
                        onetrial=cell(1,0);
                        for onefn=reshape(sess_ring,1,[])
                            ts_id=pstats.congru.(onefn{1}).ts_id();
                            tsel=ts_id(:,5)==t & ts_id(:,6)~=0;
                            time_trial={reshape(ts_id(tsel,4),1,[])};
                            onetrial(1,end+1)=time_trial;
                        end
                        ring_spikes=[ring_spikes;onetrial]; % trial x ringid cell of trial-time
                    end
                    sfn=sprintf('w%ds%d',wid,sessid);
                    stats.(sfn)=cell(0);
                    for tt=1:size(ring_spikes,1) % #trial
                        covered=zeros(1,1000);
                        for rridx=1:size(ring_spikes,2) % #loops
                            if isempty(ring_spikes{tt,rridx})
                                continue
                            end
                            ts=ring_spikes{tt,rridx}-1;
                            for ii=2:numel(ts)
                                if ts(ii)-ts(ii-1)<0.01005
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
            if ~opt.skip_save
                blame=vcs.blame();
                save(fullfile('bzdata','hebbian_ring.mat'),'stats','blame')
            end
            C=struct2cell(stats).';
            expd=[C{:}];
            expdd=[expd{:}];
            figure();
            histogram(expdd,'Normalization','probability');
            set(gca(),'YScale','log','XScale','log')
            ylim([1e-3,1])
        end


        %%
        function burst_spike_composite()
            dbfile=fullfile("bzdata","rings_wave_burst_iter_600.db");
            conn=sqlite(dbfile,"readonly");
            keys=table2array(conn.fetch("SELECT name FROM sqlite_master WHERE type='table'"));
            usess=str2double(unique(regexp(keys,'(?<=s)\d{1,3}(?=r[3-5])','match','once')));

            % statistics
            stats=struct();

            for sessid=reshape(unique(usess),1,[])
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
                            ts=ring_spikes{tt,rridx}-1+realmin;
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

            blame=vcs.blame();
            save(fullfile('bzdata','composite_burst_ring.mat'),'stats','blame')
            C=struct2cell(stats).';
            expd=[C{:}];
            expdd=[expd{:}];
            figure();
            histcompo=histcounts(expdd,[0:20:200,300:300:1200],'Normalization','pdf');
            plot([10:20:190,250,450:300:1050],histcompo,'-k')
            set(gca(),'YScale','log','XScale','log')
            ylim([1e-6,0.1])
            xlim([10,2000])
        end
    end
end

