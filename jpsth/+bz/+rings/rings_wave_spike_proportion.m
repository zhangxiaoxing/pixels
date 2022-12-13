% We found that larger proportions of spikes were associated with 
% such activity loops of same-memory neurons than that of non-memory
% neurons (fig. Sx). For all the recorded neurons, we observed xx% of
% spikes belong to some activity loops, which was much higher than
% the shuffled level (). 

function pstats=rings_wave_spike_proportion(sums_all,opt)
arguments
    sums_all
    opt.load_file (1,1) logical = true
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
    com_map=wave.get_pct_com_map(wrs_mux_meta,'curve',true,'early_smooth',false);

    rstats=cell(0,11);
    for rsize=3
        one_rsize=sums_all{rsize-2};
        curr_sess=-1;
        for ridx=1:size(one_rsize,1)
            if curr_sess~=one_rsize{ridx,1}
                curr_sess=one_rsize{ridx,1};
                sesscid=su_meta.allcid(su_meta.sess==curr_sess);
                sesswaveid=wrs_mux_meta.wave_id(su_meta.sess==curr_sess);
                sess_wave_map=containers.Map(num2cell(sesscid),num2cell(sesswaveid));
                sess_su_com_map=com_map.(['s',num2str(curr_sess)]);
            end
            curr_waveid=cell2mat(sess_wave_map.values(num2cell(one_rsize{ridx,3})));
            %         if (~all(ismember(curr_part,{'CH','BS'}),'all'))
            %             continue
            %         end
            [rwid,seltype]=bz.rings.ring_wave_type(curr_waveid);

            if ~strcmp(rwid,'congru') && ~strcmp(rwid,'nonmem')
                continue
            end
            rstats=[rstats;one_rsize(ridx,:),curr_waveid,seltype,rsize];
        end
    end

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
            
            if strcmp(rstats{ri,8},'none')
                d3=find(ismember(trials(:,5),[4 8]) & trials(:,8)==3 & all(trials(:,9:10)~=0,2));
                d6=find(ismember(trials(:,5),[4 8]) & trials(:,8)==6 & all(trials(:,9:10)~=0,2));
                tssel=(ismember(ts_id(:,5),d3) & ts_id(:,4)>=1 & ts_id(:,4)<4) ...
                    | (ismember(ts_id(:,5),d6) & ts_id(:,4)>=1 & ts_id(:,4)<7);
                rsums=ts_id(tssel,:);
                for uuid=rstats{ri,3}
                    uidtag=sprintf('s%du%d',sessid,uuid);
                    if ~isfield(pstats.nonmem,uidtag)
                        pstats.nonmem.(uidtag)=[rsums(rsums(:,2)==uuid,1),zeros(nnz(rsums(:,2)==uuid),1)];
                    end
                    pstats.nonmem.(uidtag)(rsums(:,2)==uuid & rsums(:,6)>0,2)=1;
                end
            else
                [pref3,pref6]=bz.rings.preferred_trials_rings(rstats{ri,7},trials);
                tssel=(ismember(ts_id(:,5),pref3) & ts_id(:,4)>=1 & ts_id(:,4)<4)...
                    | (ismember(ts_id(:,5),pref6) & ts_id(:,4)>=1 & ts_id(:,4)<7);
                rsums=ts_id(tssel,:);
                for uuid=rstats{ri,3}
                    uidtag=sprintf('s%du%d',sessid,uuid);
                    if ~isfield(pstats.congru,uidtag)
                        pstats.congru.(uidtag)=[rsums(rsums(:,2)==uuid,1),zeros(nnz(rsums(:,2)==uuid),1)];
                    end
                    pstats.congru.(uidtag)(rsums(:,2)==uuid & rsums(:,6)>0,2)=1;
                end
            end
        end
    end
    keyboard();
    blame=vcs.blame();
    save(fullfile('bzdata','rings_wave_spike_proportion.mat'),'pstats','blame')
else
    load(fullfile('bzdata','rings_wave_spike_proportion.mat'),'pstats')
end
%     sum_stats(pstats);
end

function sum_stats(pstats)
fns=fieldnames(pstats.congru);
ratio=[];
for fn=reshape(fns,1,[])
    ratio=[ratio;mean(pstats.congru.(fn{1})(:,6))];
end
[max(ratio),mean(ratio),std(ratio)]
end


% end

