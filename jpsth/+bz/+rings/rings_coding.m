ephys.util.dependency('ft',true,'buz',false); %data path and lib path dependency
if ~exist('ring_meta','var')
    load('ring_meta.mat','ring_meta');
end
if ~exist('sums_all','var')
    load(fullfile('bzdata','sums_ring_stats_all.mat'));
end
mtypes=["congru","nonmem"];
regtypes=["cross","within"];
rsizes=3:5;
loops_sums=struct();
for rsize=rsizes
    rstats=sums_all{rsize-2};
    rstats=rstats(cell2mat(rstats(:,6))>0.1 & cellfun(@(x) numel(unique(x)),rstats(:,3))==rsize,:);

    %     if rsize==4
    %         if ~exist('sums4','var')
    %             load(fullfile('bzdata','sums_ring_stats_4.mat'),'sums4');
    %         end
    %         rstats=sums4;
    %     else
    %         if ~exist('sums3','var')
    %             load(fullfile('bzdata','sums_ring_stats_3.mat'),'sums3');
    %         end
    %         rstats=sums3;
    %     end
    
    for mtype=mtypes
        for rtype=regtypes
            rmeta=ring_meta.(mtype).(sprintf('%s_%d',rtype,rsize)).meta;
            loops_sums.(mtype).(sprintf('%s_%d',rtype,rsize))=cell(0);
            %             keyboard()
            usess=unique(rmeta(:,1));
            pstats=struct();
            for sessid=usess.'
                %                 pstats.(mtype).(sprintf('s%d',sessid))=struct();
                [spkID,spkTS,trials,suids,folder]=ephys.getSPKID_TS(sessid);
                rid=find(rmeta(:,1)==sessid);
                for ri=reshape(rid,1,[])
                    %TODO congru vs non mem
                    disp([sessid,ri]);
                    ts_id=[];
                    cids=rmeta(ri,2:end);
                    rsids=find(cell2mat(rstats(:,1))==sessid & all(ismember(cell2mat(rstats(:,3)),cids),2));
                    if isempty(rsids)
                        continue;
                    end
                    per_cid_spk_cnt=cids;
                    for in_ring_pos=1:numel(cids) % TODO, 1:rsize
                        one_ring_sel=spkID==cids(in_ring_pos);
                        per_cid_spk_cnt(in_ring_pos)=nnz(one_ring_sel);
                        ts_id=cat(1,ts_id,[spkTS(one_ring_sel),ones(per_cid_spk_cnt(in_ring_pos),1)*in_ring_pos]);
                    end
                    ts_id=sortrows(ts_id,1);

                    for rsidi=1:numel(rsids)
                        rsid=rsids(rsidi);
                        ts_id=[ts_id,full(rstats{rsid,5}.tags)];
                        %% loops PSTH
                        ts_ring=ts_id(ts_id(:,3)~=0, 1:2);
                        SP=splitapply(@(x) {x}, ts_ring(:,1), ts_ring(:,2));
                        FT_LOOPS=struct();
                        FT_LOOPS.label=arrayfun(@(x) num2str(x),cids,'UniformOutput',false);
                        FT_LOOPS.timestamp=SP;
                        sps=30000;
                        cfg=struct();
                        cfg.trl=[trials(:,1)-3*sps,trials(:,1)+11*sps,zeros(size(trials,1),1)-3*sps,trials];
                        cfg.trlunit='timestamps';
                        cfg.timestampspersecond=sps;
                        FT_LOOPS=ft_spike_maketrials(cfg,FT_LOOPS);
                        FT_LOOPS.sessid=sessid;
                        loops_sums.(mtype).(sprintf('%s_%d',rtype,rsize)){end+1}=FT_LOOPS;
                        %                     for ni=1:numel(FT_LOOPS.label)
                        %                         nnz(FT_LOOPS.time{ni}>1 & FT_LOOPS.time{ni}<7 & ismember(FT_LOOPS.trial{ni},s1trials))
                        %                         nnz(FT_LOOPS.time{ni}>1 & FT_LOOPS.time{ni}<7 & ismember(FT_LOOPS.trial{ni},s2trials))
                        %                     end
                    end
                end
            end
        end
    end
end
if true
    save('loops_coding.mat','loops_sums');
end