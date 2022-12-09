% TODO load sums_all


function pstats=load_data(rsize,sums_all)

% TODO: also supply per_trial data

rstats=sums_all{rsize-2};
rstats=rstats(cell2mat(rstats(:,6))>0.1 & cellfun(@(x) numel(unique(x)),rstats(:,3))==rsize,:);
usess=unique(cell2mat(rstats(:,1)));
pstats=struct();
for sessid=usess.'
    pstats.(sprintf('s%d',sessid))=struct();
    [~,~,trials,suids,folder,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true);
    rid=find(cell2mat(rstats(:,1))==sessid);
    for ri=reshape(rid,1,[])
        disp([sessid,ri]);
        ts_id=[];
        cids=rstats{ri,3};
        per_cid_spk_cnt=cids;
        for in_ring_pos=1:numel(cids) % TODO, 1:rsize
            one_ring_sel=strcmp(FT_SPIKE.label,num2str(cids(in_ring_pos)));
            per_cid_spk_cnt(in_ring_pos)=numel(FT_SPIKE.timestamp{one_ring_sel});
            ts_id=cat(1,ts_id,[FT_SPIKE.timestamp{one_ring_sel};...
                repmat(in_ring_pos,1,per_cid_spk_cnt(in_ring_pos));...
                FT_SPIKE.trial{one_ring_sel};...
                FT_SPIKE.time{one_ring_sel}].'); % TS, ring_id
        end
        ts_id=sortrows(ts_id,1);
        ts_id=[ts_id,rstats{ri,5}.tags]; % join TS, ring tag

        for cidx=1:numel(rstats{ri,3})
            cid=rstats{ri,3}(cidx);
            if ~isfield(pstats.(sprintf('s%d',sessid)),sprintf('c%d',cid))
                pstats.(sprintf('s%d',sessid)).(sprintf('c%d',cid))=struct();
            end
            pstats.(sprintf('s%d',sessid)).(sprintf('c%d',cid)).(sprintf('r%d',ri))=ts_id(ts_id(:,2)==cidx & ts_id(:,3)~=0,1);
        end
    end
end
end
