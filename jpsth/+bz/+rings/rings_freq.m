function rings_freq(sessidx, rsize, opt)
arguments
    sessidx (1,1) double {mustBePositive,mustBeInteger} = 1
    rsize (1,1) double {mustBeMember(rsize,3:5)} = 3
    opt.ridx (1,:) double = []
    %%@deprecated
    %     opt.partialfile (1,:) char = []
end
load(fullfile('bzdata','rings_bz.mat'),'rings');
%loop entry
if sessidx>size(rings,1) || isempty(rings{sessidx,rsize-2})
    if isunix, quit(0);else, return;end
end
sess_rings=rings{sessidx,rsize-2};
[spkID,spkTS,trials,suids,folder]=ephys.getSPKID_TS(sessidx);
if isempty(spkID), quit(0); end

if isempty(opt.ridx)
    rids=1:size(sess_rings,1);
else
    rids=opt.ridx;
end
%% @deprecated
% if ~isempty(opt.partialfile)
%     load(opt.partialfile,'out');
%     part_sess=cell2mat(out(:,1));
%     if ~any(part_sess(:,1)==sessidx & part_sess(:,2)==rsize)
%         return
%     else
%         pfidx=find(part_sess(:,1)==sessidx & part_sess(:,2)==rsize);
%     end
% end
%%
sums=cell(0);
for ring_id=rids
    disp(ring_id);
    %% @deprecated
    %         if ~isempty(opt.partialfile) && ismember(ring_id,out{pfidx,2}(1:end-1))
    %             break;
    %         end
    %%
    cids=sess_rings(ring_id,:);
    per_cid_spk_cnt=cids;
    ts_id=[];
    for in_ring_pos=1:rsize
        one_ring_sel=spkID==cids(in_ring_pos);
        per_cid_spk_cnt(in_ring_pos)=nnz(one_ring_sel);
        ts_id=cat(1,ts_id,[spkTS(one_ring_sel),ones(per_cid_spk_cnt(in_ring_pos),1)*in_ring_pos]);
    end
    ts_id=sortrows(ts_id,1);
    ring_stats=bz.rings.relax_tag(ts_id,rsize);
    coact_count=sum(ring_stats.spk_cnt);
    
    if coact_count>(rsize+1)*10
        sums(end+1,:)={sessidx,ring_id,cids,per_cid_spk_cnt,ring_stats,coact_count./(ts_id(end,1)./30000)};
    end
end
if size(sums,1)>0
    save(sprintf('ring_stats_%d_%d_%d_%d.mat',rsize,sessidx,min(rids),max(rids)),'sums');
end
end
