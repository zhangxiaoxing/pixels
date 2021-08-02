function rings_freq(sessidx, rsize)
arguments
    sessidx (1,1) double {mustBePositive,mustBeInteger}
    rsize (1,1) double {mustBeMember(rsize,3:5)}
end
load(fullfile('bzdata','rings_bz.mat'),'rings');
sums=cell(0);
%loop entry
if sessidx>size(rings,1) || isempty(rings{sessidx,rsize-2})
    if isunix, quit(0);else, return;end
end
sess_rings=rings{sessidx,rsize-2};
[spkID,spkTS,trials,suids,folder]=ephys.getSPKID_TS(sessidx);
if isempty(spkID), quit(0); end
ts_id=[];
for ring_id=1:size(sess_rings,1)
    disp(ring_id);
%     disp('RRRRR');
    cids=sess_rings(ring_id,:);
    per_cid_spk_cnt=cids;
    for in_ring_pos=1:rsize
        one_ring_sel=spkID==cids(in_ring_pos);
        per_cid_spk_cnt(in_ring_pos)=nnz(one_ring_sel);
        ts_id=cat(1,ts_id,[spkTS(one_ring_sel),ones(per_cid_spk_cnt(in_ring_pos),1)*in_ring_pos]);
    end
    ts_id=sortrows(ts_id,1);
    ring_stats=bz.rings.relax_tag(ts_id,rsize);
    coact_count=sum(ring_stats.spk_cnt);
%%   old criteria of 0.1Hz
    if coact_count>ts_id(end,1)*0.1/30000
        sums(end+1,:)={sessidx,ring_id,cids,per_cid_spk_cnt,ring_stats};
        save(sprintf('ring_stats_%d_%d.mat',rsize,sessidx),'sums');
    end
end
end
