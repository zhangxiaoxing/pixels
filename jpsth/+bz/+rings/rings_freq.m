load(fullfile('bzdata','rings_bz.mat'));
% homedir=ephys.util.getHomedir();
% ephys.util.dependency %data path and lib path dependency
% sess=dir(fullfile(homedir,'**','spike_info.hdf5'));
% [~,idces]=sort({sess.folder});sess=sess(idces);
sums=cell(0);
%loop entry
if exist('sessidx','var') && exist('rsize','var'), disp([sessidx,rsize]); else, quit(0); end
if sessidx>size(rings,1) || isempty(rings{sessidx,rsize-2}), quit(0); end
sess_rings=rings{sessidx,rsize-2};
[spkID,spkTS,trials,suids,folder]=ephys.getSPKID_TS(sessidx);
if isempty(spkID), quit(0); end

ts_id=[];
fprintf('0000');
for ring_id=1:size(sess_rings,1)
    fprintf('\b\b\b\b%04d',ring_id);
    cids=sess_rings(ring_id,:);
    per_cid_spk_cnt=cids;
    for in_ring_pos=1:rsize
        one_ring_sel=spkID==cids(in_ring_pos);
        per_cid_spk_cnt(in_ring_pos)=nnz(one_ring_sel);
        ts_id=cat(1,ts_id,[spkTS(one_ring_sel),ones(per_cid_spk_cnt(in_ring_pos),1)*in_ring_pos]);
    end
    ts_id=sortrows(ts_id,1);
    ts_id_tagged=bz.rings.relax_tag(ts_id,rsize);
    coact_count=nnz(ts_id_tagged(:,3));
    if coact_count>ts_id(end,1)*0.1/30000
        %sess, ring_id, cids, cid spk count, cid asso spk,
        %ring_spk_count_dist <-, ring_duration_dist
        [GC,GR]=groupcounts(ts_id(ts_id_tagged(:,3)~=0,2));
        per_su_asso=arrayfun(@(x) GC(GR==x),sort(GR));
        [spk_cnt_dist,dur_dist]=get_dist(tagged);
        sums(end+1,:)={sessidx,ring_id,cids,per_cid_spk_cnt,per_su_asso,...
            spk_cnt_dist,dur_dist};
    end
end
save(sprintf('ring_stats_%d_%d.mat',rsize,sessidx),'sums');

function [spk_cnt_dist,dur_dist]=get_dist(tagged)
spk_cnt_dist=[];
dur_dist=[];
curr_ptr=1;
% fprintf('000000');
while true
    %     fprintf('\b\b\b\b\b\b%06d',curr_ptr);
    curr_entry=find(tagged(curr_ptr:end,3)~=0,1)+curr_ptr-1;
    if ~isempty(curr_entry)
        curr_sep=find(tagged(curr_entry:end,3)==0,1);
        if ~isempty(curr_sep)
            curr_ptr=curr_sep+curr_entry-1;
            spk_cnt_dist=horzcat(spk_cnt_dist,curr_sep-1);
            dur_dist=horzcat(dur_dist,tagged(curr_ptr-1,1)-tagged(curr_entry,1));
        else
            %             spk_cnt_dist=horzcat(spk_cnt_dist,size(tagged,1)-curr_entry+1);
            break
        end
    else
        break
    end
end
end



