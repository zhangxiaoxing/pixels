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

for ring_id=1:size(sess_rings,1)
    disp(ring_id);
    disp('RRRRR');
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
    if coact_count>ts_id(end,1)*0.1/30000
        %sess, ring_id, cids, cid spk count, cid asso spk,
        %ring_spk_count_dist <-, ring_duration_dist
%         [GC,GR]=groupcounts(ts_id(ring_stats(:,3)~=0,2));
%         per_su_asso=arrayfun(@(x) GC(GR==x),sort(GR));
%         [spk_cnt_dist,dur_dist]=get_dist(ring_stats);
        sums(end+1,:)={sessidx,ring_id,cids,per_cid_spk_cnt,ring_stats};
        save(sprintf('ring_stats_%d_%d.mat',rsize,sessidx),'sums');
    end
end

% 
% function [spk_cnt_dist,dur_dist]=get_dist(tagged)
% % tagged=[tagged,(1:size(tagged,1))'];
% spk_cnt_dist=[];
% dur_dist=[];
% rpt_idx=1;
% rhead=1;
% rtail=rhead;
% % fprintf('000000');
% while true
%     %     fprintf('\b\b\b\b\b\b%06d',curr_ptr);
%     rhead=find(tagged(:,3)==rpt_idx,1);
%     if ~isempty(rhead)
%         rtail=find(tagged(:,3)==rpt_idx,1,'last');
%         spk_cnt_dist=horzcat(spk_cnt_dist,nnz(tagged(rhead:rtail,3)));
%         if nnz(tagged(rhead:rtail,3))<4
%             keyboard
%         end
%         dur_dist=horzcat(dur_dist,tagged(rtail,1)-tagged(rhead,1));
%         rpt_idx=rpt_idx+1;
%     else
%         break
%     end
% end
% end



