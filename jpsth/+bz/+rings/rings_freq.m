%% enumerate single spike and burst ring activities
%% main function works, aux functions WIP

% TODO: spike screen by wave / trial

function rings_freq(opt)
arguments
    opt.burst (1,1) logical = true
    opt.burstInterval (1,1) double = 600
    opt.ccg  (1,1) logical = false
end
% bz.rings.ring_list_bz
load(fullfile('bzdata','rings_bz.mat'),'rings');

blame=vcs.blame();

%loop entry
out=cell(0);
for sessidx=1:size(rings,1)
    if all(cellfun(@(x) isempty(x), rings(sessidx,:)),'all')
        continue
    end
    [spkID,spkTS,trials,suids,folder]=ephys.getSPKID_TS(sessidx);
    for rsize=3:5
        sess_rings=rings{sessidx,rsize-2};
        if isempty(sess_rings)
            continue
        end
        rids=1:size(sess_rings,1);
        for ring_id=rids
%             disp(ring_id);
            cids=sess_rings(ring_id,:);
            per_cid_spk_cnt=cids;
            ts_id=[];
            for in_ring_pos=1:rsize
                one_ring_sel=spkID==cids(in_ring_pos);
                per_cid_spk_cnt(in_ring_pos)=nnz(one_ring_sel);
                ts_id=cat(1,ts_id,[spkTS(one_ring_sel),ones(per_cid_spk_cnt(in_ring_pos),1)*in_ring_pos]);
            end
            % TODO: match trial, screen spikes
            if opt.burst
                tsidin=cell(1,numel(cids));
                for kk=1:numel(cids)
                    tsidin{kk}=ts_id(ts_id(:,2)==kk,1);
                end
                % out=relax_tag_long(in,loopIdx,recDepth,loopCnt,perSU,opt)
                ts=bz.rings.relax_tag_long(tsidin,[],[],[],[],"burstInterval",opt.burstInterval);

                if ~isempty(ts)
                    % TODO: optional remove shorter chains for each
                    % onset-spike

                    outkey="s"+sessidx+"r"+rsize+"n"+ring_id;
%                     out.("d"+duration).(wid).(outkey).ts=ts;
%                     out.("d"+duration).(wid).(outkey).meta={cids}; % Skipped TCOM for now
%                     out.("d"+duration).(wid).(outkey).ts_id=ts_id;

                    out.(outkey).ts=ts;
                    out.(outkey).meta={cids}; % Skipped TCOM for now
                    out.(outkey).ts_id=ts_id;
                    % TODO: ccg
%                     if opt.ccg
%                         cursess=chains.sess(ring_id);
%                         sesspath=ephys.sessid2path(cursess);
%                         strippath=regexp(sesspath,'(?<=\\).*','match');
%                         sesssel=find(contains({sums_conn_str.folder},strippath));
%                         ccg=sums_conn_str(sesssel).ccg_sc;
%                         ccgid=sums_conn_str(sesssel).sig_con;
%                         chainccg=[];
%                         for ii=1:numel(cids)-1
%                             sigsel=ccgid(:,1)==cids(ii) & ccgid(:,2)==cids(ii+1);
%                             if nnz(sigsel)~=1
%                                 keyboard()
%                             end
%                             chainccg=[chainccg;ccg(sigsel,:)];
%                         end
%                         out.("d"+duration).(wid).(outkey).ccgs=chainccg;
%                     end
                end
            else
                error("Incomplete section")
%                 keyboard()
                ts_id=sortrows(ts_id,1);
                ring_stats=bz.rings.relax_tag(ts_id,rsize);
                coact_count=sum(ring_stats.spk_cnt);
                if coact_count>(rsize+1)*10
                    out(end+1,:)={sessidx,ring_id,cids,per_cid_spk_cnt,ring_stats,coact_count./(ts_id(end,1)./30000)};
                end
            end
        end
    end
    save(sprintf('rings_burst_%d.mat',opt.burstInterval),'out','blame');
end
end
