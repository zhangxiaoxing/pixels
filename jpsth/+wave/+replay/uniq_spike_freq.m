function motif_seq_activation(chain_replay,ring_replay,trials_dict)
arguments
    chain_replay = []
    ring_replay = []
    trials_dict = []
end
if isempty(chain_replay) || isempty (ring_replay)
    load(fullfile('binary','motif_replay.mat'),'ring_replay','chain_replay');
end
if isempty(trials_dict)
    load(fullfile('binary','trials_dict.mat'),'trials_dict');
end

usess=union(chain_replay.session,ring_replay.session);
for sessid=reshape(usess,1,[])
    sesscid=su_meta.allcid(su_meta.sess==sessid);
    sessreg=su_meta.reg_tree(5,su_meta.sess==sessid).';
    sess_reg_map=dictionary(sesscid(~ismissing(sessreg)),sessreg(~ismissing(sessreg)));
    trials=cell2mat(trials_dict(sessid));
    for wid="s1d6"%["s1d3","s2d3","s1d6","s2d6"]
        disp(wid+"of sess"+sessid)
        samp=str2double(regexp(wid,"(?<=s)\d(?=d)",'match','once')).*4;
        delay=str2double(regexp(wid,"(?<=d)\d$",'match','once'));
        trial_sel=find(trials(:,5)==samp & trials(:,8)==delay & all(trials(:,9:10)>0,2));
        ring_reg_sel=cellfun(@(x) all(sess_reg_map.isKey(x)),ring_replay.meta(:,2));
        chain_reg_sel=cellfun(@(x) all(sess_reg_map.isKey(x)),chain_replay.meta(:,2));
        ringsel=find(ring_replay.session==sessid & ring_replay.wave==wid & ring_reg_sel);
        chainsel=find(chain_replay.session==sessid & chain_replay.wave==wid & chain_reg_sel);

        motif_spk=cell(0,1);
        motif_id=cell(0,2);
        motif_seq=cell(0,1);
        edges=[];
        for ridx=reshape(ringsel,1,[])
            trl_align=ring_replay.trl_align{ridx};
            pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
            currts=ring_replay.ts{ridx}(pref_delay);
            if ~isempty(currts)
                motif_spk=[motif_spk;{currts}];
                motif_id=[motif_id;{"r"+ridx,ring_replay.meta{ridx,2}}];
                motif_seq=[motif_seq;{ring_replay.ts_seq{ridx}(pref_delay)}];
                edges=[edges;ring_replay.meta{ridx,2}([1:end;2:end,1]).'];
            end
        end
        for cidx=reshape(chainsel,1,[])
            trl_align=chain_replay.trl_align{cidx};
            pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
            currts=chain_replay.ts{cidx}(pref_delay,:);
            if ~isempty(currts)
                motif_spk=[motif_spk;{mat2cell(currts.',size(currts,2),ones(size(currts,1),1)).'}];
                motif_id=[motif_id;{"c"+cidx,chain_replay.meta{cidx,2}}];
                motif_seq=[motif_seq;{repmat({chain_replay.meta{cidx,2}.'},nnz(pref_delay),1)}];
                edges=[edges;chain_replay.meta{cidx,2}([1:end-1;2:end]).'];
            end
        end

        %% plot
        allspk=cell2mat(cellfun(@(x) cell2mat(x),motif_spk,'UniformOutput',false));
        allid=cell2mat(cellfun(@(x) cell2mat(x),motif_seq,'UniformOutput',false));
        uspkts=unique([allspk,allid],'rows');
    end
end
end





