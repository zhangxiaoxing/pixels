
% nonmem_chain=wave.COM_chain_nonmem(sel_meta);
% blame=vcs.blame();
% save('chains_nonmem.mat','nonmem_chains','blame')

function out=COM_chain_nonmem(sel_meta)
arguments
    sel_meta
end
% global_init
% load('sums_conn.mat','sums_conn_str');
[sig,~]=bz.load_sig_sums_conn_file('pair',false);
% meta_str=ephys.util.load_meta('skip_stats',true);
% warning('partial iteration for illustration')
% shuf_out=cell(max(shuf_idices),1);
% TODO: nonmem out
greys=ephys.getGreyRegs('range','grey');
idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));
greysId=int32(cell2mat(idmap.reg2ccfid.values(greys)));
chains=cell(0);
all_sess=reshape(unique(sig.sess),1,[]);
for sessid=all_sess
    sesssel=sig.sess==sessid;
    onecon=sig.suid(sesssel,:);
    sig=bz.join_fc_waveid(sig,sel_meta.wave_id);
    typesel=all(sig.waveid(sesssel,:)==0,2) ...
        & all(ismember(sig.reg(sesssel,5,:),greysId),3);
    subchain=cell(0);
    if nnz(typesel)<2,continue;end
    typesigcon=onecon(typesel,:);
    upre=unique(typesigcon(:,1)).';

    for i=upre
        onechain=cell(0);
        cpre=i;
        acyclic=i;
        while true % first pass-through without unfolding all chains
            newpair=typesigcon(ismember(typesigcon(:,1),cpre)...
                & ~ismember(typesigcon(:,2),acyclic),:);
            if isempty(newpair)
                if numel(onechain)>1
                    subchain=[subchain;{onechain}];
                end
                break
            else
                onechain{end+1}={newpair,zeros(size(newpair))};
                cpre=newpair(:,2);
                acyclic=[acyclic;cpre];
            end
        end
    end
    chains=extend_chain(chains,sessid,'nonmem',0,subchain);
end

% unfold chain-trees
out=struct();
[out.sess,out.wave,out.dur,out.cids,out.tcoms]=deal([]);
for ii=1:size(chains,1)
    split_chains=wave.recursive_chain(chains{ii,4},[]); %one su per order
    out.sess=[out.sess;repmat(chains{ii,1},size(split_chains.cids,1),1)];
    out.wave=[out.wave;repmat(chains{ii,2},size(split_chains.cids,1),1)];
    out.dur=[out.dur;repmat(chains{ii,3},size(split_chains.cids,1),1)];
    out.cids=[out.cids;split_chains.cids];
    out.tcoms=[out.tcoms;split_chains.tcoms];
end
end

function chains=extend_chain(chains,sessid,wvtype,dur,sub_chain)
chains=[chains;...
    num2cell(repmat(sessid,numel(sub_chain),1)),...
    repmat({wvtype},numel(sub_chain),1),...
    num2cell(repmat(dur,numel(sub_chain),1)),...%duration
    sub_chain];
end







