% Major revision around 2022.12.20
% shuf_chains=wave.COM_chain_shuf(wrs_mux_meta);

% blame=vcs.blame();
% save('chains_shuf.mat','shuf_chains','blame')

function shuf_out=COM_chain_shuf(sel_meta,shuf_idices,opt)
arguments
    sel_meta
    shuf_idices (1,:) double = 1:10
    opt.reverse (1,1) logical = false
end
% global_init
% load('sums_conn.mat','sums_conn_str');
% [sig,~]=bz.load_sig_sums_conn_file('pair',false);
% meta_str=ephys.util.load_meta('skip_stats',true);
% warning('partial iteration for illustration')
load('bz_ring_shufs.mat','shufs');
shuf_out=cell(max(shuf_idices),1);
for shufid=shuf_idices
    chains=cell(0);
    sig_str=shufs{shufid};
    all_sess=reshape(unique(sig_str.sess),1,[]);
    for sessid=all_sess
        sesssel=sig_str.sess==sessid;
        onecon=sig_str.suid(sesssel,:);

        %TODO nonmem,incongruent possible?
        onepath=ephys.sessid2path(sessid);
        onecom=wave.get_pct_com_map(sel_meta,'onepath',onepath,'curve',true); %per su center of mass, normalized FR
        skey=fieldnames(onecom);
        if isempty(skey), continue;end

        for wvtype=["s1d3","s1d6","s2d3","s2d6","olf_s1","olf_s2","dur_d3","dur_d6"]
            if ~isfield(onecom.(skey{1}),wvtype)
                continue
            end
            switch wvtype
                case "s1d3"
                    sub_chain=map2subchain(onecom,skey,'s1d3','olf_s1','com3',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,3,sub_chain);
                    sub_chain=map2subchain(onecom,skey,'s1d3','dur_d3','com3',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,3,sub_chain);

                case "s1d6"
                    sub_chain=map2subchain(onecom,skey,'s1d6','olf_s1','com6',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,6,sub_chain);
                    sub_chain=map2subchain(onecom,skey,'s1d6','dur_d6','com6',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,6,sub_chain);

                case "s2d3"
                    sub_chain=map2subchain(onecom,skey,'s2d3','olf_s2','com3',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,3,sub_chain);
                    sub_chain=map2subchain(onecom,skey,'s2d3','dur_d3','com3',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,3,sub_chain);

                case "s2d6"
                    sub_chain=map2subchain(onecom,skey,'s2d6','olf_s2','com6',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,6,sub_chain);
                    sub_chain=map2subchain(onecom,skey,'s2d6','dur_d6','com6',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,6,sub_chain);

                case "olf_s1"
                    sub_chain=map2subchain(onecom,skey,'olf_s1',[],'com3',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,3,sub_chain);
                    sub_chain=map2subchain(onecom,skey,'olf_s1',[],'com6',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,6,sub_chain);

                case "olf_s2"
                    sub_chain=map2subchain(onecom,skey,'olf_s2',[],'com3',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,3,sub_chain);
                    sub_chain=map2subchain(onecom,skey,'olf_s2',[],'com6',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,6,sub_chain);

                case "dur_d3"
                    sub_chain=map2subchain(onecom,skey,'dur_d3',[],'com3',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,3,sub_chain);

                case "dur_d6"
                    sub_chain=map2subchain(onecom,skey,'dur_d6',[],'com6',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,6,sub_chain);

                otherwise
                    keyboard()
            end
        end
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
    shuf_out{shufid}=out;
end
end

function chains=extend_chain(chains,sessid,wvtype,dur,sub_chain)
chains=[chains;...
    num2cell(repmat(sessid,numel(sub_chain),1)),...
    repmat({wvtype},numel(sub_chain),1),...
    num2cell(repmat(dur,numel(sub_chain),1)),...%duration
    sub_chain];
end

function sub_chain=map2subchain(onecom,skey,fn1,fn2,fn_sub,onecon,reverse)
if isempty(fn2) || ~isfield(onecom.(skey{1}),fn2)
    mapkeys=cell2mat(onecom.(skey{1}).(fn1).(fn_sub).keys());
    typesel=all(ismember(int32(onecon),mapkeys),2);
    curr_com_map=onecom.(skey{1}).(fn1).(fn_sub);
else
    mapkeys=[cell2mat(onecom.(skey{1}).(fn1).(fn_sub).keys()),...
        cell2mat(onecom.(skey{1}).(fn2).(fn_sub).keys())];
    typesel=all(ismember(int32(onecon),mapkeys),2);
    mapv=[cell2mat(onecom.(skey{1}).(fn1).(fn_sub).values()),...
        cell2mat(onecom.(skey{1}).(fn2).(fn_sub).values())];
    curr_com_map=containers.Map(num2cell(mapkeys),num2cell(mapv));
end
sub_chain=chain_one(onecon,typesel,curr_com_map,reverse);
end

function chains=chain_one(onecon,typesel,curr_com_map,reverse)
arguments
    onecon
    typesel
    curr_com_map
    reverse % (1,1) logical = false
end
chains=cell(0);
if nnz(typesel)<3,return;end
typesigcon=onecon(typesel,:);
con_com_prepost=cell2mat(curr_com_map.values(num2cell(int32(typesigcon))));
if reverse
    dirsel=con_com_prepost(:,2)<con_com_prepost(:,1);
else
    dirsel=con_com_prepost(:,2)>con_com_prepost(:,1); % Assuming 250ms bin
end
dirsigcon=typesigcon(dirsel,:);
upre=unique(dirsigcon(:,1)).';

for i=upre
    onechain=cell(0);
    cpre=i;
    while true % first pass-through without unfolding all chains
        newpair=dirsigcon(ismember(dirsigcon(:,1),cpre),:);
        if isempty(newpair)
            if numel(onechain)>1
                chains=[chains;{onechain}];
            end
            break
        else
            onechain{end+1}={...
                newpair,...
                cell2mat(curr_com_map.values(num2cell(newpair)))};
            cpre=newpair(:,2);
        end
    end
end
end





