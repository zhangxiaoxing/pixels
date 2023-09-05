% Pending removal as of 2023.09.05
% Mostly duplicated and fall behind updates

% Major revision around 2022.12.20
% Major update around 2023.07.16, removed duplicates

% shuf_chains=wave.COM_chain_shuf(wrs_mux_meta,1:100,'odor_only',true)
% blame=vcs.blame();save('chains_shuf.mat','shuf_chains','blame')

function shuf_chains=COM_chain_shuf(sel_meta,shuf_idices,opt)
arguments
    sel_meta
    shuf_idices (1,:) double = 1:100
    opt.reverse (1,1) logical = false
    opt.odor_only (1,1) logical = false
    opt.skip_save (1,1) logical = true
end
% global_init
% load('sums_conn.mat','sums_conn_str');
% [sig,~]=bz.load_sig_sums_conn_file('pair',false);
% meta_str=ephys.util.load_meta('skip_stats',true);
% warning('partial iteration for illustration')
load(fullfile('binary','bz_ring_shufs.mat'),'shufs');
shuf_chains=cell(max(shuf_idices),1);
for shufid=shuf_idices
    disp("SHUF"+shufid);
    chains=cell(0);
    sig_str=shufs{shufid};
    all_sess=reshape(unique(sig_str.sess),1,[]);
    for sessid=all_sess
        sesssel=sig_str.sess==sessid;
        onecon=sig_str.suid(sesssel,:);

        %TODO nonmem,incongruent possible?
        su_com_map=wave.get_pct_com_map(sel_meta,'one_sess',sessid,'curve',true); %per su center of mass, normalized FR
        skey=fieldnames(su_com_map);
        if isempty(skey), continue;end

        for delay=[3 6]
            if isfield(su_com_map.(skey{1}),"s1d"+delay) && isfield(su_com_map.(skey{1}),"olf_s1")
                sub_chain=map2subchain(su_com_map,skey,"s1d"+delay,"olf_s1","com"+delay,onecon,opt.reverse);
                chains=extend_chain(chains,sessid,"s1d"+delay,delay,sub_chain);
            elseif  isfield(su_com_map.(skey{1}),"s1d"+delay)
                sub_chain=map2subchain(su_com_map,skey,"s1d"+delay,[],"com"+delay,onecon,opt.reverse);
                chains=extend_chain(chains,sessid,"s1d"+delay,delay,sub_chain);
            elseif  isfield(su_com_map.(skey{1}),"olf_s1")
                sub_chain=map2subchain(su_com_map,skey,'olf_s1',[],"com"+delay,onecon,opt.reverse);
                chains=extend_chain(chains,sessid,"s1d"+delay,delay,sub_chain);
            end

            if isfield(su_com_map.(skey{1}),"s2d"+delay) && isfield(su_com_map.(skey{1}),"olf_s2")
                sub_chain=map2subchain(su_com_map,skey,"s2d"+delay,"olf_s2","com"+delay,onecon,opt.reverse);
                chains=extend_chain(chains,sessid,"s2d"+delay,delay,sub_chain);
            elseif  isfield(su_com_map.(skey{1}),"s2d"+delay)
                sub_chain=map2subchain(su_com_map,skey,"s2d"+delay,[],"com"+delay,onecon,opt.reverse);
                chains=extend_chain(chains,sessid,"s2d"+delay,delay,sub_chain);
            elseif  isfield(su_com_map.(skey{1}),"olf_s2")
                sub_chain=map2subchain(su_com_map,skey,'olf_s2',[],"com"+delay,onecon,opt.reverse);
                chains=extend_chain(chains,sessid,"s2d"+delay,delay,sub_chain);
            end

            if opt.odor_only
                continue
            end
            if isfield(su_com_map.(skey{1}),"s1d"+delay) && isfield(su_com_map.(skey{1}),"dur_d"+delay)
                sub_chain=map2subchain(su_com_map,skey,"s1d"+delay,"dur_d"+delay,"com"+delay,onecon,opt.reverse);
                chains=extend_chain(chains,sessid,"s1d"+delay,delay,sub_chain);
            end
            if isfield(su_com_map.(skey{1}),"s2d"+delay) && isfield(su_com_map.(skey{1}),"dur_d"+delay)
                sub_chain=map2subchain(su_com_map,skey,"s2d"+delay,"dur_d"+delay,"com"+delay,onecon,opt.reverse);
                chains=extend_chain(chains,sessid,"s2d"+delay,delay,sub_chain);
            end
            if (~isfield(su_com_map.(skey{1}),"s1d"+delay))&& (~isfield(su_com_map.(skey{1}),"s2d"+delay)) && isfield(su_com_map.(skey{1}),"dur_d"+delay)
                sub_chain=map2subchain(su_com_map,skey,"dur_d",[],"com"+delay,onecon,opt.reverse);
                chains=extend_chain(chains,sessid,"dur_"+delay,delay,sub_chain);
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
    shuf_chains{shufid}=out;
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
