% Major revision around 2022.12.20
% Major revision on 2023.07.16
% Fixed duplication

% chains_uf=wave.COM_chain(sel_meta);
% chains_uf_rev=wave.COM_chain(sel_meta,'reverse',true);
% blame=vcs.blame();
% save(fullfile('bzdata','chains_mix.mat'),'chains_uf','chains_uf_rev','blame')

function [out,chains]=COM_chain(su_meta,sel_meta,su_com_map,opt)
arguments
    su_meta
    sel_meta
    su_com_map
    opt.strict (1,1) logical = false %strict ccg criteria
    opt.reverse (1,1) logical = false
    opt.odor_only (1,1) logical = false
end
% global_init
load(fullfile("binary","sums_conn_10.mat"),'sums_conn_str');

% [sig,~]=bz.load_sig_sums_conn_file('pair',false);
% meta_str=ephys.util.load_meta('skip_stats',true);

chains=cell(0);
for fidx=1:numel(sums_conn_str)
    disp(fidx)
    % TODO update following code to reflect revised order
    ffpath=sums_conn_str(fidx).folder;
    dpath=regexp(ffpath,'(?<=SPKINFO[\\/]).*$','match','once');
    if isempty(dpath)
        dpath=ffpath;
    end
    sessid=ephys.path2sessid(dpath);
    %

    if opt.strict
        ccgqc=sums_conn_str(fidx).qc; %reference quality control parameter
        strict_sel=ccgqc(:,2)>=252 & ccgqc(:,4)>=2 & ccgqc(:,4)<=40 & ccgqc(:,5)>248;
        %1:Polarity 2:Time of peak 3:Noise peaks 4:FWHM 5:rising edge
        %6:falling edge
        oneccg=sums_conn_str(fidx).ccg_sc(strict_sel,:); %ccg
        onecon=sums_conn_str(fidx).sig_con(strict_sel,:); %jitter controlled significant functional coupling
        disp([fidx,nnz(strict_sel),numel(strict_sel)]);
    else % full input from English, Buzsaki code
        oneccg=sums_conn_str(fidx).ccg_sc;
        onecon=sums_conn_str(fidx).sig_con;
    end
    %TODO nonmem,incongruent possible if relax criteria?

    %     onecom=wave.get_pct_com_map(sel_meta,'onepath',sums_conn_str(fidx).folder,'curve',true); %per su center of mass, normalized FR
    %     skey=fieldnames(onecom);
    %     if isempty(skey), continue;end
    skey="s"+sessid;
    if ~isfield(su_com_map,skey)
        continue
    end

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
curr_sess=-1;
out.reg=cell(numel(out.sess),1);
out.cross_reg=false(numel(out.sess),1);
for ii=1:numel(out.sess)
    if out.sess(ii)~=curr_sess
        sess_cid=su_meta.allcid(su_meta.sess==out.sess(ii));
        sess_reg=su_meta.reg_tree(5,su_meta.sess==out.sess(ii));
    end
    [~,idces]=ismember(out.cids{ii},sess_cid);
    creg=sess_reg(idces);
    out.reg{ii}=creg;
    if numel(unique(creg))>1
        out.cross_reg(ii)=true;
    end
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

function chains=chain_one(onecon,typesel,curr_com_map,reverse,opt)
arguments
    onecon
    typesel
    curr_com_map
    reverse (1,1) logical = false
    opt.per_region (1,1) logical = true

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


