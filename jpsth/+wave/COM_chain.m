% Major revision around 2022.12.20
% chains_uf=wave.COM_chain(sel_meta);
% chains_uf_rev=wave.COM_chain(sel_meta,'reverse',true);
% blame=vcs.blame();
% save(fullfile('bzdata','chains_mix.mat'),'chains_uf','chains_uf_rev','blame')

function [out,chains]=COM_chain(sel_meta,su_com_map,opt)
arguments
    sel_meta
    su_com_map
    opt.strict (1,1) logical = false %strict ccg criteria
    opt.reverse (1,1) logical = false
    opt.odor_only (1,1) logical = false
end
% global_init
load('sums_conn_10.mat','sums_conn_str');
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

    if opt.odor_only
        if isfield(su_com_map.(skey),"olf_s1")
            sub_chain=map2subchain(su_com_map,skey,'olf_s1',[],'com',onecon,opt.reverse);
            chains=extend_chain(chains,sessid,"olf_s1",NaN,sub_chain);
        end
        if isfield(su_com_map.(skey),"olf_s2")
            sub_chain=map2subchain(su_com_map,skey,'olf_s2',[],'com',onecon,opt.reverse);
            chains=extend_chain(chains,sessid,"olf_s2",NaN,sub_chain);
        end
    else
        for wvtype=["s1d3","s1d6","s2d3","s2d6","olf_s1","olf_s2","dur_d3","dur_d6"]
            if ~isfield(su_com_map.(skey{1}),wvtype)
                continue
            end
            switch wvtype
                case "s1d3"
                    sub_chain=map2subchain(su_com_map,skey,'s1d3','olf_s1','com3',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,3,sub_chain);
                    sub_chain=map2subchain(su_com_map,skey,'s1d3','dur_d3','com3',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,3,sub_chain);

                case "s1d6"
                    sub_chain=map2subchain(su_com_map,skey,'s1d6','olf_s1','com6',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,6,sub_chain);
                    sub_chain=map2subchain(su_com_map,skey,'s1d6','dur_d6','com6',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,6,sub_chain);

                case "s2d3"
                    sub_chain=map2subchain(su_com_map,skey,'s2d3','olf_s2','com3',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,3,sub_chain);
                    sub_chain=map2subchain(su_com_map,skey,'s2d3','dur_d3','com3',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,3,sub_chain);

                case "s2d6"
                    sub_chain=map2subchain(su_com_map,skey,'s2d6','olf_s2','com6',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,6,sub_chain);
                    sub_chain=map2subchain(su_com_map,skey,'s2d6','dur_d6','com6',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,6,sub_chain);

                case "olf_s1"
                    sub_chain=map2subchain(su_com_map,skey,'olf_s1',[],'com3',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,3,sub_chain);
                    sub_chain=map2subchain(su_com_map,skey,'olf_s1',[],'com6',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,6,sub_chain);

                case "olf_s2"
                    sub_chain=map2subchain(su_com_map,skey,'olf_s2',[],'com3',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,3,sub_chain);
                    sub_chain=map2subchain(su_com_map,skey,'olf_s2',[],'com6',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,6,sub_chain);

                case "dur_d3"
                    sub_chain=map2subchain(su_com_map,skey,'dur_d3',[],'com3',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,3,sub_chain);

                case "dur_d6"
                    sub_chain=map2subchain(su_com_map,skey,'dur_d6',[],'com6',onecon,opt.reverse);
                    chains=extend_chain(chains,sessid,wvtype,6,sub_chain);
            end
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


