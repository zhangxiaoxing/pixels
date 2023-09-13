% reg_com_maps=cell2struct({tcom3_maps;tcom6_maps},{'tcom3_maps','tcom6_maps'});

function [out,chains]=COM_chain_reg(su_meta,sel_meta,reg_com_maps,opt)
arguments
    su_meta
    sel_meta
    reg_com_maps
    opt.strict (1,1) logical = false %strict ccg criteria
    opt.reverse (1,1) logical = false
    opt.odor_only (1,1) logical = true
    opt.non_mem (1,1) logical = false
    opt.cross_only (1,1) logical = false
    opt.shuf (1,1) logical = false
    opt.shuf_data = []
    
end
assert(opt.odor_only,"Unfinished")
if ~opt.shuf
    load(fullfile("binary","sums_conn_10.mat"),'sums_conn_str');
else
    sums_conn_str=zeros(max(opt.shuf_data.sess),1);
end

chains=cell(0);
for fidx=1:numel(sums_conn_str)
    disp(fidx)
    % TODO update following code to reflect revised order
    if opt.shuf
        sessid=fidx;
        if ~any(opt.shuf_data.sess==sessid)
            continue
        end
    else
        ffpath=sums_conn_str(fidx).folder;
        dpath=regexp(ffpath,'(?<=SPKINFO[\\/]).*$','match','once');
        if isempty(dpath)
            dpath=ffpath;
        end
        sessid=ephys.path2sessid(dpath);
    end
    sess_sel=su_meta.sess==sessid;
    %separate 3s & 6s
    for delay=[3 6]
        if delay==3
            s1_sel=ismember(sel_meta.wave_id,[1 5]);
            s2_sel=ismember(sel_meta.wave_id,[3 6]);
        else
            s1_sel=ismember(sel_meta.wave_id,[2 5]);
            s2_sel=ismember(sel_meta.wave_id,[4 6]);
        end

        grey_regs=reg_com_maps.("tcom"+delay+"_maps").odor_only.keys;
        grey_sel=ismember(su_meta.reg_tree(5,:),grey_regs).';

        if opt.non_mem
            nmm_sel=sel_meta.wave_id==0;
            nm_cids=su_meta.allcid(nmm_sel & sess_sel & grey_sel);
            nm_regs=su_meta.reg_tree(5,nmm_sel & sess_sel & grey_sel);
            reg_tcom.nm.("d"+delay)=reg_com_maps.("tcom"+delay+"_maps").odor_only.values(nm_regs);
        else
            s1_cids=su_meta.allcid(s1_sel & sess_sel & grey_sel);
            s2_cids=su_meta.allcid(s2_sel & sess_sel & grey_sel);
            s1_regs=su_meta.reg_tree(5,s1_sel & sess_sel & grey_sel);
            s2_regs=su_meta.reg_tree(5,s2_sel & sess_sel & grey_sel);

            reg_tcom.s1.("d"+delay)=reg_com_maps.("tcom"+delay+"_maps").odor_only.values(s1_regs);
            reg_tcom.s2.("d"+delay)=reg_com_maps.("tcom"+delay+"_maps").odor_only.values(s2_regs);
        end
        if false % opt.strict % TODO: deal with details later
            ccgqc=sums_conn_str(fidx).qc; %reference quality control parameter
            strict_sel=ccgqc(:,2)>=252 & ccgqc(:,4)>=2 & ccgqc(:,4)<=40 & ccgqc(:,5)>248;
            %1:Polarity 2:Time of peak 3:Noise peaks 4:FWHM 5:rising edge
            %6:falling edge
            oneccg=sums_conn_str(fidx).ccg_sc(strict_sel,:); %ccg
            sess_con=sums_conn_str(fidx).sig_con(strict_sel,:); %jitter controlled significant functional coupling
            disp([fidx,nnz(strict_sel),numel(strict_sel)]);
        else % full input from English, Buzsaki code
            if opt.shuf
                oneccg=[];
                sess_con=double(opt.shuf_data.suid(opt.shuf_data.sess==sessid,:));
            else
                oneccg=sums_conn_str(fidx).ccg_sc;
                sess_con=sums_conn_str(fidx).sig_con;
            end
        end
        % TODO: nonmem,incongruent possible if relax criteria?

        if opt.non_mem
            sub_chain=chain_one(sess_con,nm_cids,reg_tcom.nm.("d"+delay),opt.reverse);
            if ~isempty(sub_chain)
                chains=extend_chain(chains,sessid,"nmd",delay,sub_chain);
            end
        else
            % s1
            sub_chain=chain_one(sess_con,s1_cids,reg_tcom.s1.("d"+delay),opt.reverse);
            if ~isempty(sub_chain)
                chains=extend_chain(chains,sessid,"s1d",delay,sub_chain);
            end
            % s2
            sub_chain=chain_one(sess_con,s2_cids,reg_tcom.s2.("d"+delay),opt.reverse);
            if ~isempty(sub_chain)
                chains=extend_chain(chains,sessid,"s2d",delay,sub_chain);
            end
        end

        % if opt.odor_only  % TODO: .dur
        %     continue
        % end
        % if isfield(su_com_map.(skey{1}),"s1d"+delay) && isfield(su_com_map.(skey{1}),"dur_d"+delay)
        %     sub_chain=map2subchain(su_com_map,skey,"s1d"+delay,"dur_d"+delay,"com"+delay,onecon,opt.reverse);
        %     chains=extend_chain(chains,sessid,"s1d"+delay,delay,sub_chain);
        % end
        % if isfield(su_com_map.(skey{1}),"s2d"+delay) && isfield(su_com_map.(skey{1}),"dur_d"+delay)
        %     sub_chain=map2subchain(su_com_map,skey,"s2d"+delay,"dur_d"+delay,"com"+delay,onecon,opt.reverse);
        %     chains=extend_chain(chains,sessid,"s2d"+delay,delay,sub_chain);
        % end
        % if (~isfield(su_com_map.(skey{1}),"s1d"+delay))&& (~isfield(su_com_map.(skey{1}),"s2d"+delay)) && isfield(su_com_map.(skey{1}),"dur_d"+delay)
        %     sub_chain=map2subchain(su_com_map,skey,"dur_d",[],"com"+delay,onecon,opt.reverse);
        %     chains=extend_chain(chains,sessid,"dur_"+delay,delay,sub_chain);
        % end
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

    % ensure acyclic
    if numel(unique(out.cids{ii}))<numel(out.cids{ii})
        keyboard()
    end
    % cross-region check
    creg=sess_reg(idces);
    out.reg{ii}=creg;
    if numel(unique(creg))>1
        out.cross_reg(ii)=true;
    end
end

if opt.cross_only
    cross=struct();
    for fn=reshape(fieldnames(out),1,[])
        cross.(fn{1})=out.(fn{1})(out.cross_reg);
    end
    % out_all=out;
    out=cross;
end

end

function chains=extend_chain(chains,sessid,wvtype,dur,sub_chain)
chains=[chains;...
    num2cell(repmat(sessid,numel(sub_chain),1)),...
    repmat({wvtype+dur},numel(sub_chain),1),...
    num2cell(repmat(dur,numel(sub_chain),1)),...%duration
    sub_chain];
end



function chains=chain_one(sess_con,cids,tcoms,reverse)
arguments
    sess_con
    cids
    tcoms
    reverse (1,1) logical = false
end

chains=cell(0);
if nnz(cids)<3,return;end
sigcon=sess_con(all(ismember(sess_con,cids),2),:);
if numel(unique(sigcon))<3,return;end
[~,idces]=ismember(sigcon,cids);
con_com_prepost=cell2mat(tcoms(idces));
if reverse
    dirsel=con_com_prepost(:,2)<=con_com_prepost(:,1);
else
    dirsel=con_com_prepost(:,2)>=con_com_prepost(:,1); % will not rule out same-region with >= % <=
end
dirsigcon=sigcon(dirsel,:);
upre=unique(dirsigcon(:,1)).';

for i=upre
    onechain=cell(0);
    appeared=[];
    cpre=i;
    while true % first pass-through without unfolding all chains
        newpair=dirsigcon(ismember(dirsigcon(:,1),cpre) & ~ismember(dirsigcon(:,2),appeared),:);
        appeared=unique([appeared;newpair(:)]);
        if isempty(newpair)
            % cross-region check in main thread
            if numel(onechain)>1
                chains=[chains;{onechain}];
            end
            break
        else
            [~,pidces]=ismember(newpair,cids);
            onechain{end+1}={newpair,cell2mat(tcoms(pidces))};
            cpre=newpair(:,2);
        end
    end
end
end
