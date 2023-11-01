% global_init;
% su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
% % wrs_mux_meta=ephys.get_wrs_mux_meta();

%% %%%%%%%%%%%%%%%%%%%%%%%%% COM_map

load(fullfile('binary','wrs_mux_meta.mat'),'wrs_mux_meta');
load(fullfile('binary','su_meta.mat'),'su_meta');
trials_dict=behav.get_trials_dict('skip_save',true);
com_map=wave.get_pct_com_map(wrs_mux_meta,'curve',true,'odor_only',false);

tcom3_maps=struct();
tcom6_maps=struct();
grey_regs=ephys.getGreyRegs('range','grey');

for kkey=["olf","dur","mixed"]

    [fcom3.(kkey).collection,fcom3.(kkey).com_meta]=wave.per_region_COM(...
        com_map,'sel_type',kkey,'com_field','com3');
    ureg=intersect(grey_regs,fcom3.(kkey).collection(:,2));
    [~,tcidx]=ismember(ureg,fcom3.(kkey).collection(:,2));
    tcom3_maps.(kkey)=containers.Map(...
        ureg,num2cell(cellfun(@(x) x/4, fcom3.(kkey).collection(tcidx,1))));

    [fcom6.(kkey).collection,fcom6.(kkey).com_meta]=wave.per_region_COM(...
        com_map,'sel_type',kkey,'com_field','com6');
    ureg=intersect(grey_regs,fcom6.(kkey).collection(:,2));
    [~,tcidx]=ismember(ureg,fcom6.(kkey).collection(:,2));
    tcom6_maps.(kkey)=containers.Map(...
        ureg,num2cell(cellfun(@(x) x/4, fcom6.(kkey).collection(tcidx,1))));
end



%% %%%%%%%%%%%%%%%%%%%%%%%% gen chains

% reg_com_maps=cell2struct({tcom3_maps;tcom6_maps},{'tcom3_maps','tcom6_maps'});
load(fullfile("binary","sums_conn_10.mat"),'sums_conn_str');

chains=cell(0);
for fidx=1:numel(sums_conn_str)
    disp(fidx)
    ffpath=sums_conn_str(fidx).folder;
    dpath=regexp(ffpath,'(?<=SPKINFO[\\/]).*$','match','once');
    if isempty(dpath)
        dpath=ffpath;
    end
    sessid=ephys.path2sessid(dpath,'criteria','WT');
    sess_sel=su_meta.sess==sessid;
    sess_con=sums_conn_str(fidx).sig_con;
    %% 3s
    
    s1d3_cids=[su_meta.allcid(wrs_mux_meta.wave_id==1 & sess_sel ...
        & ismember(su_meta.reg_tree(5,:),tcom3_maps.mixed.keys).');...
        su_meta.allcid(wrs_mux_meta.wave_id==5 & sess_sel ...
        & ismember(su_meta.reg_tree(5,:),tcom3_maps.olf.keys).');...
        su_meta.allcid(wrs_mux_meta.wave_id==7 & sess_sel ...
        & ismember(su_meta.reg_tree(5,:),tcom3_maps.dur.keys).')];

    s1d3_reg_tcom=[tcom3_maps.mixed.values(...
        su_meta.reg_tree(5, wrs_mux_meta.wave_id==1 & sess_sel ...
        & ismember(su_meta.reg_tree(5,:),tcom3_maps.mixed.keys).')),...
        tcom3_maps.olf.values(...
        su_meta.reg_tree(5, wrs_mux_meta.wave_id==5 & sess_sel ...
        & ismember(su_meta.reg_tree(5,:),tcom3_maps.olf.keys).')),...
        tcom3_maps.dur.values(...
        su_meta.reg_tree(5, wrs_mux_meta.wave_id==7 & sess_sel ...
        & ismember(su_meta.reg_tree(5,:),tcom3_maps.dur.keys).'))];

    sub_chain=chain_one(sess_con,s1d3_cids,s1d3_reg_tcom);
    if ~isempty(sub_chain)
        chains=extend_chain(chains,sessid,"s1d",delay,sub_chain);
    end

    s2d3_cids=[su_meta.allcid(wrs_mux_meta.wave_id==3 & sess_sel ...
        & ismember(su_meta.reg_tree(5,:),tcom3_maps.mixed.keys).');...
        su_meta.allcid(wrs_mux_meta.wave_id==6 & sess_sel ...
        & ismember(su_meta.reg_tree(5,:),tcom3_maps.olf.keys).');...
        su_meta.allcid(wrs_mux_meta.wave_id==7 & sess_sel ...
        & ismember(su_meta.reg_tree(5,:),tcom3_maps.dur.keys).')];

    s2d3_reg_tcom=[tcom3_maps.mixed.values(...
        su_meta.reg_tree(5, wrs_mux_meta.wave_id==3 & sess_sel ...
        & ismember(su_meta.reg_tree(5,:),tcom3_maps.mixed.keys).')),...
        tcom3_maps.olf.values(...
        su_meta.reg_tree(5, wrs_mux_meta.wave_id==6 & sess_sel ...
        & ismember(su_meta.reg_tree(5,:),tcom3_maps.olf.keys).')),...
        tcom3_maps.dur.values(...
        su_meta.reg_tree(5, wrs_mux_meta.wave_id==7 & sess_sel ...
        & ismember(su_meta.reg_tree(5,:),tcom3_maps.dur.keys).'))];

    sub_chain=chain_one(sess_con,s2d3_cids,s2d3_reg_tcom);
    if ~isempty(sub_chain)
        chains=extend_chain(chains,sessid,"s2d",delay,sub_chain);
    end
    %% 6s
    s1d6_cids=[su_meta.allcid(wrs_mux_meta.wave_id==2 & sess_sel ...
        & ismember(su_meta.reg_tree(5,:),tcom6_maps.mixed.keys).');...
        su_meta.allcid(wrs_mux_meta.wave_id==5 & sess_sel ...
        & ismember(su_meta.reg_tree(5,:),tcom6_maps.olf.keys).');...
        su_meta.allcid(wrs_mux_meta.wave_id==8 & sess_sel ...
        & ismember(su_meta.reg_tree(5,:),tcom6_maps.dur.keys).')];

    s1d6_reg_tcom=[tcom6_maps.mixed.values(...
        su_meta.reg_tree(5, wrs_mux_meta.wave_id==2 & sess_sel ...
        & ismember(su_meta.reg_tree(5,:),tcom6_maps.mixed.keys).')),...
        tcom6_maps.olf.values(...
        su_meta.reg_tree(5, wrs_mux_meta.wave_id==5 & sess_sel ...
        & ismember(su_meta.reg_tree(5,:),tcom6_maps.olf.keys).')),...
        tcom6_maps.dur.values(...
        su_meta.reg_tree(5, wrs_mux_meta.wave_id==8 & sess_sel ...
        & ismember(su_meta.reg_tree(5,:),tcom6_maps.dur.keys).'))];

    sub_chain=chain_one(sess_con,s1d6_cids,s1d6_reg_tcom);
    if ~isempty(sub_chain)
        chains=extend_chain(chains,sessid,"s1d",delay,sub_chain);
    end

    s2d6_cids=[su_meta.allcid(wrs_mux_meta.wave_id==4 & sess_sel ...
        & ismember(su_meta.reg_tree(5,:),tcom6_maps.mixed.keys).');...
        su_meta.allcid(wrs_mux_meta.wave_id==6 & sess_sel ...
        & ismember(su_meta.reg_tree(5,:),tcom6_maps.olf.keys).');...
        su_meta.allcid(wrs_mux_meta.wave_id==8 & sess_sel ...
        & ismember(su_meta.reg_tree(5,:),tcom6_maps.dur.keys).')];

    s2d6_reg_tcom=[tcom6_maps.mixed.values(...
        su_meta.reg_tree(5, wrs_mux_meta.wave_id==4 & sess_sel ...
        & ismember(su_meta.reg_tree(5,:),tcom6_maps.mixed.keys).')),...
        tcom6_maps.olf.values(...
        su_meta.reg_tree(5, wrs_mux_meta.wave_id==6 & sess_sel ...
        & ismember(su_meta.reg_tree(5,:),tcom6_maps.olf.keys).')),...
        tcom6_maps.dur.values(...
        su_meta.reg_tree(5, wrs_mux_meta.wave_id==8 & sess_sel ...
        & ismember(su_meta.reg_tree(5,:),tcom6_maps.dur.keys).'))];

    sub_chain=chain_one(sess_con,s2d6_cids,s2d6_reg_tcom);
    if ~isempty(sub_chain)
        chains=extend_chain(chains,sessid,"s2d",delay,sub_chain);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% overlapping removal
% since 2023/10/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
overlapping=false(size(out.sess));
for sess=reshape(unique(out.sess),1,[])
    for dd=[3 6]
        sess_sel=out.sess==sess & out.dur==dd;
        cids=out.cids(sess_sel);
        chlen=cellfun(@(x) numel(x), cids);
        for ii=1:numel(cids)
            for jj=reshape(find(chlen>numel(cids{ii})),1,[])
                [in,pos]=ismember(cids{ii},cids{jj});
                if all(in) && max(diff(pos))==1
                    overlapping(subsref(find(sess_sel),substruct('()',{ii})))=true;
                end
            end
        end
    end
end
for fn=reshape(fieldnames(out),1,[])
    out.(fn{1})=out.(fn{1})(~overlapping);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dur/olf separation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff_mode=false(size(out.sess));
for sess=reshape(unique(out.sess),1,[])
    sess_sel=out.sess==sess;
    cids=out.cids(sess_sel);
    cid_wave_dict=dictionary(su_meta.allcid(su_meta.sess==sess),wrs_mux_meta.wave_id(su_meta.sess==sess));
    for ii=1:numel(cids)
        cid_wave=cid_wave_dict(cids{ii});
        if any(ismember(cid_wave,5:6)) && any(ismember(cid_wave,7:8))
            diff_mode(subsref(find(sess_sel),substruct('()',{ii})))=true;
        end
    end
end

for fn=reshape(fieldnames(out),1,[])
    out.(fn{1})=out.(fn{1})(~diff_mode);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


wave.chain_tag.tag(...
    out,3,'skip_save',false,'odor_only',false,...
    'extend_trial',true,'skip_ts_id',true,...
    'filename','chain_tag_all_trl_mux.mat','poolsize',8,'criteria','WT'); % per-spk association



%% %%%%%%% rings
load(fullfile('binary','rings_bz_vs_shuf.mat'))
for sess=1:size(rings,1)
    cid_wave_dict=dictionary(su_meta.allcid(su_meta.sess==sess),wrs_mux_meta.wave_id(su_meta.sess==sess));
    for ssize=1:3
        cid_wave=cid_wave_dict(rings{sess,ssize});
        for ii=1:numel(cids)
            cid_wave=cid_wave_dict(cids{ii});
            if any(ismember(cid_wave,5:6)) && any(ismember(cid_wave,7:8))
                diff_mode(subsref(find(sess_sel),substruct('()',{ii})))=true;
            end
        end
    end
end

% load(fullfile('binary','rings_tag_trl.mat'),'ssloop_trl')
% loopwaves=[];
% for sess=reshape(unique(ssloop_trl.session),1,[])
%     cid_wave_dict=dictionary(su_meta.allcid(su_meta.sess==sess),wrs_mux_meta.wave_id(su_meta.sess==sess));
%     sesssel=ssloop_trl.session==sess;
%     for sidx=reshape(find(sesssel),1,[])
%         cid_wave=cid_wave_dict(ssloop_trl.meta{sidx,2});
%         loopwaves=[loopwaves;cid_wave,nan(1,5-numel(cid_wave))];
%         if any(ismember(cid_wave,5:6)) && any(ismember(cid_wave,7:8))
%             keyboard();
%         end
%     end
% end
% 
%% WIP
bz.rings.rings_time_constant.stats(load_file=false,skip_save=false,compress=true,odor_only=false)
wave.replay.stats_tbl(skip_save=false,odor_only=false)




%%%%%%%%%%%%%%%%%%%%%% support function
%%%%%%%%%%%%%%%%%%%%%% gen_chain

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


