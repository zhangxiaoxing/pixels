function [sig,pair]=get_meta(sig_id,pair_id_one_dir,fpath)
arguments
    sig_id (:,2) int32
    pair_id_one_dir (:,2) int32
    fpath (1,:) char
end

idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));
meta_str=init();
reg_map=containers.Map('KeyType','int32','ValueType','any'); %reg_map(su_id)=reg
wrsp_map=containers.Map('KeyType','int32','ValueType','any'); 
selec_map=containers.Map('KeyType','int32','ValueType','any');
% pc_stem=replace(regexp(fpath,'(?<=SPKINFO/).*$','match','once'),'/','\');
sess_idx=find(startsWith(meta_str.allpath,fpath));
for suidx=reshape(sess_idx,1,[])
    suid=meta_str.allcid(suidx);
    acrontree=meta_str.reg_tree(:,suidx);
    ccfid=nan(1,6);
    for i=1:numel(acrontree), if isempty(acrontree{i}), ccfid(i)=0;else, ccfid(i)=idmap.reg2ccfid(acrontree{i});end;end
    reg_map(suid)=int32(ccfid);
    wrsp_map(suid)=meta_str.wrs_p(:,suidx);
    selec_map(suid)=meta_str.selec(:,suidx);
end
sig=struct();
sig.suid=sig_id;
sig.reg=reshape(cell2mat(arrayfun(@(x) reg_map(x),sig_id,'UniformOutput',false)),[],6,2);
sig.wrsp=reshape(cell2mat(arrayfun(@(x) wrsp_map(x),sig_id,'UniformOutput',false)),[],14,2);
sig.selec=reshape(cell2mat(arrayfun(@(x) selec_map(x),sig_id,'UniformOutput',false)),[],14,2);

pair=struct();
pair.suid=pair_id_one_dir;
pair.reg=reshape(cell2mat(arrayfun(@(x) reg_map(x),pair_id_one_dir,'UniformOutput',false)),[],6,2);
pair.wrsp=reshape(cell2mat(arrayfun(@(x) wrsp_map(x),pair_id_one_dir,'UniformOutput',false)),[],14,2);
pair.selec=reshape(cell2mat(arrayfun(@(x) selec_map(x),pair_id_one_dir,'UniformOutput',false)),[],14,2);
end

function out=init()
persistent meta_str
if isempty(meta_str)
    homedir=fullfile('K:','code','per_sec');
    
    meta_str.trial_counts=h5read(fullfile(homedir,'transient_6.hdf5'),'/trial_counts');
    meta_str.wrs_p=h5read(fullfile(homedir,'transient_6.hdf5'),'/wrs_p');
    meta_str.selec=h5read(fullfile(homedir,'transient_6.hdf5'),'/selectivity');
    meta_str.allpath=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/path'));
    meta_str.allcid=h5read(fullfile(homedir,'transient_6.hdf5'),'/cluster_id');
    meta_str.reg_tree=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/reg_tree'));
end
out=meta_str;
end


