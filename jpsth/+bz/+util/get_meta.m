%TODO move to ephys package
function [sig,pair]=get_meta(sig_id,pair_id_one_dir,fpath,opt)
arguments
    sig_id (:,2) int32
    pair_id_one_dir (:,2) int32
    fpath (1,:) char
    %     opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO','MY'})}='neupix'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'

end
if ispc
    idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
elseif isunix
    idmap=load(fullfile('~/pixels','align','reg_ccfid_map.mat'));
end
idmap.reg2ccfid('')=0;
meta=ephys.util.load_meta('criteria',opt.criteria);
% waveid=ephys.get_wave_id(meta.sess,meta.allcid);
% reg_map=containers.Map('KeyType','int32','ValueType','any'); %reg_map(su_id)=reg
wf_map=containers.Map('KeyType','int32','ValueType','any');
% wave_map=containers.Map('KeyType','int32','ValueType','int32');
sess_idx=find(startsWith(meta.allpath,replace(fpath,'/','\')));
for suidx=reshape(sess_idx,1,[])
    suid=meta.allcid(suidx);
%     ccfid=cell2mat(idmap.reg2ccfid.values(meta.reg_tree(:,suidx)));
%     reg_map(suid)=int32(ccfid);
%     wave_map(suid)=waveid(suidx);
    wf_map(suid)=meta.good_waveform(suidx);
end
sig=struct();
sig.suid=sig_id;
% tree=arrayfun(@(x) reg_map(x),sig_id,'UniformOutput',false);
% sig.reg=permute(cat(3,cell2mat(tree(:,1).'),cell2mat(tree(:,2).')),[2,1,3]);
% sig.waveid=cell2mat(arrayfun(@(x) wave_map(x),sig_id,'UniformOutput',false));
sig.wf_good=cell2mat(arrayfun(@(x) wf_map(x),sig_id,'UniformOutput',false));

pair=struct();
if numel(pair_id_one_dir)>2
    pair.suid=pair_id_one_dir;
%     tree=arrayfun(@(x) reg_map(x),pair_id_one_dir,'UniformOutput',false);
%     pair.reg=permute(cat(3,cell2mat(tree(:,1).'),cell2mat(tree(:,2).')),[2,1,3]);    % pair.wrsp=reshape(cell2mat(arrayfun(@(x) wrsp_map(x),pair_id_one_dir,'UniformOutput',false)),[],14,2);
%     pair.waveid=cell2mat(arrayfun(@(x) wave_map(x),pair_id_one_dir,'UniformOutput',false));
    pair.wf_good=cell2mat(arrayfun(@(x) wf_map(x),pair_id_one_dir,'UniformOutput',false));
end
end
