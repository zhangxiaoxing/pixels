function out=load_meta(opt)
arguments
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO','MY'})}='neupix'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
end
persistent meta_str currtype criteria

if isempty(meta_str) || ~strcmp(currtype,opt.type) || ~strcmp(criteria, opt.criteria)
    if strcmp(opt.type,'neupix') || strcmp(opt.type,'MY')
        homedir=ephys.util.getHomedir();
        if strcmp(opt.criteria,'WT'),fpath=fullfile(homedir,'transient_6.hdf5');
        elseif strcmp(opt.criteria,'Learning'),fpath=fullfile(homedir,'transient_6_complete.hdf5');end
        
        meta_str.trial_counts=h5read(fpath,'/trial_counts');
        meta_str.wrs_p=h5read(fpath,'/wrs_p');
        meta_str.selec=h5read(fpath,'/selectivity');
        meta_str.allpath=deblank(h5read(fpath,'/path'));
        meta_str.allcid=h5read(fpath,'/cluster_id');
        meta_str.reg_tree=deblank(h5read(fpath,'/reg_tree'));
%         meta_str.mem_type=h5read(fpath,'/mem_type');
        [meta_str.mem_type,meta_str.per_bin]=ephys.get_mem_type(meta_str.wrs_p,meta_str.selec);
        currtype=opt.type;
    else
        ccftree=deblank(h5read('K:\neupix\AIOPTO\META\Selectivity_AIopto_0419.hdf5','/reg'));
        meta_str.reg_tree=ccftree(3:8,:);
        meta_str.mem_type=hem2memtype(h5read('K:\neupix\AIOPTO\META\Selectivity_AIopto_0419.hdf5','/sus_trans_noPermutaion'));
        fullpath=deblank(h5read('K:\neupix\AIOPTO\META\Selectivity_AIopto_0419.hdf5','/path'));
        meta_str.allpath=regexp(fullpath,'.*(?=\\.*)','match','once');
        meta_str.allcid=h5read('K:\neupix\AIOPTO\META\Selectivity_AIopto_0419.hdf5','/cluster_id');
        currtype=opt.type;
    end
    criteria=opt.criteria;
end
out=meta_str;
end

function memtype=hem2memtype(HEM)
% 0=NM,1=S1 sust, 2=S1 trans, 3=S2 sust, 4=S2 trans,-1=switched
delay_pref=max(HEM(7:end,:));
memtype=zeros(size(HEM,2),1);
memtype(HEM(1,:)==1 & delay_pref==1)=1;
memtype(HEM(1,:)==1 & delay_pref==2)=3;

memtype(HEM(2,:)==1 & delay_pref==1)=2;
memtype(HEM(2,:)==1 & delay_pref==2)=4;

memtype(HEM(4,:)~=0)=-1;
end