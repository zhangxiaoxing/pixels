function out=load_meta(opt)
arguments
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO'})}='neupix'
end
persistent meta_str currtype

if isempty(meta_str) || ~strcmp(currtype,opt.type)
    if strcmp(opt.type,'neupix')
        homedir=ephys.util.getHomedir();
        meta_str.trial_counts=h5read(fullfile(homedir,'transient_6.hdf5'),'/trial_counts');
        meta_str.wrs_p=h5read(fullfile(homedir,'transient_6.hdf5'),'/wrs_p');
        meta_str.selec=h5read(fullfile(homedir,'transient_6.hdf5'),'/selectivity');
        meta_str.allpath=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/path'));
        meta_str.allcid=h5read(fullfile(homedir,'transient_6.hdf5'),'/cluster_id');
        meta_str.reg_tree=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/reg_tree'));
        meta_str.mem_type=h5read(fullfile(homedir,'transient_6.hdf5'),'/mem_type');
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