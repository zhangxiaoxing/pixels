%TODO Waveform based SU filter

function out=load_meta(opt)
arguments
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO','MY'})}='neupix'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
    opt.delay (1,1) double {mustBeMember(opt.delay,[3,6])} = 6
    opt.fdr (1,1) logical = true
end
persistent opt_ out_

if opt.delay==6
    warning('Delay set to default 6')
end

if isempty(out_) || ~isequaln(opt,opt_)
    if strcmp(opt.type,'neupix') || strcmp(opt.type,'MY')
        homedir=ephys.util.getHomedir();
        if strcmp(opt.criteria,'WT'),fpath=fullfile(homedir,sprintf('transient_%d.hdf5',opt.delay)); % source data from K:\code\per_sec\per_sec_stats.py
        elseif strcmp(opt.criteria,'Learning'),fpath=fullfile(homedir,sprintf('transient_%d_complete.hdf5',opt.delay));end

        meta_str.trial_counts=h5read(fpath,'/trial_counts');
        meta_str.wrs_p=h5read(fpath,'/wrs_p');
        meta_str.selec=h5read(fpath,'/selectivity');
        meta_str.allpath=cellstr(deblank(h5read(fpath,'/path')));
        meta_str.allcid=h5read(fpath,'/cluster_id');
        meta_str.reg_tree=cellstr(deblank(h5read(fpath,'/reg_tree')));
        meta_str.good_waveform=h5read(fpath,'/wf_good');
        %         meta_str.mem_type=h5read(fpath,'/mem_type');
        if ~opt.fdr
            [meta_str.mem_type,meta_str.per_bin]=ephys.get_mem_type(meta_str.wrs_p,meta_str.selec,'delay',opt.delay,'fdr',false);
        end
    else
        ccftree=deblank(h5read('K:\neupix\AIOPTO\META\Selectivity_AIopto_0419.hdf5','/reg'));
        meta_str.reg_tree=ccftree(3:8,:);
        meta_str.mem_type=hem2memtype(h5read('K:\neupix\AIOPTO\META\Selectivity_AIopto_0419.hdf5','/sus_trans_noPermutaion'));
        fullpath=deblank(h5read('K:\neupix\AIOPTO\META\Selectivity_AIopto_0419.hdf5','/path'));
        meta_str.allpath=regexp(fullpath,'.*(?=\\.*)','match','once');
        meta_str.allcid=h5read('K:\neupix\AIOPTO\META\Selectivity_AIopto_0419.hdf5','/cluster_id');
    end

    out=meta_str;
    out.sess=cellfun(@(x) ephys.path2sessid(x),out.allpath);
    if opt.fdr
        if opt.delay==3
            fdr=cell2mat(arrayfun(@(x) mafdr(out.wrs_p(4:8,x),'BHFDR',true),1:numel(out.allcid),'UniformOutput',false));
        else
            fdr=cell2mat(arrayfun(@(x) mafdr(out.wrs_p(4:11,x),'BHFDR',true),1:numel(out.allcid),'UniformOutput',false));
        end
        out.fdr=fdr;
        [out.mem_type,out.per_bin]=ephys.get_mem_type(fdr,out.selec,'delay',opt.delay,'fdr',true);
    end
    opt_=opt;
    out_=out;
end
out=out_;
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

% Temporary script for waveform sanity check
% meta=ephys.util.load_meta();
% upath=unique(meta.allpath);
%
% for ii=1:numel(upath)
%     sesssel=strcmp(meta.allpath,upath{ii});
%     wf_good_rate=nnz(meta.good_waveform(sesssel))./nnz(sesssel);
%     rate_stats(ii)=wf_good_rate;
%     fprintf('%d, %.2f, %s\n',ii,wf_good_rate,upath{ii});
% end
