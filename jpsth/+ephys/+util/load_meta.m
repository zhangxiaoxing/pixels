%TODO Waveform based SU filter

function out=load_meta(opt)
arguments
    %     opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO','MY'})}='neupix'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
    opt.n_bin (1,1) double {mustBeInteger,mustBePositive} = 3
    opt.skip_stats (1,1) logical = true
end
persistent opt_ out_

if isempty(out_) || ~isequaln(opt,opt_)
    homedir=ephys.util.getHomedir();
    fpath_6=fullfile(homedir,'transient_6.hdf5'); % source data from K:\code\per_sec\per_sec_stats.py
    out_.allpath=cellstr(deblank(h5read(fpath_6,'/path')));
    out_.allcid=h5read(fpath_6,'/cluster_id');
    out_.reg_tree=cellstr(deblank(h5read(fpath_6,'/reg_tree')));
    out_.good_waveform=h5read(fpath_6,'/wf_good');
    if ~opt.skip_stats
        fpath_3=fullfile(homedir,'transient_3.hdf5');
        out_.trial_counts=h5read(fpath_6,'/trial_counts');
        out_.wrs_p_6=h5read(fpath_6,'/wrs_p');
        out_.selec_6=h5read(fpath_6,'/selectivity');

        out_.wrs_p_3=h5read(fpath_3,'/wrs_p');
        out_.selec_3=h5read(fpath_3,'/selectivity');

        out_.fdr_3=cell2mat(arrayfun(@(x) mafdr(out_.wrs_p_3(4:8,x),'BHFDR',true),1:numel(out_.allcid),'UniformOutput',false));
        out_.fdr_6=cell2mat(arrayfun(@(x) mafdr(out_.wrs_p_6(4:11,x),'BHFDR',true),1:numel(out_.allcid),'UniformOutput',false));

        [out_.mem_type_3,out_.per_bin_3]=ephys.get_mem_type(out_.fdr_3,out_.selec_3,'delay',3,'fdr',true);
        [out_.mem_type_6,out_.per_bin_6]=ephys.get_mem_type(out_.fdr_6,out_.selec_6,'delay',6,'fdr',true);
    end
    out_.sess=cellfun(@(x) ephys.path2sessid(x),out_.allpath);
    opt_=opt;
end
out=out_;
end
