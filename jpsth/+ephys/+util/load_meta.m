%TODO Waveform based SU filter

function su_meta=load_meta(opt)
arguments
    %     opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO','MY'})}='neupix'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
    opt.n_bin (1,1) double {mustBeInteger,mustBePositive} = 3
    opt.skip_stats (1,1) logical = true
    opt.adjust_white_matter (1,1) logical
    opt.save_file (1,1) logical = false
end
persistent opt_ out_
global gather_config

if ~isempty(gather_config)
    opt.adjust_white_matter=gather_config.adjust_white_matter;
end
assert(isfield(opt,'adjust_white_matter') && ~isempty(opt.adjust_white_matter),...
    'Undefined white matter SU adjust method')

if isempty(out_) || ~isequaln(opt,opt_)
    homedir=ephys.util.getHomedir();
    fpath_6=fullfile(homedir,'transient_6.hdf5'); % source data from K:\code\per_sec\per_sec_stats.py
    out_.allpath=cellstr(deblank(h5read(fpath_6,'/path')));
    out_.allcid=h5read(fpath_6,'/cluster_id');
    out_.reg_tree=ephys.get_adjusted_reg_tree('adjust_white_matter',opt.adjust_white_matter);
    out_.good_waveform=h5read(fpath_6,'/wf_good');

    if ~opt.skip_stats
        fpath_3=fullfile(homedir,'transient_3.hdf5');
        out_.trial_counts=h5read(fpath_6,'/trial_counts');
        out_.wrs_p_6=h5read(fpath_6,'/wrs_p');
        out_.selec_6=h5read(fpath_6,'/selectivity');

        out_.wrs_p_3=h5read(fpath_3,'/wrs_p');
        out_.selec_3=h5read(fpath_3,'/selectivity');
    end
    out_.sess=cellfun(@(x) ephys.path2sessid(x),out_.allpath);
    opt_=opt;
    
end
su_meta=out_;
if opt.save_file
    blame=vcs.blame();
    save(fullfile('binary','su_meta.mat'),'su_meta','blame');
end

end

