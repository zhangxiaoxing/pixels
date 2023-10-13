function trials_dict=get_trials_dict(opt)
arguments
    opt.skip_save (1,1) logical = false
end
trials_dict=dictionary([],cell(0));
for sessid=1:116
    trials_dict(sessid)={h5read(fullfile(ephys.util.getHomedir(),replace(ephys.sessid2path(sessid),'\',filesep()),'FR_All_1000.hdf5'),'/Trials')};
end
if ~opt.skip_save
    blame=vcs.blame();
    save(fullfile('binary','trials_dict.mat'),'trials_dict','blame');
end

