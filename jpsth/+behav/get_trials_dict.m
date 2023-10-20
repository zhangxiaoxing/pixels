function trials_dict=get_trials_dict(opt)
arguments
    opt.skip_save (1,1) logical = false
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
end

trials_dict=dictionary([],cell(0));
[~,~,sesspathmap]=ephys.sessid2path(0,'criteria',opt.criteria);
for sessid=1:sesspathmap.length
    trials_dict(sessid)={h5read(fullfile(ephys.util.getHomedir(),...
        replace(ephys.sessid2path(sessid,'criteria',opt.criteria),("\"|"/"),filesep()),...
        'FR_All_1000.hdf5'),'/Trials')};
end
if ~opt.skip_save
    blame=vcs.blame();
    if strcmp(opt.criteria,'WT')
        save(fullfile('binary','trials_dict.mat'),'trials_dict','blame');
    else
        keyboard()
    end
end

