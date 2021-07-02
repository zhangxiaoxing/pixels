function [spkID_,spkTS_,trials_,SU_id_,folder_,FT_SPIKE_]=getSPKID_TS(fidx,opt)
arguments
    fidx (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(fidx,1)}
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
    opt.keep_trial (1,1) logical = false
end

persistent spkID spkTS trials SU_id folder fidx_ criteria_ FT_SPIKE

if isempty(fidx_)...
        || fidx ~= fidx_ ...
        || ~strcmp(criteria_,opt.criteria)
    homedir=ephys.util.getHomedir('type','raw');
    folder=replace(ephys.sessid2path(fidx,'criteria',opt.criteria),'\',filesep());
    trials=h5read(fullfile(homedir,folder,'FR_All_1000.hdf5'),'/Trials');
    SU_id=h5read(fullfile(homedir,folder,'FR_All_1000.hdf5'),'/SU_id');
    %     FR_All=h5read(fullfile(homedir,folder,'FR_All_1000.hdf5'),'/FR_All');
    spkID=[];spkTS=[];

    %% Behavior performance parameter controlled data retrival
    if (numel(SU_id)<2) ...
            || (strcmp(opt.criteria,'WT') && sum(trials(:,9))<40) ...  % apply well-trained criteria
            || (strcmp(opt.criteria,'Learning') && sum(trials(:,9))>=40)
        disp('Did not meet criteria');
        spkID_=[];spkTS_=[];trials_=trials;SU_id_=SU_id;folder_=folder;return;
    end
    
    cstr=h5info(fullfile(homedir,folder,'spike_info.hdf5')); % probe for available probes
    for prb=1:size(cstr.Groups,1) % concatenate same session data for cross probe function coupling
        prbName=cstr.Groups(prb).Name;
        spkID=cat(1,spkID,h5read(fullfile(homedir,folder,'spike_info.hdf5'),[prbName,'/clusters']));
        spkTS=cat(1,spkTS,h5read(fullfile(homedir,folder,'spike_info.hdf5'),[prbName,'/times']));
    end
    
    susel=ismember(spkID,SU_id); % data cleaning by FR and contam rate criteria 
    %TODO optional further cleaning by waveform
    spkID=double(spkID(susel));
    spkTS=double(spkTS(susel));
    
    if ~opt.keep_trial
        FT_SPIKE=[];
        return
    end

    ephys.util.dependency('ft',true,'buz',false); %data path and lib path dependency
    [G,ID]=findgroups(spkID);
    SP=splitapply(@(x) {x}, spkTS, G);
    FT_SPIKE=struct();
    FT_SPIKE.label=arrayfun(@(x) num2str(x),SU_id,'UniformOutput',false);
    FT_SPIKE.timestamp=SP(ismember(ID,SU_id));
    sps=30000;
    cfg=struct();
    cfg.trl=[trials(:,1)-3*sps,trials(:,1)+11*sps,zeros(size(trials,1),1)-3*sps,trials];
    cfg.trlunit='timestamps';
    cfg.timestampspersecond=sps;
    FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);
end

spkID_=spkID;
spkTS_=spkTS;
trials_=trials;
SU_id_=SU_id;
folder_=folder;
FT_SPIKE_=FT_SPIKE;
fidx_=fidx;
criteria_=opt.criteria;

end
