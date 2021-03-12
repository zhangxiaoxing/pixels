function pixFlatFR(opt)
arguments
    opt.binsize (1,1) double = 1.0
    opt.writefile (1,1) logical = false
    opt.rootdir (1,:) char = 'K:\neupix\SPKINFO'
    opt.overwrite (1,1) logical = false
end

if ispc, libroot='K:'; else, libroot='~'; end
addpath(fullfile(libroot,'Lib','fieldtrip-20200320'))
ft_defaults
sps=30000;

flist=dir(fullfile(opt.rootdir,'**','spike_info.hdf5'));
for i=1:length(flist)
    if isfile(fullfile(flist(i).folder,sprintf('FR_All_%4d.hdf5',opt.binsize*1000))) && ~opt.overwrite
        disp(strjoin({'skiped',flist(i).folder}));
        continue
    end
        
    trials=h5read(fullfile(replace(flist(i).folder,'SPKINFO','META'),'events.hdf5'),'/trials')';
    trials=behav.procPerf(trials,'mode','all');
    if isempty(trials), continue,  end
    
    cstr=h5info(fullfile(flist(i).folder,flist(i).name));
    fr_good=ephys.goodCid(replace(flist(i).folder,'SPKINFO','META')); % Good firing rate
    fr_wf_good=ephys.waveform.goodWaveform(replace(flist(i).folder,'SPKINFO','WF'),'presel',fr_good); %Good waveform
    if isempty(fr_wf_good)
        disp(flist(i).folder)
        keyboard()
        continue
    end
    spkID=[];spkTS=[];
    for prb=1:size(cstr.Groups,1)
        prbName=cstr.Groups(prb).Name;
        spkID=cat(1,spkID,h5read(fullfile(flist(i).folder,flist(i).name),[prbName,'/clusters']));
        spkTS=cat(1,spkTS,h5read(fullfile(flist(i).folder,flist(i).name),[prbName,'/times']));
    end
    susel=ismember(spkID,fr_wf_good);
    spkID=double(spkID(susel));
    spkTS=double(spkTS(susel));
    suids=unique(spkID);
    
    FT_SPIKE=struct();
    
    FT_SPIKE.label=strtrim(cellstr(num2str(suids)));
    FT_SPIKE.timestamp=cell(1,numel(suids));
    for su=1:numel(suids)
        FT_SPIKE.timestamp{su}=spkTS(spkID==suids(su))';
    end
    %  continuous format F T struct file
    cfg=struct();
    cfg.trl=[trials(:,1)-3*sps,trials(:,1)+11*sps,zeros(size(trials,1),1)-3*sps,trials];
    cfg.trlunit='timestamps';
    cfg.timestampspersecond=sps;
    
    FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);
    
    cfg=struct();
    cfg.binsize=opt.binsize;
    cfg.keeptrials='yes';
    FT_PSTH=ft_spike_psth(cfg, FT_SPIKE);
    
    if opt.writefile
        FR_File=fullfile(flist(i).folder,sprintf('FR_All_%4d.hdf5',opt.binsize*1000));
        if exist(FR_File,'file')
            delete(FR_File)
        end
        h5create(FR_File,'/FR_All',size(FT_PSTH.trial),'Datatype','double')
        h5write(FR_File,'/FR_All',FT_PSTH.trial)
        h5create(FR_File,'/Trials',size(FT_PSTH.trialinfo),'Datatype','double')
        h5write(FR_File,'/Trials',FT_PSTH.trialinfo)
        h5create(FR_File,'/SU_id',size(suids),'Datatype','double')
        h5write(FR_File,'/SU_id',suids)
    else 
        keyboard()% for devp
    end
end
end