function pixFlatFR(opt)
arguments
    opt.binsize (1,1) double = 1.0
    opt.writefile (1,1) logical = false
    opt.overwrite (1,1) logical = false
    opt.two_point_align (1,1) logical = false
end
%% external lib dependency
addpath(fullfile("..","..","Lib","npy-matlab","npy-matlab"))
ephys.util.dependency("buz",false,"ft",true)
%% constant
sps=30000; %sample per second

%% input from YCY's time-aligned spike file
[~,homedir_all]=ephys.util.getHomedir();
flist=[];
for ii=1:numel(homedir_all)
    flist=[flist;struct2table(dir(fullfile(homedir_all{ii},'**','sess_spk_t_id.npy')))];
end
% flist=dir(fullfile(rootdir,'**','spike_info.hdf5'));
unmatched=0;
for ii=1:height(flist)
    fprintf('=== %d of %d ===\n',ii,height(flist))
    if isfile(fullfile(flist.folder{ii},sprintf('FR_All_%4d.hdf5',opt.binsize*1000))) && ~opt.overwrite
        disp(strjoin({'skiped',flist.folder{ii}}));
        continue
    end
    
    %% per session behavioral data
    trials=readNPY(fullfile(flist.folder{ii},'sess_trls.npy'));
    if isempty(trials), continue,  end
    %% select SUs with low contam rate, high FR and good waveform, for all probes
    spks=readNPY(fullfile(flist.folder{ii},flist.name{ii}));
    if isempty(spks), continue; end
    if all(spks(:,1)==0) && any(spks(:,3)>9999) % aligned
        ref_prb=trials(1,1);
        aligned=all(cell2mat(transpose(accumarray(trials(:,1)+1,trials(:,2),[],@(x) {x}))),2);
        trials=trials(intersect(find(trials(:,1)==ref_prb),find(aligned)),3:10);
        trials=behav.procPerf(trials,'mode','all');
        suids=unique(spks(:,3));
        spkID=spks(:,3);
        spkTS=spks(:,2).*sps; % FIXME
        
        %% split trials with external lib
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
        %% export result as file
        if opt.writefile
            FR_File=fullfile(flist.folder{ii},sprintf('FR_All_%04d.hdf5',opt.binsize*1000));
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
    else
        unmatched=unmatched+1;
        warning("Trials do not match")
    end
end
disp("number of unmatched session")
disp(unmatched)
end