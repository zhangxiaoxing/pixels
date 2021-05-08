function [spkID_,spkTS_,trials_,SU_id_,folder_]=getSPKID_TS(fidx,opt)
arguments
    fidx (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(fidx,1)}
    %TODO EPOCH
    opt.epoch (1,:) char {mustBeMember(opt.epoch,{'delay','ITI','any'})} = 'any'
end

persistent spkID spkTS trials SU_id folder fidx_

if isempty(fidx_) || fidx ~= fidx_
    homedir=ephys.util.getHomedir('type','raw');
    folder=replace(ephys.sessid2path(fidx),'\',filesep());
    trials=h5read(fullfile(homedir,folder,'FR_All_1000.hdf5'),'/Trials');
    SU_id=h5read(fullfile(homedir,folder,'FR_All_1000.hdf5'),'/SU_id');
    %     FR_All=h5read(fullfile(homedir,folder,'FR_All_1000.hdf5'),'/FR_All');
    spkID=[];spkTS=[];
    if sum(trials(:,9))<40 || numel(SU_id)<2  % apply well-trained criteria
        trials=[];SU_id=[];return;
    end
    
    cstr=h5info(fullfile(homedir,folder,'spike_info.hdf5')); % probe for available probes
    for prb=1:size(cstr.Groups,1) % concatenate same session data for cross probe function coupling
        prbName=cstr.Groups(prb).Name;
        spkID=cat(1,spkID,h5read(fullfile(homedir,folder,'spike_info.hdf5'),[prbName,'/clusters']));
        spkTS=cat(1,spkTS,h5read(fullfile(homedir,folder,'spike_info.hdf5'),[prbName,'/times']));
    end
    
    susel=ismember(spkID,SU_id); % data cleaning by FR and contam rate criteria %TODO optional waveform cleaning
    spkID=double(spkID(susel));
    spkTS=double(spkTS(susel));
    if ~strcmp(opt.epoch,'any')
        sel=epochProcess(spkTS,trials,opt.epoch);
        spkID=spkID(sel);
        spkTS=spkTS(sel);
    end
end

spkID_=spkID;
spkTS_=spkTS;
trials_=trials;
SU_id_=SU_id;
folder_=folder;
fidx_=fidx;

end

function tssel=epochProcess(spkTS,trials,epochType)
tssel=false(size(spkTS));
switch epochType
    case 'delay'
        for i=1:size(trials,1)
            tssel(spkTS>=(trials(i,1)+45000) &...  //sample onset+1.5s, sample onset+3.5s
                spkTS<(trials(i,1)+105000))=true;
        end
    case 'ITI'
        for i=1:size(trials,1)
            tssel(spkTS>=(trials(i,1)-75000) &... //sample onset -2.5s, sample onset -0.5s
                spkTS<trials(i,1)-15000)=true;
        end
end
end
