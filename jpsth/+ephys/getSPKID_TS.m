function [spkID_,spkTS_,trials_,SU_id_,folder_]=getSPKID_TS(fidx)
arguments
    fidx (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(fidx,1)}
end

persistent spkID spkTS trials SU_id folder fidx_

if isempty(fidx_) || fidx ~= fidx_
    homedir=ephys.util.getHomedir('type','raw');
    folder=replace(ephys.sessid2path(fidx),'\',filesep());
    trials=h5read(fullfile(homedir,folder,'FR_All_1000.hdf5'),'/Trials');
    SU_id=h5read(fullfile(homedir,folder,'FR_All_1000.hdf5'),'/SU_id');

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
end

spkID_=spkID;
spkTS_=spkTS;
trials_=trials;
SU_id_=SU_id;
folder_=folder;

end
