% assuming slurm sbatch calls
% request session idx as i ,1:173 as of Mar. 2021
function x_corr_bz(i,opt)
arguments
    i (1,1) double
    opt.debug (1,1) logical = false
    opt.prefix (1,:) char = '0315'
end
if isfile(sprintf('%s_BZ_XCORR_duo_f%d.mat',opt.prefix,i))
    if isunix, quit(0); else, return; end
end
disp(i);
bz.util.dependency %data path and lib path dependency
sess=dir(fullfile(homedir,'**','spike_info.hdf5'));
bz.util.pause(i,'xcorrpause');

folder=sess(i).folder;
trials=h5read(fullfile(folder,'FR_All_1000.hdf5'),'/Trials');
SU_id=h5read(fullfile(folder,'FR_All_1000.hdf5'),'/SU_id');

if sum(trials(:,9))<40 || numel(SU_id)<2  % apply well-trained criteria
    if isunix, quit(0); else, return; end
end

cstr=h5info(fullfile(sess(i).folder,sess(i).name)); % probe for available probes
spkID=[];spkTS=[];
for prb=1:size(cstr.Groups,1) % concatenate same session data for cross probe function coupling
    prbName=cstr.Groups(prb).Name;
    spkID=cat(1,spkID,h5read(fullfile(sess(i).folder,sess(i).name),[prbName,'/clusters']));
    spkTS=cat(1,spkTS,h5read(fullfile(sess(i).folder,sess(i).name),[prbName,'/times']));
end

susel=ismember(spkID,SU_id); % data cleaning by FR and contam rate criteria %TODO optional waveform cleaning
spkID=double(spkID(susel));
spkTS=double(spkTS(susel));

mono=bz.sortSpikeIDz(spkTS,spkID); % adapted from English, Buzsaki, 2017

if opt.debug && false
    bz.util.plotCCG
end
save(sprintf('%s_BZ_XCORR_duo_f%d.mat',opt.prefix,i),'mono','-v7.3','folder')
if isunix, quit(0); else, return; end
end