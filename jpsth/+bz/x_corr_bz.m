% assuming slurm sbatch calls
% request session idx as i ,1:173 as of Mar. 2021
function x_corr_bz(fidx,opt)
arguments
    fidx (1,1) double
    opt.debug (1,1) logical = false
    opt.prefix (1,:) char = '0315'
end
if isfile(sprintf('%s_BZ_XCORR_duo_f%d.mat',opt.prefix,fidx))
    disp('File exist'); if isunix, quit(0); else, return; end
end
disp(fidx);
bz.util.pause(fidx,'xcorrpause');
[spkID,spkTS,~,~,folder]=ephys.getSPKID_TS(fidx);
if isempty(spkID)
    if isunix, quit(0); else, return; end
end
mono=bz.sortSpikeIDz(spkTS,spkID); % adapted from English, Buzsaki, 2017
if opt.debug && false
    bz.util.plotCCG
end
save(sprintf('%s_BZ_XCORR_duo_f%d.mat',opt.prefix,fidx),'mono','-v7.3','folder')
if isunix, quit(0); else, return; end
end