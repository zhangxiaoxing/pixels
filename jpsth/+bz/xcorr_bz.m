% assuming called by slurm sbatch 
% request session idx as i ,1:173 as of Mar. 2021
function xcorr_bz(fidx,opt)
arguments
    fidx (1,1) double
    opt.debug (1,1) logical = false
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
    opt.negccg (1,1) logical = false
end
ephys.util.dependency();

if opt.negccg,  suffix='_Inhibitory';
else,           suffix='';
end

if isfile(fullfile('binary','SC',sprintf('%s_BZ_XCORR_duo_f%d%s.mat',opt.criteria,fidx,suffix)))
    disp('File exist'); if isunix, quit(0); else, return; end
end
disp(fidx);
% bz.util.pause(fidx,'xcorrpause'); % TODO socket-based distributed thread control
[spkID,spkTS,~,~,folder]=ephys.getSPKID_TS(fidx,'criteria',opt.criteria);
if isempty(spkID)
    if isunix, quit(0); else, return; end
end

mono=bz.sortSpikeIDz(spkTS,spkID,'negccg',opt.negccg); % adapted from English, Buzsaki, 2017
if opt.debug && false
    bz.util.plotCCG
end

save(fullfile('binary','SC',sprintf('%s_BZ_XCORR_duo_f%d%s.mat',opt.criteria,fidx,suffix)),'mono','-v7.3','folder')
end
