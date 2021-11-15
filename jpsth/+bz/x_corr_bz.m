% assuming called by slurm sbatch 
% request session idx as i ,1:173 as of Mar. 2021
function x_corr_bz(fidx,prefix,opt)
arguments
    fidx (1,1) double
    prefix (1,:) char = 'BZWT'
    opt.debug (1,1) logical = false
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
    opt.negccg (1,1) logical = false
end
ephys.util.dependency('ft',false);
if isfile(sprintf('%s_BZ_XCORR_duo_f%d.mat',prefix,fidx))
    disp('File exist'); if isunix, quit(0); else, return; end
end
disp(fidx);
bz.util.pause(fidx,'xcorrpause');
[spkID,spkTS,~,~,folder]=ephys.getSPKID_TS(fidx,'criteria',opt.criteria);
if isempty(spkID)
    if isunix, quit(0); else, return; end
end
%% temporary validation
spkID=spkID(spkID==499 | spkID==498);
spkTS=spkTS(spkID==499 | spkID==498);
%%%%%%%%%%%%%%%%%%%%%%%%

mono=bz.sortSpikeIDz(spkTS,spkID,'negccg',opt.negccg); % adapted from English, Buzsaki, 2017
if opt.debug && false
    bz.util.plotCCG
end
if opt.negccg && strcmp(opt.criteria,'Learning'),   suffix='_Inhibitory_Learning';
elseif opt.negccg,                                  suffix='_Inhibitory';
elseif strcmp(opt.criteria,'Learning'),             suffix='_Learning';
else,                                               suffix='';
end

save(sprintf('%s_BZ_XCORR_duo_f%d%s.mat',prefix,fidx,suffix),'mono','-v7.3','folder')
if isunix, quit(0); else, return; end
end