function x_corr_delay_all(fidx,opt)
%  A peak at a negative lag for stat.xcorr(chan1,chan2,:) means that chan1 is leading
%  chan2. Thus, a negative lag represents a spike in the second dimension of
%  stat.xcorr before the channel in the third dimension of stat.stat.

arguments
    fidx (1,1) double
    opt.debug (1,1) logical = false
    opt.prefix (1,:) char = '0319'
    opt.delay (1,1) double {mustBeMember(opt.delay,[3,6])} = 6
    opt.bin_range (1,2) double = [0,1]
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
    opt.overwrite (1,1) logical =false
%     opt.overwrite (1,1) logical = false
end

if strcmp(opt.criteria,'any'), suffix='';
else, suffix=opt.criteria; end

if isfile(sprintf('%s_XCORR_duo_f%03d_delay_%d_%d_%d_%s_2msbin.mat',opt.prefix,fidx,opt.delay,opt.bin_range(1),opt.bin_range(2),suffix)) && ~opt.overwrite
    disp('File exist'); return;
end
% copied BZ logic, skipped cell selectivity type classification
disp(fidx);
if isunix,maxNumCompThreads(16);end
ephys.util.dependency %data path and lib path dependency

bz.util.pause(fidx,'xcorrpause');
[spkID,spkTS,trials,SU_id,folder]=ephys.getSPKID_TS(fidx,'criteria',opt.criteria);
if isempty(spkID)
    if isunix, quit(0); else, return; end
end
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
[xc_s1,xcshuf_s1,xc_s2,xcshuf_x2]=my.x_corr_my(FT_SPIKE,opt.delay,opt.bin_range,'criteria',opt.criteria);
sums={folder,xc_s1,xcshuf_s1,xc_s2,xcshuf_x2}; %per folder save
blame=vcs.blame();
save(sprintf('%s_XCORR_duo_f%03d_delay_%d_%d_%d_%s_2msbin.mat',...
    opt.prefix,fidx,opt.delay,opt.bin_range(1),opt.bin_range(2),suffix),...
    'sums','blame','-v7.3') %prefix
end

