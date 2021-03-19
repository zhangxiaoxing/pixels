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
end
if isfile(sprintf('%s_MY_XCORR_duo_f%d.mat',opt.prefix,fidx))
    disp('File exist'); if isunix, quit(0); else, return; end
end
%% copied BZ logic
disp(fidx);
ephys.util.dependency %data path and lib path dependency
sess=dir(fullfile(homedir,'**','spike_info.hdf5'));
[~,idces]=sort({sess.folder});sess=sess(idces);
bz.util.pause(fidx,'xcorrpause');

folder=sess(fidx).folder;
trials=h5read(fullfile(folder,'FR_All_1000.hdf5'),'/Trials');
SU_id=h5read(fullfile(folder,'FR_All_1000.hdf5'),'/SU_id');

if sum(trials(:,9))<40 || numel(SU_id)<2  % apply well-trained criteria
    if isunix, quit(0); else, return; end
end

cstr=h5info(fullfile(sess(fidx).folder,sess(fidx).name)); % probe for available probes
spkID=[];spkTS=[];
for prb=1:size(cstr.Groups,1) % concatenate same session data for cross probe function coupling
    prbName=cstr.Groups(prb).Name;
    spkID=cat(1,spkID,h5read(fullfile(sess(fidx).folder,sess(fidx).name),[prbName,'/clusters']));
    spkTS=cat(1,spkTS,h5read(fullfile(sess(fidx).folder,sess(fidx).name),[prbName,'/times']));
end

[G,ID]=findgroups(spkID);
SP=splitapply(@(x) {x}, spkTS, G);
FT_SPIKE=struct();
FT_SPIKE.label=arrayfun(@(x) num2str(x),SU_id,'UniformOutput',false);
FT_SPIKE.timestamp=SP(ismember(ID,SU_id));

%  continuous format F T struct file
sps=30000;
cfg=struct();
cfg.trl=[trials(:,1)-3*sps,trials(:,1)+11*sps,zeros(size(trials,1),1)-3*sps,trials];
cfg.trlunit='timestamps';
cfg.timestampspersecond=sps;
FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);
[xc_s1,xcshuf_s1,xc_s2,xcshuf_x2]=my.x_corr_my(FT_SPIKE,opt.delay,opt.bin_range);

sums={folder,xc_s1,xcshuf_s1,xc_s2,xcshuf_x2}; %per folder save
save(sprintf('%s_XCORR_duo_f3%d_delay_%d_%d_%d_2msbin.mat',opt.prefix,fidx,opt.delay,opt.bin_range(1),opt.bin_range(2)),'sums','-v7.3') %prefix
if isunix, quit(0); else, return; end
end

