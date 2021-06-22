function raster(fidx,suids,prefsamp)
ephys.util.dependency('buz',false);
%Fieldtrip routine
[spkID,spkTS,trials,~,~]=ephys.getSPKID_TS(fidx);
[G,ID]=findgroups(spkID);
SP=splitapply(@(x) {x}, spkTS, G);
FT_SPIKE=struct();
FT_SPIKE.label=arrayfun(@(x) num2str(x),suids,'UniformOutput',false);
FT_SPIKE.timestamp=SP(ismember(ID,suids));
sps=30000;
cfg=struct();
cfg.trl=[trials(:,1)-3*sps,trials(:,1)+11*sps,zeros(size(trials,1),1)-3*sps,trials];
cfg.trlunit='timestamps';
cfg.timestampspersecond=sps;
FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);
%%%%%%%
pref_trial_sel=find(trials(:,9)>0 & trials(:,10)>0 & trials(:,8)==6 & trials(:,5)==prefsamp);
pref_trial_sel=pref_trial_sel(1:15);
g
spk_sel=ismember(FT_SPIKE.trial{1},pref_trial_sel);
G=findgroups(FT_SPIKE.trial{1}(spk_sel));
plot(repmat(FT_SPIKE.time{1}(spk_sel),2,1),[G-0.5;G+0.3],'r-')
spk_sel=ismember(FT_SPIKE.trial{2},pref_trial_sel);
G=findgroups(FT_SPIKE.trial{2}(spk_sel));
plot(repmat(FT_SPIKE.time{2}(spk_sel),2,1),[G-0.3;G+0.5],'b-')
arrayfun(@(x) xline(x,'--k'),[0,1]);
xlim([1,7]);
ylim([0,16]);
set(gca,'XTick',[0,5])
% 
% subplot(2,1,1);
% hold on;
% 
% cfg=[];
% cfg.binsize=[0.25];
% cfg.latency=[-1,7];
% cfg.keeptrials='yes';
% cfg.trials=s1trialAll;
% psth1=ft_spike_psth(cfg,spkTrial);
% cfg.trials=s2trialAll;
% psth2=ft_spike_psth(cfg,spkTrial);
% 
% 
% ci1=bootci(1000,@(x) squeeze(mean(x)),psth1.trial);
% ci2=bootci(1000,@(x) squeeze(mean(x)),psth2.trial);
% 
% fill([psth1.time,fliplr(psth1.time)],[smooth(ci1(1,:),3);flip(smooth(ci1(2,:),3))],'b','FaceAlpha',0.2,'EdgeColor','none');
% fill([psth2.time,fliplr(psth2.time)],[smooth(ci2(1,:),3);flip(smooth(ci2(2,:),3))],'r','FaceAlpha',0.2,'EdgeColor','none');
% 
% plot(psth1.time,smooth(psth1.avg,3),'-b','LineWidth',1)
% plot(psth2.time,smooth(psth2.avg,3),'-r','LineWidth',1)
% 
% arrayfun(@(x) xline(x,'--k'),[0,1]);
% set(gca,'XTick',[0,5])
% title(sprintf('F%d,B%d,U%d',fid,bin,suids));
end

