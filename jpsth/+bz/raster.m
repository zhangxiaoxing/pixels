function psth=raster(fidx,suids,prefsamp)
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
if numel(pref_trial_sel)>=3
    center=ceil(numel(pref_trial_sel)/2);
    tsel=center-1:center+1;
else
    tsel=1:numel(pref_trial_sel);
end

%%%%%%%% TODO highlight marker for fc events %%%%%%
% spk_sel=ismember(FT_SPIKE.trial{1},pref_trial_sel(tsel));
% G=findgroups(FT_SPIKE.trial{1}(spk_sel));
% plot(repmat(FT_SPIKE.time{1}(spk_sel),2,1),[G-0.5;G+0.3],'r-')
% spk_sel=ismember(FT_SPIKE.trial{2},pref_trial_sel(tsel));
% G=findgroups(FT_SPIKE.trial{2}(spk_sel));
% plot(repmat(FT_SPIKE.time{2}(spk_sel),2,1),[G-0.3;G+0.5],'b-')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yidx=1;
evtstats=zeros(1,3);
for onetrial=pref_trial_sel(tsel).'
    spk1=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==onetrial);
    spk2=FT_SPIKE.time{2}(FT_SPIKE.trial{2}==onetrial);
    deltaT=spk2-spk1.';
    evts=deltaT>=0.0008 & deltaT<0.01;
    evtstats(1)=evtstats(1)+nnz(evts);
    evtstats(2)=evtstats(2)+numel(spk1);
    evtstats(3)=evtstats(3)+numel(spk2);
    [evt1,evt2]=find(evts);
    if nnz(evts)>0
        plot(repmat(spk1(evt1),2,1),[yidx-0.4;yidx],'r-')
        plot(repmat(spk2(evt2),2,1),[yidx;yidx+0.4],'b-')
    end
    if numel(spk1)>nnz(evts)
        plot(repmat(setdiff(spk1,spk1(evt1)),2,1),[yidx-0.4;yidx],'-','Color',[1,0.7,0.7])
    end
    if numel(spk2)>nnz(evts)
        plot(repmat(setdiff(spk2,spk2(evt2)),2,1),[yidx;yidx+0.4],'-','Color',[0.7,0.7,1])
    end
    yidx=yidx+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


arrayfun(@(x) xline(x,'--k'),[0,1]);
xlim([1,7]);
ylim([0,yidx-0.6]);
set(gca,'XTick',1:2:7,'XTickLabel',0:2:6)
xlabel('Delay time (s)')
ylabel('Trial #')
title(sprintf('%.3f, %.3f',evtstats(1)./evtstats(2),evtstats(1)./evtstats(3)));
cfg=struct();
cfg.binsize=[0.25];
cfg.trials=pref_trial_sel;
cfg.latency=[1,7];
psth=ft_spike_psth(cfg,FT_SPIKE);
% keyboard
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

