%%
if ~exist('sc','var')
    dursel=wave.delay_dur_selective('window',1:2);
    sc=dursel(dursel(:,10)>0.7 & min(dursel(:,7:8),[],2)>2 & dursel(:,6)<0.05,:);
end

for pi=89

ii=sc(pi,1);
se=sc(pi,2);
su=sc(pi,3);



window=5:20;%17:26;

[spkID,spkTS,trials,suids,folder]=ephys.getSPKID_TS(se);

meta=ephys.util.load_meta('type','neupix','delay',6);

% homedir=ephys.util.getHomedir('type','raw');
% fl=dir(fullfile(homedir,'**','FR_All_ 250.hdf5'));
% 
% dpath=regexp(fl(ii).folder,'(?<=SPKINFO[\\/]).*$','match','once');
% fpath=fullfile(fl(ii).folder,fl(ii).name);
% fpath=replace(fpath,'\',filesep());
% pc_stem=replace(dpath,'/','\');

% fr=h5read(fpath,'/FR_All');
% trial=h5read(fpath,'/Trials');
% suid=h5read(fpath,'/SU_id');

%SPKTS -> PSTH
sps=30000;

FT_SPIKE=struct();

[G,ID]=findgroups(spkID);
SP=splitapply(@(x) {x}, spkTS, G);
FT_SPIKE=struct();
FT_SPIKE.label=arrayfun(@(x) num2str(x),su,'UniformOutput',false);
FT_SPIKE.timestamp=SP(ismember(ID,su));
if ispc, libroot='K:'; else, libroot='~'; end
addpath(fullfile(libroot,'Lib','fieldtrip-20200320'))
ft_defaults

cfg=struct();
cfg.trl=[trials(:,1)-6*sps,trials(:,1)+8*sps,zeros(size(trials,1),1)-6*sps,trials];
cfg.trlunit='timestamps';
cfg.timestampspersecond=sps;

FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);

cfg=struct();
cfg.binsize=0.25;
cfg.keeptrials='yes';
FT_PSTH=ft_spike_psth(cfg, FT_SPIKE);
fr=FT_PSTH.trial;

s3c=(trials(:,8)==3 & trials(:,9)>0 & trials(:,10)>0);
s6c=(trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0);
s3e=(trials(:,8)==3 & trials(:,9)>0 & trials(:,10)==0);
s6e=(trials(:,8)==6 & trials(:,9)>0 & trials(:,10)==0);

ss=1;

ci3=bootci(500,@(x) squeeze(mean(x)),fr(s3c,ss,:));
ci6=bootci(500,@(x) squeeze(mean(x)),fr(s6c,ss,:));

fh=figure('Color','w','Position',[32,32,465,220]);
subplot(1,2,1)
hold on
fill([1:56,56:-1:1],[smooth(ci3(1,:));smooth(fliplr(ci3(2,:)))],'r','FaceAlpha',0.1,'EdgeColor','none')
fill([1:56,56:-1:1],[smooth(ci6(1,:));smooth(fliplr(ci6(2,:)))],'b','FaceAlpha',0.1,'EdgeColor','none')

ph3=plot(smooth(squeeze(mean(fr(s3c,ss,:)))),'-r');
ph6=plot(smooth(squeeze(mean(fr(s6c,ss,:)))),'-b');

arrayfun(@(x) xline(x,':k'),[3,4,7,8,10,11]*4+12.5)
% xlim([8.5,44.5])
% ylim([2,25])
% set(gca(),'XTick',16.5:20:44.5,'XTickLabel',0:5:5)
xlabel('Time (s)')
ylabel('FR (Hz)')
set(gca,'XTick',(8:20:56)+0.5,'XTickLabel',-5:5:15)

subplot(1,2,2)
hold on
plot(find(s3c),mean(fr(s3c,ss,window),3),'ro','MarkerFaceColor','r','MarkerSize',4);
plot(find(s6c),mean(fr(s6c,ss,window),3),'bo','MarkerFaceColor','b','MarkerSize',4);
plot(find(s3e),mean(fr(s3e,ss,window),3),'ro','MarkerSize',4);
plot(find(s6e),mean(fr(s6e,ss,window),3),'bo','MarkerSize',4);
disp(meta.reg_tree(:,meta.sess==se & meta.allcid==su))
% xlim([176.5,200.5])
% ylim([0,25])
% set(gca,'XTick',181:5:204,'XTickLabel',5:5:25)
xlabel('Trial #')
ylabel('Baseline FR (Hz)')
set(gca,'XTick',192:10:240,'XTickLabel',0:10:60)
sgtitle(num2str([ii,se,su]))

exportgraphics(fh,'dur_sel_SC_ITI.pdf','ContentType','vector');
% exportgraphics(fh,fullfile('SC',sprintf('dur_sel_SC_ITI_%d.png',pi)),'ContentType','image');
% close(fh)
end