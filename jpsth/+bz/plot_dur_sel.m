window=29:36;%5:20
%%

if ~exist('sc','var')
    dursel=wave.delay_dur_selective('window',5:6);
    %1:fileid,2:sessid,3:suid,4:s1p,5:s2p,6:sp,7:fr3,8:fr6,9:si,10:sie,11:auc
    sc=dursel(dursel(:,11)>0.7 & min(dursel(:,7:8),[],2)>2 & dursel(:,6)<0.05,:);
end

for pi=39%1:size(sc,1)

ii=sc(pi,1);
se=sc(pi,2);
su=sc(pi,3);


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

s3c=(trials(:,8)==3 & trials(:,9)>0 & trials(:,10)>0);
s6c=(trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0);
s3e=(trials(:,8)==3 & trials(:,9)>0 & trials(:,10)==0);
s6e=(trials(:,8)==6 & trials(:,9)>0 & trials(:,10)==0);

%% raster
s3ci=find(s3c);
s6ci=find(s6c);
s3ci=s3ci(1:20);
frh=figure('Color','w','Position',[32,32,225,215]);
hold on
for yy=1:10
    ts=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==s3ci(yy));
    plot(repmat(ts,2,1),repmat([yy-0.3;yy+0.3],1,numel(ts)),'-b');
end
for yy=1:10
    ts=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==s6ci(yy));
    plot(repmat(ts,2,1),repmat([yy-0.3;yy+0.3],1,numel(ts))+10,'-r');
end
ylim([0.5,20.5])
xlim([-1,8])
set(gca(),'XTick',1:3:8,'XTickLabel',0:3:6)
xlabel('Time (s)')
ylabel('Trial #')
arrayfun(@(x) xline(x,'--k'),[0,1,4,5,7,8])
exportgraphics(frh,'dur_sel_SC_raster.pdf','ContentType','vector');
%% psth
fr=FT_PSTH.trial;


ss=1;

% ci3=bootci(50,@(x) squeeze(mean(x)),fr(s3c,ss,:));
% ci6=bootci(50,@(x) squeeze(mean(x)),fr(s6c,ss,:));

mm3=squeeze(mean(fr(s3c,ss,:)));
mm6=squeeze(mean(fr(s6c,ss,:)));

sem3=squeeze(std(fr(s3c,ss,:),0,1)./sqrt(nnz(s3c)));
sem6=squeeze(std(fr(s6c,ss,:),0,1)./sqrt(nnz(s6c)));


fh=figure('Color','w','Position',[32,32,700,220]);
subplot(1,3,1)
hold on
fill([1:56,56:-1:1],[smooth(mm3+sem3);smooth(flip(mm3-sem3))],'b','FaceAlpha',0.1,'EdgeColor','none')
fill([1:56,56:-1:1],[smooth(sem6+mm6);smooth(flip(mm6-sem6))],'r','FaceAlpha',0.1,'EdgeColor','none')

ph3=plot(smooth(mm3),'-b');
ph6=plot(smooth(mm6),'-r');

arrayfun(@(x) xline(x,':k'),[3,4,7,8,10,11]*4+12.5)
xlim([19.5,56.6])
ylim([0,20])
% set(gca(),'XTick',16.5:20:44.5,'XTickLabel',0:5:5)
xlabel('Time (s)')
ylabel('FR (Hz)')
set(gca,'XTick',(8:20:56)+0.5,'XTickLabel',-5:5:15)

subplot(1,3,2)
hold on
mm=[mean(fr(s6c,ss,window),'all'),mean(fr(s3c,ss,window),'all'),mean(fr(s6e,ss,window),'all'),mean(fr(s3e,ss,window),'all')];
sems=[std(mean(fr(s6c,ss,window),3)),std(mean(fr(s3c,ss,window),3)),std(mean(fr(s6e,ss,window),3)),std(mean(fr(s3e,ss,window),3))]...
    ./sqrt([nnz(s6c),nnz(s3c),nnz(s6e),nnz(s3e)]);
bh=bar(1:4,diag(mm),'stacked','FaceColor','r','EdgeColor','k');
bh(2).FaceColor='b';
bh(3).FaceColor=[1,0.5,0.5];
bh(4).FaceColor=[0.5,0.5,1];
errorbar(1:4,mm,sems,'k.')
set(gca(),'XTick',1:4,'XTickLabel',{'6s correct','3s correct','6s error','3s error'})
ylabel('FR (Hz)')

pc=ranksum(reshape(fr(s6c,ss,window),1,[]),reshape(fr(s3c,ss,window),1,[]))
pe3=ranksum(reshape(fr(s6c,ss,window),1,[]),reshape(fr(s3e,ss,window),1,[]))
pe6=ranksum(reshape(fr(s6c,ss,window),1,[]),reshape(fr(s6e,ss,window),1,[]))

subplot(1,3,3)
hold on
plot(find(s3c),mean(fr(s3c,ss,window),3),'bo','MarkerFaceColor','b','MarkerSize',4);
plot(find(s6c),mean(fr(s6c,ss,window),3),'ro','MarkerFaceColor','r','MarkerSize',4);
plot(find(s3e),mean(fr(s3e,ss,window),3),'bo','MarkerSize',4);
plot(find(s6e),mean(fr(s6e,ss,window),3),'ro','MarkerSize',4);
disp(meta.reg_tree(:,meta.sess==se & meta.allcid==su))
% xlim([176.5,200.5])
% ylim([0,25])
% set(gca,'XTick',181:5:204,'XTickLabel',5:5:25)
xlabel('Trial #')
ylabel('Baseline FR (Hz)')
% set(gca,'XTick',192:10:240,'XTickLabel',0:10:60)
sgtitle(num2str([pi,ii,se,su]))
keyboard()
% exportgraphics(fh,'dur_sel_SC.pdf','ContentType','vector');
% exportgraphics(fh,fullfile('SC',sprintf('dur_sel_SC_ITI_%d.png',pi)),'ContentType','image');
%  exportgraphics(fh,fullfile('SC',sprintf('dur_sel_SC_%d.png',pi)),'ContentType','image');
% waitfor(fh)
close(fh)
end