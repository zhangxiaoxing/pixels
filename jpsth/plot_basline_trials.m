function plot_basline_trials(state)
if ~exist('state','var')
    state=false;
end

try
    evalin('base','ft_ready');
catch me    
    init()
end
su_list=importdata('..\baseline_stats.csv',',',0);
accu=0;
if state
    for idx=find(abs(su_list.data(:,1))>0.9)'
        accu=accu+1;
        disp(accu)
        fstr=su_list.textdata{idx,1};
        suid=su_list.textdata{idx,2};
        plot_su(fstr,suid,su_list.textdata{idx,3});
    end
else
    for idx=find(abs(su_list.data(:,1))<0.2 & su_list.data(:,3)<0.001)'
        accu=accu+1;
        disp(accu)
        fstr=su_list.textdata{idx,1};
        suid=su_list.textdata{idx,2};
        plot_su_pair(replace(fstr,'K:','D:'),suid,su_list.textdata{idx,3});
    end
end
end



function FT_SPIKE=plot_su(f_path,suid,currReg)

tic

[avail,FT_SPIKE]=pre_process(f_path,suid);
if ~avail
    return
end

% psth=getPsth(FT_SPIKE,num2str(suid));
trl=FT_SPIKE.trialinfo;
ts=arrayfun(@(x) FT_SPIKE.time{1}(FT_SPIKE.trial{1}==x),1:size(trl,1),'UniformOutput',false);

fig=figure('Color','w','Position',[100,100,1200,800]);
p1=[0.15,0.5,0.7,0.40];
p2=[0.15,0.1,0.7,0.35];

rsps=zeros(size(trl,1),1);
rsps(trl(:,7)>0 & trl(:,10)==1)=1; %hit
rsps(trl(:,7)>0 & trl(:,10)==0)=2; %false
rsps(trl(:,7)<0 & trl(:,10)==1)=3; %CR
rsps(trl(:,7)<0 & trl(:,10)==0)=4; %miss

[axu,axd]=plotOne(ts,11,p1,p2,rsps);

[~,fp,~]=fileparts(f_path);
sgtitle(strjoin({currReg,replace(fp,'_cleaned',''),'SU#',suid}),'Interpreter','none')
print(fig,strjoin({'bs_',currReg,replace(fp,'_cleaned',''),suid,'.png'},'_'),'-dpng','-painters');
close(fig)
toc
end

function init()
addpath('npy-matlab-master\npy-matlab')
addpath('fieldtrip-20200320')
ft_defaults
assignin('base','ft_ready',true)
end


function [avail,out]=pre_process(f_path,suid)
sps=30000;
rootpath=f_path;
disp(rootpath)
trials=markWTPerf(h5read(fullfile(rootpath,'events.hdf5'),'/trials')');
if isempty(trials)
    avail=false;
    out=[];
    return
end
spkTS=readNPY(fullfile(rootpath,'spike_times.npy'));
spkId=readNPY(fullfile(rootpath,'spike_clusters.npy'));

FT_SPIKE=struct();
FT_SPIKE.label{1}=suid;
FT_SPIKE.timestamp{1,1}=spkTS(spkId==str2double(suid));
%  continuous format F T struct file
cfg=struct();
cfg.trl=[trials(:,1)-3*sps,trials(:,1)+11*sps,zeros(size(trials,1),1)-3*sps,trials];
cfg.trlunit='timestamps';
cfg.timestampspersecond=sps;

FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);

avail=true;
out=FT_SPIKE;

end


function psth=getPsth(spikeTrials,suid)
cfg         = [];
cfg.binsize =  [0.25];
cfg.latency = [-2 12];
cfg.spikechannel = spikeTrials.label(strcmp(spikeTrials.label,suid));
cfg.keeptrials = 'yes';
psth = ft_spike_psth(cfg,spikeTrials);
end


function [out,WTIdx]=markWTPerf(facSeq)
facSeq(:,9)=0;
i=40;
while i<=length(facSeq)
    goodOff=nnz(xor(facSeq(i-39:i,5)==facSeq(i-39:i,6) , facSeq(i-39:i,7)>0));
    if goodOff>=30 %.75 correct rate
        facSeq(i-39:i,9)=1;
    end
    i=i+1;
end
out=[facSeq,xor(facSeq(:,5)==facSeq(:,6) , facSeq(:,7)>0)];
WTIdx=9;
end



function [axu,axd]=plotOne(ts,rLim,p1,p2,trial_sel)
%         pf=nan(0,0);
%         bn=nan(0,0);
axu=subplot('Position',p1);
hold on;
ts=ts';
counts=cellfun(@(x) numel(x),ts);
skips=ceil(max(counts)/200);

posRange=1:size(ts,1);
cellfun(@(x) plot([x{1}(1:skips:end);x{1}(1:skips:end)],repmat([x{2};x{2}+0.8],1,ceil(length(x{1})/skips)),'-r'),cellfun(@(x,y) {x,y},ts(posRange),num2cell(1:length(posRange),1)','UniformOutput',false));
xlim([-2.5,rLim]);
ylim([0,size(ts,1)+1]);
set(gca,'XTick',[0,5,10]);
xlabel(sprintf('time (s), 1:%d downsample',skips));
ylabel('trial #');
arrayfun(@(x) xline(x,'--k'),[0,1,4,7]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
axd=subplot('Position',p2);
binSize=0.25;
hold on;
perTrial=cellfun(@(x) nnz(x<-0.5 & x>=-2.5)/2,ts);
hh=plot(find(trial_sel==1),perTrial(trial_sel==1),'.r'); %hit
hf=plot(find(trial_sel==2),perTrial(trial_sel==2),'.m'); %false
hc=plot(find(trial_sel==3),perTrial(trial_sel==3),'.b'); %cr
hm=plot(find(trial_sel==4),perTrial(trial_sel==4),'.c'); %miss

legend([hh,hf,hc,hm],{'hit','false','reject','miss'})
xlim([0,size(ts,1)+1]);
xlabel('trial #')
ylabel('baseline FR (Hz)')

end


function [FT_SPIKE,avail]=plot_su_pair(f_path,suid,currReg)


tic

[avail,FT_SPIKE]=pre_process(f_path,suid);
if ~avail
    return
end

% psth=getPsth(FT_SPIKE,num2str(suid));
trl=FT_SPIKE.trialinfo;
correct_sel=find(trl(:,9)==1 & trl(:,10)==1);
err_sel=find(trl(:,10)==0);
if nnz(correct_sel)<20 || nnz(err_sel)<20
    return
end

ts{1}=correct_sel;
ts{2}=arrayfun(@(x) FT_SPIKE.time{1}(FT_SPIKE.trial{1}==x),correct_sel,'UniformOutput',false);

ts{3}=err_sel;
ts{4}=arrayfun(@(x) FT_SPIKE.time{1}(FT_SPIKE.trial{1}==x),err_sel,'UniformOutput',false);

fig=figure('Color','w','Position',[100,100,600,800]);
p1=[0.15,0.5,0.7,0.40];
p2=[0.15,0.1,0.7,0.35];
[axu,axd]=plotOnePair(ts,11,p1,p2);


[~,fp,~]=fileparts(f_path);
sgtitle(strjoin({currReg,fp,'SU#',suid}),'Interpreter','none')
% print(fig,strjoin({currReg,fp,suid,'.png'},'_'),'-dpng','-painters');
close(fig)
toc
end


function [axu,axd]=plotOnePair(ts,rLim,p1,p2)
%         pf=nan(0,0);
%         bn=nan(0,0);
axu=subplot('Position',p1);
hold on;
ts=ts';
counts=cellfun(@(x) numel(x),[ts{1};ts{2}]);
skips=ceil(max(counts)/200);

count1=size(ts{2},1);
count2=size(ts{4},1);

arrayfun(@(x) plot([ts{2}{x,1:skips:end};ts{2}{x,1:skips:end}],repmat([ts{1}(x);ts{1}(x)+0.8],1,ceil(length(ts{2}{x})/skips)),'-r'),1:length(ts{1}));

arrayfun(@(x) plot([ts{4}{x,1:skips:end};ts{4}{x,1:skips:end}],repmat([ts{3}(x);ts{3}(x)+0.8],1,ceil(length(ts{4}{x})/skips)),'-b'),1:length(ts{3}));

xlim([-1,rLim]);

set(gca,'XTick',[]);

ylim([0,count1+count2]);
ylabel('trial #');
arrayfun(@(x) xline(x,'--k'),[0,1,4,7]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
axd=subplot('Position',p2);
binSize=0.25;
hold on;

pfHist=cell2mat(cellfun(@(x) histcounts(x,-4:binSize:rLim)./binSize,ts{2},'UniformOutput',false));
cia=bootci(1000,@(x) mean(x), pfHist);
fill([-4+binSize/2:binSize:rLim,rLim-binSize/2:-binSize:-4],[cia(1,:),fliplr(cia(2,:))],[1,0.8,0.8],'EdgeColor','none');


bnHist=cell2mat(cellfun(@(x) histcounts(x,-4:binSize:rLim)./binSize,ts{4},'UniformOutput',false));
cib=bootci(1000,@(x) mean(x), bnHist);
fill([-4+binSize/2:binSize:rLim,rLim-binSize/2:-binSize:-4],[cib(1,:),fliplr(cib(2,:))],[0.8,0.8,1],'EdgeColor','none');

plot(-4+binSize/2:binSize:rLim,(mean(pfHist,1))','-r');
plot(-4+binSize/2:binSize:rLim,(mean(bnHist,1))','-b');

xlim([-1,rLim]);

ylim([min(ylim()),max([cia(:);cib(:)])]);
ylabel('avearged FR (Hz)');
arrayfun(@(x) xline(x,'--k'),[0,1,4,7]);


set(gca,'XTick',[0,5,10]);
if skips > 1
    xlabel(sprintf('time (s), raster 1:%d downsample',skips))
else
    xlabel('time (s)')
end

end