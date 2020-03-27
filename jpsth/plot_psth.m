function plot_psth
su_list=importdata('su_list.csv',',',0);
sus=importdata('transient_6.csv',',',1);
accu=0;
for idx=find(sus.data(1,:)==1)
    accu=accu+1;
        disp(accu)

        fstr=su_list.textdata{idx};
        suid=su_list.data(idx);
        plot_su(replace(fstr,'D:','K:'),suid);


end
end

function FT_SPIKE=plot_su(f_path,suid)
if ~evalin('base','ft_ready')
    init()
end
% rootlist=dir('K:\neupix\DataSum\**\spike_times.npy');
% 
% idx=contains({rootlist.folder},fstr);
% if nnz(idx)~=1
%     return
% end
% f=rootlist(idx);
tic

f=fopen(fullfile(f_path,'su_id2reg.csv'),'r');
su_reg=textscan(f,'%d %s %s','Delimiter',',','HeaderLines',1);
fclose(f);
currReg=su_reg{2}{su_reg{1}==suid};


[avail,FT_SPIKE]=pre_process(f_path,suid);
if ~avail
    return
end

% psth=getPsth(FT_SPIKE,num2str(suid));
trl=FT_SPIKE.trialinfo;
S1_6_sel=find(trl(:,5)==4 & trl(:,8)==6 & trl(:,9)==1 & trl(:,10)==1);
S2_6_sel=find(trl(:,5)==8 & trl(:,8)==6 & trl(:,9)==1 & trl(:,10)==1);
ts{1}=arrayfun(@(x) FT_SPIKE.time{1}(FT_SPIKE.trial{1}==x),S1_6_sel,'UniformOutput',false);
ts{2}=arrayfun(@(x) FT_SPIKE.time{1}(FT_SPIKE.trial{1}==x),S2_6_sel,'UniformOutput',false);
fig=figure('Color','w','Position',[100,100,1200,400]);
p1=[0.02,0.5,0.2,0.40];
p2=[0.02,0.1,0.2,0.35];
[axu,axd]=plotOne(ts,6,11,p1,p2);
title(axu,'6s correct');


S1_6_sel=find(trl(:,5)==4 & trl(:,8)==6 & trl(:,10)==0);
S2_6_sel=find(trl(:,5)==8 & trl(:,8)==6 & trl(:,10)==0);
ts{1}=arrayfun(@(x) FT_SPIKE.time{1}(FT_SPIKE.trial{1}==x),S1_6_sel,'UniformOutput',false);
ts{2}=arrayfun(@(x) FT_SPIKE.time{1}(FT_SPIKE.trial{1}==x),S2_6_sel,'UniformOutput',false);

p1=[0.27,0.5,0.2,0.40];
p2=[0.27,0.1,0.2,0.35];
[axu,axd]=plotOne(ts,6,11,p1,p2);
title(axu,'6s error');

S1_6_sel=find(trl(:,5)==4 & trl(:,8)==3 & trl(:,9)==1 & trl(:,10)==1);
S2_6_sel=find(trl(:,5)==8 & trl(:,8)==3 & trl(:,9)==1 & trl(:,10)==1);
ts{1}=arrayfun(@(x) FT_SPIKE.time{1}(FT_SPIKE.trial{1}==x),S1_6_sel,'UniformOutput',false);
ts{2}=arrayfun(@(x) FT_SPIKE.time{1}(FT_SPIKE.trial{1}==x),S2_6_sel,'UniformOutput',false);

p1=[0.52,0.5,0.2,0.40];
p2=[0.52,0.1,0.2,0.35];
[axu,axd]=plotOne(ts,3,11,p1,p2);
title(axu,'3s correct');

S1_6_sel=find(trl(:,5)==4 & trl(:,8)==3 & trl(:,10)==0);
S2_6_sel=find(trl(:,5)==8 & trl(:,8)==3 & trl(:,10)==0);
ts{1}=arrayfun(@(x) FT_SPIKE.time{1}(FT_SPIKE.trial{1}==x),S1_6_sel,'UniformOutput',false);
ts{2}=arrayfun(@(x) FT_SPIKE.time{1}(FT_SPIKE.trial{1}==x),S2_6_sel,'UniformOutput',false);

p1=[0.77,0.5,0.2,0.40];
p2=[0.77,0.1,0.2,0.35];
[axu,axd]=plotOne(ts,3,11,p1,p2);
title(axu,'3s error');

[~,fp,~]=fileparts(f_path);
sgtitle(strjoin({currReg,fp,num2str(suid)}),'Interpreter','none')
print(fig,strjoin({currReg,fp,num2str(suid),'.png'},'_'),'-dpng','-painters');
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
FT_SPIKE.label{1}=num2str(suid);
FT_SPIKE.timestamp{1,1}=spkTS(spkId==suid);
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



function [axu,axd]=plotOne(ts,delayLen,rLim,p1,p2)
%         pf=nan(0,0);
%         bn=nan(0,0);
axu=subplot('Position',p1);
hold on;

if size(ts{1},1)==0
elseif size(ts{1},1)==1
    plot([ts{1}{1};ts{1}{1}],repmat([21;21.8]+1,1,length(ts{1}{1})),'-r');
elseif size(ts{1},1)>20
    pos=round((length(ts{1})-20)/2);
    posRange=pos:pos+19;
    cellfun(@(x) plot([x{1};x{1}],repmat([x{2};x{2}+0.8],1,length(x{1})),'-r'),cellfun(@(x,y) {x,y},ts{1}(posRange),num2cell(1:length(posRange),1)','UniformOutput',false));
else
    posRange=1:size(ts{1},1);
    cellfun(@(x) plot([x{1};x{1}],repmat([x{2};x{2}+0.8],1,length(x{1})),'-r'),cellfun(@(x,y) {x,y},ts{1}(posRange),num2cell(1:length(posRange),1)','UniformOutput',false));
end



if size(ts{2},1)==0
elseif size(ts{2},1)==1
    plot([ts{2}{1};ts{2}{1}],repmat([21;21.8]+1,1,length(ts{2}{1})),'-b');
elseif size(ts{2},1)>20
    pos=round((length(ts{2})-20)/2);
    posRange=pos:pos+19;
    cellfun(@(x) plot([x{1};x{1}],repmat([x{2};x{2}+0.8],1,length(x{1})),'-b'),cellfun(@(x,y) {x,y},ts{2}(posRange),num2cell([1:length(posRange)]+24,1)','UniformOutput',false));
else
    posRange=1:size(ts{2},1);
    cellfun(@(x) plot([x{1};x{1}],repmat([x{2};x{2}+0.8],1,length(x{1})),'-b'),cellfun(@(x,y) {x,y},ts{2}(posRange),num2cell([1:length(posRange)]+24,1)','UniformOutput',false));
end


xlim([-1,rLim]);
count1=min(size(ts{1},1),20);
if count1<2
    count1=2;
end

count2=min(size(ts{2},1),20);
if count2<2
    count2=2;
end
set(gca,'XTick',[],'YTick',[1,count1,25,count2+24],'YTickLabel',[1,count1,1,count2]);

ylim([0,45]);
assignin('base','axT',gca());
plotSegs(delayLen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
axd=subplot('Position',p2);
binSize=0.25;
hold on;
if size(ts{1},1)==0
    pfHist=zeros(size(-4+binSize/2:binSize:rLim));
elseif size(ts{1},1)<=2
    pfHist=histcounts(ts{1}{1},-4:binSize:rLim)./binSize;
else
    pfHist=cell2mat(cellfun(@(x) histcounts(x,-4:binSize:rLim)./binSize,ts{1},'UniformOutput',false));
    cia=bootci(1000,@(x) mean(x), pfHist);
    fill([-4+binSize/2:binSize:rLim,rLim-binSize/2:-binSize:-4],[cia(1,:),fliplr(cia(2,:))],[1,0.8,0.8],'EdgeColor','none');
end
if size(ts{2},1)==0
    bnHist=zeros(size(-4+binSize/2:binSize:rLim));
elseif size(ts{2},1)<=2
    bnHist=histcounts(ts{2}{1},-4:binSize:rLim)./binSize;
else
    bnHist=cell2mat(cellfun(@(x) histcounts(x,-4:binSize:rLim)./binSize,ts{2},'UniformOutput',false));
    cib=bootci(1000,@(x) mean(x), bnHist);
    fill([-4+binSize/2:binSize:rLim,rLim-binSize/2:-binSize:-4],[cib(1,:),fliplr(cib(2,:))],[0.8,0.8,1],'EdgeColor','none');
end


plot(-4+binSize/2:binSize:rLim,(mean(pfHist,1))','-r');
plot(-4+binSize/2:binSize:rLim,(mean(bnHist,1))','-b');

xlim([-1,rLim]);
if min(size(ts{1},1),size(ts{2},1))>2
    ylim([min(ylim()),max([cia(:);cib(:)])]);
end
assignin('base','axB',gca());
plotSegs(delayLen);
set(gca,'XTick',[0,5,10]);
xlabel('time(s)')
% [~,lpwd,~]=fileparts(rootpath);
% lpwd=replace(lpwd,'_cleaned','');
% if ~exist(lpwd,'dir')
%     mkdir(lpwd)
% end
% print('-dpng','-painters',sprintf('%s\\%s_%03d_%d.png',lpwd,lpwd,uid,delayLen));

end


function plotSegs(delayLen)
    vertLine=[0,1,1,2]+[0,0,ones(1,2).*delayLen];
    plot(repmat(vertLine,2,1),repmat(ylim()',1,length(vertLine)),'--k');
end