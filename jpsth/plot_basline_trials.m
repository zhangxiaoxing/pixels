function plot_basline_trials
try
    evalin('base','ft_ready');
catch me    
    init()
end
su_list=importdata('..\baseline_stats.csv',',',0);
accu=0;
for idx=find(abs(su_list.data(:,1))>0.9)'
    accu=accu+1;
    disp(accu)
    fstr=su_list.textdata{idx,1};
    suid=su_list.textdata{idx,2};
    plot_su(fstr,suid,su_list.textdata{idx,3});
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

fig=figure('Color','w','Position',[100,100,600,800]);
p1=[0.15,0.5,0.7,0.40];
p2=[0.15,0.1,0.7,0.35];

rsps=zeros(size(trl,1),1);
rsps(trl(:,7)>0 & trl(:,10)==1)=1; %hit
rsps(trl(:,7)>0 & trl(:,10)==0)=2; %false
rsps(trl(:,7)<0 & trl(:,10)==1)=3; %CR
rsps(trl(:,7)<0 & trl(:,10)==0)=4; %miss

[axu,axd]=plotOne(ts,11,p1,p2,rsps);

[~,fp,~]=fileparts(f_path);
sgtitle(strjoin({currReg,replace(fp,'_cleaned',''),'SU#',num2str(suid)}),'Interpreter','none')
print(fig,strjoin({'bs_',currReg,replace(fp,'_cleaned',''),num2str(suid),'.png'},'_'),'-dpng','-painters');
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
FT_SPIKE.timestamp{1,1}=spkTS(spkId==str2num(suid));
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


function plotSegs(delayLen)
    vertLine=[0,1,1,2]+[0,0,ones(1,2).*delayLen];
    plot(repmat(vertLine,2,1),repmat(ylim()',1,length(vertLine)),'--k');
end