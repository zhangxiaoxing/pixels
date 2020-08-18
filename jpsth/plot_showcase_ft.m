target='K:\code\showcase\0621_xcorr_showcase_24_218_219_b1.png';
tok=regexp(target,'(?<=showcase_)(\d+)_(\d+)_(\d+)_b(\d+).png','tokens');
bin=str2double(tok{1}{4});
fid=str2double(tok{1}{1});
uid1=str2double(tok{1}{2});
uid2=str2double(tok{1}{3});
close all
fh=figure('Color','w','Position',[100,100,200,200]);
plotOne(bin,fid,uid1,0)
plotOne(bin,fid,uid2,1)
% print(fh,replace(target,'.png','_psth.png'),'-r300','-dpng');
exportgraphics(fh,replace(target,'.png','_psth.pdf'),'ContentType','vector')

function plotOne(bin,fid,uid1,pos)
% cd('~/pixels/jpsth')
if isunix
    homedir=fullfile('home','zx','neupix','wyt','Datasum');
elseif ispc
    homedir=fullfile('k:','neupix','wyt','Datasum');
end


addpath(fullfile('npy-matlab-master','npy-matlab'))
addpath('fieldtrip-20200320')
ft_defaults
load('bin_file_list.mat');

% udi2=10225;
fpath=file_bins{bin}{fid};

if ~isfolder(fullfile(homedir,fpath))
    homedir=replace(homedir,'Datasum',['Datasum',filesep,'singleProbe']);
    spkFolder=fullfile(homedir,fpath);
    metaFolder=fullfile(homedir,fpath);
    folderType=1;
else
    spkFolder=replace(fullfile(homedir,fpath),'imec1','imec0');
    if uid1>10000
        metaFolder=replace(fullfile(homedir,fpath),'imec0','imec1');
    else
        metaFolder=replace(fullfile(homedir,fpath),'imec1','imec0');
    end
    folderType=2;
end
[avail,spkTrial]=process(spkFolder,metaFolder,uid1,folderType);

s1trial=find(spkTrial.trialinfo(:,5)==4 & spkTrial.trialinfo(:,8)==6 & spkTrial.trialinfo(:,10));
s1trialAll=s1trial;
if length(s1trial)>20
    s1trial=s1trial(1:20);
end
s1sel=ismember(spkTrial.trial{1},s1trial);
s1trial_list=unique(spkTrial.trial{1}(s1sel));
rearr_trial_s1=zeros(size(spkTrial.trial{1}));
for i=1:length(s1trial_list)
    rearr_trial_s1(spkTrial.trial{1}==s1trial_list(i))=i;
end

s2trial=find(spkTrial.trialinfo(:,5)==8 & spkTrial.trialinfo(:,8)==6 & spkTrial.trialinfo(:,10));
s2trialAll=s2trial;
if length(s2trial)>20
    s2trial=s2trial(1:20);
end
s2sel=ismember(spkTrial.trial{1},s2trial);
s2trial_list=unique(spkTrial.trial{1}(s2sel));
rearr_trial_s2=zeros(size(spkTrial.trial{1}));
for i=1:length(s2trial_list)
    rearr_trial_s2(spkTrial.trial{1}==s2trial_list(i))=i;
end
rearr_trial_s2=rearr_trial_s2+2+length(s1trial_list);


cfg=[];
cfg.binsize=[0.25];
cfg.latency=[-1,7];
cfg.keeptrials='yes';
cfg.trials=s1trialAll;
psth1=ft_spike_psth(cfg,spkTrial);
cfg.trials=s2trialAll;
psth2=ft_spike_psth(cfg,spkTrial);

subplot(2,2,3+pos);
hold on;
plot(repmat(spkTrial.time{1}(s2sel),2,1),[rearr_trial_s2(s2sel)-0.5;rearr_trial_s2(s2sel)+0.5],'m-')
plot(repmat(spkTrial.time{1}(s1sel),2,1),[rearr_trial_s1(s1sel)-0.5;rearr_trial_s1(s1sel)+0.5],'c-')
arrayfun(@(x) xline(x,'--k'),[0,1]);
xlim([-1,7]);
ylim([0,length(s1trial)+length(s2trial)+3]);
set(gca,'XTick',[0,5])

subplot(2,2,1+pos);
hold on;

ci1=bootci(1000,@(x) squeeze(mean(x)),psth1.trial);
ci2=bootci(1000,@(x) squeeze(mean(x)),psth2.trial);

fill([psth1.time,fliplr(psth1.time)],[smooth(ci1(1,:),5);flip(smooth(ci1(2,:),5))],'b','FaceAlpha',0.2,'EdgeColor','none');
fill([psth2.time,fliplr(psth2.time)],[smooth(ci2(1,:),5);flip(smooth(ci2(2,:),5))],'r','FaceAlpha',0.2,'EdgeColor','none');

plot(psth1.time,smooth(psth1.avg,5),'-b','LineWidth',1)
plot(psth2.time,smooth(psth2.avg,5),'-r','LineWidth',1)

arrayfun(@(x) xline(x,'--k'),[0,1]);
set(gca,'XTick',[0,5])
title(sprintf('F%d,B%d,U%d',fid,bin,uid1));

end



function [avail,out]=process(spkFolder,metaFolder,ids,folderType)
sps=30000;
trials=clearBadPerf(h5read(fullfile(metaFolder,'events.hdf5'),'/trials')','corr');

if isempty(trials)
    avail=false;
    out=[];
    return
end
cluster_ids=ids;

%  single-unit candidate

if folderType==1
    spkTS=readNPY(fullfile(spkFolder,'spike_times.npy'));
    spkId=readNPY(fullfile(spkFolder,'spike_clusters.npy'));
elseif folderType==2
    fstr=load(fullfile(spkFolder,'spike_info.mat'));
    spkId=double([fstr.spike_info{1}{1};fstr.spike_info{1}{2}]);
    spkTS=double([fstr.spike_info{2}{1};fstr.spike_info{2}{2}]);
end
FT_SPIKE=struct();

FT_SPIKE.label=strtrim(cellstr(num2str(cluster_ids)));
FT_SPIKE.timestamp=cell(1,numel(cluster_ids));
for i=1:numel(cluster_ids)
    FT_SPIKE.timestamp{i}=spkTS(spkId==cluster_ids(i))';
end
%  continuous format F T struct file
cfg=struct();
cfg.trl=[trials(:,1)-3*sps,trials(:,1)+11*sps,zeros(size(trials,1),1)-3*sps,trials];
cfg.trlunit='timestamps';
cfg.timestampspersecond=sps;

FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);

out=FT_SPIKE;
avail=true;
end


function out=clearBadPerf(facSeq, model)
if strcmp(model, 'error')
    if length(facSeq)>=40
        errorsel=~xor(facSeq(:,5)==facSeq(:,6) , facSeq(:,7)>0);
        out=facSeq(errorsel,:);
    else
        out=[];
    end
else
    correct=(xor(facSeq(:,5)==facSeq(:,6) , facSeq(:,7)>0));
    facSeq(:,10)=correct;
    if length(facSeq)>=40
        facSeq(:,9)=0;
        i=40;
        while i<=length(facSeq)
            goodOff=nnz(xor(facSeq(i-39:i,5)==facSeq(i-39:i,6) , facSeq(i-39:i,7)>0));
            if goodOff>=30 %.75 correct rate
                facSeq(i-39:i,9)=1;
            end
            i=i+1;
        end
        out=facSeq(facSeq(:,9)==1,:);
    else
        out=[];
    end
end
end

