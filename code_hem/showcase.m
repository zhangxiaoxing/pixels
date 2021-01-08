clc
clear
close all
addpath('D:\code\fieldtrip-20200320')
ft_defaults
%%
homedir='I:\WT';

SU_id=216;
folder='191108-DPA_31_g0\191108-DPA_31_g0_imec1_cleaned';

[folderType,~,spkFolder,metaFolder,~]=jointFolder(folder,cell(0),homedir);
if folderType>1 && contains(folder,'imec1')
    SU_id=SU_id+10000;
end
[~,spktrial]=pre_process(folderType,spkFolder,metaFolder,SU_id,0,0,'single');
trials=spktrial.trialinfo;
for j=1:size(spktrial.trialtime,1)
    FR(j,:)=histcounts(spktrial.time{1}(1,spktrial.trial{1,1}==j),spktrial.trialtime(1,1):0.25:spktrial.trialtime(1,2))./0.25;
    trial_spike{j,:}=spktrial.time{1}(1,spktrial.trial{1,1}==j);
end
%% Transient 3s single unit(Correct)
FA_6s=FR(trials(:,8)==6&trials(:,5)==4&trials(:,9)==1&trials(:,10)==1,:);
FB_6s=FR(trials(:,8)==6&trials(:,5)==8&trials(:,9)==1&trials(:,10)==1,:);

RasterA_6s=trial_spike(trials(:,8)==6&trials(:,5)==4&trials(:,9)==1&trials(:,10)==1,1);
RasterB_6s=trial_spike(trials(:,8)==6&trials(:,5)==8&trials(:,9)==1&trials(:,10)==1,1);

%% Plot
close all
fh=figure('Color','w','Position',[100,100,200,200]);

num=10; % trial plot num
b=10;
subplot(2,1,2)
for itr =1:num % trial
    if isempty(RasterA_6s{itr+b,:})
        RasterA_6s{itr+b,:}=-3;
    end
    if isempty(RasterB_6s{itr+b,:})
        RasterB_6s{itr+b,:}=-3;
    end
    plot([RasterA_6s{itr+b,:};RasterA_6s{itr+b,:}], [num+1-itr-0.25 num+1-itr+0.25],'r');
    hold on
    plot([RasterB_6s{itr+b,:};RasterB_6s{itr+b,:}], [num*2+1-itr-0.25 num*2+1-itr+0.25],'b');
    hold on
end
xlim([-1 7])
ylim([0 num*2])
plot([0 0],[0 40],'k--')
hold on
plot([1 1],[0 40],'k--')
hold on

set(gca,'XTick',-1:1:10,'XTickLabel',{' ','0','1','','','','','','7','8','',''})
ylabel('Trial ID','FontSize',10);
box off

subplot(2,1,1)
time=size(FA_6s,2);
Highervalue = (smooth(mean(FA_6s,1)))' + std(FA_6s,0,1)/sqrt(size(FA_6s,1));
Lowervalue = (smooth(mean(FA_6s,1)))' - std(FA_6s,0,1)/sqrt(size(FA_6s,1));
Time = [1:time, fliplr(1:time)];
value = [Highervalue, fliplr(Lowervalue)];
a = fill(Time, value, [1 0 0], 'edgecolor','none');
alpha(a,0.2);
hold on
plot(1:time,smooth(mean(FA_6s,1))','r-')
hold on

Highervalue = (smooth(mean(FB_6s,1)))' + std(FB_6s,0,1)/sqrt(size(FB_6s,1));
Lowervalue = (smooth(mean(FB_6s,1)))' - std(FB_6s,0,1)/sqrt(size(FB_6s,1));
Time = [1:time, fliplr(1:time)];
value = [Highervalue, fliplr(Lowervalue)];
a = fill(Time, value, [0 0 1], 'edgecolor','none');
alpha(a,0.2);
hold on
plot(1:time,smooth(mean(FB_6s,1))','b-')
hold on
plot([13 13],[0 20],'k--')
hold on
plot([16 16],[0 20],'k--')
hold on
xlim([9 40])
ylim([0 25])
set(gca,'XTick',9:4:40,'XTickLabel',{' ','0',' ','','','','5','','7','8','',''})
xlabel('Time(s)','FontSize',10);
ylabel('Firing Rate(Hz)','FontSize',10);
box off

%%
name=strcat('F:\Npxl_paper\WT-showcase\',num2str(SU_id),'_',reg_list{sust(SU_id),1},'_',file{end,1},'.pdf');
exportgraphics(fh,'219.pdf','ContentType','vector')

%% Function
function [folderType,file,spkFolder,metaFolder,error_list]=jointFolder(folder,error_list,homedir)
metaFolder=replace(folder,'\','/');
metaFolder=fullfile(fullfile(homedir,'DataSum'),metaFolder);
if isfolder(metaFolder)
    spkFolder=replace(metaFolder,'imec1','imec0');
    file=dir(fullfile(spkFolder,'spike_info.mat'));
    if isempty(file)
        folderType=-1;
        file=[];
        spkFolder=[];
        disp('Error processing file 2-tracks');
        disp(metaFolder);
        error_list(end+1,:)={folderType,metaFolder};
        %             pause;
        return
    end
    folderType=2;
else
    metaFolder=replace(metaFolder,'DataSum','DataSum/singleProbe');
    spkFolder=metaFolder;
    file=dir(fullfile(spkFolder,'spike_times.npy'));
    if isempty(file)
        folderType=-1;
        file=[];
        spkFolder=[];
        disp('Error processing file 1-track');
        disp(metaFolder);
        error_list(end+1,:)={folderType,metaFolder};
        return
    end
    folderType=1;
end
end

function [avail,out]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,model)
sps=30000;

trials=clearBadPerf(h5read(fullfile(metaFolder,'events.hdf5'),'/trials')',model);
if isempty(trials)
    avail=false;
    out=[];
    return
end

%     trials=double(trials);
%     info=[trials(:,1)/s1s,trials(:,2)/s1s,trials(:,5),trials(:,6),trials(:,7),trials(:,8)];
if strcmp(model, 'full')
    cluster_ids=[sustIds;transIds;nonselIds];
elseif startsWith(model,'selec')
    cluster_ids=[sustIds;transIds];
elseif startsWith(model,'nonsel')
    cluster_ids=nonselIds;
else
    cluster_ids=sustIds;
end

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
cfg.trl=[trials(:,1)-3*sps,trials(:,1)+14*sps,zeros(size(trials,1),1)-3*sps,trials];
cfg.trlunit='timestamps';
cfg.timestampspersecond=sps;

FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);

out=FT_SPIKE;
avail=true;
end


function out=clearBadPerf(facSeq, mode)
if strcmp(mode, 'error')
    if length(facSeq)>=40
        errorsel=~xor(facSeq(:,5)==facSeq(:,6) , facSeq(:,7)>0);
        out=facSeq(errorsel,:);
    else
        out=[];
    end
else
    if length(facSeq)>=40
        facSeq(:,9)=0;
        i=40;
        while i<=length(facSeq)
            good=xor(facSeq(i-39:i,5)==facSeq(i-39:i,6) , facSeq(i-39:i,7)>0);
            facSeq(i-39:i,10)=good;
            if nnz(good)>=30 %.75 correct rate
                facSeq(i-39:i,9)=1;
            end
            i=i+1;
        end
        out=facSeq(all(facSeq(:,9:10),2),:);
    else
        out=[];
    end
end
end


