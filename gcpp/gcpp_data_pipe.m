%% this is basically a copy from x_corr_delay.m to reuse parts of the data pipeline
% key change marked with "gcpp" tag.
%%


% cd('~/pixels/jpsth')
% homedir='/home/zx/neupix/wyt';

cd('K:\code\jpsth')
homedir='K:\neupix\wyt\DataSum';

currmodel='selec';
prefix='0807';
delay=6;
% bin_range=[4,5];
% addpath(fullfile('npy-matlab-master','npy-matlab'))
addpath('fieldtrip-20200320')
ft_defaults
sus_trans=h5read('../transient_6.hdf5','/sus_trans'); %export_arr = np.vstack((sust, transient, switched, unclassified, early_in_6s, late_in_6s, prefer_s))
reg_list=h5read('../transient_6.hdf5','/reg');
cid_list=h5read('../transient_6.hdf5','/cluster_id');
path_list=h5read('../transient_6.hdf5','/path');

if startsWith(currmodel,'selec')
    sust=find(sus_trans(:,1));
    trans=find(sus_trans(:,2));
    supool=[sust;trans]';
elseif startsWith(currmodel,'nonsel')
    sust=[];
    trans=[];
    nonsel_logic=~(sus_trans(:,1) | sus_trans(:,2) | sus_trans(:,3)|sus_trans(:,4));
    supool=find(nonsel_logic);
end
counter=[];
done=[];


error_list=cell(0);

for i=1:length(supool)
    if ismember(supool(i),done)
        continue
    end

    folder=regexp(path_list{supool(i)},'(\w|\\|-)*','match','once');
    [folderType,file,spkFolder,metaFolder,error_list]=jointFolder(folder,error_list);
    if folderType<0
        continue
    end
    wffile=fullfile(metaFolder,'wf_stats.hdf5');% posix
    if isfile(wffile)
        if startsWith(currmodel,'selec')
            if folderType==1
                sustIds=cid_list(startsWith(path_list,folder) & sus_trans(:,1));
                transIds=cid_list(startsWith(path_list,folder) & sus_trans(:,2));
                sameFolder=find(startsWith(path_list,folder) & (sus_trans(:,1)| sus_trans(:,2)));
                done=[done;sameFolder];
                sustCount=numel(sustIds);
                transCount=numel(transIds);
            elseif folderType==2
                fimec0=replace(folder,'imec1','imec0');
                fimec1=replace(folder,'imec0','imec1');

                sustIds0=cid_list(startsWith(path_list,fimec0) & sus_trans(:,1));
                transIds0=cid_list(startsWith(path_list,fimec0) & sus_trans(:,2));

                sustIds1=cid_list(startsWith(path_list,fimec1) & sus_trans(:,1))+10000;
                transIds1=cid_list(startsWith(path_list,fimec1) & sus_trans(:,2))+10000;

                sameFolder=find((startsWith(path_list,fimec0)|startsWith(path_list,fimec1)) & (sus_trans(:,1)| sus_trans(:,2)));
                done=[done;sameFolder];

                sustIds=[sustIds0(:);sustIds1(:)];
                transIds=[transIds0(:);transIds1(:)];

                sustCount=numel(sustIds);
                transCount=numel(transIds);            

            end

            if transCount<1
                continue
            end
            nonselIds=[];
            [avail,spktrial]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,currmodel); % posix
        elseif startsWith(currmodel,'nonsel')
            if folderType==1
                sameFolder=find(startsWith(path_list,folder) & nonsel_logic);
                nonselIds=cid_list(startsWith(path_list,folder) & nonsel_logic);
                done=[done;sameFolder];
            elseif folderType==2
                fimec0=replace(folder,'imec1','imec0');
                fimec1=replace(folder,'imec0','imec1');

                nonselIds0=cid_list(startsWith(path_list,fimec0) & nonsel_logic);

                nonselIds1=cid_list(startsWith(path_list,fimec1) & nonsel_logic)+10000;


                sameFolder=find((startsWith(path_list,fimec0)|startsWith(path_list,fimec1)) & nonsel_logic);
                done=[done;sameFolder];

                nonselIds=[nonselIds0(:);nonselIds1(:)];
                sustIds=[];
                transIds=[];

            end

            if numel(nonselIds)<1
                continue
            end
            [avail,spktrial]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,currmodel); % posix
        
        end
        if avail
            %% create mat X for futher processing
            
            %todo for s in samples, for bin in bins, (trial*1000) x SU
            trialsel=find(spktrial.trialinfo(:,5)==4 & spktrial.trialinfo(:,8)==6)';
            for bin=1
                Xraw=zeros(1000*numel(trialsel),numel(spktrial.label));
                for suid=1:numel(spktrial.label)
                    for rearrIdx=1:numel(trialsel)
                        tIdx=trialsel(rearrIdx);
                        tt=spktrial.time{suid}(spktrial.trial{suid}==tIdx);
                        xIdx=int32((tt(tt>bin & tt<bin+1)-bin+rearrIdx-1)*1000+1);
%                         disp(xIdx);
                        Xraw(xIdx,suid)=1;
                    end
                end
                
                keyboard
            end
        end
    else
        continue
    end
    
    %% gcpp
%     sums={i,folder,sustIds,transIds,xc_s1,xcshuf_s1,xc_s2,xcshuf_x2}; %per folder save
%     save(sprintf('%s_%s_XCORR_duo_f%d_delay_%d_%d_%d_2msbin.mat',prefix,currmodel,i,delay,bin_range(1),bin_range(2)),'sums','-v7.3','sust','trans','supool','counter','done') %prefix
%  	fprintf('%d of %d\n',i,length(supool))
end

return 
function [folderType,file,spkFolder,metaFolder,error_list]=jointFolder(folder,error_list)
    metaFolder=replace(folder,'\','/');
    homedir=evalin('base','homedir');
    metaFolder=fullfile(homedir,metaFolder);
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
    s1s=30000;
    FR_Th=1.0;

    metaf=strtrim(ls(fullfile(metaFolder,'*.meta')));
    fh=fopen(metaf);
    ts=textscan(fh,'%s','Delimiter',{'\n'});
    nSample=str2double(replace(ts{1}{startsWith(ts{1},'fileSizeBytes')},'fileSizeBytes=',''));
    spkNThresh=nSample/385/s1s/2*FR_Th;
    clusterInfo = readtable(fullfile(metaFolder,'cluster_info.tsv'),'FileType','text','Delimiter','tab');
    waveformGood=strcmp(clusterInfo{:,4},'good');
    freqGood=clusterInfo{:,10}>spkNThresh;
    cluster_ids = table2array(clusterInfo(waveformGood & freqGood,1));
elseif startsWith(model,'selec')
    cluster_ids=[sustIds;transIds];
elseif startsWith(model,'nonsel')
    cluster_ids=nonselIds;
else
    keyboard
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

function postHocStats()
mm=mean(Xraw);
[~,idx]=sort(mm,'descend');
% X=XX(:,idx(1:32));
lblstr=spktrial.label(idx(1:32));
lbl=cellfun(@(x) str2double(x), lblstr);
metastat=zeros(113,1);
for f=1:113
    disp(f);
stats=[];
    for i=1:numel(lbl)
        for j=1:numel(lbl)
            if i==j
                continue
            end
            xcorrType=nnz(conn_chain_S1(:,1)==lbl(j)+f*100000 & conn_chain_S1(:,2)==lbl(i)+f*100000);
            stats=[stats;j,i,lbl(j),lbl(i),Phi(i,j),Psi2(i,j),xcorrType];
        end
    end
    metastat(f)=nnz(stats(:,6)>0 & stats(:,7)>0);
end
end



