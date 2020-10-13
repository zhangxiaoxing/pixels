keyboard
if to_load
    addpath(fullfile('k:','code','jpsth','npy-matlab-master','npy-matlab'))
    addpath('k:\code\jpsth\fieldtrip-20200320')
    ft_defaults
    
    load k:\code\jpsth\114_sorted_file_path.mat
    if isunix
        cd('~/pixels/jpsth')
        homedir='/home/zx/neupix/wyt';
    elseif ispc
        homedir='k:\neupix\wyt';
    end
    sus_trans=h5read('../transient_6.hdf5','/sus_trans'); %export_arr = np.vstack((sust, transient, switched, unclassified, early_in_6s, late_in_6s, prefer_s))
    reg_list=h5read('../transient_6.hdf5','/reg');
    cid_list=h5read('../transient_6.hdf5','/cluster_id');
    path_list=h5read('../transient_6.hdf5','/path');
end


sessIdx=23;
folder=regexp(sorted_fpath{sessIdx},'(\w|\\|-)*','match','once');
stub=regexp(folder,'.*(?=\\)','match','once');
sesssel=startsWith(path_list,stub);
p1sel=contains(path_list,'imec0');
p2sel=contains(path_list,'imec1');
memorysel=any(sus_trans(:,1:2)>0,2);
p1id=cid_list(sesssel & p1sel & memorysel);
p2id=cid_list(sesssel & p2sel & memorysel)+10000;
[folderType,file,spkFolder,metaFolder,~]=jointFolder(folder,cell(0),homedir);
currmodel='selec';
sustIds=[];nonselIds=[];
transIds=int32([p1id;p2id]);

%% %%%%%% localize %%%%%%%%%%%%%
% freg=load('0831_selec_conn_chain_duo_6s_1_2.mat','pair_reg','pair_chain');
% load reg_keep.mat
% suid=ringsm(ssidx,:);
% reg=[];
% for i1=suid(:)'
%     reg=[reg;{i1,reg_set{freg.pair_reg(find(freg.pair_chain==i1,1))}}];
% end
%% end of localize

% msize=numel(transIds);
[avail,spktrial]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,currmodel); % posix
delay=6;
TT1 = find(spktrial.trialinfo(:,5)==4 & spktrial.trialinfo(:,8)==delay);
TT2 = find(spktrial.trialinfo(:,5)==8 & spktrial.trialinfo(:,8)==delay);

TT = [TT1(1:5);TT2(1:5)];
%% reduced number of trials for debug purpose only
%     TT = TT(1);
%%
% TT=32
NN=size(spktrial.label,1);
bin=[1,7];
spiketrains=cell(1,NN);
spiketrials=cell(1,NN);
for i=1:NN
    sptr=struct();
    sptr.times=[];
    sptr.trials=[];
    for t=1:length(TT)
    ttimes=spktrial.time{i}...
        (spktrial.trial{i}==TT(t) ...
        & spktrial.time{i}<bin(2) ...
        & spktrial.time{i}>=bin(1)).*1000+(t-1)*7000;
    ttrials=ones(size(ttimes))*TT(t);
    sptr.times=[sptr.times,ttimes];
    sptr.trials=[sptr.trials,ttrials];
    end
    sptr.times_units='ms';
    sptr.t_start=1000;
    sptr.t_start_units='ms';
    sptr.t_stop=7*length(TT)*1000;
    sptr.t_stop_units='ms';
    spiketrains{i}=sptr;
end
seg.spiketrains=spiketrains;
block=struct();
block.name=num2str(sessIdx);
block.segments={seg};
save(sprintf('spkt_%d_%d.mat',sessIdx,TT(1)),'block','-v7');



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

function out=tag_ring(in,msize)
out=in;
out(:,3)=0;
if msize==3
    for i=1:length(in)-3
        if (isequal(in(i:i+3,2),[1 2 3 1]') || isequal(in(i:i+3,2),[2 3 1 2]') || isequal(in(i:i+3,2),[3 1 2 3]')) ...
                && max(diff(in(i:i+3,1)))<0.01 && min(diff(in(i:i+3,1)))>0.001
            out(i:i+3,3)=1;
        end
    end
elseif msize==4
    for i=1:length(in)-4
        if (isequal(in(i:i+4,2),[1 2 3 4 1]')...
                || isequal(in(i:i+4,2),[2 3 4 1 2]')...
                || isequal(in(i:i+4,2),[3 4 1 2 3]')...
                || isequal(in(i:i+4,2),[4 1 2 3 4]'))...
                && max(diff(in(i:i+4,1)))<0.01 ...
                && min(diff(in(i:i+4,1)))>0.001
            out(i:i+3,3)=1;
        end
    end
end
end

function out=relax_tag(in,msize)
out=in;
out(:,3)=0;
skiptag=0;
for i=1:length(in)
    if i<skiptag
        continue
    end
    curr_su=in(i,2);
    targets=[(curr_su+1):(curr_su+msize-1),curr_su];
    targets(targets>msize)=targets(targets>msize)-msize;
    tsseq=[i,in(i,1:2)];
    for t=targets
        rows=tsseq(end,1)+(1:msize*10);
        rows(rows>length(in))=[];
        if isempty(rows)
            break
        end
        didx=find( ...
            in(rows,2)==t ... %post unit
            & in(rows,1)<tsseq(end,2)+0.01 ...
            & in(rows,1)>tsseq(end,2)+0.0005 ... %matching time window, assuming 1kHz
            ,1);
        if isempty(didx)
            break
        else
            tsseq=[tsseq;tsseq(end,1)+didx,in(tsseq(end,1)+didx,1:2)];
        end
    end
    if length(tsseq)<msize+1
        continue
    else
        out(tsseq(:,1),3)=1;
        skiptag=tsseq(2,1);
    end
end
end

function plotseq(in,msize,tidx,ids,sessIdx,ssidx,cocount)
fh=figure('Color','w','Position',[100,100,400,400]);
hold on
color={'r','b','c','m','g'};
ph=matlab.graphics.chart.primitive.Line.empty(0,5);
for bin=1:6
    for i=1:msize
        j=i;
        selp=in(:,2)==i & in(:,3)==1 & in(:,1)>=bin & in(:,1)<bin+1;
        seln=in(:,2)==i & in(:,3)==0 & in(:,1)>=bin & in(:,1)<bin+1;
        tph=plot(repmat(in(selp,1)',2,1)-bin,(bin-1)*(msize+1)+repmat([j-0.6;j+0.6],1,nnz(selp)),'-','Color',color{i});
        %         ph(i)=tph(1);
        plot(repmat(in(seln,1)',2,1)-bin,(bin-1)*(msize+1)+repmat([j-0.6;j+0.6],1,nnz(seln)),'-','Color',[0.8,0.8,0.8]);
        
    end
end
xlabel('Time within time-bin(s)')
ylabel('Time-bin')
set(gca,'XTick',0:0.5:1,'YTick',(msize/2):(msize+1):((msize+1)*6),'YTickLabel',1:6)
if msize==3
    title(sprintf('Sess#%d, Trial#%d, SU%d, %d, %d',sessIdx,tidx,ids(1),ids(2),ids(3)));
    %     print(fh,sprintf('spike_seq_%d_%d_%d.png',ids(1),ids(2),ids(3)),'-dpng');
elseif msize==4
    title(sprintf('Sess#%d, Trial#%d, SU%d, %d, %d, %d',sessIdx,tidx,ids(1),ids(2),ids(3),ids(4)));
    %     print(fh,sprintf('spike_seq_%d_%d_%d_%d.png',ids(1),ids(2),ids(3),ids(4)),'-dpng');
elseif msize==4
    title(sprintf('Sess#%d, Trial#%d',sessIdx,tidx));
end
print(fh,sprintf('r%d_spike_seq_%d_%d.png',msize,cocount,ssidx),'-dpng');
if false
    %     legend(ph,arrayfun(@(x) ['cell ',num2str(x)],1:5,'UniformOutput',false),'NumColumns',3,'Location','northoutside')
    ylim([0,6*(msize+1)]);
    fh.Position(3:4)=[225,225];
    set(gca,'FontSize',10);
    exportgraphics(fh,sprintf('r%d_spike_seq_%d_%d.pdf',msize,cocount,ssidx),'ContentType','vector');
end

close(fh);
end

