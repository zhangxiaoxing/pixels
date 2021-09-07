% keyboard()
%scatter(max(scores_S1(:,[1:4,7:10]),[],2),totalcount_S1',2,'Marker','o','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none','MarkerFaceColor','k')
load 0115_path
if isunix
    addpath(fullfile('npy-matlab-master','npy-matlab'))
    addpath(fullfile('fieldtrip-20200320'))
    homedir=fullfile('/','home','zx','neupix','wyt','DataSum');
else
    addpath(fullfile('K:','Lib','npy-matlab-master','npy-matlab'))
    addpath(fullfile('K:','Lib','fieldtrip-20200320'))
    homedir=fullfile('K:','neupix','wyt','DataSum');
end
ft_defaults
sums=cell(0);
for bin=binin
if denovo
    load(sprintf('0115_nonsel_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
    scoresel=max(scores_S1(:,[1:4,7:10]),[],2)>5;
    countsel=totalcount_S1'>=2000;
%     prefsel=pref_chain_S1(:,6)==0 & pref_chain_S1(:,12)==0;
    regsel=reg_chain_S1(:,1)~=reg_chain_S1(:,2);
    sel_chain=conn_chain_S1(scoresel & countsel & regsel,:);
else
    load FC_SC_Sample.mat
    sel_chain=selected;    
end

for idx=9%1:size(sel_chain,1)
% for idx=53
sessidx=idivide(sel_chain(idx,1),int32(100000));
path=sorted_path_0115{sessidx};
if contains(path,'imec1') && exist(replace(path,'imec1','imec0'),'dir')
    path=replace(path,'imec1','imec0');
end
if isunix
    path=replace(path,'\',filesep);
end
if ~exist(fullfile(homedir,path,'spike_info.mat'),'file') || ~exist(fullfile(homedir,path,'events.hdf5'),'file')
    continue
end


trials=h5read(fullfile(homedir,path,'events.hdf5'),'/trials')';
ctrials=clearBadPerf(trials);
fstr=load(fullfile(homedir,path,'spike_info.mat'));
spkId=double([fstr.spike_info{1}{1};fstr.spike_info{1}{2}]);
spkTS=double([fstr.spike_info{2}{1};fstr.spike_info{2}{2}]);

FT_SPIKE=struct();
cluster_ids=rem(sel_chain(idx,:),100000);
FT_SPIKE.label=strtrim(cellstr(num2str(cluster_ids')));
FT_SPIKE.timestamp=cell(1,2);
for i=1:numel(cluster_ids)
    FT_SPIKE.timestamp{i}=spkTS(spkId==cluster_ids(i))';
end
%  continuous format F T struct file
sps=30000;
cfg=struct();
cfg.trl=[trials(:,1)-3*sps,trials(:,1)+11*sps,zeros(size(trials,1),1)-3*sps,trials];
cfg.trlunit='timestamps';
cfg.timestampspersecond=sps;
FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);


%dimord={spkN1,spkN2,fc_event,reverse fc_event,+50ms fc_event} x trl(240) x bin(-3:11)

stats=nan(5,size(trials,1),14);
for trlIdx=1:size(trials,1)
    ts1=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==trlIdx)';
    ts2=FT_SPIKE.time{2}(FT_SPIKE.trial{2}==trlIdx)';

    ts1(:,2)=1;
    ts2(:,2)=2;
    ts_id=[ts1;ts2];

    [~,s]=sort(ts_id(:,1));
    ts_id=ts_id(s,:);
    ts_id_tagged=fc_tag(ts_id);
    stats(1,trlIdx,:)=histcounts(ts1(:,1),-3:11);
    stats(2,trlIdx,:)=histcounts(ts2(:,1),-3:11);
    stats(3,trlIdx,:)=histcounts(ts_id_tagged(ts_id_tagged(:,3)==1 & ts_id_tagged(:,2)==1,1),-3:11);
    stats(4,trlIdx,:)=histcounts(ts_id_tagged(ts_id_tagged(:,4)==1 & ts_id_tagged(:,2)==1,1),-3:11);
    stats(5,trlIdx,:)=histcounts(ts_id_tagged(ts_id_tagged(:,5)==1 & ts_id_tagged(:,2)==1,1),-3:11);
end

% if nnz(isnan(stats(1:3,:,:)))>0
%     keyboard()
% end

sel_6s_S1=trials(:,5)==4 & trials(:,8)== 6;
sel_6s_S2=trials(:,5)==8 & trials(:,8)== 6;

sel_3s_S1=trials(:,5)==4 & trials(:,8)== 3;
sel_3s_S2=trials(:,5)==8 & trials(:,8)== 3;
p3=ones(1,3);
if nnz(sel_3s_S1)>5 && nnz(sel_3s_S2)>5
    for sbin=1:3
        p3(sbin)=ranksum(squeeze(stats(3,sel_3s_S1,sbin+4)),squeeze(stats(3,sel_3s_S2,sbin+4)));
    end
end
if nnz(sel_6s_S1)>5 && nnz(sel_6s_S2)>5
p6=ones(1,6);
    for sbin=1:6
        p6(sbin)=ranksum(squeeze(stats(3,sel_6s_S1,sbin+4)),squeeze(stats(3,sel_6s_S2,sbin+4)));

    end
end
% 
% 
% if any([p3,p6]<0.01)
%     disp([bin,idx,p3,p6])
% %     keyboard()
% end
sums(end+1,:)={bin,sel_chain(idx,:),trials,stats,[p3,p6]};
% find(any(sums(:,3:5)<0.01,2) & any(sums(:,6:11)<0.01,2))

end
save(sprintf('nonsel_silent_fc_%d.mat',binin),'sums')
end
return

function out=clearBadPerf(facSeq, mode)
if exist('mode','var') && strcmp(mode, 'error')
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




function out=fc_tag(in)
out=in;
out(:,3:5)=0;

for i=1:size(in,1)-1
    if in(i,2)==1
        j=i+1;
        while in(j,2)==2 && j<size(in,1)
            if in(j,1)<=in(i,1)+0.01 && in(j,1)>in(i,1)+0.002
                out(i,3)=1;
                out(j,3)=1;
            end
            j=j+1;
        end
    end
end
for i=1:size(in,1)-1
    if in(i,2)==2
        j=i+1;
        while in(j,2)==1 && j<size(in,1)
            if in(j,1)<=in(i,1)+0.01 && in(j,1)>in(i,1)+0.002
                out(i,4)=1;
                out(j,4)=1;
            end
            j=j+1;
        end
    end
end
for i=1:size(in,1)-1
    if in(i,2)==1
        j=i+1;
        while in(j,2)==2 && j<size(in,1)
            if in(j,1)<=in(i,1)+0.06 && in(j,1)>in(i,1)+0.052
                out(i,5)=1;
                out(j,5)=1;
            end
            j=j+1;
        end
    end
end
end
