
addpath(fullfile('k:','Lib','npy-matlab-master','npy-matlab'))
addpath('K:\Lib\fieldtrip-20200320')
ft_defaults
load 114_sorted_file_path.mat
delay=6;
if isunix
    cd('~/pixels/jpsth')
    homedir='/home/zx/neupix/wyt';
elseif ispc
    homedir='k:\neupix\wyt';
end


all_tagged=[];
patt_suids=[57,82,85,10070,10097,10099,10101,913,339,414,445,451,452,475,485];
fc_list=get_func_conn(patt_suids);
sessIdx=4;

trials=[88916040,89126095,2963,2970,4,8,1,6,1,1;163631315,163841370,5454,5461,8,8,-1,6,0,1];

folder=regexp(sorted_fpath{sessIdx},'(\w|\\|-)*','match','once');
[folderType,file,spkFolder,metaFolder,~]=jointFolder(folder,cell(0),homedir);

fstr=load(fullfile(spkFolder,'spike_info.mat'));
spkId=double([fstr.spike_info{1}{1};fstr.spike_info{1}{2}]);
spkTS=double([fstr.spike_info{2}{1};fstr.spike_info{2}{2}]);

for ssidx=1:size(fc_list,1)
    sample=fc_list(ssidx,3);    
    transIds=int32(rem(fc_list(ssidx,1:2),100000))';
    %%
    sps=30000;
    cluster_ids=transIds;

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
    spktrial=FT_SPIKE;
    %%

    examples=cell(0);
    tidx=sample;
    ts_id=[];
    for seqid=1:2
        tsel=spktrial.trial{seqid}==tidx;
        ts=spktrial.time{seqid}(tsel);
        ts=ts(ts>1 & ts<7);
        ts(2,:)=seqid;
        ts_id=[ts_id;ts'];
    end
    [~,s]=sort(ts_id(:,1));
    ts_id=ts_id(s,:);
    ts_id_tagged=relax_tag(ts_id,msize);
    ts_id_tagged(:,4:7)=repmat([transIds',sample,tidx],size(ts_id_tagged,1),1);
    all_tagged=[all_tagged;ts_id_tagged];

end

%%TS,FC_SU_ID,FC_TAG,SU1_CID,SU2_CID,FC_SAMPLE,TRIAL_ID(1->22,2->112)

tag_ts=[];
for samp=1:2
    for pid=patt_suids
        %pre
        ts_pre=unique(all_tagged(all_tagged(:,4)==pid ...
            & all_tagged(:,2)==1 ...
            & all_tagged(:,3)==1 ...
            & all_tagged(:,7)==samp,1));

        %post
        ts_post=unique(all_tagged(all_tagged(:,5)==pid ...
            & all_tagged(:,2)==2 ...
            & all_tagged(:,3)==1 ...
            & all_tagged(:,7)==samp,1));

        ts=[ts_pre;ts_post];
        ts(:,2)=pid;
        ts(:,3)=samp;    
        ts=ts(:,[2,3,1]);
        tag_ts=[tag_ts;ts];
    end
end
tag_ts=unique(tag_ts,'rows');
save('func_conn_showcase_tagged_spikes_sess4_trial_22_122.mat','all_tagged','tag_ts')

return

function out=get_func_conn(suids)
patt_su=suids+400000;

out=[];
fstr=cell(1,6);
for bin=1:6
    fstr{bin}=load(sprintf('0831_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
end
I=4;
lbound=100000*I;
ubound=100000*(I+1);
for bin=1:6
    sel11=fstr{bin}.conn_chain_S1(:,1)>=lbound & fstr{bin}.conn_chain_S1(:,1)<ubound & diff(fstr{bin}.reg_chain_S1,1,2);
    sel12=fstr{bin}.conn_chain_S2(:,1)>=lbound & fstr{bin}.conn_chain_S2(:,1)<ubound & diff(fstr{bin}.reg_chain_S2,1,2);
    out1=fstr{bin}.conn_chain_S1(sel11,:);
    out2=fstr{bin}.conn_chain_S2(sel12,:);
    out1(:,3)=1;
    out2(:,3)=2;
    out=[out;out1;out2];
end
out=unique(out,'rows');
patt_sel=ismember(out(:,1),patt_su) & ismember(out(:,2),patt_su);
out=out(patt_sel,:);

end

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

function out=relax_tag(in,msize)
out=in;
out(:,3)=0;

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
   
%             & in(rows,1)<tsseq(end,2)+0.01 ...
%             & in(rows,1)>tsseq(end,2)+0.0005 ... %matching time window, assuming 1kHz
end
end
