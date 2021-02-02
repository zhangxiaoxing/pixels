disk='K:\neupix\DataSum';
s1s=30000;
FR_Th=1.0;
fl=dir([disk,'\**\cluster_info.tsv']);
missing_disk=cell(0);
all_wfstats=cell(0);
for onefile=fl'
    fclose('all');
    rootpath=onefile.folder;
    disp(rootpath);
    if ~isfile(fullfile(rootpath,'wf_stats.hdf5'))
        metaf=ls(fullfile(rootpath,'*.ap.meta'));
        fh=fopen(fullfile(rootpath,metaf));
        ts=textscan(fh,'%s','Delimiter',{'\n'});
        nSample=str2double(replace(ts{1}{startsWith(ts{1},'fileSizeBytes')},'fileSizeBytes=',''));
        spkNThresh=nSample/385/s1s/2*FR_Th;
        clusterInfo = readtable(fullfile(rootpath,'cluster_info.tsv'),'FileType','text','Delimiter','tab');
        waveformGood=strcmp(clusterInfo{:,4},'good');
        freqGood=clusterInfo{:,10}>spkNThresh;
        freq=clusterInfo{:,10}/spkNThresh; %thresh happens to be 1.0 Hz
        cluster_ids = table2array(clusterInfo(waveformGood & freqGood,1));
        if isempty(cluster_ids)
            continue
        end
        wfpath=replace(rootpath,disk,'K:\neupix\WF\neuropixel');
        if isfile(fullfile(wfpath,'waveform.mat'))
%             continue
            wf_fstr=load(fullfile(wfpath,'waveform.mat'));
            wf_all=wf_fstr.waveform;
            if numel(wf_all)<4
                continue
            end

            folder_wf_stats=[];

            for cid=cluster_ids'
                wfidx=find([wf_all{:,2}]==cid);
                if isempty(wfidx)
                    continue
                end
                wfStat=process_wf(wf_all{wfidx,4});
                folder_wf_stats(end+1,:)=[cid,freq(clusterInfo{:,1}==cid),wfStat];
                all_wfstats(end+1,:)={rootpath,folder_wf_stats(end,:)};
            end
            save(fullfile(rootpath,'wf_stats.mat'),'folder_wf_stats');
            syncH5=fullfile(rootpath,'wf_stats.hdf5');
            h5create(syncH5,'/wf',size(folder_wf_stats),'Datatype','double')
            h5write(syncH5,'/wf',folder_wf_stats)
        else
            missing_disk{end+1}=rootpath;
        end
    end
end

keyboard()

all_wf_mat=reshape([all_wfstats{:,2}],5,[])';
% save('assign_type.mat','all_wfstats','fl','missing_disk','all_wf_mat')
fh=figure();
sel=all_wf_mat(:,3)==0;
scatter(all_wf_mat(sel,2),all_wf_mat(sel,4),5,'k','o','MarkerFaceColor','k','MarkerFaceAlpha',0.2,'MarkerEdgeColor','none')
xlim([0,60])
ylim([100,1100])
xlabel('firing rate (Hz)')
ylabel('trough-peak distance')

fh=figure();
histogram(all_wf_mat(sel,2),0:1:30)

fh=figure();
histogram(all_wf_mat(sel,4),0:20:1100)


function out=process_wf(wf)
if max(wf)>-min(wf)
    out=[-1,0,0];
    return
end
%criteria 2
[lc_pk,~]=findpeaks(wf,'MinPeakProminence',-0.05*min(wf));
if numel(lc_pk)>=6
    findpeaks(wf,'MinPeakProminence',-0.05*min(wf));
    out=[-2,0,0];
    return
end
%criteria 3
[~,t_ts]=min(wf);
[~,p_ts]=max(wf(t_ts+1:end));
if p_ts>3
    [lc_pk,~]=findpeaks(wf(t_ts:(t_ts+p_ts-1)),'MinPeakProminence',-0.05*min(wf));
end
if isempty(p_ts) || (p_ts<=3) || (~isempty(lc_pk))
    out=[-3,0,0];
    return
end
%trough_peak dist
wf=spline(1:91,wf,1:0.03:91);
scale=max(abs(wf));
wf=wf./scale;
[~,troughTS]=min(wf);
[~,deltaTS]=max(wf((troughTS+1):end));%to late peak

%fwhm
lcross=find(wf<-0.5,1);
rcross=find(wf(troughTS:end)>-0.5,1)+troughTS;
if numel(lcross)==1 && numel(rcross)==1
    fwhm=rcross-lcross;
    out=[0,deltaTS,fwhm];
    return
else
    out=[-4,0,0];
    return
end

end