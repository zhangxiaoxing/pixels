function errors=extract_waveform(disk)
% to_plot=false;
s1s=30000;
FR_Th=1.0;
disk='k';
addpath('K:\code\neuropixel-utils')
% addpath(genpath('D:\code\Kilosort2'))
channelMapFile='K:\code\neuropixel-utils\map_files\neuropixPhase3B2_kilosortChanMap.mat';

fl=dir([disk,':\neuropixel\**\cluster_info.tsv']);
errors=cell(0);
for onefile=fl'
    rootpath=onefile.folder;
    try
        if isfile(fullfile(replace(rootpath,[disk,':'],'\neupix\WF\'),'waveform.mat'))
            fstr=load(fullfile(replace(rootpath,[disk,':'],'\neupix\WF\'),'waveform.mat'));
            if size(fstr.waveform,2)==4
                fprintf('Exist file %s, skipped\n\n',fl.folder)
                continue
            end
        end
        disp(rootpath)

        tic
        
        if isfolder(replace(rootpath,[disk,':'],'\neupix\WF\'))
            disp('Folder exists, please double check! Press Ctrl-C to break')
        else
%             mkdir(replace(rootpath,[disk,':'],'\neupix\WF\'));
              mkdir(replace(rootpath,[disk,':'],'K:\neupix\WF\'));
        end

        metaf=ls(fullfile(rootpath,'*.ap.meta'));
        fh=fopen(fullfile(rootpath,metaf));
        ts=textscan(fh,'%s','Delimiter',{'\n'});
        nSample=str2double(replace(ts{1}{startsWith(ts{1},'fileSizeBytes')},'fileSizeBytes=',''));
        spkNThresh=nSample/385/s1s/2*FR_Th;
        clusterInfo = readtable(fullfile(rootpath,'cluster_info.tsv'),'FileType','text','Delimiter','tab');
        waveformGood=strcmp(clusterInfo{:,4},'good');
        freqGood=clusterInfo{:,10}>spkNThresh;
        cluster_ids = table2array(clusterInfo(waveformGood & freqGood,1));

        if numel(cluster_ids)>0
            waveform=cell(0,4);
            ks=Neuropixel.KiloSortDataset(rootpath,'channelMap',channelMapFile);
            ks.load();

            for cidx=1:numel(cluster_ids)
                try
                    snippetSetTet = ks.getWaveformsFromRawData('cluster_ids',cluster_ids(cidx),'num_waveforms', 100, 'best_n_channels', 4, 'car', true, ...
                        'subtractOtherClusters', false,'window', [-30 60]);

                    snippetSetBest = ks.getWaveformsFromRawData('cluster_ids',cluster_ids(cidx),'num_waveforms', 500, 'best_n_channels', 1, 'car', true, ...
                        'subtractOtherClusters', false,'window', [-30 60]);

                    waveform{end+1,1}=rootpath;
                    waveform{end,2}=cluster_ids(cidx);
                    waveform{end,3}=snippetSetTet.data;
                    waveform{end,4}=mean(snippetSetBest.data,3);
                catch ME
                    errors{end+1}=sprintf('%s %d waveform error',onefile.folder,cidx);
                end
%                 if to_plot
%                     fh=figure();
%                     plot(cell2mat(arrayfun(@(x) mean(squeeze(snippetSet.data(x,:,:))'), 1:4,'UniformOutput',false)));
%                     pause
%                     close(fh);
%                 end
            end
        end
        save(fullfile(replace(rootpath,[disk,':'],'K:\neupix\WF\'),'waveform.mat'),'waveform');
        toc
    catch ME
        errors{end+1}=onefile.folder;
    end
end
end
