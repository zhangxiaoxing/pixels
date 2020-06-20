function pixFlatFR()
addpath('R:\ZX\neupix\npy-matlab\npy-matlab');
binSize=0.25;
binsOnset=-3;
FR_Th=1.0;
s1s=30000;
% 
javaaddpath('R:\ZX\java\spk2fr\build\classes');
% javaaddpath('R:\ZX\java\spk2fr\lib\jmatio.jar');
% javaaddpath('R:\ZX\java\spk2fr\lib\commons-math3-3.5.jar');
% javaaddpath('R:\ZX\java\DualEvtParser\build\classes');


bf=spk2fr.multiplesample.BuildFileList();
fs=string(bf.buildPixels('R:\ZX\neupix\DataSum\').toArray);
clear('bf');

% i=325
% 
% [toPath,~,~]=fileparts(fs(i));
% plotOneDir(char(toPath))
% return

% parpool(12);
futures=parallel.FevalFuture.empty(0,length(fs));
for i=1:length(fs)
    [toPath,~,~]=fileparts(fs(i));
    futures(i)=parfeval(@plotOneDir,0,char(toPath));
end

for i=1:length(fs)
    fetchOutputs(futures(i));
    fprintf('%d of %d\n',i,length(fs));
end


% for i=1:length(fs)
%     [toPath,~,~]=fileparts(fs(i));
%     plotOneDir(char(toPath));
%     
% end


    function plotOneDir(rootpath)
        %rootpath=char(toPath);
        if exist(fullfile(rootpath,'FR_All.hdf5'),'file') || exist(fullfile(rootpath,'NOSU'),'file')
            return
        end
        spkTS=readNPY(fullfile(rootpath,'spike_times.npy'));
        spkId=readNPY(fullfile(rootpath,'spike_clusters.npy'));
        
        metaf=ls(fullfile(rootpath,'*.meta'));
        fh=fopen(fullfile(rootpath,metaf));
        ts=textscan(fh,'%s','Delimiter',{'\n'});
        nSample=str2double(replace(ts{1}{startsWith(ts{1},'fileSizeBytes')},'fileSizeBytes=',''));
        spkNThresh=nSample/385/s1s/2*FR_Th;
        clusterInfo = readtable(fullfile(rootpath,'cluster_info.tsv'),'FileType','text','Delimiter','tab');
        trials=h5read(fullfile(rootpath,'events.hdf5'),'/trials')';
        trials=double(trials);
        spk=[double(ones(size(spkId))),double(spkId),double(spkTS)/s1s];
        clear spkTS
        waveformGood=strcmp(clusterInfo{:,4},'good');
        freqGood=clusterInfo{:,10}>spkNThresh;
        cluster_ids = table2array(clusterInfo(waveformGood & freqGood,1));
        spk=spk(ismember(spk(:,2),cluster_ids),:);
        
        info=[trials(:,1)/s1s,trials(:,2)/s1s,trials(:,5),trials(:,6),trials(:,7),trials(:,8)];
        
        genAll(spk,info,rootpath);
        
%         delayLen=6;
%         filter=(trials(:,8)==delayLen);
%         rLim=delayLen+7;% 9 for 4s delay, 13 for 8s delay
%         genAll(spk,info,filter,isCorrect,rootpath,delayLen,rLim);
    end




    function genAll(spk,info,rootpath)
       
        if(numel(spk)>1000)
            javaaddpath('R:\ZX\java\spk2fr\build\classes');
            javaaddpath('R:\ZX\java\spk2fr\lib\jmatio.jar');
            javaaddpath('R:\ZX\java\spk2fr\lib\commons-math3-3.5.jar');
            javaaddpath('R:\ZX\java\DualEvtParser\build\classes');
            
            s2f=spk2fr.Spk2fr;
            s2f.setRefracRatio(1);
            s2f.setLeastFR('all');
            
 
            cb=s2f.getAllFiringRate(s2f.buildData(info,spk,'pixels'), ...
                'everytrialall', ...
                s2f.setBin(binsOnset,binSize,14.0));
            FR_All=cb.getFRDataA();
            keyIdx=cb.getKeyIdx();
            keyIdx=uint16(keyIdx(:,2));
            
            FR_File=fullfile(rootpath,'FR_All.hdf5');
            if exist(FR_File,'file')
                delete(FR_File)
            end
            h5create(FR_File,'/FR_All',size(FR_All),'Datatype','double')
            h5write(FR_File,'/FR_All',FR_All)
            h5create(FR_File,'/Trials',size(info),'Datatype','double')
            h5write(FR_File,'/Trials',info)
            h5create(FR_File,'/SU_id',size(keyIdx),'Datatype','uint16')
            h5write(FR_File,'/SU_id',keyIdx)

        else
            fid=fopen(fullfile(rootpath,'NOSU'),'w');
            fclose(fid);
        end
    end


end