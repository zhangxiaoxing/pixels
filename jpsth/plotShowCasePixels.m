function plotShowCasePixels()
addpath('R:\ZX\neupix\npy-matlab\npy-matlab');
byLick=false;
binSize=0.25;
FR_Th=1.0;
javaaddpath('R:\ZX\java\spk2fr\build\classes');
javaaddpath('R:\ZX\java\spk2fr\lib\jmatio.jar');
javaaddpath('R:\ZX\java\spk2fr\lib\commons-math3-3.5.jar');
javaaddpath('R:\ZX\java\DualEvtParser\build\classes');
s2f=spk2fr.Spk2fr;
% s2f.setRefracRatio(0.0015);
% s2f.setLeastFR('Average2Hz');
s2f.setRefracRatio(1);
s2f.setLeastFR('all');
s1s=30000;

bf=spk2fr.multiplesample.BuildFileList();
fs=string(bf.buildPixels('R:\ZX\neupix\DataSum\').toArray);
for i=1:length(fs)
    [toPath,~,~]=fileparts(fs(i));
    plotOneDir(char(toPath));
    
end


    function plotOneDir(rootpath)
       
        spkTS=readNPY(fullfile(rootpath,'spike_times.npy'));
        spkId=readNPY(fullfile(rootpath,'spike_clusters.npy'));
        
        spkNThresh=double(spkTS(end))/s1s*FR_Th;
        clusterInfo = readtable(fullfile(rootpath,'cluster_info.tsv'),'FileType','text','Delimiter','tab');
        trials=h5read(fullfile(rootpath,'events.hdf5'),'/trials')';
        trials=double(trials);
        s1s=30000;
        spk=[double(ones(size(spkId))),double(spkId),double(spkTS)/s1s];
        
        waveformGood=strcmp(clusterInfo{:,4},'good');
        freqGood=clusterInfo{:,10}>spkNThresh;
        cluster_ids = table2array(clusterInfo(waveformGood & freqGood,1));
        spk=spk(ismember(spk(:,2),cluster_ids),:);
        
        info=[trials(:,1)/s1s,trials(:,2)/s1s,trials(:,5),trials(:,6),trials(:,7)];
        isCorrect=true;
        
        delayLen=3;
        filter=(trials(:,8)==delayLen);
        rLim=delayLen+7;% 9 for 4s delay, 13 for 8s delay
        genAll(spk,info,filter,isCorrect,rootpath,delayLen,rLim);
        
        delayLen=6;
        filter=(trials(:,8)==delayLen);
        rLim=delayLen+7;% 9 for 4s delay, 13 for 8s delay
        genAll(spk,info,filter,isCorrect,rootpath,delayLen,rLim);
    end




    function genAll(spk,info,filter,isCorrect,rootpath,delayLen,rLim)
        if(numel(spk)>1000)
            ts=s2f.getTS(info(filter,:),spk,'pixels',byLick,isCorrect);
            uKeys=s2f.getKeyIdx();
            for k=1:size(uKeys,1)
                
                try
                    plotOne(ts{k},uKeys(k,2),rootpath,delayLen,rLim);
                catch ME
                    disp(rootpath);
                    disp(k);
                    disp(ME.identifier); 
                end
%                 pause();
                close all;
            end
        end
    end



    function plotOne(ts,uid,rootpath,delayLen,rLim)
        %         pf=nan(0,0);
        %         bn=nan(0,0);
        figure('Color','w','Position',[100,100,280,280]);
        
        subplot('Position',[0.1,0.5,0.85,0.40]);
        hold on;
        
        if size(ts{1},1)==0
        elseif size(ts{1},1)==1
            plot([ts{1};ts{1}],repmat([21;21.8]+1,1,length(ts{1})),'-r');
        elseif size(ts{1},1)>20
            pos=round((length(ts{1})-20)/2);
            posRange=pos:pos+19;
            cellfun(@(x) plot([x{1}';x{1}'],repmat([x{2};x{2}+0.8]+1,1,length(x{1})),'-r'),cellfun(@(x,y) {x,y},ts{1}(posRange),num2cell(1:length(posRange),1)','UniformOutput',false));
        else
            posRange=1:size(ts{1},1);
            cellfun(@(x) plot([x{1}';x{1}'],repmat([x{2};x{2}+0.8]+1,1,length(x{1})),'-r'),cellfun(@(x,y) {x,y},ts{1}(posRange),num2cell(1:length(posRange),1)','UniformOutput',false));
        end
        
        
        
        if size(ts{2},1)==0
        elseif size(ts{2},1)==1
            plot([ts{2};ts{2}],repmat([21;21.8]+1,1,length(ts{2})),'-b');
        elseif size(ts{2},1)>20
            pos=round((length(ts{2})-20)/2);
            posRange=pos:pos+19;
            cellfun(@(x) plot([x{1}';x{1}'],repmat([x{2};x{2}+0.8]+1,1,length(x{1})),'-b'),cellfun(@(x,y) {x,y},ts{2}(posRange),num2cell([1:length(posRange)]+20,1)','UniformOutput',false));
        else
            posRange=1:size(ts{2},1);
            cellfun(@(x) plot([x{1}';x{1}'],repmat([x{2};x{2}+0.8]+1,1,length(x{1})),'-b'),cellfun(@(x,y) {x,y},ts{2}(posRange),num2cell([1:length(posRange)]+20,1)','UniformOutput',false));
        end
        
        
        xlim([-1,rLim]);
        ylim([0,41]);
        set(gca,'XTick',[],'YTick',[0,18,23,40],'YTickLabel',[0,20,0,20]);
        assignin('base','axT',gca());
        plotSegs(delayLen);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot('Position',[0.1,0.1,0.85,0.35]);
        hold on;
        if size(ts{1},1)==0
            pfHist=zeros(size(-4+binSize/2:binSize:rLim));
        elseif size(ts{1},1)==1
            pfHist=histcounts(ts{1},-4:binSize:rLim)./binSize;
        else
            pfHist=cell2mat(cellfun(@(x) histcounts(x,-4:binSize:rLim)./binSize,ts{1},'UniformOutput',false));
            cia=bootci(1000,@(x) mean(x), pfHist);
            fill([-4+binSize/2:binSize:rLim,rLim-binSize/2:-binSize:-4],[cia(1,:),fliplr(cia(2,:))],[1,0.8,0.8],'EdgeColor','none');
        end
        if size(ts{2},1)==0
            bnHist=zeros(size(-4+binSize/2:binSize:rLim));
        elseif size(ts{2},1)==1
            bnHist=histcounts(ts{2},-4:binSize:rLim)./binSize;
        else
            bnHist=cell2mat(cellfun(@(x) histcounts(x,-4:binSize:rLim)./binSize,ts{2},'UniformOutput',false));
            cib=bootci(1000,@(x) mean(x), bnHist);
            fill([-4+binSize/2:binSize:rLim,rLim-binSize/2:-binSize:-4],[cib(1,:),fliplr(cib(2,:))],[0.8,0.8,1],'EdgeColor','none');
        end
        
        
        plot(-4+binSize/2:binSize:rLim,(mean(pfHist,1))','-r');
        plot(-4+binSize/2:binSize:rLim,(mean(bnHist,1))','-b');
        
        xlim([-1,rLim]);
        ylim([min(ylim()),max([cia(:);cib(:)])]);
        assignin('base','axB',gca());
        plotSegs(delayLen);
        set(gca,'XTick',[0,5,10]);
        
        [~,lpwd,~]=fileparts(rootpath);
        lpwd=replace(lpwd,'_cleaned','');
        if ~exist(lpwd,'dir')
            mkdir(lpwd)
        end
        print('-dpng','-painters',sprintf('%s\\%s_%03d_%d.png',lpwd,lpwd,uid,delayLen));
        
    end




    function plotSegs(delayLen)
        vertLine=[0,1,1,2]+[0,0,ones(1,2).*delayLen];
        plot(repmat(vertLine,2,1),repmat(ylim()',1,length(vertLine)),'--k');
    end



end