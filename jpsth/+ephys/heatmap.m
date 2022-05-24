ephys.util.dependency
dualprobe=dir(fullfile(homedir,'**','spike_info.mat'));
% accu=0;
sps=30000;
for i=1:length(dualprobe)
    disp(i)
    folder=dualprobe(i).folder;
    if isempty(behav.procPerf(h5read(fullfile(folder,'events.hdf5'),'/trials')',''))
        continue
    end
%     if accu>8
%         return
%     end
%     accu=accu+1;
    
    fstr=load(fullfile(dualprobe(i).folder,dualprobe(i).name));
    spkID=[fstr.spike_info{1}{1};fstr.spike_info{1}{2}];
    spkTS=[fstr.spike_info{2}{1};fstr.spike_info{2}{2}];
    suids=ephys.goodCid(dualprobe(i).folder);
    susel=ismember(spkID,suids);
    spkID=double(spkID(susel));
    spkTS=double(spkTS(susel));
    
    trials=h5read(fullfile(folder,'events.hdf5'),'/trials')';
    trials=behav.procPerf(trials);
    
    sel_6s_S1=trials(:,5)==4 & trials(:,8)== 6;
    sel_6s_S2=trials(:,5)==8 & trials(:,8)== 6;
    
    sel_3s_S1=trials(:,5)==4 & trials(:,8)== 3;
    sel_3s_S2=trials(:,5)==8 & trials(:,8)== 3;
    
    FT_SPIKE=struct();

    FT_SPIKE.label=strtrim(cellstr(num2str(suids)));
    FT_SPIKE.timestamp=cell(1,numel(suids));
    for su=1:numel(suids)
        FT_SPIKE.timestamp{su}=spkTS(spkID==suids(su))';
    end
    %  continuous format F T struct file
    cfg=struct();
    cfg.trl=[trials(:,1)-3*sps,trials(:,1)+11*sps,zeros(size(trials,1),1)-3*sps,trials];
    cfg.trlunit='timestamps';
    cfg.timestampspersecond=sps;
    
    FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);
    
    heatS1=zeros(size(suids,1),56);
    heatS2=zeros(size(suids,1),56);

    for su=1:numel(suids)
        S1=cell2mat(arrayfun(@(x) histcounts(FT_SPIKE.time{su}(FT_SPIKE.trial{su}==x),-3:0.25:11),find(sel_6s_S1),'UniformOutput',false));
        S2=cell2mat(arrayfun(@(x) histcounts(FT_SPIKE.time{su}(FT_SPIKE.trial{su}==x),-3:0.25:11),find(sel_6s_S2),'UniformOutput',false));
        base=cell2mat(arrayfun(@(x) histcounts(FT_SPIKE.time{su}(FT_SPIKE.trial{su}==x),-3:0.25:0),1:size(trials,1)','UniformOutput',false));
        bmm=mean(base);
        bst=std(base);
        heatS1(su,:)=(mean(S1)-bmm)./bst;
        heatS2(su,:)=(mean(S2)-bmm)./bst;
    end
    [~,idx]=sort(sum(heatS2(:,13:40),2)-sum(heatS1(:,13:40),2));
    fh=figure('Color','w','Position',[100,100,800,600]);
    for subidx=1:2
        subplot(1,2,subidx);
        imagesc(eval(sprintf('heatS%d(idx,:)',subidx)),[-2,2]);
        colormap('jet');
        arrayfun(@(x) xline(x,'--w','LineWidth',1,'Alpha',1),[12.5,16.5,40.5,44.5]);
        colorbar();
        set(gca,'XTick',[12.5,32.5,52.5],'XTickLabel',[0,5,10]);
        xlabel('Time (s)');
        ylabel('Neuron #');
    end
    exportgraphics(fh,sprintf('heatmap_dual_f%d.pdf',i));
    close all;
    
end