function fh=sens_dur_SC(metaidx)
set(groot,'defaultTextFontSize',10);
meta=ephys.util.load_meta();
homedir=ephys.util.getHomedir('type','raw');
for ii=reshape(metaidx,1,[])
    fpath=fullfile(homedir,meta.allpath{ii},'FR_All_ 250.hdf5');
%     trials=h5read(fpath,'/Trials');
    suid=h5read(fpath,'/SU_id');
    fr=h5read(fpath,'/FR_All');
    
    sesspath=regexp(fpath,'(?<=\\)[0-9a-zA-Z_\-]*\\[0-9a-zA-Z_\-]*(?=\\FR_All)','match','once');
    fidx=ephys.path2sessid(sesspath);
    ephys.util.dependency('buz',false);
    [spkID,spkTS,trials,~,~]=ephys.getSPKID_TS(fidx);
    oneuid=suid(suid==meta.allcid(ii));
    FT_SPIKE=struct();
    FT_SPIKE.label={num2str(oneuid)};
    FT_SPIKE.timestamp={spkTS(spkID==oneuid)};
    sps=30000;
    cfg=struct();
    cfg.trl=[trials(:,1)-3*sps,trials(:,1)+11*sps,zeros(size(trials,1),1)-3*sps,trials];
    cfg.trlunit='timestamps';
    cfg.timestampspersecond=sps;
    FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);
    s1d3t=find(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==4 & trials(:,8)==3);
    s2d3t=find(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==8 & trials(:,8)==3);
    s1d6t=find(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==4 & trials(:,8)==6);
    s2d6t=find(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==8 & trials(:,8)==6);
    if min([numel(s1d3t),numel(s2d3t),numel(s1d6t),numel(s2d6t)])<2
        continue
    end

    fh=figure('Color','w');
    subplot(2,2,1)
    hold on
    s1d3p=s1d3t(max(floor(numel(s1d3t)./2)-4,1):min(floor(numel(s1d3t)./2)+5,numel(s1d3t)));
    s2d3p=s2d3t(max(floor(numel(s2d3t)./2)-4,1):min(floor(numel(s2d3t)./2)+5,numel(s2d3t)));
    pidx=0;
    for ti=reshape(s2d3p,1,[])
        ts=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==ti);
        plot(repmat(ts,2,1),repmat([pidx+0.1;pidx+0.9],1,numel(ts)),'b-')
        pidx=pidx+1;
    end
    for ti=reshape(s1d3p,1,[])
        ts=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==ti);
        plot(repmat(ts,2,1),repmat([pidx+0.1;pidx+0.9],1,numel(ts)),'r-')
        pidx=pidx+1;
    end
    arrayfun(@(x) xline(x,'--k'),[0 1 4 5]);
    xlim([-1,7])
    set(gca,'XTick',-1:7,'XTickLabel',-2:6)

    subplot(2,2,2)
    hold on
    s2d6p=s2d6t(max(floor(numel(s2d6t)./2)-4,1):min(floor(numel(s2d6t)./2)+5,numel(s2d6t)));
    s1d6p=s1d6t(max(floor(numel(s1d6t)./2)-4,1):min(floor(numel(s1d6t)./2)+5,numel(s1d6t)));
    pidx=0;
    for ti=reshape(s2d6p,1,[])
        ts=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==ti);
        plot(repmat(ts,2,1),repmat([pidx+0.1;pidx+0.9],1,numel(ts)),'b-')
        pidx=pidx+1;
    end
    for ti=reshape(s1d6p,1,[])
        ts=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==ti);
        plot(repmat(ts,2,1),repmat([pidx+0.1;pidx+0.9],1,numel(ts)),'r-')
        pidx=pidx+1;
    end

    xline(0,'--k')
    xline(1,'--k')
    xlim([-1,7])
    set(gca,'XTick',-1:7,'XTickLabel',-2:6)

    ax3=subplot(2,2,3);
    hold on
    s1d3c=squeeze(mean(fr(s1d3t,suid==meta.allcid(ii),5:44),2));
    s2d3c=squeeze(mean(fr(s2d3t,suid==meta.allcid(ii),5:44),2));
    
    s1d3ci=[-1;1]*std(s1d3c)./sqrt(size(s1d3c,1)).'+mean(s1d3c);
    s2d3ci=[-1;1]*std(s2d3c)./sqrt(size(s2d3c,1)).'+mean(s2d3c);
    
    fill([1:40,fliplr(1:40)],[smooth(s1d3ci(1,:),5);flip(smooth(s1d3ci(2,:),5))],'r','EdgeColor','none','FaceAlpha',0.1);
    fill([1:40,fliplr(1:40)],[smooth(s2d3ci(1,:),5);flip(smooth(s2d3ci(2,:),5))],'b','EdgeColor','none','FaceAlpha',0.1);
    
    plot(smooth(mean(s1d3c),5),'-r')
    plot(smooth(mean(s2d3c),5),'-b')
    arrayfun(@(x) xline(x,'--k'),[8.5 12.5 24.5 28.5]);
    title(ii)
    set(gca,'XTick',12.5:4:40.5,'XTickLabel',0:6)
    xlim([4.5,36.5]);
    xlabel('Time (s)');
    ylabel('F.R.');


    ax6=subplot(2,2,4);
    hold on
    s1d6c=squeeze(mean(fr(s1d6t,suid==meta.allcid(ii),5:44),2));
    s2d6c=squeeze(mean(fr(s2d6t,suid==meta.allcid(ii),5:44),2));
    
    s1d6ci=[-1;1]*std(s1d6c)./sqrt(size(s1d6c,1)).'+mean(s1d6c);
    s2d6ci=[-1;1]*std(s2d6c)./sqrt(size(s2d6c,1)).'+mean(s2d6c);
    
    fill([1:40,fliplr(1:40)],[smooth(s1d6ci(1,:),5);flip(smooth(s1d6ci(2,:),5))],'r','EdgeColor','none','FaceAlpha',0.1);
    fill([1:40,fliplr(1:40)],[smooth(s2d6ci(1,:),5);flip(smooth(s2d6ci(2,:),5))],'b','EdgeColor','none','FaceAlpha',0.1);

    plot(smooth(mean(s1d6c),5),'-r')
    plot(smooth(mean(s2d6c),5),'-b')
    
    xline(8.5,'--k')
    xline(12.5,'--k')
    title(ii)
    set(gca,'XTick',12.5:4:40.5,'XTickLabel',0:6)
    xlim([4.5,36.5]);
    xlabel('Time (s)');
    ylabel('F.R.');

    arrayfun(@(x) set(x,'YLim',[min([ax3.YLim,ax6.YLim].'),max([ax3.YLim,ax6.YLim]).']),[ax3,ax6]);
end
% exportgraphics(fh,'correct_error_sc.pdf','ContentType','vector')
