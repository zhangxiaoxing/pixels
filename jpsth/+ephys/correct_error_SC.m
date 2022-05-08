type ='Transient';
set(groot,'defaultTextFontSize',10);
meta=ephys.util.load_meta();
% cadidate=find(meta.per_bin_6(6,:)==1 & cellfun(@(x) ~isempty(x),meta.reg_tree(5,:)) & strcmp(meta.reg_tree(2,:),'CTX') & ismember(meta.mem_type_6,2:2:4) & meta.good_waveform.' & abs(meta.selec(10,:))>0.2 & meta.wrs_p(10,:)<0.01);
homedir=ephys.util.getHomedir('type','raw');
for ii=642%26328%[3993,16071,16763,26328]%18030,14337,26179
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
    cs1t=find(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==4 & trials(:,8)==6);
    cs2t=find(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==8 & trials(:,8)==6);
    es1t=find(trials(:,10)==0 &trials(:,5)==4 & trials(:,8)==6);
    es2t=find(trials(:,10)==0 &trials(:,5)==8 & trials(:,8)==6);
    if min([numel(cs1t),numel(cs2t),numel(es1t),numel(es2t)])<2
        continue
    end

    fh=figure('Color','w')
    subplot(2,1,1)
    hold on
    cs1p=cs1t(max(floor(numel(cs1t)./2)-4,1):min(floor(numel(cs1t)./2)+5,numel(cs1t)));
    cs2p=cs2t(max(floor(numel(cs2t)./2)-4,1):min(floor(numel(cs2t)./2)+5,numel(cs2t)));
    es2p=es2t(max(floor(numel(es2t)./2)-4,1):min(floor(numel(es2t)./2)+5,numel(es2t)));
    es1p=es1t(max(floor(numel(es1t)./2)-4,1):min(floor(numel(es1t)./2)+5,numel(es1t)));
    pidx=0;
    for ti=reshape(cs2p,1,[])
        ts=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==ti);
        plot(repmat(ts,2,1),repmat([pidx+0.1;pidx+0.9],1,numel(ts)),'b-')
        pidx=pidx+1;
    end
%     for ti=reshape(es1p,1,[])
%         ts=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==ti);
%         plot(repmat(ts,2,1),repmat([pidx+0.1;pidx+0.9],1,numel(ts)),'k-')
%         pidx=pidx+1;
%     end
    for ti=reshape(es2p,1,[])
        ts=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==ti);
        plot(repmat(ts,2,1),repmat([pidx+0.1;pidx+0.9],1,numel(ts)),'k-')
        pidx=pidx+1;
    end
    for ti=reshape(cs1p,1,[])
        ts=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==ti);
        plot(repmat(ts,2,1),repmat([pidx+0.1;pidx+0.9],1,numel(ts)),'r-')
        pidx=pidx+1;
    end
    xline(0,'--k')
    xline(1,'--k')
    xlim([-1,7])
    set(gca,'XTick',-1:7,'XTickLabel',-2:6)
    subplot(2,1,2)
    hold on
    cs1c=squeeze(mean(fr(cs1t,suid==meta.allcid(ii),5:44),2));
    cs2c=squeeze(mean(fr(cs2t,suid==meta.allcid(ii),5:44),2));
    es1c=squeeze(mean(fr(es1t,suid==meta.allcid(ii),5:44),2));
    es2c=squeeze(mean(fr(es2t,suid==meta.allcid(ii),5:44),2));
    
    cs1ci=bootci(500,@(x) mean(x), cs1c);
    cs2ci=bootci(500,@(x) mean(x), cs2c);
    es1ci=bootci(500,@(x) mean(x), es1c);
    es2ci=bootci(500,@(x) mean(x), es2c);
    
    fill([1:40,fliplr(1:40)],[smooth(cs1ci(1,:),5);flip(smooth(cs1ci(2,:),5))],'r','EdgeColor','none','FaceAlpha',0.1);
    fill([1:40,fliplr(1:40)],[smooth(cs2ci(1,:),5);flip(smooth(cs2ci(2,:),5))],'b','EdgeColor','none','FaceAlpha',0.1);
    fill([1:40,fliplr(1:40)],[smooth(es2ci(1,:),5);flip(smooth(es2ci(2,:),5))],'k','EdgeColor','none','FaceAlpha',0.1);
%     fill([1:40,fliplr(1:40)],[smooth(es1ci(1,:),5);flip(smooth(es1ci(2,:),5))],'k','EdgeColor','none','FaceAlpha',0.1);
    
    plot(smooth(mean(cs1c),5),'-r')
    plot(smooth(mean(cs2c),5),'-b')
    plot(smooth(mean(es2c),5),'-k')
%     plot(smooth(mean(es1c),5),'-k')
    xline(8.5,'--k')
    xline(12.5,'--k')
    title(ii)
    set(gca,'XTick',12.5:4:40.5,'XTickLabel',0:6)
    xlim([4.5,36.5]);
    xlabel('Time (s)');
    ylabel('F.R.');
    keyboard();
end
exportgraphics(fh,'correct_error_sc.pdf','ContentType','vector')
