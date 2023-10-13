meta3=ephys.util.load_meta('delay',3);
meta6=ephys.util.load_meta('delay',6);

load('cross_ep_anovameta.mat','cross_ep_anovameta'); % K:\code\jpsth\+ephys\selectivity_anova.m


bothsel=any(meta3.wrs_p(5:7,:)<0.01 & meta6.wrs_p(5:7,:)<0.01);
selidxsel=any(abs(meta3.selec(5:7,:))>0.4 & abs(meta6.selec(5:7,:))>0.4);
regsel=strcmp(meta6.reg_tree(1,:),'CH') & cellfun(@(x) ~isempty(x)>0,meta6.reg_tree(5,:));

anysel=any(meta3.wrs_p(5:7,:)<0.01 | meta6.wrs_p(5:7,:)<0.01);
anyselidxsel=any(abs(meta3.selec(5:7,:))>0.4 | abs(meta6.selec(5:7,:))>0.4);

bothscsel=bothsel & regsel & selidxsel;
eitherscsel=anysel & ~bothsel & regsel & anyselidxsel;

homedir=ephys.util.getHomedir();

sem=@(x) std(x)./sqrt(size(x,1));
semfill=@(x) [smooth(mean(x)-sem(x),3);smooth(fliplr(mean(x)+sem(x)),3)];

get10=@(x) x(ceil(numel(x)/2)+(-4:5));

%% either
usess=unique(meta6.sess(eitherscsel | bothscsel));
for sess=reshape(usess,1,[])
    folder=replace(ephys.sessid2path(sess,'criteria','WT'),'\',filesep());
    trials=h5read(fullfile(homedir,folder,'FR_All_ 250.hdf5'),'/Trials');

    d3s1=trials(:,5)==4 & trials(:,8)==3 & all(trials(:,9:10)~=0,2);
    d6s1=trials(:,5)==4 & trials(:,8)==6 & all(trials(:,9:10)~=0,2);
    d3s2=trials(:,5)==8 & trials(:,8)==3 & all(trials(:,9:10)~=0,2);
    d6s2=trials(:,5)==8 & trials(:,8)==6 & all(trials(:,9:10)~=0,2);

    if min([nnz(d3s1),nnz(d3s2),nnz(d6s1),nnz(d6s2)])<10
        continue
    end
    
    SU_id=h5read(fullfile(homedir,folder,'FR_All_ 250.hdf5'),'/SU_id');
    FR=h5read(fullfile(homedir,folder,'FR_All_ 250.hdf5'),'/FR_All');

    % both will be temporarily negative 
%     sessSU=[int32(meta6.allcid(meta6.sess==sess & eitherscsel.'));...
%         -int32(meta6.allcid(meta6.sess==sess & bothscsel.'))];
    
    for tsu=reshape(sessSU,1,[])
        su=abs(tsu);
        
        if mean(FRsu(d3s1|d3s2|d6s1|d6s2,17:28),'all')<4
            continue
        end
        
        FRsu=squeeze(FR(:,SU_id==su,:));
        [~,~,~,~,~,FT_SPIKE]=ephys.getSPKID_TS(sess,'suids',[su],'keep_trial',true);
        
        fh=figure('Color','w');
        subplot(2,2,1);hold on;
        yidx=1;
        for tt=reshape(get10(find(d3s1)),1,[])
            xx=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==tt);
            yy=[yidx-0.4;yidx+0.4]*ones(size(xx));
            plot([1;1]*xx,yy,'-b');
            yidx=yidx+1;
        end
        for tt=reshape(get10(find(d3s2)),1,[])
            xx=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==tt);
            yy=[yidx-0.4;yidx+0.4]*ones(size(xx));
            plot([1;1]*xx,yy,'-r');
            yidx=yidx+1;
        end
        ylim([0,21])
        xlim([-1,4])
        arrayfun(@(x) xline(x,'k:'),[0,1]);
        set(gca,'XTick',-2:3:4,'XTickLabel',-3:3:3)

        
        subplot(2,2,2);hold on;
        yidx=1;
        for tt=reshape(get10(find(d6s1)),1,[])
            xx=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==tt);
            yy=[yidx-0.4;yidx+0.4]*ones(size(xx));
            plot([1;1]*xx,yy,'-b');
            yidx=yidx+1;
        end
        for tt=reshape(get10(find(d6s2)),1,[])
            xx=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==tt);
            yy=[yidx-0.4;yidx+0.4]*ones(size(xx));
            plot([1;1]*xx,yy,'-r');
            yidx=yidx+1;
        end
        ylim([0,21])
        xlim([-1,7])
        arrayfun(@(x) xline(x,'k:'),[0,1]);
        set(gca,'XTick',-2:3:7,'XTickLabel',-3:3:6)

        xx=[1:56,56:-1:1]*0.25-4+0.125;

        subplot(2,2,3);hold on;
        fill(xx,semfill(FRsu(d3s1,:)),'b','EdgeColor','none','FaceAlpha',0.2)
        fill(xx,semfill(FRsu(d3s2,:)),'r','EdgeColor','none','FaceAlpha',0.2)
        plot(xx(1:56),smooth(mean(FRsu(d3s1,:)),3),'b-');
        plot(xx(1:56),smooth(mean(FRsu(d3s2,:)),3),'r-');

        arrayfun(@(x) xline(x,'k:'),[-1,0])
        set(gca,'XTick',-3:3:3)
        xlabel('Time (s)')
        ylabel('Firing rate (Hz)')
        xlim([-2,3])
        text(-1.5:1:2.5,ones(1,5)*max(ylim()),num2str(meta3.wrs_p(3:7,meta3.allcid==su & meta3.sess==sess),'%.3f'),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',5)
        

        subplot(2,2,4);hold on;
        fill(xx,semfill(FRsu(d6s1,:)),'b','EdgeColor','none','FaceAlpha',0.2)
        fill(xx,semfill(FRsu(d6s2,:)),'r','EdgeColor','none','FaceAlpha',0.2)
        plot(xx(1:56),smooth(mean(FRsu(d6s1,:)),3),'b--');
        plot(xx(1:56),smooth(mean(FRsu(d6s2,:)),3),'r--');

        arrayfun(@(x) xline(x,'k:'),[-1,0])
        set(gca,'XTick',-3:3:6)
        xlabel('Time (s)')
        ylabel('Firing rate (Hz)')
        xlim([-2,6])
        text(-1.5:1:5.5,ones(1,8)*max(ylim()),num2str(meta6.wrs_p(3:10,meta6.allcid==su & meta6.sess==sess),'%.3f'),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',6)
        if tsu<0
            sgtitle(sprintf('S%d#%d, both',sess,su));
        else
            sgtitle(sprintf('S%d#%d, either',sess,su));
        end
        
%         exportgraphics(fh,fullfile('SC',sprintf('Both_either_S%dC%d.png',sess,su)),'ContentType','image');
%         close(fh)
        exportgraphics(fh,sprintf('Both_either_S%dC%d.pdf',sess,su),'ContentType','vector');
    end
end




