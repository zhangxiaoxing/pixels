function fh=sens_dur_SC(metaidx,su_meta,opt)
arguments
    metaidx
    su_meta
    opt.sens_meta
    opt.dur_meta   % backward incompatible change as of AUG-16-2022
    opt.anova_meta
    opt.skip_raster (1,1) logical = false
    opt.dim (1,1) double = 2
    opt.xlim_sec (1,2) double = [-1,7]
    opt.xlim_bin (1,2) double = [4.5,36.5]
    opt.skip_single_modal (1,1) logical = true
    opt.skip_fill (1,1) logical = false % do NOT plot variance (i.e. S.E.M) in psth
end
% close all
set(groot,'defaultTextFontSize',10);
homedir=ephys.util.getHomedir('type','raw');
for ii=reshape(metaidx,1,[])

    fpath=fullfile(homedir,su_meta.allpath{ii},'FR_All_ 250.hdf5');
    sessid=su_meta.sess(ii);
    trials=h5read(fpath,'/Trials');
    suid=h5read(fpath,'/SU_id');
    fr=h5read(fpath,'/FR_All');


    s1d3t=find(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==4 & trials(:,8)==3);
    s2d3t=find(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==8 & trials(:,8)==3);
    s1d6t=find(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==4 & trials(:,8)==6);
    s2d6t=find(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==8 & trials(:,8)==6);
    [mintrls,idx]=min([numel(s1d3t),numel(s2d3t),numel(s1d6t),numel(s2d6t)]);
    if mintrls<8
        disp("Minimal number of trials <8 "+num2str([mintrls,idx]))
        fh=[];
        continue
    end

    if opt.skip_single_modal
        fh=figure('Color','w','Position',[32,32,800,300]);
        tiledlayout(1,3)
    else
        fh=figure('Color','w','Position',[32,32,1280,720]);
        tiledlayout(2,5)
    end
    

    %% raster
    if ~opt.skip_raster
        ephys.util.dependency('buz',false);
        [spkID,spkTS,~,~,~]=ephys.getSPKID_TS(sessid);
        oneuid=suid(suid==su_meta.allcid(ii));
        FT_SPIKE=struct();
        FT_SPIKE.label={num2str(oneuid)};
        FT_SPIKE.timestamp={spkTS(spkID==oneuid)};
        sps=30000;
        cfg=struct();
        cfg.trl=[trials(:,1)-3*sps,trials(:,1)+11*sps,zeros(size(trials,1),1)-3*sps,trials];
        cfg.trlunit='timestamps';
        cfg.timestampspersecond=sps;
        FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);
        s1d3p=s1d3t(max(floor(numel(s1d3t)./2)-3,1):min(floor(numel(s1d3t)./2)+4,numel(s1d3t)));
        s2d3p=s2d3t(max(floor(numel(s2d3t)./2)-3,1):min(floor(numel(s2d3t)./2)+4,numel(s2d3t)));
        s2d6p=s2d6t(max(floor(numel(s2d6t)./2)-3,1):min(floor(numel(s2d6t)./2)+4,numel(s2d6t)));
        s1d6p=s1d6t(max(floor(numel(s1d6t)./2)-3,1):min(floor(numel(s1d6t)./2)+4,numel(s1d6t)));

        %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        if opt.skip_single_modal
            nexttile(1) %d3, s1 vs s2
        else
            nexttile(5) %d3, s1 vs s2
        end
        hold on

        pidx=1;
        for ti=reshape(s2d3p,1,[])
            ts=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==ti & FT_SPIKE.time{1}<4);
            plot(repmat(ts,2,1),repmat([pidx+0.1;pidx+0.9],1,numel(ts)),'-','Color',[0.5,0.5,0.48]) %previous color 'c'
            pidx=pidx+1;
        end
        pidx=pidx+2;
        for ti=reshape(s1d3p,1,[])
            ts=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==ti & FT_SPIKE.time{1}<4);
            plot(repmat(ts,2,1),repmat([pidx+0.1;pidx+0.9],1,numel(ts)),'r','Color',[0.5,0.48,0.48])
            pidx=pidx+1;
        end
        if opt.skip_single_modal
            line([0,1,4;0,1,4],[0,0,0;20,20,20],'LineStyle','--','Color','k')
            xlim(opt.xlim_sec)
            set(gca,'XTick',-1:8,'XTickLabel',-2:7)

            nexttile(2) %d3, s1 vs s2
            pidx=1;
            hold on
        else
            pidx=pidx+2;
        end
        
        for ti=reshape(s2d6p,1,[])
            ts=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==ti);
            plot(repmat(ts,2,1),repmat([pidx+0.1;pidx+0.9],1,numel(ts)),'b','Color',[0.48,0.48,0.5])
            pidx=pidx+1;
        end
        pidx=pidx+2;
        for ti=reshape(s1d6p,1,[])
            ts=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==ti);
            plot(repmat(ts,2,1),repmat([pidx+0.1;pidx+0.9],1,numel(ts)),'m','Color',[0.5,0.48,0.5])
            pidx=pidx+1;
        end
        if opt.skip_single_modal
            line([0,1,7;0,1,7],[0,0,0;20,20,20],'LineStyle','--','Color','k')
        else
            line([0,0,1,1,4,7;0,0,1,1,4,7],[0,20,0,20,0,20;20,40,20,40,20,40],'LineStyle','--','Color','k')
        end
        xlim(opt.xlim_sec)
        set(gca,'XTick',-1:8,'XTickLabel',-2:7)

        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        if ~opt.skip_single_modal
            nexttile(1) %d3, s1 vs s2
            hold on

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

            arrayfun(@(x) xline(x,'--k'),[0 1 4 5 7 8]);
            xlim(opt.xlim_sec)
            set(gca,'XTick',-1:8,'XTickLabel',-2:7)

            %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            nexttile(2) %s2 d3 v d6
            hold on

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

            arrayfun(@(x) xline(x,'--k'),[0 1 4 5 7 8]);
            xlim(opt.xlim_sec)
            set(gca,'XTick',-1:8,'XTickLabel',-2:7)
            % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


            nexttile(3) %s1,d3 v d6
            hold on

            pidx=0;
            for ti=reshape(s1d3p,1,[])
                ts=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==ti);
                plot(repmat(ts,2,1),repmat([pidx+0.1;pidx+0.9],1,numel(ts)),'-','Color',[0,0,0.5].*opt.dim)
                pidx=pidx+1;
            end

            for ti=reshape(s1d6p,1,[])
                ts=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==ti);
                plot(repmat(ts,2,1),repmat([pidx+0.1;pidx+0.9],1,numel(ts)),'-','Color',[0.5,0,0].*opt.dim)
                pidx=pidx+1;
            end

            arrayfun(@(x) xline(x,'--k'),[0 1 4 5 7 8]);
            xlim(opt.xlim_sec)
            set(gca,'XTick',-1:8,'XTickLabel',-2:7)



            nexttile(4) %s2 d3 v d6
            hold on

            pidx=0;
            for ti=reshape(s2d3p,1,[])
                ts=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==ti);
                plot(repmat(ts,2,1),repmat([pidx+0.1;pidx+0.9],1,numel(ts)),'-','Color',[0,0,0.5].*opt.dim)
                pidx=pidx+1;
            end

            for ti=reshape(s2d6p,1,[])
                ts=FT_SPIKE.time{1}(FT_SPIKE.trial{1}==ti);
                plot(repmat(ts,2,1),repmat([pidx+0.1;pidx+0.9],1,numel(ts)),'r-','Color',[0.5,0,0].*opt.dim)
                pidx=pidx+1;
            end

            arrayfun(@(x) xline(x,'--k'),[0 1 4 5 7 8]);
            xlim(opt.xlim_sec)
            set(gca,'XTick',-1:8,'XTickLabel',-2:7)

        end
    end
    %%
    if ~opt.skip_single_modal
        ax3=nexttile(6);
        hold on
        s1d3c=squeeze(fr(s1d3t,suid==su_meta.allcid(ii),5:44));
        s2d3c=squeeze(fr(s2d3t,suid==su_meta.allcid(ii),5:44));

        %     mean(s2d3c)
        s1d3ci=[-1;1]*std(s1d3c)./sqrt(size(s1d3c,1)).'+mean(s1d3c);
        s2d3ci=[-1;1]*std(s2d3c)./sqrt(size(s2d3c,1)).'+mean(s2d3c);

        fill([1:40,fliplr(1:40)],[smooth(s1d3ci(1,:),5);flip(smooth(s1d3ci(2,:),5))],'r','EdgeColor','none','FaceAlpha',0.1);
        fill([1:40,fliplr(1:40)],[smooth(s2d3ci(1,:),5);flip(smooth(s2d3ci(2,:),5))],'b','EdgeColor','none','FaceAlpha',0.1);

        plot(smooth(mean(s1d3c),5),'-r')
        plot(smooth(mean(s2d3c),5),'-b')

        title('s1 v. s2, d3')
        arrayfun(@(x) xline(x,'--k'),[8.5 12.5 24.5 28.5 36.5 40.5]);
        set(gca,'XTick',12.5:4:40.5,'XTickLabel',0:7,'XLim',opt.xlim_bin)
        xlabel('Time (s)');
        ylabel('F.R.');

        if isfield(opt,'sens_meta')
            arrayfun(@(x) text(10.5+x*4,max(ylim()),formatp(opt.sens_meta.fdr_d3(ii,x)),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9),1:3)
        end
        %%
        ax6=nexttile(7);
        hold on
        s1d6c=squeeze(fr(s1d6t,suid==su_meta.allcid(ii),5:44));
        s2d6c=squeeze(fr(s2d6t,suid==su_meta.allcid(ii),5:44));

        s1d6ci=[-1;1]*std(s1d6c)./sqrt(size(s1d6c,1)).'+mean(s1d6c);
        s2d6ci=[-1;1]*std(s2d6c)./sqrt(size(s2d6c,1)).'+mean(s2d6c);

        fill([1:40,fliplr(1:40)],[smooth(s1d6ci(1,:),5);flip(smooth(s1d6ci(2,:),5))],'r','EdgeColor','none','FaceAlpha',0.1);
        fill([1:40,fliplr(1:40)],[smooth(s2d6ci(1,:),5);flip(smooth(s2d6ci(2,:),5))],'b','EdgeColor','none','FaceAlpha',0.1);

        plot(smooth(mean(s1d6c),5),'-r')
        plot(smooth(mean(s2d6c),5),'-b')

        title('s1 v. s2, d6')
        arrayfun(@(x) xline(x,'--k'),[8.5 12.5 24.5 28.5 36.5 40.5]);
        set(gca,'XTick',12.5:4:40.5,'XTickLabel',0:7,'XLim',opt.xlim_bin)
        xlabel('Time (s)');
        ylabel('F.R.');

        if isfield(opt,'sens_meta')
            arrayfun(@(x) text(10.5+x*4,max(ylim()),formatp(opt.sens_meta.fdr_d6(ii,x)),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9),1:6)
        end

        %%
        axs1=nexttile(8);
        hold on
        s1d3c=squeeze(fr(s1d3t,suid==su_meta.allcid(ii),5:44));
        s1d6c=squeeze(fr(s1d6t,suid==su_meta.allcid(ii),5:44));

        s1d3ci=[-1;1]*std(s1d3c)./sqrt(size(s1d3c,1)).'+mean(s1d3c);
        s1d6ci=[-1;1]*std(s1d6c)./sqrt(size(s1d6c,1)).'+mean(s1d6c);

        fill([1:40,fliplr(1:40)],[smooth(s1d3ci(1,:),5);flip(smooth(s1d3ci(2,:),5))],[0,0,0.5].*opt.dim,'EdgeColor','none','FaceAlpha',0.1);
        fill([1:40,fliplr(1:40)],[smooth(s1d6ci(1,:),5);flip(smooth(s1d6ci(2,:),5))],[0.5,0,0].*opt.dim,'EdgeColor','none','FaceAlpha',0.1);

        plot(smooth(mean(s1d3c),5),'-','Color',[0,0,0.5].*opt.dim)
        plot(smooth(mean(s1d6c),5),'-','Color',[0.5,0,0].*opt.dim)
        title('d3 v. d6 ,s1')
        arrayfun(@(x) xline(x,'--k'),[8.5 12.5 24.5 28.5 36.5 40.5]);
        set(gca,'XTick',12.5:4:40.5,'XTickLabel',0:7,'XLim',opt.xlim_bin)
        xlabel('Time (s)');
        ylabel('F.R.');
        if isfield(opt,'dur_meta')
            arrayfun(@(x) text(10.5+x*4,max(ylim()),formatp(opt.dur_meta.fdr_s1(ii,x+1)),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9),1:3)
        end
        %%

        axs2=nexttile(9);
        hold on
        s2d3c=squeeze(fr(s2d3t,suid==su_meta.allcid(ii),5:44));
        s2d6c=squeeze(fr(s2d6t,suid==su_meta.allcid(ii),5:44));

        s2d3ci=[-1;1]*std(s2d3c)./sqrt(size(s2d3c,1)).'+mean(s2d3c);
        s2d6ci=[-1;1]*std(s2d6c)./sqrt(size(s2d6c,1)).'+mean(s2d6c);

        fill([1:40,fliplr(1:40)],[smooth(s2d3ci(1,:),5);flip(smooth(s2d3ci(2,:),5))],[0,0,0.5].*opt.dim,'EdgeColor','none','FaceAlpha',0.1);
        fill([1:40,fliplr(1:40)],[smooth(s2d6ci(1,:),5);flip(smooth(s2d6ci(2,:),5))],[0.5,0,0].*opt.dim,'EdgeColor','none','FaceAlpha',0.1);

        plot(smooth(mean(s2d3c),5),'-','Color',[0,0,0.5].*opt.dim)
        plot(smooth(mean(s2d6c),5),'-','Color',[0.5,0,0].*opt.dim)


        title('d3 v. d6 ,s2')
        arrayfun(@(x) xline(x,'--k'),[8.5 12.5 24.5 28.5 36.5 40.5]);
        set(gca,'XTick',12.5:4:40.5,'XTickLabel',0:7,'XLim',opt.xlim_bin)
        xlabel('Time (s)');
        ylabel('F.R.');
        if isfield(opt,'sens_meta')
            arrayfun(@(x) text(10.5+x*4,max(ylim()),formatp(opt.dur_meta.fdr_s2(ii,x+1)),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9),1:3)
        end

        arrayfun(@(x) set(x,'YLim',[min([ax3.YLim,ax6.YLim].'),max([ax3.YLim,ax6.YLim]).']),[ax3,ax6,axs1,axs2]);

        %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        arrayfun(@(x) set(x,'YLim',[min([ax3.YLim,ax6.YLim].'),max([ax3.YLim,ax6.YLim]).']),[ax3,ax6,axs1,axs2]);
    end
    if opt.skip_single_modal
        axs2=nexttile(3);
    else
        axs2=nexttile(10);
    end

    hold on
    s1d3c=squeeze(fr(s1d3t,suid==su_meta.allcid(ii),5:44));
    s1d6c=squeeze(fr(s1d6t,suid==su_meta.allcid(ii),5:44));
    s2d3c=squeeze(fr(s2d3t,suid==su_meta.allcid(ii),5:44));
    s2d6c=squeeze(fr(s2d6t,suid==su_meta.allcid(ii),5:44));


    s1d3ci=[-1;1]*std(s1d3c)./sqrt(size(s1d3c,1)).'+mean(s1d3c);
    s1d6ci=[-1;1]*std(s1d6c)./sqrt(size(s1d6c,1)).'+mean(s1d6c);
    s2d3ci=[-1;1]*std(s2d3c)./sqrt(size(s2d3c,1)).'+mean(s2d3c);
    s2d6ci=[-1;1]*std(s2d6c)./sqrt(size(s2d6c,1)).'+mean(s2d6c);
    
%     pd3=@(x) x(1:24);
    if ~opt.skip_fill
        fill([1:24,fliplr(1:24)],[smoothdata(s1d3ci(1,1:24),'sgolay',5),flip(smoothdata(s1d3ci(2,1:24),'sgolay',5))],'r','EdgeColor','none','FaceAlpha',0.1);
        fill([1:40,fliplr(1:40)],[smoothdata(s1d6ci(1,:),'sgolay',5),flip(smoothdata(s1d6ci(2,:),'sgolay',5))],'m','EdgeColor','none','FaceAlpha',0.1);
    
        fill([1:24,fliplr(1:24)],[smoothdata(s2d3ci(1,1:24),'sgolay',5),flip(smoothdata(s2d3ci(2,1:24),'sgolay',5))],'k','EdgeColor','none','FaceAlpha',0.1);
        fill([1:40,fliplr(1:40)],[smoothdata(s2d6ci(1,:),'sgolay',5),flip(smoothdata(s2d6ci(2,:),'sgolay',5))],'b','EdgeColor','none','FaceAlpha',0.1);
    end
    plot(smoothdata(mean(s1d3c(:,1:24)),'sgolay',5),'-','Color','r','LineWidth',1)
    plot(smoothdata(mean(s1d6c),'sgolay',5),'-','Color','m','LineWidth',1)

    plot(smoothdata(mean(s2d3c(:,1:24)),'sgolay',5),'-','Color','k','LineWidth',1)
    plot(smoothdata(mean(s2d6c),'sgolay',5),'-','Color','b','LineWidth',1)


    title('4 conditions')
    xline([8.5 12.5 24.5 36.5],'--k','Alpha',1);

%     line([8.5 12.5 24.5 36.5],
    set(gca,'XTick',12.5:4:40.5,'XTickLabel',0:7,'XLim',opt.xlim_bin)
    xlabel('Time (s)');
    ylabel('F.R.');
    if isfield(opt,'sens_meta')
        arrayfun(@(x) text(10.5+x*4,max(ylim()),formatp(opt.dur_meta.fdr_s2(ii,x+1)),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9),1:3)
    end
    
    %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    %TODO stats

    if exist('opt','var') && isfield(opt,'anova_meta')
        ttstr=sprintf('#%d,ANOVA SEN%.3f,DUR%.3f,INT%.3f',ii,opt.anova_meta.anovap(ii,1),opt.anova_meta.anovap(ii,2),opt.anova_meta.anovap(ii,3));
        sgtitle(ttstr);
    else
        sgtitle(ii);
    end
    % exportgraphics(fh,'correct_error_sc.pdf','ContentType','vector')
end
end
function s=formatp(v)
if v>=0.1
    s='ns';
elseif v<0.05
    s=sprintf('%0.3f',v);
    s=s(2:end);
else
    s=sprintf('%0.3f',v);
    s=['#',s(3:end)];
end
end