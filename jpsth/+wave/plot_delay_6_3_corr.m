function plot_delay_6_3_corr(opt)
arguments
    opt.plot_3_6early (1,1) logical = false
    opt.plot_3_6late (1,1) logical = false
    opt.plot_3_6full (1,1) logical = true
    opt.plot_6early_6late (1,1) logical = false
    opt.plot_6_6early (1,1) logical = false
end

persistent com_str_early3 com_str_late3 com_str_3 com_str_6

%% 3 and 6early3
if opt.plot_3_6early
    if isempty(com_str_early3) || isempty(com_str_3)
        com_str_early3=wave.get_com_map('partial','early3in6');
        com_str_3=wave.get_com_map('delay',3);
    end
    allsess=intersect(fieldnames(com_str_3),fieldnames(com_str_early3));
    corrmat=[];
    for sess=reshape(allsess,1,[])
        for pref=["s1","s2"]
            allcid=intersect(cell2mat(com_str_early3.(sess{1}).(pref).keys),cell2mat(com_str_3.(sess{1}).(pref).keys));
            for onecid=allcid(:).'
                corrmat=[corrmat;...
                    com_str_early3.(sess{1}).(pref)(onecid),com_str_3.(sess{1}).(pref)(onecid)];
            end
        end
    end
    [r3e,p3e]=corr(corrmat(:,1),corrmat(:,2));
    % rci3e=bootci(500,@(x) corr(x(:,1),x(:,2)),corrmat);
    binavg=[];
    [~,~,kbin]=histcounts(corrmat(:,2),0:2:12);
    for bb=reshape(unique(kbin),1,[])
        ci=bootci(500,@(x) median(x),corrmat(kbin==bb,1));% ci of 6early3
        binavg=[binavg;median(corrmat(kbin==bb,1)),median(corrmat(kbin==bb,2)),ci(1),ci(2)];%[6early3,3]
    end

    fh=figure('Color','w','Position',[100,100,230,260]);
    hold on
    pha=plot(binavg(:,2)./4,binavg(:,1)./4,'-k','LineWidth',1);
    errorbar(binavg(:,2)./4,binavg(:,1)./4,diff(binavg(:,2:3),1,2)./4,diff(binavg(:,2:2:4),1,2)./4,'k.','LineWidth',1)
    scatter(corrmat(:,2)./4,corrmat(:,1)./4,4,'k','filled','o','MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
    text(max(xlim()),max(ylim()),sprintf('r=%.3f,p=%.3f,n=%d',r3e,p3e,size(corrmat,1)),'HorizontalAlignment','right','VerticalAlignment','top');
    phi=plot([0,3],[0,3],'--r','LineWidth',1);
    ylabel('TCOM in 6s delay trials early half (s)');
    xlabel('TCOM in 3s delay trials (s)');
    legend([pha,phi],{'data-median','identical'},'Location','northoutside')
    exportgraphics(fh,'TCOM_corr_3s_6early3.pdf','ContentType','vector');
end

%% 3 and 6late3
if opt.plot_3_6late
    if isempty(com_str_late3) || isempty(com_str_3)
        com_str_late3=wave.get_com_map('partial','late3in6');
        com_str_3=wave.get_com_map('delay',3);
    end
    allsess=intersect(fieldnames(com_str_3),fieldnames(com_str_late3));
    corrmat=[];
    for sess=reshape(allsess,1,[])
        for pref=["s1","s2"]
            allcid=intersect(cell2mat(com_str_late3.(sess{1}).(pref).keys),cell2mat(com_str_3.(sess{1}).(pref).keys));
            for onecid=allcid(:).'
                corrmat=[corrmat;...
                    com_str_late3.(sess{1}).(pref)(onecid),com_str_3.(sess{1}).(pref)(onecid)];
            end
        end
    end
    [r3e,p3e]=corr(corrmat(:,1),corrmat(:,2));
    % rci3e=bootci(500,@(x) corr(x(:,1),x(:,2)),corrmat);
    binavg=[];
    [~,~,kbin]=histcounts(corrmat(:,2),0:2:12);
    for bb=reshape(unique(kbin),1,[])
        ci=bootci(500,@(x) median(x),corrmat(kbin==bb,1));% ci of 6late3
        binavg=[binavg;median(corrmat(kbin==bb,1)),median(corrmat(kbin==bb,2)),ci(1),ci(2)];%[6late3,3]
    end

    fh=figure('Color','w','Position',[100,100,230,260]);
    hold on
    pha=plot(binavg(:,2)./4,binavg(:,1)./4,'-k','LineWidth',1);
    errorbar(binavg(:,2)./4,binavg(:,1)./4,diff(binavg(:,[1 3]),1,2)./4,diff(binavg(:,[1 4]),1,2)./4,'k.','LineWidth',1)
    scatter(corrmat(:,2)./4,corrmat(:,1)./4,4,'k','filled','o','MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
    text(max(xlim()),max(ylim()),sprintf('r=%.3f,p=%.3f,n=%d',r3e,p3e,size(corrmat,1)),'HorizontalAlignment','right','VerticalAlignment','top');
%     phi=plot([0,3],[0,3],'--r','LineWidth',1);
    ylabel('TCOM in 6s delay trials late half (s)');
    xlabel('TCOM in 3s delay trials (s)');
    legend([pha],{'Time-bin median'},'Location','northoutside')
    keyboard()
    exportgraphics(fh,'TCOM_corr_3s_6late3.pdf','ContentType','vector');
end

%% 3 and 6
if opt.plot_3_6full
    if isempty(com_str_6) || isempty(com_str_3)
        com_str_6=wave.get_com_map('delay',6);
        com_str_3=wave.get_com_map('delay',3);
    end
    allsess=intersect(fieldnames(com_str_3),fieldnames(com_str_6));
    corrmat=[];
    for sess=reshape(allsess,1,[])
        for pref=["s1","s2"]
            allcid=intersect(cell2mat(com_str_6.(sess{1}).(pref).keys),cell2mat(com_str_3.(sess{1}).(pref).keys));
            for onecid=allcid(:).'
                corrmat=[corrmat;...
                    com_str_6.(sess{1}).(pref)(onecid),com_str_3.(sess{1}).(pref)(onecid)]; %[6,3]
            end
        end
    end
    [r36,p36]=corr(corrmat(:,1),corrmat(:,2));
    % rci36=bootci(500,@(x) corr(x(:,1),x(:,2)),corrmat);

    binavg=[];
    [~,~,kbin]=histcounts(corrmat(:,2),0:2:12);
    for bb=reshape(unique(kbin),1,[])
        ci=bootci(500,@(x) median(x),corrmat(kbin==bb,1));% ci in 6s data
        binavg=[binavg;median(corrmat(kbin==bb,1)),median(corrmat(kbin==bb,2)),ci(1),ci(2)];
    end

    fh=figure('Color','w','Position',[100,100,460,230]);
    hold on
    pha=plot(binavg(:,1)./4,binavg(:,2)./4,'-k','LineWidth',1);
    errorbar(binavg(:,1)./4,binavg(:,2)./4,diff(binavg(:,[1 3]),1,2)./4,diff(binavg(:,[1 4]),1,2)./4,'k.','LineWidth',1)
    scatter(corrmat(:,1)./4,corrmat(:,2)./4,4,'k','filled','o','MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
    phi=plot([0,3],[0,3],'--r','LineWidth',1);
    phs=plot([0,6],[0,3],':r','LineWidth',1);
    text(max(xlim()),max(ylim()+0.5),sprintf('r=%.3f,p=%.3f, n=%d',r36,p36,size(corrmat,1)),'HorizontalAlignment','right','VerticalAlignment','top');
    xlabel('TCOM in 6s delay trials (s)');
    ylabel('TCOM in 3s delay trials (s)');
    legend([pha,phi,phs],{'data-median','y = x','y = 1/2x'},'Location','northoutside')
    keyboard()
    exportgraphics(fh,'TCOM_corr_3s_6s.pdf','ContentType','vector');
end
%% 6 and 6early3
if opt.plot_6_6early
    if isempty(com_str_early3) || isempty(com_str_6)
        com_str_early3=wave.get_com_map('partial','early3in6');
        com_str_6=wave.get_com_map('delay',6);
    end
    allsess=intersect(fieldnames(com_str_early3),fieldnames(com_str_6));
    corrmat=[];
    for sess=reshape(allsess,1,[])
        for pref=["s1","s2"]
            allcid=intersect(cell2mat(com_str_6.(sess{1}).(pref).keys),cell2mat(com_str_early3.(sess{1}).(pref).keys));
            for onecid=allcid(:).'
                corrmat=[corrmat;...
                    com_str_6.(sess{1}).(pref)(onecid),com_str_early3.(sess{1}).(pref)(onecid)];
            end
        end
    end

    [r6e,p6e]=corr(corrmat(:,1),corrmat(:,2));
    fh=figure('Color','w','Position',[32,32,600,230]);
    hold on
    scatter(corrmat(:,1)./4,corrmat(:,2)./4,4,'red','filled','o','MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
    text(max(xlim()),max(ylim()),sprintf('r=%.3f,p=%.3f',r6e,p6e),'HorizontalAlignment','right','VerticalAlignment','top');
    plot([0,3],[0,3],'--k')
end


%% 6 early 3 and late 3
if opt.plot_6early_6late
    if isempty(com_str_early3) || isempty(com_str_late3)
        com_str_early3=wave.get_com_map('partial','early3in6');
        com_str_late3=wave.get_com_map('partial','late3in6');
    end
    allsess=intersect(fieldnames(com_str_early3),fieldnames(com_str_late3));
    corrmat=[];
    for sess=reshape(allsess,1,[])
        for pref=["s1","s2"]
            allcid=intersect(cell2mat(com_str_late3.(sess{1}).(pref).keys),cell2mat(com_str_early3.(sess{1}).(pref).keys));
            for onecid=allcid(:).'
                corrmat=[corrmat;...
                    com_str_early3.(sess{1}).(pref)(onecid),com_str_late3.(sess{1}).(pref)(onecid)];%[3,6]
            end
        end
    end

    [r6el,p6el]=corr(corrmat(:,1),corrmat(:,2));

    binavg=[];
    [~,~,kbin]=histcounts(corrmat(:,1),0:2:12);
    for bb=reshape(unique(kbin),1,[])
        ci=bootci(500,@(x) median(x),corrmat(kbin==bb,2));% ci of 6late3
        binavg=[binavg;median(corrmat(kbin==bb,1)),median(corrmat(kbin==bb,2)),ci(1),ci(2)];%[6early3,3]
    end

    fh=figure('Color','w','Position',[32,32,230,230]);
    hold on
    scatter(corrmat(:,1)./4,corrmat(:,2)./4,4,'k','filled','o','MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
    ph=plot(binavg(:,1)./4,binavg(:,2)./4,'-k','LineWidth',1);
    errorbar(binavg(:,1)./4,binavg(:,2)./4,diff(binavg(:,[2 3]),1,2)./4,diff(binavg(:,[2 4]),1,2)./4,'k.','LineWidth',1)
    text(max(xlim()),max(ylim()),sprintf('r=%.3f,p=%.3f',r6el,p6el),'HorizontalAlignment','right','VerticalAlignment','top');
%     plot([0,3],[0,3],'--k')
    ylabel('Late-3s TCOM of 6s delay');
    xlabel('Early-3s TCOM of 6s delay');
    keyboard()
    exportgraphics(fh,'TCOM_corr_6searly_6slate.pdf','ContentType','vector');
end
end
function misc_func()
%%
counter=0;
for ss=reshape(fieldnames(com_str_3),1,[])
    counter=counter+numel(com_str_3.(ss{1}).s1.keys);
    counter=counter+numel(com_str_3.(ss{1}).s2.keys);
end
%%
fh=figure('Color','w','Position',[100,100,120,150]);
hold on
bar([r36,r3e],'w')
errorbar(1:2,[r36,r3e],[rci36(1)-r36,rci3e(1)-r3e],[rci36(2)-r36,rci3e(2)-r3e],'k.','CapSize',15)
xlim([0.5,2.5])
set(gca,'XTick',1:2,'XTickLabel',{'In 6s delay','In early 3s'},'XTickLabelRotation',90)
ylabel('Pearson''s r with 3s delay' )
exportgraphics(fh,'corr_3s_6s.pdf')
end