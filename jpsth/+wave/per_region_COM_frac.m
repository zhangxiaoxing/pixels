function [fcom,ffrac]=per_region_COM_frac(opt)
arguments
    opt.frac_COM (1,1) logical = false
    opt.frac_PVSST (1,1) logical = false
    opt.COM_PVSST (1,1) logical = false
    opt.sust_type (1,:) char {mustBeMember(opt.sust_type,{'any','sust','trans'})} = 'any'
    opt.extent (1,:) char {mustBeMember(opt.extent,{'CH','CTX','CTXpl'})} = 'CH'
    opt.corr (1,:) char {mustBeMember(opt.corr,{'Pearson','Spearman'})} = 'Pearson'
    opt.export (1,1) logical = false
end


%data generated from wave.per_region_COM
%data of interest, region,branch level, count
[fcom.collection,fcom.com_meta]=wave.per_region_COM('stats_method','mean');
ffrac.collection=ephys.per_region_fraction('memtype',opt.sust_type);
[~,~,ratiomap]=ref.get_pv_sst();
idmap=load(fullfile('..','align','reg_ccfid_map.mat'));

uregall=fcom.collection(cell2mat(fcom.collection(:,4))>20 & cell2mat(fcom.collection(:,3))==5,2);
ureg=[];
for ri=1:numel(uregall)
    if idmap.reg2tree.isKey(uregall{ri}) && any(ismember(idmap.reg2tree(uregall{ri}),{opt.extent}))
        ureg=[ureg;uregall(ri)];
    end
end

%% com vs frac
if opt.frac_COM
    fh=figure('Color','w','Position',[100,100,245,235]);
    hold on;
    coord=[];
    regs=[];
    
    for ri=1:numel(ureg)
        fridx=find(strcmp(ffrac.collection(:,2),ureg(ri)));
        if ~isempty(fridx) && ffrac.collection{fridx,4}>40
            comidx=find(strcmp(fcom.collection(:,2),ureg(ri)));
            xx=fcom.collection{comidx,1}./4;
            yy=ffrac.collection{fridx,1}.*100;
            coord=[coord;xx,yy];
            regs=[regs,fcom.collection{comidx,2}];
            %Notice, FRP still miss reg-color-group
            scatter(xx,yy,9,'o','MarkerFaceColor',ephys.getRegColor(ureg{ri}),'MarkerEdgeColor','none');
            text(xx,yy,fcom.collection{comidx,2},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',7,'Color',ephys.getRegColor(ureg{ri}));
        end
    end
    if strcmp(opt.sust_type,'sust')
        ylim([0,4]);
        set(gca(),'YTick',0:2:4);
    else
        ylim([10,60]);
        set(gca(),'YTick',0:20:60);
    end
    xlim([2.4,3.5]);
    set(gca(),'XTick',2.5:0.5:3.5);
    ylabel('Fraction of delay selective neuron')
    xlabel('F.R. center of mass (s)')
    [r,p]=corr(coord(:,1),coord(:,2),'type',opt.corr);
    text(max(xlim()),max(ylim()),sprintf('r = %.3f, p = %.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');
    if opt.export
        exportgraphics(fh,'per_region_FC_COM_PVSST.pdf');
    end

end
%% com vs pv/sst
if opt.COM_PVSST
    fh=figure('Color','w','Position',[100,100,245,235]);
    hold on;
    coord=[];
    regs=[];
    
    for ri=1:numel(ureg)
        if  ratiomap.isKey(ureg{ri})
            comidx=find(strcmp(fcom.collection(:,2),ureg(ri)));
            xx=fcom.collection{comidx,1}./4;
            yy=ratiomap(ureg{ri});
            coord=[coord;xx,yy];
            regs=[regs,fcom.collection{comidx,2}];
            scatter(xx,yy,9,'o','MarkerFaceColor',ephys.getRegColor(ureg{ri}),'MarkerEdgeColor','none');
            text(xx,yy,fcom.collection{comidx,2},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',7,'Color',ephys.getRegColor(ureg{ri}));
        else
            warning('Missing PVSST map key');
            disp(ureg{ri})
        end
    end
    ylim([0,0.6])
    xlim([2.4,3.5]);
    set(gca(),'XTick',2.5:0.5:3.5);
    ylabel('Hierarchy index (Low->High)')
    xlabel('F.R. center of mass (s)')
    [r,p]=corr(coord(:,1),coord(:,2),'type',opt.corr);
    text(max(xlim()),max(ylim()),sprintf('r = %.3f, p = %.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');
    if opt.export
        exportgraphics(fh,'per_region_COM_pv_sst.pdf');
    end


end

%% frac vs pv/sst
if opt.frac_PVSST
    uregall=ffrac.collection(cell2mat(ffrac.collection(:,4))>40 & cell2mat(ffrac.collection(:,3))==5,2);
    ureg=[];
    for ri=1:numel(uregall)
        if idmap.reg2tree.isKey(uregall{ri}) && any(ismember(idmap.reg2tree(uregall{ri}),{opt.extent}))
            ureg=[ureg;uregall(ri)];
        end
    end
    
    fh=figure('Color','w','Position',[100,100,245,235]);
    hold on;
    coord=[];
    regs=[];
    for ri=1:numel(ureg)
        fridx=find(strcmp(ffrac.collection(:,2),ureg(ri)));
        if ~isempty(fridx) && ratiomap.isKey(ureg{ri})
            xx=ffrac.collection{fridx,1}.*100;
            yy=ratiomap(ureg{ri});
            coord=[coord;xx,yy];
            regs=[regs,ffrac.collection{fridx,2}];
            scatter(xx,yy,9,'o','MarkerFaceColor',ephys.getRegColor(ureg{ri}),'MarkerEdgeColor','none');
            text(xx,yy,ffrac.collection{fridx,2},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',7,'Color',ephys.getRegColor(ureg{ri}));
        end
    end
    ylim([0,1])
    if strcmp(opt.sust_type,'sust')
        xlim([0,4]);
        set(gca(),'XTick',0:2:4,'YTick',0:0.2:1);
    else
        xlim([0,60]);
        set(gca(),'XTick',0:20:60,'YTick',0:0.2:1);
    end
    ylabel('Hierarchy index (Low->High)')
    xlabel('Selective fraction')
    ylim([0,0.6]);
    [r,p]=corr(coord(:,1),coord(:,2),'type',opt.corr);
    text(max(xlim()),max(ylim()),sprintf('r = %.3f, p = %.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','bottom');
    if opt.export
        exportgraphics(fh,'per_region_frac_pv_sst.pdf');
    end
end
end