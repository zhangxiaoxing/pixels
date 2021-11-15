function [fcom,ffrac]=per_region_COM_frac(opt)
arguments
    opt.frac_COM (1,1) logical = false
    opt.frac_PVSST (1,1) logical = true
    opt.COM_PVSST (1,1) logical = false
    opt.frac_sensemotor (1,1) logical = false
    opt.COM_sensemotor (1,1) logical = false
    opt.sust_type (1,:) char {mustBeMember(opt.sust_type,{'any','sust','trans'})} = 'any'
    opt.extent (1,:) char {mustBeMember(opt.extent,{'CH','CTX','CTXpl'})} = 'CTX'
    opt.corr (1,:) char {mustBeMember(opt.corr,{'Pearson','Spearman'})} = 'Pearson'
    opt.export (1,1) logical = true
    opt.selidx (1,1) logical = false % calculate COM of selectivity index
    opt.delay (1,1) double {mustBeMember(opt.delay,[3,6])} = 6 % DPA delay duration
end


%data generated from wave.per_region_COM
%data of interest, region,branch level, count
[fcom.collection,fcom.com_meta]=wave.per_region_COM('stats_method','mean','selidx',opt.selidx,'delay',opt.delay);
ffrac.collection=ephys.per_region_fraction('memtype',opt.sust_type,'delay',opt.delay);
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
            yy=fcom.collection{comidx,1}./4;
            xx=ffrac.collection{fridx,1}.*100;
            coord=[coord;xx,yy];
            regs=[regs,fcom.collection{comidx,2}];
            %Notice, FRP still miss reg-color-group
            scatter(xx,yy,9,'o','MarkerFaceColor',ephys.getRegColor(ureg{ri}),'MarkerEdgeColor','none');
            text(xx,yy,fcom.collection{comidx,2},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',7,'Color',ephys.getRegColor(ureg{ri}));
        end
    end
    coord(:,3)=1;
    regres=coord(:,[1,3])\coord(:,2);
    if strcmp(opt.sust_type,'sust')
        xlim([0,4]);
        set(gca(),'XTick',0:2:4);
    else
        xlim([10,60]);
        set(gca(),'XTick',0:20:60);
    end
    if opt.delay==6
        ylim([2.4,3.5]);
        set(gca(),'YTick',2.5:0.5:3.5);
    else
        ylim([1.4,2.0]);
        set(gca(),'YTick',1.2:0.2:2.0);
    end
    plot(xlim(),xlim().*regres(1)+regres(2),'--k');
    xlabel('Fraction of delay selective neuron')
    ylabel('F.R. center of mass (s)')
    [r,p]=corr(coord(:,1),coord(:,2),'type',opt.corr);
    text(max(xlim()),max(ylim()),sprintf('r = %.3f, p = %.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');
    if opt.export
        exportgraphics(fh,sprintf('per_region_TCOM_FRAC_%d.pdf',opt.delay));
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
    coord(:,3)=1;
    regres=coord(:,[1,3])\coord(:,2);
    if opt.delay==6
        xlim([2.4,3.5]);
        set(gca(),'XTick',2.5:0.5:3.5);
    else
        xlim([1.4,2.0]);
        set(gca(),'XTick',1.2:0.2:2.0);
    end
    plot(xlim(),xlim().*regres(1)+regres(2),'--k');
    
    ylabel('Hierarchy index (Low->High)')
    xlabel('F.R. center of mass (s)')
    [r,p]=corr(coord(:,1),coord(:,2),'type',opt.corr);
    text(max(xlim()),max(ylim()),sprintf('r = %.3f, p = %.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');
    if opt.export
        exportgraphics(fh,sprintf('per_region_COM_pv_sst_%d.pdf',opt.delay));
    end


end

%% com vs SMI
if opt.COM_sensemotor
    load('OBM1Map.mat','OBM1map')
    fh=figure('Color','w','Position',[100,100,245,235]);
    hold on;
    coord=[];
    regs=[];
    
    for ri=1:numel(ureg)
        if  OBM1map.isKey(ureg{ri})
            comidx=find(strcmp(fcom.collection(:,2),ureg(ri)));
            xx=OBM1map(ureg{ri});
            yy=fcom.collection{comidx,1}./4;
            coord=[coord;xx,yy];
            regs=[regs,fcom.collection{comidx,2}];
            scatter(xx,yy,9,'o','MarkerFaceColor',ephys.getRegColor(ureg{ri}),'MarkerEdgeColor','none');
            text(xx,yy,fcom.collection{comidx,2},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',7,'Color',ephys.getRegColor(ureg{ri}));
        else
            warning('Missing PVSST map key');
            disp(ureg{ri})
        end
    end
    
    xlabel('Sensory Motor Index (A.U.)')
    ylabel('F.R. center of mass (s)')

    coord(:,3)=1;
    regres=coord(:,[1,3])\coord(:,2);
    if opt.delay==6
        ylim([2.4,3.3]);
        set(gca(),'YTick',2.5:0.5:3.5);
    else
        ylim([1.4,2.0]);
        set(gca(),'YTick',1.2:0.2:2.0);
    end
    plot(xlim(),xlim().*regres(1)+regres(2),'--k');
    xlim([-7,7])
    [r,p]=corr(coord(:,1),coord(:,2),'type',opt.corr);
    text(max(xlim()),max(ylim()),sprintf('r = %.3f, p = %.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');
    if opt.export
        exportgraphics(fh,sprintf('per_region_COM_SMI_%d.pdf',opt.delay));
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
            yy=ffrac.collection{fridx,1}.*100;
            xx=ratiomap(ureg{ri});
            coord=[coord;xx,yy];
            regs=[regs,ffrac.collection{fridx,2}];
            scatter(xx,yy,9,'o','MarkerFaceColor',ephys.getRegColor(ureg{ri}),'MarkerEdgeColor','none');
            text(xx,yy,ffrac.collection{fridx,2},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',7,'Color',ephys.getRegColor(ureg{ri}));
        end
    end
    coord(:,3)=1;
    regres=coord(:,[1,3])\coord(:,2);
    plot(xlim(),xlim().*regres(1)+regres(2),'--k');
    if strcmp(opt.sust_type,'sust')
        ylim([0,4]);
        set(gca(),'YTick',0:2:4,'XTick',0:0.2:1);
    else
        ylim([10,60]);
        set(gca(),'YTick',0:20:60,'XTick',0:0.2:1);
    end
    xlabel('PV / (PV + SST) density')
    ylabel('Proportion of memory selective neuron')
    xlim([0,0.6]);
    [r,p]=corr(coord(:,1),coord(:,2),'type',opt.corr);
    text(max(xlim()),max(ylim()),sprintf('r = %.3f, p = %.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','bottom');
    if opt.export
        keyboard()
        exportgraphics(fh,sprintf('per_region_frac_pv_sst_%d.pdf',opt.delay));
    end
end
%% frac vs sensory motor transfer index

if opt.frac_sensemotor
    load('OBM1Map.mat','OBM1map')
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
        if ~isempty(fridx) && OBM1map.isKey(ureg{ri})
            xx=ffrac.collection{fridx,1}.*100;
            yy=OBM1map(ureg{ri});
            coord=[coord;xx,yy];
            regs=[regs,ffrac.collection{fridx,2}];
            scatter(xx,yy,9,'o','MarkerFaceColor',ephys.getRegColor(ureg{ri}),'MarkerEdgeColor','none');
            text(xx,yy,ffrac.collection{fridx,2},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',7,'Color',ephys.getRegColor(ureg{ri}));
        end
    end
    coord(:,3)=1;
    regres=coord(:,[1,3])\coord(:,2);
    plot(xlim(),xlim().*regres(1)+regres(2),'--k');
    ylim([-7,7])
    if strcmp(opt.sust_type,'sust')
        xlim([0,4]);
        set(gca(),'XTick',0:2:4,'YTick',0:0.2:1);
    else
        xlim([10,60]);
        set(gca(),'XTick',0:20:60,'YTick',-6:3:6);
    end
    ylabel('Sensory-motor transfer index')
    xlabel('Selective fraction')
    [r,p]=corr(coord(:,1),coord(:,2),'type',opt.corr);
    text(max(xlim()),max(ylim()),sprintf('r = %.3f, p = %.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','bottom');
    if opt.export
        exportgraphics(fh,sprintf('per_region_frac_sense_motor_%d.pdf',opt.delay));
    end
end
end