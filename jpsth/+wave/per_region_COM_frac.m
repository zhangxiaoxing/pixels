function per_region_COM_frac(opt)
arguments
    opt.frac_COM (1,1) logical = false
    opt.frac_PVSST (1,1) logical = false
    opt.COM_PVSST (1,1) logical = false
    opt.memtype (1,:) char {mustBeMember(opt.memtype,{'any','sust','trans'})} = 'any'
end


%data generated from wave.per_region_COM
fcom=load('per_region_com_collection.mat','collection');
ffrac.collection=ephys.per_region_fraction('memtype',opt.memtype);
[~,~,ratiomap]=ref.get_pv_sst();
idmap=load(fullfile('..','align','reg_ccfid_map.mat'));

uregall=fcom.collection(cell2mat(fcom.collection(:,4))>20 & cell2mat(fcom.collection(:,3))==5,2);
ureg=[];
for ri=1:numel(uregall)
    if idmap.reg2tree.isKey(uregall{ri}) && any(ismember(idmap.reg2tree(uregall{ri}),{'CH'}))
        ureg=[ureg;uregall(ri)];
    end
end

%% com vs frac
if opt.frac_COM
    fh=figure('Color','w');
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
            plot(xx,yy,'wo','MarkerFaceColor',[255,127,127]./255);
            text(xx,yy,fcom.collection{comidx,2},'HorizontalAlignment','center','VerticalAlignment','middle');
        end
    end
    if strcmp(opt.memtype,'sust')
        ylim([0,4]);
        set(gca(),'YTick',0:2:4);
    else
        ylim([0,60]);
        set(gca(),'YTick',0:20:60);
    end
    xlim([2.4,3.5]);
    set(gca(),'XTick',2.5:0.5:3.5);
    ylabel('Fraction of delay selective neuron')
    xlabel('F.R. center of mass (s)')
    [r,p]=corr(coord(:,1),coord(:,2));
    text(max(xlim()),max(ylim()),sprintf('r = %.3f, p = %.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');
end
%% com vs pv/sst
if opt.COM_PVSST
    fh=figure('Color','w');
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
            plot(xx,yy,'wo','MarkerFaceColor',[255,127,127]./255);
            text(xx,yy,fcom.collection{comidx,2},'HorizontalAlignment','center','VerticalAlignment','middle');
        end
    end
    ylim([0.0,1])
    xlim([2.4,3.5]);
    set(gca(),'XTick',2.5:0.5:3.5);
    ylabel('PV/PV+SST density')
    xlabel('F.R. center of mass (s)')
    [r,p]=corr(coord(:,1),coord(:,2),'type','Pearson');
    text(max(xlim()),max(ylim()),sprintf('r = %.3f, p = %.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');
end

%% frac vs pv/sst
if opt.frac_PVSST
    uregall=ffrac.collection(cell2mat(ffrac.collection(:,4))>40 & cell2mat(ffrac.collection(:,3))==5,2);
    ureg=[];
    for ri=1:numel(uregall)
        if idmap.reg2tree.isKey(uregall{ri}) && any(ismember(idmap.reg2tree(uregall{ri}),{'CH'}))
            ureg=[ureg;uregall(ri)];
        end
    end
    
    fh=figure('Color','w');
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
            plot(xx,yy,'wo','MarkerFaceColor',[255,127,127]./255);
            text(xx,yy,ffrac.collection{fridx,2},'HorizontalAlignment','center','VerticalAlignment','middle');
        end
    end
    ylim([0,1])
    if strcmp(opt.memtype,'sust')
        xlim([0,4]);
        set(gca(),'XTick',0:2:4,'YTick',0:0.2:1);
    else
        xlim([0,60]);
        set(gca(),'XTick',0:20:60,'YTick',0:0.2:1);
    end
    ylabel('PV/PV+SST density')
    xlabel('Selective fraction')
    [r,p]=corr(coord(:,1),coord(:,2),'type','Pearson');
    text(max(xlim()),min(ylim()),sprintf('r = %.3f, p = %.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','bottom');
end
end