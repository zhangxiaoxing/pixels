function [fh,t]=per_region_COM_frac(fcom,ffrac,opt)
arguments
    fcom
    ffrac
    opt.frac_COM (1,1) logical = true
    opt.frac_PVSST (1,1) logical = true
    opt.COM_PVSST (1,1) logical = true
    opt.COM_hieridx (1,1) logical = true
    opt.range (1,:) char {mustBeMember(opt.range,{'CH','CTX','grey'})} = 'grey'
    opt.corr (1,:) char {mustBeMember(opt.corr,{'Pearson','Spearman'})} = 'Spearman'
    opt.export (1,1) logical = false
    opt.selidx (1,1) logical = false % calculate COM of selectivity index
    opt.delay (1,1) double {mustBeMember(opt.delay,[3,6])} = 6 % DPA delay duration
    opt.hier_reg (1,:) char = 'AON'
    opt.density_scale (1,:) char = 'log'
end


%data generated from wave.per_region_COM
%data of interest, region,branch level, count

fh=figure('Color','w','Position',[100,100,600,600]);
t=tiledlayout('flow');

%% hiermap
idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
sink_ccfid=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/grey_targets');
src_ccfid=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/grey_srcs');
sink_src_mat=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/src_target_matrix');
anov=sink_src_mat(:,src_ccfid==idmap.reg2ccfid(opt.hier_reg));
hiermap=containers.Map(cellfun(@(x) x,idmap.ccfid2reg.values(num2cell(sink_ccfid))),...
    anov);
%%
[~,~,pvsst_map]=ref.get_pv_sst();
reg_range=ephys.getGreyRegs('range',opt.range);


%% com vs frac
if opt.frac_COM
%     fh.frac_COM=figure('Color','w','Position',[100,100,245,235]);
    nexttile
    hold on;
    coord=[];
    regs=[];
    inter_reg=fcom.collection(:,2); %TODO refactor call stack
    for ri=1:numel(inter_reg)
        fridx=find(strcmp(ffrac.collection(:,2),inter_reg(ri)));
        if ~isempty(fridx)
            comidx=find(strcmp(fcom.collection(:,2),inter_reg(ri)));
            yy=fcom.collection{comidx,1}./4;
            xx=ffrac.collection{fridx,1};
            coord=[coord;xx,yy];
            regs=[regs,fcom.collection{comidx,2}];
            %Notice, FRP still miss reg-color-group
            scatter(xx,yy,9,'o','MarkerFaceColor',ephys.getRegColor(inter_reg{ri},'large_area',true),'MarkerEdgeColor','none');
            text(xx,yy,fcom.collection{comidx,2},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',7,'Color',ephys.getRegColor(inter_reg{ri},'large_area',true));
        end
    end
    coord(:,3)=1;
    regres=coord(:,[1,3])\coord(:,2);
    plot(xlim(),xlim().*regres(1)+regres(2),'--k');
    xlabel('Regional proportion of sensory neuron')
    ylabel('Center of FR modulation (s)')
    if strcmp(opt.corr,'Pearson') && strcmp(opt.density_scale,'log')
        vsel=coord(:,1)>0;
        [r,p]=corr(log10(coord(vsel,1)),coord(vsel,2),'type',opt.corr);
    else
        [r,p]=corr(coord(:,1),coord(:,2),'type',opt.corr);
    end
    set(gca(),'XScale',opt.density_scale);
    text(max(xlim()),max(ylim()),sprintf('r = %.3f, p = %.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');
    if opt.export
        exportgraphics(fh,sprintf('per_region_TCOM_FRAC_%d.pdf',opt.delay));
    end

end
%% com vs pv/sst
if opt.COM_PVSST
%     fh.COM_PVSST=figure('Color','w','Position',[100,100,245,235]);
    nexttile
    hold on;
    coord=[];
    regs=[];
    ureg=intersect(fcom.collection(:,2),ffrac.collection(:,2));
    for ri=1:numel(ureg)
        if  pvsst_map.isKey(ureg{ri})
            comidx=find(strcmp(fcom.collection(:,2),ureg(ri)));
            yy=fcom.collection{comidx,1}./4;
            xx=pvsst_map(ureg{ri});
            coord=[coord;xx,yy];
            regs=[regs,fcom.collection{comidx,2}];
            scatter(xx,yy,9,'o','MarkerFaceColor',ephys.getRegColor(ureg{ri},'large_area',true),'MarkerEdgeColor','none');
            text(xx,yy,fcom.collection{comidx,2},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',7,'Color',ephys.getRegColor(ureg{ri},'large_area',true));
        else
            warning('Missing PVSST map key');
            disp(ureg{ri})
        end
    end
%     ylim([0,0.6])
    coord(:,3)=1;
    regres=coord(:,[1,3])\coord(:,2);
%     if opt.delay==6
%         xlim([2.4,3.5]);
%         set(gca(),'XTick',2.5:0.5:3.5);
%     else
%         xlim([1.4,2.0]);
%         set(gca(),'XTick',1.2:0.2:2.0);
%     end
%     plot(xlim(),xlim().*regres(1)+regres(2),'--k');
    
    xlabel('PV/(PV+SST) interneuron ratio')
    ylabel('Center of FR modulation')
    [r,p]=corr(coord(:,1),coord(:,2),'type',opt.corr);
    text(max(xlim()),max(ylim()),sprintf('r = %.3f, p = %.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');
    if opt.export
        exportgraphics(fh,sprintf('per_region_COM_pv_sst_%d.pdf',opt.delay));
    end


end

%% com vs SMI
if opt.COM_hieridx
%     load('OBM1Map.mat','OBM1map')
%     fh.COM_hieridx=figure('Color','w','Position',[100,100,245,235]);
    nexttile
    hold on;
    coord=[];
    regs=[];
    inter_reg=intersect(fcom.collection(:,2),ffrac.collection(:,2));

    for ri=1:numel(inter_reg)
        if  hiermap.isKey(inter_reg{ri})
            comidx=find(strcmp(fcom.collection(:,2),inter_reg(ri)));
            xx=hiermap(inter_reg{ri});
            yy=fcom.collection{comidx,1}./4;
            coord=[coord;xx,yy];
            regs=[regs,fcom.collection{comidx,2}];
            scatter(xx,yy,9,'o','MarkerFaceColor',ephys.getRegColor(inter_reg{ri},'large_area',true),'MarkerEdgeColor','none');
            text(xx,yy,fcom.collection{comidx,2},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',7,'Color',ephys.getRegColor(inter_reg{ri},'large_area',true));
        else
            warning('Missing PVSST map key');
            disp(inter_reg{ri})
        end
    end
    
    xlabel(['Projection density from ',opt.hier_reg])
    ylabel('Center of FR modulation')

    coord(:,3)=1;
    regres=coord(:,[1,3])\coord(:,2);
%     if opt.delay==6
%         ylim([2.4,3.3]);
%         set(gca(),'YTick',2.5:0.5:3.5);
%     else
%         ylim([1.4,2.0]);
%         set(gca(),'YTick',1.2:0.2:2.0);
%     end
%     plot(xlim(),xlim().*regres(1)+regres(2),'--k');

    set(gca(),'XScale','log')
%     xlim([-7,7])
    if strcmp(opt.corr,'Pearson')
        vsel=coord(:,1)>0;
        [r,p]=corr(log10(coord(vsel,1)),coord(vsel,2),'type',opt.corr);
    else
        [r,p]=corr(coord(:,1),coord(:,2),'type',opt.corr);
    end
    text(max(xlim()),max(ylim()),sprintf('r = %.3f, p = %.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');
    if opt.export
        exportgraphics(fh,sprintf('per_region_COM_SMI_%d.pdf',opt.delay));
    end


end



%% frac vs pv/sst
if opt.frac_PVSST
    ureg=intersect(fcom.collection(:,2),ffrac.collection(:,2));
    
%     fh.frac_PVSST=figure('Color','w','Position',[100,100,245,235]);
    nexttile
    hold on;
    coord=[];
    regs=[];
    for ri=1:numel(ureg)
        fridx=find(strcmp(ffrac.collection(:,2),ureg(ri)));
        if ~isempty(fridx) && pvsst_map.isKey(ureg{ri})
            yy=ffrac.collection{fridx,1}.*100;
            xx=pvsst_map(ureg{ri});
            coord=[coord;xx,yy];
            regs=[regs,ffrac.collection{fridx,2}];
            scatter(xx,yy,9,'o','MarkerFaceColor',ephys.getRegColor(ureg{ri},'large_area',true),'MarkerEdgeColor','none');
            text(xx,yy,ffrac.collection{fridx,2},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',7,'Color',ephys.getRegColor(ureg{ri},'large_area',true));
        end
    end
    coord(:,3)=1;
    regres=coord(:,[1,3])\coord(:,2);
    plot(xlim(),xlim().*regres(1)+regres(2),'--k');
%     if strcmp(opt.sust_type,'sust')
%         ylim([0,4]);
%         set(gca(),'YTick',0:2:4,'XTick',0:0.2:1);
%     else
%         ylim([10,60]);
%         set(gca(),'YTick',0:20:60,'XTick',0:0.2:1);
%     end
    set(gca(),'YScale','log');
    xlabel('PV / (PV + SST) density')
    ylabel('Regional proportion of sensory neuron')
%     xlim([0,0.6]);
    if strcmp(opt.corr,'Pearson')
        vsel=coord(:,2)>0;
        [r,p]=corr(coord(vsel,1),log10(coord(vsel,2)),'type',opt.corr);
    else
        [r,p]=corr(coord(:,1),coord(:,2),'type',opt.corr);
    end
    text(max(xlim()),max(ylim()),sprintf('r = %.3f, p = %.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','bottom');
    if opt.export
        keyboard()
        exportgraphics(fh,sprintf('per_region_frac_pv_sst_%d.pdf',opt.delay));
    end
end

end