function [map_cells,fh]=pct_reg_bars(su_meta,sel_meta,opt) %TODO refactoring function name
arguments
    su_meta
    sel_meta
    opt.skip_plot (1,1) logical = false % return map_cell data for GLM etc.
    % opt.range (1,:) char {mustBeMember(opt.range,{'grey','CH','CTX'})} = 'grey'
    % opt.stats_type (1,:) char = "percentile"
    % opt.data_type (1,:) char = "mixed-olf-dur"
    opt.skip_error_bar (1,1) logical = true
    opt.xyscale (1,2) cell = {'Linear','Linear'}
    opt.skip_export (1,1) logical = true
    opt.skip_dur (1,1) logical = false
    opt.only_odor (1,1) logical = true
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','Naive','any'})} = 'WT'
    opt.min_reg_count (1,1) double =100
    opt.area_set (1,1) logical = true
end

idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
if opt.area_set
    areatbl=readtable('../allensdk/mouse_areas_3.csv');
    ureg=areatbl.acronym;
else
    uregsel=ismember(su_meta.reg_tree(1,:),{'CH','BS'}) & ~ismissing(su_meta.reg_tree(5,:));
    ureg=unique(su_meta.reg_tree(5,uregsel));
end
sums=[];

for reg=reshape(ureg,1,[])

    regsel=any(strcmp(su_meta.reg_tree,reg),1).';
    cnt=nnz(regsel);
    if cnt<opt.min_reg_count
        continue
    end
    mixed_cnt=nnz(regsel & ismember(sel_meta.wave_id,1:4));
    olf_cnt=nnz(regsel & ismember(sel_meta.wave_id,5:6));
    dur_cnt=nnz(regsel & ismember(sel_meta.wave_id,7:8));
    legends={'Mixed modality','Olfactory only','Duration only'};

    grp=idmap.reg2tree(reg{1});
    sums=[sums;idmap.reg2ccfid(grp{end-1}),idmap.reg2ccfid(reg{1}),cnt,mixed_cnt,olf_cnt,dur_cnt];
    %=================1LVL6======================2LVL7==========3=========4======5======6=
end

sums(:,7:15)=0;

for ii=1:size(sums,1)
    if opt.only_odor
        [olf_hat,olf_ci]=binofit(sums(ii,5)+sums(ii,4),sums(ii,3));
        sums(ii,10:12)=[olf_hat,olf_ci];
    else
        [mixed_hat,mixed_ci]=binofit(sums(ii,4),sums(ii,3));
        [olf_hat,olf_ci]=binofit(sums(ii,5),sums(ii,3));
        [dur_hat,dur_ci]=binofit(sums(ii,6),sums(ii,3));
        sums(ii,7:15)=[mixed_hat,mixed_ci,olf_hat,olf_ci,dur_hat,dur_ci];
        %=================7=======8|9======10======11|12====13=====14:15
    end
end

%map_cells

bardata=sortrows(sums,10,'descend');
regstr=idmap.ccfid2reg.values(num2cell(bardata(:,2)));
if ~opt.skip_export
    exp=input('export for 3d render? yes/no\n','s');
    if strcmpi(exp,'yes') % export for brain renderer
        for rr=1:size(bardata,1)
            disp("['"+string(idmap.ccfid2reg(bardata(rr,2))) ...
                +"',"...
                + num2str(bardata(rr,[10 13 7]),"%.3f,")...
                + "]");
        end
        keyboard()
    end
end
if opt.only_odor
    map_cells.olf=containers.Map([regstr{:}],num2cell([bardata(:,10),bardata(:,5)+bardata(:,4),bardata(:,3)],2));
else
    map_cells.mixed=containers.Map([regstr{:}],num2cell(bardata(:,[7,4,3]),2));
    map_cells.olf=containers.Map([regstr{:}],num2cell(bardata(:,[10,5,3]),2));
    map_cells.dur=containers.Map([regstr{:}],num2cell(bardata(:,[13,5,3]),2));
end

if opt.skip_plot
    fh=[];
    return
end

%===============================================================
fh=figure('Color','w','Position',[32,32,1280,640]);
th=tiledlayout(2,4);

if ~opt.only_odor
    nexttile(5);
    hold on
    for ll=1:size(bardata,1)
        c=ephys.getRegColor(regstr{ll},'large_area',true);
        scatter(bardata(ll,10),bardata(ll,7),4,c,'filled','o')
        text(bardata(ll,10),bardata(ll,7),regstr{ll},'HorizontalAlignment','center','VerticalAlignment','top','Color',c);
    end
    xlabel(legends{2})
    ylabel(legends{1})
    set(gca,'XScale',opt.xyscale{1},'YScale',opt.xyscale{2})
    xxdata=bardata(:,10);
    yydata=bardata(:,7);
    if strcmp(opt.xyscale{1},'log')
        xxdata=log10(xxdata);
    end
    if strcmp(opt.xyscale{2},'log')
        yydata=log10(yydata);
    end

    [r,p]=corr(xxdata,yydata);
    title(sprintf(' r = %.3f, p = %.3f',r,p));

    coord=[xxdata,yydata];
    coord(:,3)=1;
    regres=coord(:,[1,3])\coord(:,2);
    xx=minmax(xxdata.');

    if all(strcmp(opt.xyscale,'log'),'all')
        yy=10.^(xx.*regres(1)+regres(2));
        plot(10.^xx,yy,'--k');
    else
        yy=xx.*regres(1)+regres(2);
        plot(xx,yy,'--k');
    end
    %==============================================================
    nexttile(6);
    hold on

    for ll=1:size(bardata,1)
        c=ephys.getRegColor(regstr{ll},'large_area',true);
        scatter(bardata(ll,7),bardata(ll,13),4,c,'filled','o')
        text(bardata(ll,7),bardata(ll,13),regstr{ll},'HorizontalAlignment','center','VerticalAlignment','top','Color',c);
    end
    xlabel(legends{1})
    ylabel(legends{3})
    set(gca,'XScale',opt.xyscale{1},'YScale',opt.xyscale{2})
    xxdata=bardata(:,7);
    yydata=bardata(:,13);
    if strcmp(opt.xyscale{1},'log')
        xxdata=log10(xxdata);
    end
    if strcmp(opt.xyscale{2},'log')
        yydata=log10(yydata);
    end
    [r,p]=corr(xxdata,yydata);
    title(sprintf(' r = %.3f, p = %.3f',r,p));

    coord=[xxdata,yydata];
    coord(:,3)=1;
    regres=coord(:,[1,3])\coord(:,2);
    xx=minmax(xxdata.');
    if all(strcmp(opt.xyscale,'log'))
        yy=10.^(xx.*regres(1)+regres(2));
        plot(10.^xx,yy,'--k');
    else
        yy=xx.*regres(1)+regres(2);
        plot(xx,yy,'--k');
    end


    %==============================================================
    nexttile(7);
    hold on

    for ll=1:size(bardata,1)
        c=ephys.getRegColor(regstr{ll},'large_area',true);
        scatter(bardata(ll,10),bardata(ll,13),4,c,'filled','o')
        text(bardata(ll,10),bardata(ll,13),regstr{ll},'HorizontalAlignment','center','VerticalAlignment','top','Color',c);
    end
    xlabel(legends{2})
    ylabel(legends{3})
    set(gca,'XScale',opt.xyscale{1},'YScale',opt.xyscale{2})
    xxdata=bardata(:,10);
    yydata=bardata(:,13);
    if strcmp(opt.xyscale{1},'log')
        xxdata=log10(xxdata);
    end
    if strcmp(opt.xyscale{2},'log')
        yydata=log10(yydata);
    end
    [r,p]=corr(xxdata,yydata);
    title(sprintf(' r = %.3f, p = %.3f',r,p));
    coord=[xxdata,yydata];
    coord(:,3)=1;
    regres=coord(:,[1,3])\coord(:,2);
    xx=minmax(xxdata.');
    if all(strcmp(opt.xyscale,'log'))
        yy=10.^(xx.*regres(1)+regres(2));
        plot(10.^xx,yy,'--k');
    else
        yy=xx.*regres(1)+regres(2);
        plot(xx,yy,'--k');
    end


    %%======================
    th=nexttile(4);
    tbl=cell(0);
    % if isfield(opt,'stats_type') && ~isempty(opt.stats_type)
    %     tbl=[tbl;'Stats';opt.stats_type];
    % end
    % if isfield(opt,'data_type') && ~isempty(opt.data_type)
    %     tbl=[tbl;'Data';opt.data_type];
    % end
    % tbl=[tbl;'Range';opt.range];
    ephys.util.figtable(fh,th,tbl)
else
    th=nexttile(4);
    tbl={'Odor only';'grey'};
    ephys.util.figtable(fh,th,tbl)
end

%======================
nexttile(1,[1,3]);
hold on
if opt.only_odor
    bh=bar(bardata(:,10),0.75);
    if ~opt.skip_error_bar
        errorbar(bh(1).XEndPoints,bh(1).YEndPoints,diff(bardata(:,10:11),1,2),diff(bardata(:,[10,12]),1,2),'k.');
    end
    bh(1).FaceColor='k';% olf
    fid=fopen(fullfile('binary','upload','WM_neuron_per_region.json'),'w');
    fprintf(fid, jsonencode(table(bardata(:,5)+bardata(:,4),bardata(:,3),bardata(:,10),regstr,'VariableNames',{'WM_neurons','All_neurons','Proportion','Region_abbreviation'}))); 
    fclose(fid);
elseif opt.skip_dur
    bh=bar(bardata(:,[10,7]),1,'grouped');
    if ~opt.skip_error_bar
        errorbar(bh(2).XEndPoints,bh(2).YEndPoints,diff(bardata(:,7:8),1,2),diff(bardata(:,[7,9]),1,2),'k.');
        errorbar(bh(1).XEndPoints,bh(1).YEndPoints,diff(bardata(:,10:11),1,2),diff(bardata(:,[10,12]),1,2),'k.');
    end
    bh(2).FaceColor='k';% mixed
    bh(1).FaceColor='w';% olf
else
    bh=bar(bardata(:,[10,13,7]),1,'grouped');
    if ~opt.skip_error_bar
        errorbar(bh(3).XEndPoints,bh(3).YEndPoints,diff(bardata(:,7:8),1,2),diff(bardata(:,[7,9]),1,2),'k.');
        errorbar(bh(1).XEndPoints,bh(1).YEndPoints,diff(bardata(:,10:11),1,2),diff(bardata(:,[10,12]),1,2),'k.');
        errorbar(bh(2).XEndPoints,bh(2).YEndPoints,diff(bardata(:,13:14),1,2),diff(bardata(:,[13,15]),1,2),'k.');
    end
    bh(3).FaceColor='w';% mixed
    bh(1).FaceColor='r';% olf
    bh(2).FaceColor='b';% dur
end
set(gca(),'YScale',opt.xyscale{2})
% if any(strcmp(opt.xyscale,'log'),'all')
% 
% else
%     ylim([0,0.6])
% end
set(gca(),'XTick',1:size(bardata,1),'XTickLabel',[regstr{:}],'XTickLabelRotation',90)
% colorize region
ymax=max(ylim());
for rr=1:numel(regstr)
    c=ephys.getRegColor(regstr{rr}{1},'large_area',true);
    text(rr,ymax,regstr{rr}{1},'HorizontalAlignment','center','VerticalAlignment','middle','Color',c,'Rotation',90)
end
% exportgraphics(fh.reg_bar,'Both_either_proportion_bars.pdf','ContentType','vector');
legend(bh,legends([2,1]),'Location','northoutside','Orientation','horizontal');


%% individual duration distribution

if opt.skip_dur
    [~,duridx]=sort(bardata(:,13),'descend');
    figure()
    bh=bar(bardata(duridx,13),1,'grouped');
    if ~opt.skip_error_bar
        errorbar(bh(1).XEndPoints,bh(1).YEndPoints,diff(bardata(duridx,10:11),1,2),diff(bardata(duridx,[13,15]),1,2),'k.');
    end
    bh(1).FaceColor='w';% olf

    set(gca(),'YScale',opt.xyscale{2})
    if any(strcmp(opt.xyscale,'log'),'all')

    else
        ylim([0,0.2])
    end
    set(gca(),'XTick',1:size(bardata,1),'XTickLabel',[regstr{duridx}],'XTickLabelRotation',90)
    %% colorize region
    ymax=max(ylim());
    for rr=1:numel(regstr)
        c=ephys.getRegColor(regstr{duridx(rr)}{1},'large_area',true);
        text(rr,ymax,regstr{duridx(rr)}{1},'HorizontalAlignment','center','VerticalAlignment','middle','Color',c,'Rotation',90)
    end
    % exportgraphics(fh.reg_bar,'Both_either_proportion_bars.pdf','ContentType','vector');
    legend(bh,legends([3]),'Location','northoutside','Orientation','horizontal');
end


