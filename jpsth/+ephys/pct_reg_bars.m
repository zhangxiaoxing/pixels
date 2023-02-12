function [map_cells,fh]=pct_reg_bars(pct_meta,opt) %TODO refactoring function name
arguments
    pct_meta
    opt.skip_plot (1,1) logical = false % return map_cell data for GLM etc.
    opt.range (1,:) char {mustBeMember(opt.range,{'grey','CH','CTX'})} = 'grey'
    opt.stats_type (1,:) char = "percentile"
    opt.data_type (1,:) char = "mixed-olf-dur"
    opt.skip_error_bar (1,1) logical = true
    opt.xyscale (1,2) cell = {'Linear','Linear'}
    opt.skip_export (1,1) logical = true
end

idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
meta=ephys.util.load_meta();
ureg=ephys.getGreyRegs('range',opt.range);
sums=[];


for reg=reshape(ureg,1,[])
    regsel=strcmp(meta.reg_tree(5,:),reg).';
    cnt=nnz(regsel);

    mixed_cnt=nnz(regsel & ismember(pct_meta.wave_id,1:4));
    olf_cnt=nnz(regsel & ismember(pct_meta.wave_id,5:6));
    dur_cnt=nnz(regsel & ismember(pct_meta.wave_id,7:8));
    legends={'Mixed modality','Olfactory only','Duration only'};
            
    grp=idmap.reg2tree(reg{1});
    sums=[sums;idmap.reg2ccfid(grp{6}),idmap.reg2ccfid(reg{1}),cnt,mixed_cnt,olf_cnt,dur_cnt];
    %====================1=========================2============3=========4==================================5========
end

sums(:,7:15)=0;

for ii=1:size(sums,1)
    [mixed_hat,mixed_ci]=binofit(sums(ii,4),sums(ii,3));
    [olf_hat,olf_ci]=binofit(sums(ii,5),sums(ii,3));
    [dur_hat,dur_ci]=binofit(sums(ii,6),sums(ii,3));
    sums(ii,7:15)=[mixed_hat,mixed_ci,olf_hat,olf_ci,dur_hat,dur_ci];
    %=================6=======7|8======9======10|11====12=====13:14
end

%map_cells

flatten=@(y) cellfun(@(x) x,y);
bardata=sortrows(sums,10,'descend');
regstr=flatten(idmap.ccfid2reg.values(num2cell(bardata(:,2))));
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
mixed_map=containers.Map(regstr,num2cell(bardata(:,[7,4,3]),2));
olf_map=containers.Map(regstr,num2cell(bardata(:,[10,5,3]),2));
dur_map=containers.Map(regstr,num2cell(bardata(:,[13,5,3]),2));

map_cells={mixed_map,olf_map,dur_map};

if opt.skip_plot
    fh=[];
    return
end

%===============================================================
fh=figure('Color','w','Position',[32,32,1280,640]);
th=tiledlayout(2,4);
nexttile(5);
hold on

for ll=1:size(bardata,1)
    c=ephys.getRegColor(regstr{ll},'large_area',true);
    scatter(bardata(ll,7),bardata(ll,10),4,c,'filled','o')
    text(bardata(ll,7),bardata(ll,10),regstr{ll},'HorizontalAlignment','center','VerticalAlignment','top','Color',c);
end
xlabel(legends{1})
ylabel(legends{2})
set(gca,'XScale',opt.xyscale{1},'YScale',opt.xyscale{2})
xxdata=bardata(:,7);
yydata=bardata(:,10);
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
yy=10.^(xx.*regres(1)+regres(2));
plot(10.^xx,yy,'--k');

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
yy=10.^(xx.*regres(1)+regres(2));
plot(10.^xx,yy,'--k');


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
yy=10.^(xx.*regres(1)+regres(2));
plot(10.^xx,yy,'--k');


%             exportgraphics(fh.corr1,sprintf('frac_allen_scatter_selec%d_epoch%d.pdf',featii,epochii),'ContentType','vector')

%======================


nexttile(1,[1,3]);
hold on
bh=bar(bardata(:,[10,13,7]),1,'grouped');
if ~opt.skip_error_bar
    errorbar(bh(3).XEndPoints,bh(3).YEndPoints,diff(bardata(:,7:8),1,2),diff(bardata(:,[7,9]),1,2),'k.');
    errorbar(bh(1).XEndPoints,bh(1).YEndPoints,diff(bardata(:,10:11),1,2),diff(bardata(:,[10,12]),1,2),'k.');
    errorbar(bh(2).XEndPoints,bh(2).YEndPoints,diff(bardata(:,13:14),1,2),diff(bardata(:,[13,15]),1,2),'k.');
end
bh(3).FaceColor='w';% mixed
bh(1).FaceColor='r';% olf
bh(2).FaceColor='b';% dur

set(gca(),'YScale',opt.xyscale{2})
ylim([0.02,1])
set(gca(),'XTick',1:size(bardata,1),'XTickLabel',regstr,'XTickLabelRotation',90)
%% colorize region
ymax=max(ylim());
for rr=1:numel(regstr)
    c=ephys.getRegColor(regstr{rr},'large_area',true);
    text(rr,ymax,regstr{rr},'HorizontalAlignment','center','VerticalAlignment','bottom','Color',c,'Rotation',90)
end


% exportgraphics(fh.reg_bar,'Both_either_proportion_bars.pdf','ContentType','vector');
legend(bh,legends([2,3,1]),'Location','northoutside','Orientation','horizontal');
% cellfun(@(x) [x{6},' ',x{7}],idmap.reg2tree.values(regstr),'UniformOutput',false);
%%------------------------
% nexttile(6,[1,2]);
% hold on
% [~,sidx]=sort(summat(:,1),'descend');
% bh=bar(summat(sidx,1),'FaceColor',ones(1,3)/2);
% if ~opt.skip_error_bar
%     errorbar(bh.XEndPoints,bh.YEndPoints,diff(summat(sidx,[1,4]),1,2),diff(summat(sidx,[1,5]),1,2),'k.');
% end
% set(gca(),'XTick',1:size(bardata,1),'XTickLabel',regstr(sidx),'XTickLabelRotation',90,'Yscale',opt.xyscale{2})
% legend(bh,{'Total'},'Location','northeast','Orientation','horizontal');


%%======================
th=nexttile(4);
tbl=cell(0);
if isfield(opt,'stats_type') && ~isempty(opt.stats_type)
    tbl=[tbl;'Stats';opt.stats_type];
end
if isfield(opt,'data_type') && ~isempty(opt.data_type)
    tbl=[tbl;'Data';opt.data_type];
end
tbl=[tbl;'Range';opt.range];
ephys.util.figtable(fh,th,tbl)