function [map_cells,fh]=Both_either_reg_bars(opt)
arguments
    opt.skip_plot (1,1) logical = false % return map_cell data for GLM etc.
    opt.stats_model (1,:) char {mustBeMember(opt.stats_model,{'ANOVA2','ANOVA3','RANKSUM','RANKSUM2'})} = 'ANOVA2'  % 2-way anova, 3-way anova (with time bin), s1-s2 rannksum, 2-way ranksum (w/ duration)
    opt.waveid
    opt.meta
end

idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
meta=ephys.util.load_meta();
ureg=ephys.getGreyRegs();
sums=[];
% switch opt.stats_model
%     case 'RANKSUM'
%         waveid=ephys.get_wave_id(meta.sess,meta.allcid);
%     case 'ANOVA2'
%         anovameta=ephys.selectivity_anova('merge_time_bin',true);
% end

match_ranksum=(contains(opt.stats_model,'RANKSUM') && isfield(opt,'waveid') && ~isempty(opt.waveid));
match_ANOVA2=(strcmp(opt.stats_model,'ANOVA2') && isfield(opt,'meta') && ~isempty(opt.meta));
assert(match_ranksum || match_ANOVA2, 'Unmatch stats and data');


for reg=reshape(ureg,1,[])
    regsel=strcmp(meta.reg_tree(5,:),reg).';
    cnt=nnz(regsel);

    switch opt.stats_model
        case 'RANKSUM'
            waveid=opt.waveid;
            context_independent_cnt=nnz(regsel & waveid>4);
            context_depdent_cnt=nnz(regsel &waveid>0 & waveid<5);

        case 'RANKSUM2'
            waveid=opt.waveid;
            context_independent_cnt=nnz(regsel & waveid==1);
            context_depdent_cnt=nnz(regsel & waveid==2);

        case 'ANOVA2'
            anovameta=opt.meta;
            context_depdent_cnt=nnz(regsel & anovameta.mixed);
            context_independent_cnt=nnz(regsel & anovameta.sens_only);
    end
    grp=idmap.reg2tree(reg{1});
    sums=[sums;idmap.reg2ccfid(grp{6}),idmap.reg2ccfid(reg{1}),cnt,context_independent_cnt,context_depdent_cnt];
    %====================1=========================2============3=========4=======================5========
end

sums(:,6:11)=0;

for ii=1:size(sums,1)
    [phat_indi,pci_indi]=binofit(sums(ii,4),sums(ii,3));
    [phat_dep,pci_dep]=binofit(sums(ii,5),sums(ii,3));
    sums(ii,6:11)=[phat_indi,pci_indi,phat_dep,pci_dep];
    %=================6=========7|8======9======10|11====
end

%map_cells

flatten=@(y) cellfun(@(x) x,y);
bardata=sortrows(sums,6,'descend');
regstr=flatten(idmap.ccfid2reg.values(num2cell(bardata(:,2))));

context_indi_map=containers.Map(regstr,num2cell(bardata(:,[6,4,3]),2));
context_dep_map=containers.Map(regstr,num2cell(bardata(:,[9,5,3]),2));
map_cells={context_indi_map,context_dep_map};
if false % one-time export for 3d heatmap in brainrender
    for ii=1:numel(regstr)
        fprintf('[''%s'',%.3f,%.3f],\n',regstr{ii},bardata(ii,6),bardata(ii,9))
    end
end

if opt.skip_plot
    return
end

%===============================================================
fh.corr=figure('Color','w','Position',[32,32,320,320]);
hold on

for ll=1:size(bardata,1)
    c=ephys.getRegColor(regstr{ll},'large_area',true);
    scatter(bardata(ll,6),bardata(ll,9),4,c,'filled','o')
    text(bardata(ll,6),bardata(ll,9),regstr{ll},'HorizontalAlignment','center','VerticalAlignment','top','Color',c);
end
xlabel('Cross context')
ylabel('Within context')
set(gca,'XScale','log','YScale','log')
[r,p]=corr(bardata(:,6),bardata(:,9));
title(sprintf(' r = %.3f, p = %.3f',r,p));
%             exportgraphics(fh.corr1,sprintf('frac_allen_scatter_selec%d_epoch%d.pdf',featii,epochii),'ContentType','vector')

%======================


fh.reg_bar=figure('Color','w');
hold on
bh=bar(bardata(:,[6,9]),1,'grouped');
errorbar(bh(1).XEndPoints,bh(1).YEndPoints,diff(bardata(:,6:7),1,2),diff(bardata(:,[6,8]),1,2),'k.');
errorbar(bh(2).XEndPoints,bh(2).YEndPoints,diff(bardata(:,9:10),1,2),diff(bardata(:,[9,11]),1,2),'k.');
bh(1).FaceColor='k';% independent
bh(2).FaceColor='w';% dependent
set(gca(),'YScale','Log')
ylim([0.005,0.5])
set(gca(),'XTick',1:size(bardata,1),'XTickLabel',regstr,'XTickLabelRotation',90)
exportgraphics(fh.reg_bar,'Both_either_proportion_bars.pdf','ContentType','vector');
cellfun(@(x) [x{6},' ',x{7}],idmap.reg2tree.values(regstr),'UniformOutput',false)
end