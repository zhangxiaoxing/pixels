function [map_cells,fh]=Both_either_reg_bars(opt) %TODO refactoring function name
arguments
    opt.skip_plot (1,1) logical = false % return map_cell data for GLM etc.
    opt.stats_model (1,:) char {mustBeMember(opt.stats_model,{'Combinatorial','ANOVA2','ANOVA3','RANKSUM','RANKSUM2','SINGLE_MIX'})} = 'ANOVA2'  % 2-way anova, 3-way anova (with time bin), s1-s2 rannksum, 2-way ranksum (w/ duration)
    opt.waveid
    opt.meta
    opt.single_field (1,:) char {mustBeMember(opt.single_field,{'sens','dur','time_bin'})}='sens'
    opt.range (1,:) char {mustBeMember(opt.range,{'grey','CH','CTX'})} = 'grey'
    opt.stats_type (1,:) char;
    opt.data_type (1,:) char;
    opt.skip_error_bar (1,1) logical = true
    opt.yscale (1,2) cell = {'Log','Linear'}
end

idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
meta=ephys.util.load_meta();
ureg=ephys.getGreyRegs('range',opt.range);
sums=[];
% switch opt.stats_model
%     case 'RANKSUM'
%         waveid=ephys.get_wave_id(meta.sess,meta.allcid);
%     case 'ANOVA2'
%         anovameta=ephys.selectivity_anova('merge_time_bin',true);
% end

match_ranksum=(contains(opt.stats_model,'RANKSUM') && isfield(opt,'waveid') && ~isempty(opt.waveid));
match_ANOVA2=(strcmp(opt.stats_model,'ANOVA2') && isfield(opt,'meta') && ~isempty(opt.meta));
match_singlemix=(strcmp(opt.stats_model,'SINGLE_MIX') && isfield(opt,'meta') && ~isempty(opt.meta));
match_comb=(strcmp(opt.stats_model,'Combinatorial') && isfield(opt,'meta') && ~isempty(opt.meta));
assert(match_ranksum || match_ANOVA2 || match_singlemix || match_comb, 'Unmatch stats and data');


for reg=reshape(ureg,1,[])
    regsel=strcmp(meta.reg_tree(5,:),reg).';
    cnt=nnz(regsel);

    switch opt.stats_model
        case 'RANKSUM'
            waveid=opt.waveid;
            ctxt_indi_maineffect_mix_cnt=nnz(regsel & waveid>4);
            ctxt_depd_interact_single_cnt=nnz(regsel &waveid>0 & waveid<5);
            legends={'Context indepedent','Context dependent'};
        case 'RANKSUM2'
            waveid=opt.waveid;
            ctxt_indi_maineffect_mix_cnt=nnz(regsel & waveid==1);
            ctxt_depd_interact_single_cnt=nnz(regsel & waveid==2);
            legends={'Context indepedent','Context dependent'};
        case 'ANOVA2'
            anovameta=opt.meta;
            ctxt_indi_maineffect_mix_cnt=nnz(regsel & anovameta.(opt.single_field));
            ctxt_depd_interact_single_cnt=nnz(regsel & anovameta.interact);
            legends={'Maineffect','Interaction'};
        case 'SINGLE_MIX'
            singlemix_meta=opt.meta;
            ctxt_indi_maineffect_mix_cnt=nnz(regsel & (singlemix_meta.single1 & singlemix_meta.single2) );
            ctxt_depd_interact_single_cnt=nnz(regsel & (singlemix_meta.single1 | singlemix_meta.single2) & ~(singlemix_meta.single1 & singlemix_meta.single2));
            legends={'Mixed modality','Single modality'};
        case 'Combinatorial'
            comb_meta=opt.meta;
            ctxt_indi_maineffect_mix_cnt=nnz(regsel & comb_meta.typeAsel);
            ctxt_depd_interact_single_cnt=nnz(regsel & comb_meta.typeBsel);
            legends={'Mixed modality','Single modality'};   

            
    end
    grp=idmap.reg2tree(reg{1});
    sums=[sums;idmap.reg2ccfid(grp{6}),idmap.reg2ccfid(reg{1}),cnt,ctxt_indi_maineffect_mix_cnt,ctxt_depd_interact_single_cnt];
    %====================1=========================2============3=========4==================================5========
end

sums(:,6:11)=0;

for ii=1:size(sums,1)
    [phat_indi_single,pci_indi]=binofit(sums(ii,4),sums(ii,3));
    [phat_dep_mix,pci_dep]=binofit(sums(ii,5),sums(ii,3));
    sums(ii,6:11)=[phat_indi_single,pci_indi,phat_dep_mix,pci_dep];
    %=================6=============7|8======9==============10|11====
end

%map_cells

flatten=@(y) cellfun(@(x) x,y);
bardata=sortrows(sums,6,'descend');
regstr=flatten(idmap.ccfid2reg.values(num2cell(bardata(:,2))));

context_indi_single_map=containers.Map(regstr,num2cell(bardata(:,[6,4,3]),2));
context_dep_mix_map=containers.Map(regstr,num2cell(bardata(:,[9,5,3]),2));

summat=[bardata(:,3),bardata(:,4)+bardata(:,5),bardata(:,[3,3,3])];
for ii=1:size(summat,1)
    [summat(ii,1),summat(ii,4:5)]=binofit(summat(ii,2),summat(ii,3));
end
sum_map=containers.Map(regstr,num2cell(summat(:,1:3),2));
map_cells={context_indi_single_map,context_dep_mix_map,sum_map};
if false % one-time export for 3d heatmap in brainrender
    for ii=1:numel(regstr)
        fprintf('[''%s'',%.3f,%.3f,%.3f],\n',regstr{ii},bardata(ii,6),bardata(ii,9),summat(ii,1))
    end
end

if opt.skip_plot
    return
end

%===============================================================
fh=figure('Color','w','Position',[32,32,1280,640]);
th=tiledlayout(2,4);
nexttile(1);
hold on

for ll=1:size(bardata,1)
    c=ephys.getRegColor(regstr{ll},'large_area',true);
    scatter(bardata(ll,6),bardata(ll,9),4,c,'filled','o')
    text(bardata(ll,6),bardata(ll,9),regstr{ll},'HorizontalAlignment','center','VerticalAlignment','top','Color',c);
end
xlabel(legends{1})
ylabel(legends{2})
set(gca,'XScale','log','YScale','log')
[r,p]=corr(bardata(:,6),bardata(:,9));
title(sprintf(' r = %.3f, p = %.3f',r,p));
%             exportgraphics(fh.corr1,sprintf('frac_allen_scatter_selec%d_epoch%d.pdf',featii,epochii),'ContentType','vector')

%======================


nexttile(2,[1,2]);
hold on
bh=bar(bardata(:,[6,9]),1,'grouped');
if ~opt.skip_error_bar
    errorbar(bh(1).XEndPoints,bh(1).YEndPoints,diff(bardata(:,6:7),1,2),diff(bardata(:,[6,8]),1,2),'k.');
    errorbar(bh(2).XEndPoints,bh(2).YEndPoints,diff(bardata(:,9:10),1,2),diff(bardata(:,[9,11]),1,2),'k.');
end
bh(1).FaceColor='k';% independent
bh(2).FaceColor='w';% dependent
set(gca(),'YScale',opt.yscale{1})
% ylim([0.005,0.5])
set(gca(),'XTick',1:size(bardata,1),'XTickLabel',regstr,'XTickLabelRotation',90)
% exportgraphics(fh.reg_bar,'Both_either_proportion_bars.pdf','ContentType','vector');
legend(bh,legends,'Location','northoutside','Orientation','horizontal');
% cellfun(@(x) [x{6},' ',x{7}],idmap.reg2tree.values(regstr),'UniformOutput',false);
%%------------------------
nexttile(6,[1,2]);
hold on
[~,sidx]=sort(summat(:,1),'descend');
bh=bar(summat(sidx,1),'FaceColor',ones(1,3)/2);
if ~opt.skip_error_bar
    errorbar(bh.XEndPoints,bh.YEndPoints,diff(summat(sidx,[1,4]),1,2),diff(summat(sidx,[1,5]),1,2),'k.');
end
set(gca(),'XTick',1:size(bardata,1),'XTickLabel',regstr(sidx),'XTickLabelRotation',90,'Yscale',opt.yscale{2})
legend(bh,{'Total'},'Location','northeast','Orientation','horizontal');


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