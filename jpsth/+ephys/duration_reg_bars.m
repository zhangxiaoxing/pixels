function [map_cells,fh]=duration_reg_bars(opt)
arguments
    opt.skip_plot (1,1) logical = false % return map_cell data for GLM etc.
    opt.stats_model (1,:) char {mustBeMember(opt.stats_model,{'ANOVA2','ANOVA3','RANKSUM','RANKSUM2'})} = 'ANOVA2'  % 2-way anova, 3-way anova (with time bin), s1-s2 rannksum, 2-way ranksum (w/ duration)
end

ureg=ephys.getGreyRegs();
meta=ephys.util.load_meta();
idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
sums=[];

switch opt.stats_model
    case 'RANKSUM'
        assert(false,'Unimplemented')
        waveid=ephys.get_wave_id(meta.sess,meta.allcid);
    case 'ANOVA2'
        anovameta=ephys.selectivity_anova('merge_time_bin',true);
end

for reg=reshape(ureg,1,[])
    regsel=strcmp(meta.reg_tree(5,:),reg).';
    cnt=nnz(regsel);
    switch opt.stats_model
        case 'RANKSUM'
            assert(false,'Unimplemented')
%             sens_only_cnt=nnz(regsel & waveid>4);
%             mixed_cnt=nnz(regsel &waveid>0 & waveid<5);
        case 'ANOVA2'
            mixed_cnt=nnz(regsel & anovameta.mixed);
            dur_only_cnt=nnz(regsel & anovameta.dur_only);
    end
    grp=idmap.reg2tree(reg{1});
    sums=[sums;idmap.reg2ccfid(grp{6}),idmap.reg2ccfid(reg{1}),cnt,dur_only_cnt,mixed_cnt];
    %====================1=========================2============3======4=======5========
end
sums(:,6:11)=0;

for ii=1:size(sums,1)
    [phatb,pcib]=binofit(sums(ii,4),sums(ii,3));
    [phate,pcie]=binofit(sums(ii,5),sums(ii,3));
    sums(ii,6:11)=[phatb,pcib,phate,pcie];
    %================6===7|8====9===10|11====
end

%map_cells
flatten=@(y) cellfun(@(x) x,y);
bardata=sortrows(sums,9,'descend');
regstr=flatten(idmap.ccfid2reg.values(num2cell(bardata(:,2))));

mixed_map=containers.Map(regstr,num2cell(bardata(:,[6,4,3]),2));
exclu_map=containers.Map(regstr,num2cell(bardata(:,[9,5,3]),2));
map_cells={mixed_map,exclu_map};
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
xlabel('duration selective only')
ylabel('Mixed modality')
set(gca,'XScale','log','YScale','log')
[r,p]=corr(bardata(:,6),bardata(:,9));
title(sprintf(' r = %.3f, p = %.3f',r,p));
%             exportgraphics(fh.corr1,sprintf('frac_allen_scatter_selec%d_epoch%d.pdf',featii,epochii),'ContentType','vector')

%======================

fh.reg_bar=figure('Color','w');
hold on
bh=bar(bardata(:,[9,6]),1,'grouped');
errorbar(bh(2).XEndPoints,bh(2).YEndPoints,diff(bardata(:,6:7),1,2),diff(bardata(:,[6,8]),1,2),'k.');
errorbar(bh(1).XEndPoints,bh(1).YEndPoints,diff(bardata(:,9:10),1,2),diff(bardata(:,[9,11]),1,2),'k.');
bh(1).FaceColor='k';
bh(2).FaceColor='w';
set(gca(),'YScale','log')
ylim([0.02,0.5])
set(gca(),'XTick',1:size(bardata,1),'XTickLabel',regstr,'XTickLabelRotation',90,'YTick',[0.02,0.1,0.2],'YTickLabel',[0.02,0.1,0.2]*100)
exportgraphics(fh.reg_bar,'Dur_proportion_bars.pdf','ContentType','vector');
cellfun(@(x) [x{6},' ',x{7}],idmap.reg2tree.values(regstr),'UniformOutput',false)
end