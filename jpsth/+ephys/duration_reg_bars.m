idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
meta=ephys.util.load_meta();
load('anovameta.mat','anovameta');%K:\code\jpsth\+ephys\selectivity_anova.m
%1:Samp,2:Dur,3:Bin,4:Samp*Dur,5:Samp*Bin,6:Dur*Bin,7:Samp*Dur*Bin
dur_sel_mix=any(anovameta.anovap(:,[2 4 6 7])<0.05,2);
sens_sel_mix=any(anovameta.anovap(:,[1 4 5 7])<0.05,2);
dur_sel_exclu=dur_sel_mix & ~sens_sel_mix;
ureg=ephys.getGreyRegs();
sums=[];
for reg=reshape(ureg,1,[])
    regsel=strcmp(meta.reg_tree(5,:),reg).';
    cnt=nnz(regsel);
    
    mix_cnt=nnz(regsel & dur_sel_mix);
    exclu_cnt=nnz(regsel & dur_sel_exclu);
    grp=idmap.reg2tree(reg{1});
    sums=[sums;idmap.reg2ccfid(grp{6}),idmap.reg2ccfid(reg{1}),cnt,mix_cnt,exclu_cnt];
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
keyboard(); % egress point for alternative data processing

figure('Color','w');
scatter(bardata(:,6),bardata(:,9))
[r,p]=corr(bardata(:,6),bardata(:,9));
text(bardata(:,6),bardata(:,9),regstr)

set(gca(),'XScale','linear','YScale','linear')
xlim([0.01,0.5])
ylim([0.01,0.5])
xlabel('Odor selective both duration')
ylabel('Odor selective either duration')

fh=figure('Color','w')
hold on
bh=bar(bardata(:,[9,6]),1,'grouped');
errorbar(bh(2).XEndPoints,bh(2).YEndPoints,diff(bardata(:,6:7),1,2),diff(bardata(:,[6,8]),1,2),'k.');
errorbar(bh(1).XEndPoints,bh(1).YEndPoints,diff(bardata(:,9:10),1,2),diff(bardata(:,[9,11]),1,2),'k.');
bh(1).FaceColor='k';
bh(2).FaceColor='w';
set(gca(),'YScale','log')
ylim([0.02,0.5])
set(gca(),'XTick',1:size(bardata,1),'XTickLabel',regstr,'XTickLabelRotation',90,'YTick',[0.02,0.1,0.2],'YTickLabel',[0.02,0.1,0.2]*100)
exportgraphics(fh,'Dur_proportion_bars.pdf','ContentType','vector');
cellfun(@(x) [x{6},' ',x{7}],idmap.reg2tree.values(regstr),'UniformOutput',false)
