idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
meta=ephys.util.load_meta();
waveid=ephys.get_wave_id(meta.sess,meta.allcid);
% ureg=unique(meta6.reg_tree(5,strcmp(meta6.reg_tree(1,:),'CH')));
ureg=ephys.getGreyRegs();
sums=[];
for reg=reshape(ureg,1,[])
    regsel=strcmp(meta.reg_tree(5,:),reg).';
    cnt=nnz(regsel);
    bothcnt=nnz(regsel & waveid>4);
    eithercnt=nnz(regsel &waveid>0 & waveid<5);
    grp=idmap.reg2tree(reg{1});
    sums=[sums;idmap.reg2ccfid(grp{6}),idmap.reg2ccfid(reg{1}),cnt,bothcnt,eithercnt];
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
bardata=sortrows(sums,6,'descend');
regstr=flatten(idmap.ccfid2reg.values(num2cell(bardata(:,2))));

figure('Color','w');
scatter(bardata(:,6),bardata(:,9))
flatten=@(y) cellfun(@(x) x,y);
text(bardata(:,6),bardata(:,9),regstr)

set(gca(),'XScale','log','YScale','log')
xlim([0.004,0.5])
ylim([0.004,0.5])
xlabel('Odor selective both duration')
ylabel('Odor selective either duration')

both_map=containers.Map(regstr,num2cell(bardata(:,[6,4,3]),2));
either_map=containers.Map(regstr,num2cell(bardata(:,[9,5,3]),2));
map_cells={both_map,either_map};

fh=figure('Color','w')
hold on
bh=bar(bardata(:,[6,9]),1,'grouped');
errorbar(bh(1).XEndPoints,bh(1).YEndPoints,diff(bardata(:,6:7),1,2),diff(bardata(:,[6,8]),1,2),'k.');
errorbar(bh(2).XEndPoints,bh(2).YEndPoints,diff(bardata(:,9:10),1,2),diff(bardata(:,[9,11]),1,2),'k.');
bh(1).FaceColor='k';
bh(2).FaceColor='w';
set(gca(),'YScale','Log')
ylim([0.005,0.5])
set(gca(),'XTick',1:size(bardata,1),'XTickLabel',regstr,'XTickLabelRotation',90)
exportgraphics(fh,'Both_either_proportion_bars.pdf','ContentType','vector');
cellfun(@(x) [x{6},' ',x{7}],idmap.reg2tree.values(regstr),'UniformOutput',false)
