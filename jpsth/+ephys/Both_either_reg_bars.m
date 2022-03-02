idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));
meta6=ephys.util.load_meta('delay',6);
waveid=ephys.get_wave_id(meta6.sess,meta6.allcid);
ureg=unique(meta6.reg_tree(5,strcmp(meta6.reg_tree(1,:),'CH')));
sums=[];
for reg=reshape(ureg,1,[])
    if isempty(reg{1}),continue;end
    regsel=strcmp(meta6.reg_tree(5,:),reg).';
    cnt=nnz(regsel);
    bothcnt=nnz(regsel & waveid>4);
    eithercnt=nnz(regsel &waveid>0 & waveid<5);
    grp=idmap.reg2tree(reg{1});
    sums=[sums;idmap.reg2ccfid(grp{6}),idmap.reg2ccfid(reg{1}),cnt,bothcnt,eithercnt];
end

sums(:,6:11)=0;

for ii=1:size(sums,1)
    [phatb,pcib]=binofit(sums(ii,4),sums(ii,3));
    [phate,pcie]=binofit(sums(ii,5),sums(ii,3));
    sums(ii,6:11)=[phatb,pcib,phate,pcie];
end


pctsel=sums(:,3)>100;

figure('Color','w');
scatter(sums(pctsel,4)./sums(pctsel,3),sums(pctsel,5)./sums(pctsel,3))
flatten=@(y) cellfun(@(x) x,y);

text(sums(pctsel,4)./sums(pctsel,3),sums(pctsel,5)./sums(pctsel,3),flatten(idmap.ccfid2reg.values(num2cell(sums(pctsel,2)))))

set(gca(),'XScale','log','YScale','log')
xlim([0.004,0.5])
ylim([0.004,0.5])
xlabel('Unspecific selective proportion')
ylabel('Specific selective proportion')

bardata=sortrows(sums(pctsel,:),6,'descend');

fh=figure('Color','w')
hold on
bh=bar(bardata(:,[6,9]),'grouped');
errorbar(bh(1).XEndPoints,bh(1).YEndPoints,diff(bardata(:,6:7),1,2),diff(bardata(:,[6,8]),1,2),'k.');
errorbar(bh(2).XEndPoints,bh(2).YEndPoints,diff(bardata(:,9:10),1,2),diff(bardata(:,[9,11]),1,2),'k.');
bh(1).FaceColor='k';
bh(2).FaceColor='w';
set(gca(),'YScale','Log')
ylim([0.002,0.5])
set(gca(),'XTick',1:size(bardata,1),'XTickLabel',flatten(idmap.ccfid2reg.values(num2cell(bardata(:,2)))),'XTickLabelRotation',90)
exportgraphics(fh,'Both_either_proportion_bars.pdf','ContentType','vector');
cellfun(@(x) [x{6},' ',x{7}],idmap.reg2tree.values(flatten(idmap.ccfid2reg.values(num2cell(bardata(:,2))))),'UniformOutput',false)