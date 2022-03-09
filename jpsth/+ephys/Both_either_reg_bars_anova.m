load('anovameta.mat','anovameta'); % K:\code\jpsth\+ephys\selectivity_anova.m
meta=ephys.util.load_meta();
% waveid=ephys.get_wave_id(meta.sess,meta.allcid);

dur_indep_sel=(anovameta.anovap(:,1)<0.05 | anovameta.anovap(:,5)<0.05);
dur_dep_sel=(anovameta.anovap(:,4)<0.05 | anovameta.anovap(:,7)<0.05) &~dur_indep_sel;
dur_only_sel=(anovameta.anovap(:,2)<0.05 | anovameta.anovap(:,6)<0.05) & ~dur_dep_sel;
idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));
ureg=unique(meta.reg_tree(5,strcmp(meta.reg_tree(1,:),'CH') | strcmp(meta.reg_tree(1,:),'BS')));
sums=[];
for reg=reshape(ureg,1,[])
    if isempty(reg{1}),continue;end
    regsel=strcmp(meta.reg_tree(5,:),reg).';
    cnt=nnz(regsel);
    bothcnt=nnz(regsel & dur_indep_sel);
    eithercnt=nnz(regsel & dur_dep_sel);
    durcnt=nnz(regsel & dur_only_sel);
    grp=idmap.reg2tree(reg{1});
    sums=[sums;idmap.reg2ccfid(grp{6}),idmap.reg2ccfid(reg{1}),cnt,bothcnt,eithercnt,durcnt];
end

sums(:,7:15)=0;

for ii=1:size(sums,1)
    [phatb,pcib]=binofit(sums(ii,4),sums(ii,3));
    [phate,pcie]=binofit(sums(ii,5),sums(ii,3));
    [phatd,pcid]=binofit(sums(ii,6),sums(ii,3));
    sums(ii,7:15)=[phatb,pcib,phate,pcie,phatd,pcid];
end


%% dur-dep|indep sense-selective, scatter
pctsel=sums(:,3)>100;
figure('Color','w');
scatter(sums(pctsel,4)./sums(pctsel,3),sums(pctsel,5)./sums(pctsel,3))
flatten=@(y) cellfun(@(x) x,y);

text(sums(pctsel,4)./sums(pctsel,3),sums(pctsel,5)./sums(pctsel,3),flatten(idmap.ccfid2reg.values(num2cell(sums(pctsel,2)))),'HorizontalAlignment','center','VerticalAlignment','top')

[r,p]=corr(sums(pctsel,4)./sums(pctsel,3),sums(pctsel,5)./sums(pctsel,3));

set(gca(),'XScale','linear','YScale','linear')
xlim([0,0.8])
ylim([0,0.3])
xlabel('Unspecific selective proportion')
ylabel('Specific selective proportion')
text(min(xlim()),min(ylim()),sprintf('r=%.2f,p=%.2f',r,p),'HorizontalAlignment','left','VerticalAlignment','bottom')

%% dur-only(sense-exclueded) scatter
figure('Color','w');
scatter(sums(pctsel,4)./sums(pctsel,3),sums(pctsel,6)./sums(pctsel,3))
flatten=@(y) cellfun(@(x) x,y);
text(sums(pctsel,4)./sums(pctsel,3),sums(pctsel,6)./sums(pctsel,3),flatten(idmap.ccfid2reg.values(num2cell(sums(pctsel,2)))),'HorizontalAlignment','center','VerticalAlignment','top')

[r,p]=corr(sums(pctsel,4)./sums(pctsel,3),sums(pctsel,6)./sums(pctsel,3));

set(gca(),'XScale','linear','YScale','linear')
xlim([0,0.8])
ylim([0,0.3])
xlabel('Unspecific selective proportion')
ylabel('Duration-only proportion')
text(min(xlim()),min(ylim()),sprintf('r=%.2f,p=%.2f',r,p),'HorizontalAlignment','left','VerticalAlignment','bottom')
%%%%%%%%%%%
figure('Color','w');
scatter(sums(pctsel,5)./sums(pctsel,3),sums(pctsel,6)./sums(pctsel,3))
flatten=@(y) cellfun(@(x) x,y);
text(sums(pctsel,5)./sums(pctsel,3),sums(pctsel,6)./sums(pctsel,3),flatten(idmap.ccfid2reg.values(num2cell(sums(pctsel,2)))),'HorizontalAlignment','center','VerticalAlignment','top')

[r,p]=corr(sums(pctsel,5)./sums(pctsel,3),sums(pctsel,6)./sums(pctsel,3));

set(gca(),'XScale','linear','YScale','linear')
xlim([0,0.3])
ylim([0,0.3])
xlabel('Duration-dependent olfaction selective proportion')
ylabel('Duration-only proportion')
text(min(xlim()),min(ylim()),sprintf('r=%.2f,p=%.2f',r,p),'HorizontalAlignment','left','VerticalAlignment','bottom')




%% dur-dep|indep sense-selective, bar
bardata=sortrows(sums(pctsel,:),7,'descend');

fh=figure('Color','w');
hold on
bh=bar(bardata(:,[7,10]),'grouped');
errorbar(bh(1).XEndPoints,bh(1).YEndPoints,diff(bardata(:,7:8),1,2),diff(bardata(:,[7,9]),1,2),'k.');
errorbar(bh(2).XEndPoints,bh(2).YEndPoints,diff(bardata(:,10:11),1,2),diff(bardata(:,[10,12]),1,2),'k.');
bh(1).FaceColor='r';
bh(2).FaceColor='b';
set(gca(),'YScale','Linear')
ylim([0,0.8])
set(gca(),'XTick',1:size(bardata,1),'XTickLabel',flatten(idmap.ccfid2reg.values(num2cell(bardata(:,2)))),'XTickLabelRotation',90)
exportgraphics(fh,'Both_either_proportion_bars.pdf','ContentType','vector');
cellfun(@(x) [x{6},' ',x{7}],idmap.reg2tree.values(flatten(idmap.ccfid2reg.values(num2cell(bardata(:,2))))),'UniformOutput',false)
ylabel('Proportion of selective neuron')
set(gca(),'YTick',0:0.1:0.8,'YTickLabel',0:10:80)
legend(bh,{'Duration-independent','Duration-dependent'})


%% dur-dep|indep sense-selective, bar
bardata=sortrows(sums(pctsel,:),7,'descend');

fh=figure('Color','w');
hold on
bh=bar(bardata(:,[7,10]),'grouped');
errorbar(bh(1).XEndPoints,bh(1).YEndPoints,diff(bardata(:,7:8),1,2),diff(bardata(:,[7,9]),1,2),'k.');
errorbar(bh(2).XEndPoints,bh(2).YEndPoints,diff(bardata(:,10:11),1,2),diff(bardata(:,[10,12]),1,2),'k.');
bh(1).FaceColor='r';
bh(2).FaceColor='b';
set(gca(),'YScale','Linear')
ylim([0,0.8])
set(gca(),'XTick',1:size(bardata,1),'XTickLabel',flatten(idmap.ccfid2reg.values(num2cell(bardata(:,2)))),'XTickLabelRotation',90)
exportgraphics(fh,'Both_either_proportion_bars.pdf','ContentType','vector');
cellfun(@(x) [x{6},' ',x{7}],idmap.reg2tree.values(flatten(idmap.ccfid2reg.values(num2cell(bardata(:,2))))),'UniformOutput',false)
ylabel('Proportion of selective neuron')
set(gca(),'YTick',0:0.1:0.8,'YTickLabel',0:10:80)
legend(bh,{'Duration-independent','Duration-dependent'})




%% dur-dep|indep sense-selective, bar
bardata=sortrows(sums(pctsel,:),13,'descend'); %7:dur-indep sense, 10:dur-dep sense, 13:dur-only(no-sense)

fh=figure('Color','w');
hold on
bh=bar(bardata(:,[13,10]),'grouped');
errorbar(bh(1).XEndPoints,bh(1).YEndPoints,diff(bardata(:,13:14),1,2),diff(bardata(:,[13,15]),1,2),'k.');
errorbar(bh(2).XEndPoints,bh(2).YEndPoints,diff(bardata(:,10:11),1,2),diff(bardata(:,[10,12]),1,2),'k.');
bh(1).FaceColor='r';
bh(2).FaceColor='b';
set(gca(),'YScale','Linear')
ylim([0,0.4])
set(gca(),'XTick',1:size(bardata,1),'XTickLabel',flatten(idmap.ccfid2reg.values(num2cell(bardata(:,2)))),'XTickLabelRotation',90)
exportgraphics(fh,'Both_either_proportion_bars.pdf','ContentType','vector');
cellfun(@(x) [x{6},' ',x{7}],idmap.reg2tree.values(flatten(idmap.ccfid2reg.values(num2cell(bardata(:,2))))),'UniformOutput',false)
ylabel('Proportion of selective neuron (%)')
set(gca(),'YTick',0:0.1:0.4,'YTickLabel',0:10:40)
legend(bh,{'Duration-only','Duration-dependent sensory'})
