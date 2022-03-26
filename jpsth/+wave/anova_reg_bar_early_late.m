
load('map_cells.mat','map_cells');
keys=map_cells{3,1}.keys();
late_samp=cell2mat(map_cells{3,1}.values().');
[late_samp_srt,sidx]=sort(late_samp(:,1),'descend');
skeys=keys(sidx);
late_samp=late_samp(sidx,:);
early_samp=cell2mat(map_cells{2,1}.values(skeys).');
late_dur=cell2mat(map_cells{3,3}.values(skeys).');
early_dur=cell2mat(map_cells{2,3}.values(skeys).');

barmat=[early_samp(:,1),late_samp_srt,early_dur(:,1),late_dur(:,1)];

%%

fh=figure('Color','w','Position',[32,32,1500,900]);
subplot(2,1,1)
hold on
bh=bar(1:size(barmat,1),barmat(:,1:2).*100,1,'grouped');
bh(1).FaceColor='#4DBEEE';
bh(2).FaceColor='#D95319';
[~,pcies]=binofit(early_samp(:,2),early_samp(:,3));
[~,pcils]=binofit(late_samp(:,2),late_samp(:,3));
errorbar(bh(1).XEndPoints,bh(1).YEndPoints,pcies(:,1),pcies(:,2),'k.')
errorbar(bh(2).XEndPoints,bh(2).YEndPoints,pcils(:,1),pcils(:,2),'k.')
set(gca(),'XTick',1:numel(skeys),'XTickLabel',skeys,'XTickLabelRotation',90)
ylabel('Sample selective neuron (%)')
ylim([0,30])
subplot(2,1,2)
hold on
bh=bar(barmat(:,3:4).*100,1,'grouped');
bh(1).FaceColor='#0072BD';
bh(2).FaceColor='#A2142F';
[~,pcied]=binofit(early_dur(:,2),early_dur(:,3));
[~,pcild]=binofit(late_dur(:,2),late_dur(:,3));
errorbar(bh(1).XEndPoints,bh(1).YEndPoints,pcied(:,1),pcied(:,2),'k.')
errorbar(bh(2).XEndPoints,bh(2).YEndPoints,pcild(:,1),pcild(:,2),'k.')
set(gca(),'XTick',1:numel(skeys),'XTickLabel',skeys,'XTickLabelRotation',90)
ylabel('Duration selective neuron (%)')
ylim([0,30])
exportgraphics(fh,'anova_reg_bar_early_late.pdf','ContentType','vector')
%%
fh=figure('Color','w','Position',[32,32,215,225]);
subplot(1,2,1)
hold on
bh=plot(barmat(:,1:2).'.*100,'-k');
plot(ones(size(barmat(:,1))),barmat(:,1).'.*100,'o','MarkerEdgeColor','#4DBEEE');
plot(ones(size(barmat(:,2))).*2,barmat(:,2).'.*100,'o','MarkerEdgeColor','#D95319');
plot(0.5,median(barmat(:,1)).*100,'o','MarkerFaceColor','#4DBEEE','MarkerEdgeColor','#4DBEEE');
plot(2.5,median(barmat(:,2)).*100,'o','MarkerFaceColor','#D95319','MarkerEdgeColor','#D95319');
iqres=(prctile(barmat(:,1),[25,75])-median(barmat(:,1))).*100;
iqrls=(prctile(barmat(:,2),[25,75])-median(barmat(:,2))).*100;
errorbar(0.5,median(barmat(:,1)).*100,iqres(1),iqres(2),'k.')
errorbar(2.5,median(barmat(:,2)).*100,iqrls(1),iqrls(2),'k.')
xlim([0,3])
ylim([0,30])
ylabel('Sample selective neuron (%)')
set(gca(),'XTick',[])
subplot(1,2,2)
hold on
bh=plot(barmat(:,3:4).'.*100,'-k');
plot(ones(size(barmat(:,3))),barmat(:,3).'.*100,'o','MarkerEdgeColor','#0072BD');
plot(ones(size(barmat(:,4))).*2,barmat(:,4).'.*100,'o','MarkerEdgeColor','#A2142F');
plot(0.5,median(barmat(:,3)).*100,'o','MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD');
plot(2.5,median(barmat(:,4)).*100,'o','MarkerFaceColor','#A2142F','MarkerEdgeColor','#A2142F');
iqred=(prctile(barmat(:,3),[25,75])-median(barmat(:,3))).*100;
iqrld=(prctile(barmat(:,4),[25,75])-median(barmat(:,4))).*100;
errorbar(0.5,median(barmat(:,3)).*100,iqred(1),iqred(2),'k.')
errorbar(2.5,median(barmat(:,4)).*100,iqrld(1),iqrld(2),'k.')
xlim([0,3])
ylim([0,30])
ylabel('Duration selective neuron (%)')
set(gca(),'XTick',[])
exportgraphics(fh,'anova_reg_bar_early_late.pdf','ContentType','vector','Append',true)
format long
disp([signrank(barmat(:,1),barmat(:,2)),signrank(barmat(:,3),barmat(:,4)),signrank(diff(barmat(:,1:2),1,2),diff(barmat(:,3:4),1,2))])
format short

%%
fh=figure('Color','w','Position',[32,32,600,300]);
mpath=mfilename('fullpath');
if contains(mpath,'LiveEditorEvaluationHelper'),mpath='';end
pathstr=strjoin({mpath,matlab.desktop.editor.getActiveFilename},'|');
text(0.5,0.5,pathstr,'HorizontalAlignment','center','VerticalAlignment','middle')
exportgraphics(fh,'anova_reg_bar_early_late.pdf','ContentType','vector','Append',true)


