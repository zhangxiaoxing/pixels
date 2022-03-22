
load('map_cells.mat','map_cells');
keys=map_cells{3,1}.keys();
late_samp=cell2mat(map_cells{3,1}.values().');
[late_samp_srt,sidx]=sort(late_samp(:,1),'descend');
skeys=keys(sidx);
early_samp=cell2mat(map_cells{2,1}.values(skeys).');
late_dur=cell2mat(map_cells{3,3}.values(skeys).');
early_dur=cell2mat(map_cells{2,3}.values(skeys).');

barmat=[early_samp(:,1),late_samp_srt,early_dur(:,1),late_dur(:,1)];



figure('Color','w')
subplot(2,1,1)
bh=bar(barmat(:,1:2).*100,'grouped');
bh(1).FaceColor='b';
bh(2).FaceColor='r';
set(gca(),'XTick',1:numel(skeys),'XTickLabel',skeys)
ylabel('Sample selective neuron (%)')
subplot(2,1,2)
bh=bar(barmat(:,3:4).*100,'grouped');
bh(1).FaceColor='b';
bh(2).FaceColor='r';
set(gca(),'XTick',1:numel(skeys),'XTickLabel',skeys)
ylabel('Duration selective neuron (%)')

%%
figure('Color','w')
subplot(1,2,1)
hold on
bh=plot(barmat(:,1:2).'.*100,'-k');
plot(ones(size(barmat(:,1))),barmat(:,1).'.*100,'ro');
plot(ones(size(barmat(:,2))).*2,barmat(:,2).'.*100,'bo');
plot(0.5,median(barmat(:,1)).*100,'ro','MarkerFaceColor','r')
plot(2.5,median(barmat(:,2)).*100,'bo','MarkerFaceColor','b')
xlim([0,3])
ylabel('Sample selective neuron (%)')
subplot(1,2,2)
hold on
bh=plot(barmat(:,3:4).'.*100,'-k');
plot(ones(size(barmat(:,3))),barmat(:,3).'.*100,'ro');
plot(ones(size(barmat(:,4))).*2,barmat(:,4).'.*100,'bo');
plot(0.5,median(barmat(:,3)).*100,'ro','MarkerFaceColor','r')
plot(2.5,median(barmat(:,4)).*100,'bo','MarkerFaceColor','b')
xlim([0,3])

