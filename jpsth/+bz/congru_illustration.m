congru=[1,1;2,2;3,3;4,4;5,5;6,6;7,7;8,8;1,5;2,5;3,6;4,6;1,7;2,8;3,7;4,8];
congru=unique([congru;fliplr(congru)],'rows');
incongru=[1,2;1,3;1,4;2,3;2,4;3,4;5,6;7,8;1,6;2,6;3,5;4,5;1,8;2,7;3,8;4,7];
incongru=unique([incongru;fliplr(incongru)],'rows');
nonmem=[9,9];
figure();
hold on;
scatter(congru(:,1),congru(:,2),400,'ro','filled')
scatter(incongru(:,1),incongru(:,2),400,'bo','filled')
scatter(nonmem(:,1),nonmem(:,2),400,'ko','filled')
xlim([0.5,9.5]);
ylim([0.5,9.5]);
lbl={'S1,3s','S1,6s','S2,3s','S2,6s','S1','S2','3s','6s','Non-memory'}
set(gca,'XTick',1:9,'XTickLabel',lbl,'YTick',1:9,'YTickLabel',lbl,'XTickLabelRotation',90)
xlabel('Leading neuron preferred condition')
ylabel('following neuron preferred condition')