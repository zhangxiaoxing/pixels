[cross_3,within_3]=bz.rings.rings_span('ring_size',3,'memtype','congru');
[cross_4,within_4]=bz.rings.rings_span('ring_size',4,'memtype','congru');
[cross_5,within_5]=bz.rings.rings_span('ring_size',5,'memtype','congru');
plotone(cross_3,cross_4,cross_5,within_3,within_4,within_5,'Congruent')


[cross_3,within_3]=bz.rings.rings_span('ring_size',3,'memtype','nonmem');
[cross_4,within_4]=bz.rings.rings_span('ring_size',4,'memtype','nonmem');
[cross_5,within_5]=bz.rings.rings_span('ring_size',5,'memtype','nonmem');
plotone(cross_3,cross_4,cross_5,within_3,within_4,within_5,'Non-mem')


function plotone(cross_3,cross_4,cross_5,within_3,within_4,within_5,ftitle)
[~,~,ratiomap]=ref.get_pv_sst();

figure('Color','w')
subplot(1,3,1)
[B,BG]=groupcounts(within_3.reg);
[B,BI]=sort(B,'descend');
BG=BG(BI);
B=B./sum(B);
bar(B);
set(gca(),'XTick',1:numel(B),'XTickLabel',BG,'XTickLabelRotation',90);
% set(gca(),'YScale','log');
ylabel('Fraction of all in-loop neuron')
title('3-neuron loops');

subplot(1,3,2)
[B,BG]=groupcounts(within_4.reg);
[B,BI]=sort(B,'descend');
BG=BG(BI);
B=B./sum(B);
bar(B);
set(gca(),'XTick',1:numel(B),'XTickLabel',BG,'XTickLabelRotation',90);
% set(gca(),'YScale','log');
ylabel('Fraction of all in-loop neuron')
title('4-neuron loops');

subplot(1,3,3)
[B,BG]=groupcounts(within_5.reg);
[B,BI]=sort(B,'descend');
BG=BG(BI);
B=B./sum(B);
bar(B);
set(gca(),'XTick',1:numel(B),'XTickLabel',BG,'XTickLabelRotation',90);
% set(gca(),'YScale','log');
ylabel('Fraction of all in-loop neuron')
title('5-neuron loops');
sgtitle(strjoin({ftitle,'within region'}));

[B,BG]=groupcounts([cross_3.reg{:},cross_4.reg{:},cross_5.reg{:}].');
[B,BI]=sort(B,'descend');
BG=BG(BI);
BX=cellfun(@(x) ratiomap(x),BG);
maxbx=min(numel(BX),8);

[~,maxI]=sort(min(cross_3.pv_ratio,[],2),'descend');
figure('Color','w','Position',[32,32,1800,600]);
subplot(1,3,1);
hold on
for i=1:size(cross_3.pv_ratio,1)
    plot(cross_3.pv_ratio(maxI(i),:),i,'k.')
end
ylabel('Loops No.');
xlabel('PV/(PV+SST) ratio');
text(BX(1:maxbx),mean(ylim())*ones(maxbx,1),BG(1:maxbx),'Rotation',90,'Color','r','FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom');
set(gca(),'XDir','reverse')
title('3-neuron loops')

[~,maxI]=sort(min(cross_4.pv_ratio,[],2),'descend');
subplot(1,3,2);
hold on
for i=1:size(cross_4.pv_ratio,1)
    plot(cross_4.pv_ratio(maxI(i),:),i,'k.')
end
ylabel('Loops No.');
xlabel('PV/(PV+SST) ratio');
text(BX(1:maxbx),mean(ylim())*ones(maxbx,1),BG(1:maxbx),'Rotation',90,'Color','r','FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom');
set(gca(),'XDir','reverse')
title('4-neuron loops')

[~,maxI]=sort(min(cross_5.pv_ratio,[],2),'descend');
subplot(1,3,3);
hold on
for i=1:size(cross_5.pv_ratio,1)
    plot(cross_5.pv_ratio(maxI(i),:),i,'k.')
end
xlabel('PV/(PV+SST) ratio');
ylabel('Loops No.');
text(BX(1:maxbx),mean(ylim())*ones(maxbx,1),BG(1:maxbx),'Rotation',90,'Color','r','FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom');
set(gca(),'XDir','reverse')
title('5-neuron loops')
sgtitle(strjoin({ftitle,'cross region'}));
end