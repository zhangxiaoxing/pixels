%% TODO merge with normalized_within_loops.m
keyboard()
%% count number
ring_meta=bz.rings.get_ring_meta('loadfile',true);
if false
[ratio3,z3]=get_ratio(ring_meta.congru.cross_3,ring_meta.congru.within_3_shuf);
[ratio4,z4]=get_ratio(ring_meta.congru.cross_4,ring_meta.congru.within_4_shuf);
[ratio5,z5]=get_ratio(ring_meta.congru.cross_5,ring_meta.congru.within_5_shuf);
fn='Relative_loops_number_congru.pdf';
else
[ratio3,z3]=get_ratio(ring_meta.nonmem.cross_3,ring_meta.nonmem.within_3_shuf);
[ratio4,z4]=get_ratio(ring_meta.nonmem.cross_4,ring_meta.nonmem.within_4_shuf);
[ratio5,z5]=get_ratio(ring_meta.nonmem.cross_5,ring_meta.nonmem.within_5_shuf);
fn='Relative_loops_number_nonmem.pdf';
end


rdata=[ratio3;ratio4;ratio5];
fh=figure('Color','w','Position',[32,32,155,235]);
hold on
bar(rdata(:,1),'FaceColor','w','EdgeColor','k')
errorbar(1:3,rdata(:,1),diff(rdata(:,1:2),1,2),diff(rdata(:,1:2:3),1,2),'k.')
set(gca(),'XTick',1:3,'XTickLabel',{'3-Neuron','4-Neuron','5-Neuron'},'XTickLabelRotation',90)
ylabel('Relative number of loops')
xlim([0.5,3.5])
ylim([0,50])
exportgraphics(fh,fn)






plotone(cross_3,cross_4,cross_5,within_3,within_4,within_5,'Congruent')
% 
% [cross_3,within_3]=bz.rings.rings_span('ring_size',3,'memtype','nonmem');
% [cross_4,within_4]=bz.rings.rings_span('ring_size',4,'memtype','nonmem');
% [cross_5,within_5]=bz.rings.rings_span('ring_size',5,'memtype','nonmem');
% plotone(cross_3,cross_4,cross_5,within_3,within_4,within_5,'Non-mem')

%% span
fh=figure('Color','w','Position',[32,32,400,235]);
plotonespan(cross_3,cross_3_shuf,'ffrac',1,'Selective fraction span (%)')
plotonespan(cross_4,cross_4_shuf,'ffrac',2,'Selective fraction span (%)')
plotonespan(cross_5,cross_5_shuf,'ffrac',3,'Selective fraction span (%)')
exportgraphics(fh,'loops_span_frac.pdf')

fh=figure('Color','w','Position',[32,32,400,235]);
plotonespan(cross_3,cross_3_shuf,'pv_ratio',1,'Hierarchy index')
plotonespan(cross_4,cross_4_shuf,'pv_ratio',2,'Hierarchy index')
plotonespan(cross_5,cross_5_shuf,'pv_ratio',3,'Hierarchy index')
exportgraphics(fh,'loops_span_PVSST.pdf')

fh=figure('Color','w','Position',[32,32,400,235]);
plotonespan(cross_3,cross_3_shuf,'frcom',1,'Firing rate center of mass (sec)')
plotonespan(cross_4,cross_4_shuf,'frcom',2,'Firing rate center of mass (sec)')
plotonespan(cross_5,cross_5_shuf,'frcom',3,'Firing rate center of mass (sec)')
exportgraphics(fh,'loops_span_COM.pdf')

%% dependency
function plotone(cross_3,cross_4,cross_5,within_3,within_4,within_5,ftitle)
[~,~,ratiomap]=ref.get_pv_sst();
%% within
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

%% cross
[B,BG]=groupcounts([cross_3.reg{:},cross_4.reg{:},cross_5.reg{:}].');
[B,BI]=sort(B,'descend');
BG=BG(BI);
BX=cellfun(@(x) ratiomap(x).*100,BG);
maxbx=min(numel(BX),8);
[~,maxI]=sort(min(cross_3.pv_ratio,[],2));

figure('Color','w','Position',[32,32,1800,600]);
subplot(1,3,1);
hold on
for i=1:size(cross_3.pv_ratio,1)
    plot(cross_3.pv_ratio(maxI(i),:),i,'k.')
end
ylabel('Loops No.');
xlabel('Hierarchy Index');
text(BX(1:maxbx),mean(ylim())*ones(maxbx,1),BG(1:maxbx),'Rotation',90,'Color','r','FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom');
title('3-neuron loops')

[~,maxI]=sort(min(cross_4.pv_ratio,[],2));
subplot(1,3,2);
hold on
for i=1:size(cross_4.pv_ratio,1)
    plot(cross_4.pv_ratio(maxI(i),:),i,'k.')
end
ylabel('Loops No.');
xlabel('Hierarchy Index');
text(BX(1:maxbx),mean(ylim())*ones(maxbx,1),BG(1:maxbx),'Rotation',90,'Color','r','FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom');
title('4-neuron loops')

[~,maxI]=sort(min(cross_5.pv_ratio,[],2));
subplot(1,3,3);
hold on
for i=1:size(cross_5.pv_ratio,1)
    plot(cross_5.pv_ratio(maxI(i),:),i,'k.')
end
xlabel('Hierarchy Index');
ylabel('Loops No.');
text(BX(1:maxbx),mean(ylim())*ones(maxbx,1),BG(1:maxbx),'Rotation',90,'Color','r','FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom');
title('5-neuron loops')
sgtitle(strjoin({ftitle,'cross region'}));
end

function [ratio,zscore]=get_ratio(rreal,shuf)
    shufn=cellfun(@(x) size(x.reg,1),shuf);
    realn=size(rreal.reg,1);
    shufmm=mean(shufn);
    shufstd=std(shufn);
    shufci=bootci(1000,@(x) mean(x),shufn);
    ratio=realn./([shufmm,reshape(shufci,1,[])]);
    zscore=(realn-shufmm)./shufstd;
end

function plotonespan(rreal,shuf,field,subidx,ylbl)
rspan=max(rreal.(field),[],2)-min(rreal.(field),[],2);
sspan=cell2mat(cellfun(@(x) max(x.(field),[],2)-min(x.(field),[],2),shuf,'UniformOutput',false));
mm=[median(rspan),median(sspan)];
rci=bootci(500,@(x) median(x),rspan);
sci=bootci(500,@(x) median(x),sspan);
p=ranksum(rspan,sspan);
% sspan=randsample(sspan,ceil(numel(sspan)./100));
subplot(1,3,subidx);
hold on
% swarmchart([ones(size(rspan));2*ones(size(sspan))],[rspan;sspan],4,[0.5,0.5,.5],'o','MarkerFaceAlpha',0.5,'MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor','none')
% boxplot([rspan;sspan],[ones(size(rspan));2*ones(size(sspan))],'Symbol','','Widths',0.4)
bh=bar(mm,'FaceColor','w','EdgeColor','k');
errorbar(bh.XEndPoints,mm,[rci(1),sci(1)]-mm,[rci(2),sci(2)]-mm,'k.')

text(1.5,max(ylim()),sprintf('%.3f',p),'HorizontalAlignment','center','VerticalAlignment','bottom');
set(gca(),'XTick',1:2,'XTickLabel',{'Real','Shuffled'},'XTickLabelRotation',45);
ylabel(ylbl);

end