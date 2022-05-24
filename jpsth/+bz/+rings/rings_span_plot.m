%% TODO merge with normalized_within_loops.m
% TODO waveid
keyboard()
%% count number
if ~exist('ring_meta','var')
    ring_meta=bz.rings.get_ring_meta('loadfile',true);
end
%congruent within
[ratio3cw,z3cw]=get_ratio(ring_meta.congru.within_3,ring_meta.congru.within_3_shuf);
[ratio4cw,z4cw]=get_ratio(ring_meta.congru.within_4,ring_meta.congru.within_4_shuf);
[ratio5cw,z5cw]=get_ratio(ring_meta.congru.within_5,ring_meta.congru.within_5_shuf);
%congruent cross
[ratio3cc,z3cc]=get_ratio(ring_meta.congru.cross_3,ring_meta.congru.cross_3_shuf);
[ratio4cc,z4cc]=get_ratio(ring_meta.congru.cross_4,ring_meta.congru.cross_4_shuf);
[ratio5cc,z5cc]=get_ratio(ring_meta.congru.cross_5,ring_meta.congru.cross_5_shuf);
% fn='Relative_loops_number_congru.pdf';
%nonmem within
[ratio3nw,z3nw]=get_ratio(ring_meta.nonmem.within_3,ring_meta.nonmem.within_3_shuf);
[ratio4nw,z4nw]=get_ratio(ring_meta.nonmem.within_4,ring_meta.nonmem.within_4_shuf);
[ratio5nw,z5nw]=get_ratio(ring_meta.nonmem.within_5,ring_meta.nonmem.within_5_shuf);
%nonmem cross
[ratio3nc,z3nc]=get_ratio(ring_meta.nonmem.cross_3,ring_meta.nonmem.cross_3_shuf);
[ratio4nc,z4nc]=get_ratio(ring_meta.nonmem.cross_4,ring_meta.nonmem.cross_4_shuf);
[ratio5nc,z5nc]=get_ratio(ring_meta.nonmem.cross_5,ring_meta.nonmem.cross_5_shuf);
% fn='Relative_loops_number_nonmem.pdf';
fn='Relative_loops_number.pdf';


% rdata=[ratio3c,ratio3n;ratio4c,ratio4n;ratio5c,ratio5n];
rdata=[z4nw,z4nc,z4cw,z4cc,;z5nw,z5nc,z5cw,z5cc];

fh=figure('Color','w','Position',[32,32,155,235]);
hold on
bh=bar(rdata(:,[1,4,7,10]));
[bh(1).FaceColor,bh(2).FaceColor,bh(3).FaceColor,bh(4).FaceColor]=deal('k','k','r','r');
[bh(1).FaceAlpha,bh(3).FaceAlpha]=deal(0.5);


errorbar([bh(1).XEndPoints,bh(2).XEndPoints,bh(3).XEndPoints,bh(4).XEndPoints],...
    [rdata(:,1);rdata(:,4);rdata(:,7);rdata(:,10)],...
    [diff(rdata(:,1:2),1,2);diff(rdata(:,4:5),1,2);diff(rdata(:,7:8),1,2);diff(rdata(:,10:11),1,2)],...
    [diff(rdata(:,1:2:3),1,2);diff(rdata(:,4:2:6),1,2);diff(rdata(:,7:2:9),1,2);diff(rdata(:,10:2:12),1,2)],'k.');
set(gca(),'XTick',1:3,'XTickLabel',{'4-Neuron','5-Neuron'},'XTickLabelRotation',90)
ylabel('Relative number of loops')
% xlim([0.5,3.5])
ylim([-15,150])
legend([bh(1),bh(2),bh(3),bh(4)],{'Same-memory within','Non-memory within','Same-memory cross','Non-memory cross'},'Location','northoutside');
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

function [ratio,zscore]=get_ratio(mydata,shuf,opt)
arguments
    mydata
    shuf
end
persistent OBM1map
if isempty(OBM1map)
    fstr=load('OBM1map.mat');
    OBM1map=fstr.OBM1map;
end
if iscell(mydata.reg{1,1})
    shufn=cellfun(@(x) nnz(all(OBM1map.isKey(cellfun(@(x) x,x.reg)),2)),shuf);
    realn=nnz(all(OBM1map.isKey(cellfun(@(x) x,mydata.reg)),2));
else
    shufn=cellfun(@(x) nnz(all(OBM1map.isKey(x.reg),2)),shuf);
    realn=nnz(all(OBM1map.isKey(mydata.reg),2));
end
shufmm=mean(shufn);
shufstd=std(shufn);
shufci=bootci(1000,@(x) mean(x),shufn);
ratio=realn./([shufmm,reshape(shufci,1,[])]);
zscore=(realn-[shufmm,reshape(shufci,1,[])])./shufstd;
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