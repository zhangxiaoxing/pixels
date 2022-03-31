
any6=nnz(any(meta.fdr_6(2:7,:)<0.05));
persist6=nnz(all(meta.fdr_6(2:7,:)<0.05));
persist3=nnz(all(meta.fdr_3(2:4,:)<0.05));
any3=nnz(any(meta.fdr_3(2:4,:)<0.05));
tot=numel(meta.sess);

[p3hat,p3ci]=binofit(persist3,tot);
[p6hat,p6ci]=binofit(persist6,tot);


[t3hat,t3ci]=binofit(any3-persist3,tot);
[t6hat,t6ci]=binofit(any6-persist6,tot);

fh=figure('Color','w');
hold on
bh=bar([p3hat,t3hat;p6hat,t6hat],1,'grouped');
bh(1).FaceColor='k';
bh(2).FaceColor='w';

errorbar(bh(1).XEndPoints,bh(1).YEndPoints,[p3ci(1),p6ci(1)]-[p3hat,p6hat],[p3ci(2),p6ci(2)]-[p3hat,p6hat],'k.');
errorbar(bh(2).XEndPoints,bh(2).YEndPoints,[t3ci(1),t6ci(1)]-[t3hat,t6hat],[t3ci(2),t6ci(2)]-[t3hat,t6hat],'k.');
ylabel('Fraction of all neurons (%)')
set(gca(),'XTick',1:2,'XTickLabel',{'3s','6s'},'FontSize',10,'YTickLabel',get(gca(),'YTick').*100)

% text(1,mean(ylim()),num2str(nnz(ismember(meta.mem_type,[1 3]))),'Rotation',90,'FontSize',10)
% text(2,mean(ylim()),num2str(nnz(ismember(meta.mem_type,[2 4]))),'Rotation',90,'FontSize',10)
% text(max(xlim()),max(ylim()),num2str(numel(meta.mem_type)),'HorizontalAlignment','right','VerticalAlignment','top')

exportgraphics(fh,'sust_trans_fraction.pdf');