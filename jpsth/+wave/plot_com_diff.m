samedist=wave.com_diff('level',5,'diff',false,'to_plot',false,'per_sec_stats',false);
diffdist=wave.com_diff('level',5,'diff',true,'to_plot',false,'per_sec_stats',false);
bin_edge=-2000:100:2000;
sc=histcounts(samedist,bin_edge,'Normalization','cdf');
dc=histcounts(diffdist,bin_edge,'Normalization','cdf');

boots=bootstrp(1000,@(x) histcounts(x,bin_edge,'Normalization','cdf'),samedist);
bootd=bootstrp(1000,@(x) histcounts(x,bin_edge,'Normalization','cdf'),diffdist);

cis=[prctile(boots,2.5);prctile(boots,97.5)];
cid=[prctile(bootd,2.5);prctile(bootd,97.5)];

close all
fh=figure('Color','w','Position',[100,100,215,215]);
hold onsamedist=wave.com_diff('level',5,'diff',false,'to_plot',false,'per_sec_stats',false);
diffdist=wave.com_diff('level',5,'diff',true,'to_plot',false,'per_sec_stats',false);
bin_edge=-2000:100:2000;
sc=histcounts(samedist,bin_edge,'Normalization','cdf');
dc=histcounts(diffdist,bin_edge,'Normalization','cdf');

boots=bootstrp(1000,@(x) histcounts(x,bin_edge,'Normalization','cdf'),samedist);
bootd=bootstrp(1000,@(x) histcounts(x,bin_edge,'Normalization','cdf'),diffdist);

cis=[prctile(boots,2.5);prctile(boots,97.5)];
cid=[prctile(bootd,2.5);prctile(bootd,97.5)];

close all
fh=figure('Color','w','Position',[100,100,215,215]);
hold on
fill([bin_edge(1:end-1)+50,fliplr(bin_edge(1:end-1)+50)],...
    [cis(1,:),fliplr(cis(2,:))],'k','EdgeColor','none','FaceAlpha',0.2)
fill([bin_edge(1:end-1)+50,fliplr(bin_edge(1:end-1)+50)],...
    [cid(1,:),fliplr(cid(2,:))],'r','EdgeColor','none','FaceAlpha',0.2)
hs=plot(bin_edge(1:end-1)+50,sc,'-k');
hd=plot(bin_edge(1:end-1)+50,dc,'-r');

xline(mean(samedist),'--k')
xline(mean(diffdist),'--r')
fprintf('same mean %.3f, diff mean %.3f\n',mean(samedist),mean(diffdist));
[h,p]=ttest2(samedist,diffdist);
% legend([hs,hd],{sprintf('Within region,\nmean=%.1fms',mean(samedist)),...
%     sprintf('Between region,\nmean=%.1fms',mean(diffdist))},...
%     'Location','northoutside');

legend([hs,hd],{'Within','Between'},...
    'Location','southeast');

xlim([-1000,1000])
ylim([0,1])
text(min(xlim()),max(ylim()),sprintf('p = %.3f',p),'HorizontalAlignment','left','VerticalAlignment','top')
xlabel('COM lag (ms)');
ylabel('Cumulated probability')

exportgraphics(fh,'FC_COM_bias_lvl_5.pdf');

fill([bin_edge(1:end-1)+50,fliplr(bin_edge(1:end-1)+50)],...
    [cis(1,:),fliplr(cis(2,:))],'k','EdgeColor','none','FaceAlpha',0.2)
fill([bin_edge(1:end-1)+50,fliplr(bin_edge(1:end-1)+50)],...
    [cid(1,:),fliplr(cid(2,:))],'r','EdgeColor','none','FaceAlpha',0.2)
hs=plot(bin_edge(1:end-1)+50,sc,'-k');
hd=plot(bin_edge(1:end-1)+50,dc,'-r');

xline(mean(samedist),'--k')
xline(mean(diffdist),'--r')
fprintf('same mean %.3f, diff mean %.3f\n',mean(samedist),mean(diffdist));
[h,p]=ttest2(samedist,diffdist);
% legend([hs,hd],{sprintf('Within region,\nmean=%.1fms',mean(samedist)),...
%     sprintf('Between region,\nmean=%.1fms',mean(diffdist))},...
%     'Location','northoutside');

legend([hs,hd],{'Within','Between'},...
    'Location','southeast');

xlim([-1000,1000])
ylim([0,1])
text(min(xlim()),max(ylim()),sprintf('p = %.3f',p),'HorizontalAlignment','left','VerticalAlignment','top')
xlabel('COM lag (ms)');
ylabel('Cumulated probability')

exportgraphics(fh,'FC_COM_bias_lvl_5.pdf');
