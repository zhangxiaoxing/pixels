sus_trans=h5read('../transient_6.hdf5','/sus_trans');
reg_list=h5read('../transient_6.hdf5','/reg');
wrs_p_list=h5read('../transient_6.hdf5','/wrs_p');
auc_list=h5read('../transient_6.hdf5','/auc');
fr_list=h5read('../transient_6.hdf5','/fr');
sel_list=h5read('../transient_6.hdf5','/raw_selectivity');
path_list=h5read('../transient_6.hdf5','/path');



fh=figure('Color','w');
hold on;
scatter(sel_list(wrs_p_list<0.05),log10(fr_list(wrs_p_list<0.05)),3,...
    'Marker','o','MarkerFaceAlpha',0.15,'MarkerEdgeColor','none','MarkerFaceColor','r')
scatter(sel_list(wrs_p_list>=0.05),log10(fr_list(wrs_p_list>=0.05)),3,...
    'Marker','o','MarkerFaceAlpha',0.15,'MarkerEdgeColor','none','MarkerFaceColor','k')
plot(xlim(),[0,0],'k:')
ylim([-0.5,2])
ylabel('log10 firing rate');
xlabel('Selectivity index');
% exportgraphics(fh,'memory_stats_WRS_sel.pdf')

fh=figure('Color','w');
hold on;
scatter(auc_list(wrs_p_list<0.05),log10(fr_list(wrs_p_list<0.05)),3,...
    'Marker','o','MarkerFaceAlpha',0.15,'MarkerEdgeColor','none','MarkerFaceColor','r')
scatter(auc_list(wrs_p_list>=0.05),log10(fr_list(wrs_p_list>=0.05)),3,...
    'Marker','o','MarkerFaceAlpha',0.15,'MarkerEdgeColor','none','MarkerFaceColor','k')
ylim([-0.5,2])
ylabel('log10 firing rate');
xlabel('auc index');
% exportgraphics(fh,'memory_stats_WRS_sel.pdf')

fh=figure('Color','w');
hold on;
scatter(sel_list(wrs_p_list<0.05),log10(wrs_p_list(wrs_p_list<0.05)),3,...
    'Marker','o','MarkerFaceAlpha',0.15,'MarkerEdgeColor','none','MarkerFaceColor','r')
scatter(sel_list(wrs_p_list>=0.05),log10(wrs_p_list(wrs_p_list>=0.05)),3,...
    'Marker','o','MarkerFaceAlpha',0.15,'MarkerEdgeColor','none','MarkerFaceColor','k')
ylim([-6,0])
ylabel('log WRS p-value');
xlabel('selectivity index');
exportgraphics(fh,'memory_stats_WRS_sel.pdf')


fh=figure('Color','w');
hold on;
scatter(sel_list(wrs_p_list<0.05),auc_list(wrs_p_list<0.05),3,...
    'Marker','o','MarkerFaceAlpha',0.15,'MarkerEdgeColor','none','MarkerFaceColor','r')
scatter(sel_list(wrs_p_list>=0.05),auc_list(wrs_p_list>=0.05),3,...
    'Marker','o','MarkerFaceAlpha',0.15,'MarkerEdgeColor','none','MarkerFaceColor','k')
ylabel('AUROC');
xlabel('selectivity index');
exportgraphics(fh,'memory_stats_AUC_sel.pdf')


fh=figure('Color','w');
hold on;
yyaxis left
histogram(log10(wrs_p_list(:)),-6:0.1:0,'Normalization','probability');
ylabel('Probability')
yyaxis right
cdfplot(log10(wrs_p_list(:)))
ylabel('CDF')
xlim([-6,0])
set(gca,'XDir','reverse')
plot([log10(0.05),log10(0.05)],ylim(),'k:');
text(log10(0.05),1,'p = 0.05','HorizontalAlignment','left','VerticalAlignment','top')
exportgraphics(fh,'memory_stats_WRS_hist.pdf')

fh=figure('Color','w');
hold on;
yyaxis left
histogram(log10(fr_list(:)),-1:0.1:2,'Normalization','probability');
ylabel('Probability')
yyaxis right
cdfplot(log10(fr_list(:)))
ylabel('CDF')
xlim([-1,2])
plot([0,0],ylim(),'k:');
exportgraphics(fh,'memory_stats_FR_hist.pdf')

