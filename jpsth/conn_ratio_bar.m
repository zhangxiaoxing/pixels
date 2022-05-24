close all;
load('0604_nonsel__conn_mat_duo_6s_1_2.mat')
load('0604_nonsel__pair_mat_duo_6s_1_2.mat')

non_pair=sum(pair_mat(:));
non_conn=sum(conn_mat(:));
[non_ratio,non_pci]=binofit(non_conn,non_pair);


load('0604_sel__conn_mat_duo_6s_1_2.mat')
load('0604_sel__pair_mat_duo_6s_1_2.mat')

sel_pair=sum(pair_mat(:));
sel_conn=sum(conn_mat(:));
[sel_ratio,sel_pci]=binofit(sel_conn,sel_pair);

bin_pair=sum(pair_sel_mat(:));
bin_conn=sum(conn_sel_mat(:));
[bin_ratio,bin_pci]=binofit(bin_conn,bin_pair);

close all
fh=figure('Color','w','Position',[100,100,110,200]);
hold on
errorbar([bin_ratio,sel_ratio,non_ratio],[diff(bin_pci)/2,diff(sel_pci)/2,diff(non_pci)/2],'k.');
bar([bin_ratio,sel_ratio,non_ratio],0.7,'FaceColor','k');
ax=gca();
ax.XTick=1:3;
ax.XTickLabel={'bin selective','delay selective','non-selective'};
ax.XTickLabelRotation=90;
ylabel('neuron connectivity ratio')
exportgraphics(fh,'conn_ratio.eps','ContentType','vector')


[~,~,pbin_sel]=crosstab([zeros(bin_pair,1);ones(sel_pair,1)],...
    [zeros(bin_conn,1);ones(bin_pair-bin_conn,1);...
    zeros(sel_conn,1);ones(sel_pair-sel_conn,1)]);

[~,~,pnon_sel]=crosstab([zeros(non_pair,1);ones(sel_pair,1)],...
    [zeros(non_conn,1);ones(non_pair-non_conn,1);...
    zeros(sel_conn,1);ones(sel_pair-sel_conn,1)]);
