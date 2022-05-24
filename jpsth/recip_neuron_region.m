load('reg_keep.mat');
load('su_region_recip.mat')
load('io_sel.mat')

sus_trans=h5read('../transient_6.hdf5','/sus_trans');%4+2+7
reg_all=deblank(h5read('../transient_6.hdf5','/reg'));


wing=diff(ioselstats{1}(:,[3 7]),1,2);
per_reg_recip_idx=nanmean(recip_count./input_to_count,2);
per_reg_shuf=nanmean(recip_count_shuf./input_to_count_shuf,2);


per_reg_recip_idx=[per_reg_recip_idx,(1:140)',wing];
recip_idx=per_reg_recip_idx(sum(input_to_count,2)>120,:);
recip_shuf=squeeze(per_reg_shuf(sum(input_to_count,2)>120,:,:));

mmshuf=nanmean(recip_shuf,2);
stdshuf=nanstd(recip_shuf,0,2);
z=(recip_idx(:,1)-mmshuf)./stdshuf;
%%figure 3
fh=figure('Color','w','Position',[100,100,175,235]);
hold on;
sh=scatter(mmshuf(z>1.96),recip_idx(z>1.96,1),'MarkerFaceAlpha',0.5,'MarkerFaceColor','r','MarkerEdgeColor','none');
sh=scatter(mmshuf(z<-1.96),recip_idx(z<-1.96,1),'MarkerFaceAlpha',0.5,'MarkerFaceColor','b','MarkerEdgeColor','none');
sh=scatter(mmshuf(z>=-1.96 & z<=1.96),recip_idx(z>=-1.96 & z<=1.96,1),'MarkerFaceAlpha',0.5,'MarkerFaceColor','k','MarkerEdgeColor','none');

plot([0,1],[0,1],'k:')
xlim([0,1]);
ylim([0,1]);
set(gca,'YTick',0:0.5:1)
xlabel('shuffled reciprocal fraction')
ylabel('recorded reciprocal fraction')
exportgraphics(fh,'recip_neuron_region.pdf')

figure('Color','w','Position',[100,100,50,235])
hold on
cic=bootci(1000,@(x) mean(x), mmshuf);
cid=bootci(1000,@(x) mean(x), recip_idx(:,1));

mm=mean([mmshuf,recip_idx(:,1)]);
hc=bar(1,mm(1),'FaceColor','k');
hd=bar(2,mm(2),'FaceColor','w');

errorbar(1:2,mean([mmshuf,recip_idx(:,1)]),...
[cic(1),cid(1)]-mm,[cic(2),cid(2)]-mm,'k.')

xlim([0.5,2.5]);
ylim([0,0.6]);

set(gca,'XTick',1:2,'XTickLabel',{'shuffled','recorded'},'XTickLabelRotation',90,'YTick',0:0.3:0.6)
ylabel('reciprocal fraction')
exportgraphics(fh,'recip_neuron_region_bar.pdf')

signrank(mmshuf,recip_idx(:,1))
