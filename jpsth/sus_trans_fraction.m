sus_trans=h5read('../transient_6.hdf5','/sus_trans');
reg_list=h5read('../transient_6.hdf5','/reg');
cid_list=h5read('../transient_6.hdf5','/cluster_id');
path_list=h5read('../transient_6.hdf5','/path');
wrsp=h5read('../transient_6.hdf5','/wrs_p');

upath=unique(path_list);
spath=upath;
for i=1:length(upath)
    spath{i}=regexp(upath{i},'^.*?(?=\\)','match','once');
end
spath=unique(spath);

load reg_keep.mat

reg_sel=ismember(deblank(reg_list),reg_set(1:115));
sums=[];
for i=1:length(spath)
    sesssel=startsWith(path_list,spath{i});
    sus=nnz(sus_trans(:,1) & sesssel & reg_sel);
    trans=nnz((sus_trans(:,2) | sus_trans(:,4))& sesssel & reg_sel);
    cnt=nnz(sesssel & reg_sel);
    sums=[sums;cnt,sus,trans];
end

sums=sums(sums(:,1)>30,:);
susr=sums(:,2)./sums(:,1);
transr=sums(:,3)./sums(:,1);
susci=bootci(1000,@(x) mean(x),susr);
transci=bootci(1000,@(x) mean(x),transr);
fh=figure('Color','w','Position',[100,100,150,150]);
hold on
bar(1,mean(susr),'FaceColor','k','EdgeColor','k');
bar(2,mean(transr),'FaceColor','w','EdgeColor','k');
errorbar(1,mean(susr),susci(1)-mean(susr),susci(2)-mean(susr),'k.');
errorbar(2,mean(transr),transci(1)-mean(transr),transci(2)-mean(transr),'k.');
set(gca,'XTick',1:2,'XTickLabel',{'Sust.','Trans.'},'XTickLabelRotation',45);
xlim([0.5,2.5])
ylabel('Fraction out of 24444 neurons')
exportgraphics(fh,'sust_trans_ratio.pdf')
% 
% 
%     sesssel=startsWith(path_list,'M48_20191205_g0');
%     sus=nnz(sus_trans(:,1) & sesssel & reg_sel);
%     trans=nnz((sus_trans(:,2) | sus_trans(:,4))& sesssel & reg_sel);
%     cnt=nnz(sesssel & reg_sel);
%     sums=[sums;cnt,sus,trans];
% 
%     [GC,GR]=groupcounts(reg_list(sesssel & reg_sel))
