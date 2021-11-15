fpath=fullfile('..','allensdk','proj_mat_hier.hdf5');
mob_ccfids=h5read(fpath,'/mob_targets');
mob_mat=h5read(fpath,'/mob_matrix');
mop_ccfids=h5read(fpath,'/mop_src');
mop_mat=h5read(fpath,'/mop_matrix');


idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));

mob_reg=arrayfun(@(x) idmap.ccfid2reg(x),mob_ccfids);
mop_reg=arrayfun(@(x) idmap.ccfid2reg(x),mop_ccfids);

[~,~,ratiomap]=ref.get_pv_sst();
ffrac.collection=ephys.per_region_fraction('memtype','any');

corr_reg_mob=cell(0);
corr_mat_mob=[]; % MOB density, pv ratio, selective frac
for ri=1:numel(mob_reg)
    if ratiomap.isKey(mob_reg{ri}) && any(strcmp(mob_reg{ri},ffrac.collection(:,2))) && ffrac.collection{strcmp(mob_reg{ri},ffrac.collection(:,2)),4}>50
        corr_reg_mob=[corr_reg_mob;mob_reg{ri}];
        corr_mat_mob=[corr_mat_mob;mob_mat(ri),ratiomap(mob_reg{ri}),ffrac.collection{strcmp(mob_reg{ri},ffrac.collection(:,2)),1}];
    end
end

fh=figure('Color','w','Position',[100,100,800,400]);
subplot(1,2,1);
scatter(corr_mat_mob(:,1),corr_mat_mob(:,2),4,'ro','MarkerFaceColor','r')
text(corr_mat_mob(:,1),corr_mat_mob(:,2),corr_reg_mob,'HorizontalAlignment','center','VerticalAlignment','top');
set(gca(),'XScale','log')
xlabel('Projection density from OB');
ylabel('PV/(PV+SST)');
[r,p]=corr(corr_mat_mob(:,1),corr_mat_mob(:,2),'type','Spearman');
text(max(xlim()),max(ylim()),sprintf('r=%.2f,p=%.2f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');

subplot(1,2,2);
scatter(corr_mat_mob(:,1),corr_mat_mob(:,3),4,'ro','MarkerFaceColor','r')
text(corr_mat_mob(:,1),corr_mat_mob(:,3),corr_reg_mob,'HorizontalAlignment','center','VerticalAlignment','top');
set(gca(),'XScale','log')
xlabel('Projection density from OB');
ylabel('Memory neuron proportion');
[r,p]=corr(corr_mat_mob(:,1),corr_mat_mob(:,3),'type','Spearman');
text(max(xlim()),max(ylim()),sprintf('r=%.2f,p=%.2f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');


corr_reg_mop=cell(0);
corr_mat_mop=[]; % MOB density, pv ratio, selective frac
for ri=1:numel(mop_reg)
    if ratiomap.isKey(mop_reg{ri}) && any(strcmp(mop_reg{ri},ffrac.collection(:,2))) && ffrac.collection{strcmp(mop_reg{ri},ffrac.collection(:,2)),4}>50
        corr_reg_mop=[corr_reg_mop;mop_reg{ri}];
        corr_mat_mop=[corr_mat_mop;mop_mat(ri),mob_mat(strcmp(mop_reg{ri},mob_reg)),ratiomap(mop_reg{ri}),ffrac.collection{strcmp(mop_reg{ri},ffrac.collection(:,2)),1}];
    end
end


fh=figure('Color','w','Position',[100,100,800,400]);
subplot(1,2,1);
hold on
logmat=log(corr_mat_mop(:,1:2));
logmat(:,3)=1;
for ii=1:size(logmat,1)
scatter(logmat(ii,1),logmat(ii,2),9,'o','MarkerFaceColor',ephys.getRegColor(corr_reg_mop{ii}),'MarkerEdgeColor','none')
text(logmat(ii,1),logmat(ii,2),corr_reg_mop{ii},'HorizontalAlignment','center','VerticalAlignment','top','Color',ephys.getRegColor(corr_reg_mop{ii}));
% set(gca(),'XScale','log','YScale','log')
end
xlabel('Log projection density to M1');
ylabel('Log projection density from OB');
[r,p]=corr(logmat(:,1),logmat(:,2),'type','Pearson');
text(max(xlim()),max(ylim()),sprintf('r=%.2f,p=%.2f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');
regres=logmat(:,[1,3])\logmat(:,2);
plot(xlim(),xlim().*regres(1)+regres(2),'--k');

theta=atan(regres(1));
xlim([-14,-2])
ylim([-14,-2])
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

rotmat=logmat(:,1:2)*R;
rotmat=rotmat-mean(rotmat);
subplot(1,2,2)
hold on
for ii=1:size(rotmat,1)
scatter(rotmat(ii,1),rotmat(ii,2),4,'o','MarkerFaceColor',ephys.getRegColor(corr_reg_mop{ii}),'MarkerEdgeColor','none');
text(rotmat(ii,1),rotmat(ii,2),corr_reg_mop{ii},'HorizontalAlignment','center','VerticalAlignment','top','Color',ephys.getRegColor(corr_reg_mop{ii}));
end
yline(0,'--k');
xlabel('OB-M1 hierarchy index');
ylabel('Regression residual')
xlim([-7,7])
ylim([-7,7])

OBM1map=containers.Map('KeyType','char','ValueType','double');
for ri=1:numel(corr_reg_mop)
    reg=corr_reg_mop{ri};
    OBM1map(reg)=rotmat(ri,1);
end
save('OBM1Map.mat','OBM1map')
exportgraphics(fh,'OBM1Map.pdf','ContentType','vector')

fh=figure('Color','w','Position',[100,100,800,400]);
subplot(1,2,1)
scatter(corr_mat_mop(:,1),corr_mat_mop(:,3),9,'ro','MarkerFaceColor','r')
text(corr_mat_mop(:,1),corr_mat_mop(:,3),corr_reg_mop,'HorizontalAlignment','center','VerticalAlignment','top');
set(gca(),'XScale','log')
xlabel('Projection density to M1');
ylabel('PV/(PV+SST)');
[r,p]=corr(corr_mat_mop(:,1),corr_mat_mop(:,3),'type','Spearman');
text(max(xlim()),max(ylim()),sprintf('r=%.2f,p=%.2f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');
subplot(1,2,2)
scatter(corr_mat_mop(:,1),corr_mat_mop(:,4),4,'ro','MarkerFaceColor','r')
text(corr_mat_mop(:,1),corr_mat_mop(:,4),corr_reg_mop,'HorizontalAlignment','center','VerticalAlignment','top');
set(gca(),'XScale','log')
xlabel('Projection density to M1');
ylabel('Memory neuron proportion');
[r,p]=corr(corr_mat_mop(:,1),corr_mat_mop(:,4),'type','Spearman');
text(max(xlim()),max(ylim()),sprintf('r=%.2f,p=%.2f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');



fh=figure('Color','w','Position',[100,100,800,400]);
subplot(1,2,1)
scatter(rotmat(:,1),corr_mat_mop(:,3),4,'ro','MarkerFaceColor','r')
text(rotmat(:,1),corr_mat_mop(:,3),corr_reg_mop,'HorizontalAlignment','center','VerticalAlignment','top');
% set(gca(),'XScale','log')
xlabel('OB-M1 hierarchy index');
ylabel('PV/(PV+SST)');
[r,p]=corr(rotmat(:,1),corr_mat_mop(:,3),'type','Spearman');
text(max(xlim()),max(ylim()),sprintf('r=%.2f,p=%.2f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');

subplot(1,2,2)
scatter(rotmat(:,1),corr_mat_mop(:,4),4,'ro','MarkerFaceColor','r')
text(rotmat(:,1),corr_mat_mop(:,4),corr_reg_mop,'HorizontalAlignment','center','VerticalAlignment','top');
% set(gca(),'XScale','log')
xlabel('OB-M1 hierarchy index');
ylabel('Memory neuron proportion');
[r,p]=corr(rotmat(:,1),corr_mat_mop(:,4),'type','Spearman');
text(max(xlim()),max(ylim()),sprintf('r=%.2f,p=%.2f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');


