function [fh,mdl]=conn_prob_spatial_dist(sig,pair,opt)
arguments
    sig (1,1) struct
    pair (1,1) struct
    opt.dist (1,1) double = 5
    opt.pair_count = 1000
end

addpath('k:\code\align\')
sess_cnt=max(sig.sess);
same_stats=struct();
[same_stats.nm_nm,same_stats.congr,same_stats.incon,same_stats.mem_nm,same_stats.nm_mem]...
    =deal(nan(sess_cnt,1));
l2h_stats=same_stats;
h2l_stats=same_stats;
 [sig_same, sig_h2l, sig_l2h,pair_same, pair_h2l, pair_l2h]=hier.get_reg_hier_relation();
idmap=load(fullfile('..','align','reg_ccfid_map.mat'));

diff_reg_pair=squeeze(pair.reg(pair_h2l(:,opt.dist) | pair_l2h(:,opt.dist),opt.dist,:));
ureg=unique(diff_reg_pair(:,1));
same_stats=[];
for ridx=1:numel(ureg)
    pair_count=nnz(all(pair.reg(:,opt.dist,:)==ureg(ridx),3));
    if pair_count<opt.pair_count
        continue;
    end
    sig_count=nnz(all(sig.reg(:,opt.dist,:)==ureg(ridx),3));
    same_stats=[same_stats;double(ureg(ridx)),double(ureg(ridx)),0,sig_count,pair_count];
end

reg_comb=nchoosek(ureg,2);
dist_stats=[];
for ridx=1:size(reg_comb,1)
    pair_count=nnz(pair.reg(:,opt.dist,1)==reg_comb(ridx,1) & pair.reg(:,opt.dist,2)==reg_comb(ridx,2));
    if pair_count<opt.pair_count
        continue;
    end
    [avail,dist]=bz.get_spatial_dist(idmap.ccfid2reg(reg_comb(ridx,1)),idmap.ccfid2reg(reg_comb(ridx,2)));
    if ~avail
        continue
    end
    sig_count=nnz(sig.reg(:,opt.dist,1)==reg_comb(ridx,1) & sig.reg(:,opt.dist,2)==reg_comb(ridx,2));
    dist_stats=[dist_stats;double(reg_comb(ridx,1)),double(reg_comb(ridx,2)),dist,sig_count,pair_count];
end
[samemm,sameci]=binofit(sum(same_stats(:,4)),sum(same_stats(:,5)));
dist_sums=[zeros(size(samemm)),samemm,sameci-samemm];
Y=discretize(dist_stats(:,3)./100,0:9);
for di=reshape(unique(Y),1,[])
    disp(di)
    xx=di-0.5;
    dsel=Y==di;
    [mm,ci]=binofit(sum(dist_stats(dsel,4)),sum(dist_stats(dsel,5)));
    dist_sums=[dist_sums;xx*ones(size(mm)),mm,ci-mm];
end

dist_sums(:,2:end)=dist_sums(:,2:end).*100;
% [r,p]=corr([same_stats(:,3);dist_stats(:,3)]./100,[same_stats(:,4);dist_stats(:,4)]./[same_stats(:,5);dist_stats(:,5)].*100);

xx=[same_stats(:,3);dist_stats(:,3)].*10;% micro-meter unit
yy=[same_stats(:,4);dist_stats(:,4)]./[same_stats(:,5);dist_stats(:,5)].*100;
tbl=table(xx,yy);
save('fcrate_distance.mat','tbl')

% [fxy,gof] = fit(tbl.xx,tbl.yy,'exp1');
modelfun=@(b,x) b(1)*(x(:,1).^b(2))+b(3);
mdl=fitnlm(tbl,modelfun,[-0.5,0.5,2]);
disp(mdl)

fh=figure('Color','w','Position',[32,32,135,220]);
hold on;
mh=plot(dist_sums(:,1),dist_sums(:,2),'-k','LineWidth',1);
% fplot(@(x) mdl.Coefficients.Estimate(1).*(x.^mdl.Coefficients.Estimate(2))+mdl.Coefficients.Estimate(3),[0,7],'-k')
pltxy=sortrows([mdl.Variables.xx,mdl.Fitted]);
fith=plot(pltxy(:,1)./1000,pltxy(:,2),'-r','LineWidth',1); %micro-meter unit
errorbar(dist_sums(:,1),dist_sums(:,2),dist_sums(:,3),dist_sums(:,4),'k.');
scatter(same_stats(:,3)./100,same_stats(:,4)./same_stats(:,5).*100,4,...
    'o','MarkerFaceColor','k','MarkerFaceAlpha',0.4,'MarkerEdgeColor','none');
scatter(dist_stats(:,3)./100,dist_stats(:,4)./dist_stats(:,5).*100,4,...
    'o','MarkerFaceColor','k','MarkerFaceAlpha',0.4,'MarkerEdgeColor','none');
xlabel('Region distance (mm)');
ylabel('Coupling rate (%)');
xlim([-0.5,7])
ylim([0,3])
legend([mh,fith],{'Mean','Power law fit'},'Location','northoutside')
text(max(xlim()),max(ylim()),sprintf('%.2f,%.2f',sqrt(mdl.Rsquared.Ordinary),0),'HorizontalAlignment','right','VerticalAlignment','top');
keyboard()
exportgraphics(fh,'fc_rate_vs_spatial_dist.pdf')
