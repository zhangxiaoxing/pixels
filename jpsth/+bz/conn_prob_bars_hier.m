% input data from
% [sig,pair]=bz.load_sig_pair('pair',true)

function [hier_stats,fh,bh]=conn_prob_bars_hier(sig,pair,opt)
arguments
    sig (1,1) struct
    pair (1,1) struct
    opt.dist (1,1) double = 5
    opt.per_region_within (1,1) logical = false
end

addpath('k:\code\align\')
sess_cnt=max(sig.sess);
same_stats=struct();
[same_stats.nm_nm,same_stats.congr,same_stats.incon,same_stats.mem_nm,same_stats.nm_mem]...
    =deal(nan(sess_cnt,1));
l2h_stats=same_stats;
h2l_stats=same_stats;
[~,sig_same,sig_h2l,sig_l2h]=bz.util.diff_at_level(sig.reg,'hierarchy',true);
[~,pair_same,pair_h2l,pair_l2h]=bz.util.diff_at_level(pair.reg,'hierarchy',true);

sess_sig_type=sig.mem_type(sig_same(:,opt.dist),:);
sess_pair_type=pair.mem_type(pair_same(:,opt.dist),:);
same_stats=bz.bars_util.get_ratio(sess_sig_type,sess_pair_type);

%% %%%%%%%%%per region local FC%%%%%%%%%%
if opt.per_region_within
    same_reg_pair=squeeze(pair.reg(pair_same(:,opt.dist),opt.dist,:));
    ureg=unique(same_reg_pair(:,1));
    for ridx=1:numel(ureg)
        onereg=ureg(ridx);
        reg_sig_type=sig.mem_type(sig_same(:,opt.dist) & sig.reg(:,opt.dist,1)==onereg,:);
        reg_pair_type=pair.mem_type(pair_same(:,opt.dist) & pair.reg(:,opt.dist,1)==onereg,:);
        reg_stats(ridx)=bz.bars_util.get_ratio(reg_sig_type,reg_pair_type);
    end
    [~,~,ratiomap]=ref.get_pv_sst();
    idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));
    
    congru_corr=cell2mat(arrayfun(@(x) [ratiomap(char(idmap.ccfid2reg(ureg(x)))),reg_stats(x).congr(1)],(1:numel(ureg)).','UniformOutput',false));
    incong_corr=cell2mat(arrayfun(@(x) [ratiomap(char(idmap.ccfid2reg(ureg(x)))),reg_stats(x).incon(1)],(1:numel(ureg)).','UniformOutput',false));
    nmnm_corr=cell2mat(arrayfun(@(x) [ratiomap(char(idmap.ccfid2reg(ureg(x)))),reg_stats(x).nm_nm(1)],(1:numel(ureg)).','UniformOutput',false));
    
    figure('Color','w')
    hold on;
    ch=scatter(congru_corr(:,1),congru_corr(:,2),100,'.','r');
    ih=scatter(incong_corr(:,1),incong_corr(:,2),100,'.','b');
    nh=scatter(nmnm_corr(:,1),nmnm_corr(:,2),100,'.','k');
    xlabel('Hierarchy Index');
    ylabel('Coupling rate');
    legend([ch,ih,nh],{'Congruent','Incongruent','Nonmemory'});
end

%h2l
sess_sig_type=sig.mem_type(sig_h2l(:,opt.dist),:);
sess_pair_type=pair.mem_type(pair_h2l(:,opt.dist),:);
h2l_stats=bz.bars_util.get_ratio(sess_sig_type,sess_pair_type);

%l2h
sess_sig_type=sig.mem_type(sig_l2h(:,opt.dist),:);
sess_pair_type=pair.mem_type(pair_l2h(:,opt.dist),:);
l2h_stats=bz.bars_util.get_ratio(sess_sig_type,sess_pair_type);

hier_stats=struct('same_stats',same_stats,'l2h_stats',l2h_stats,'h2l_stats',h2l_stats);
assignin('base','hier_stats',hier_stats);

fh=figure('Color','w','Position',[32,32,235,235]);
hold on
bh=bar([same_stats.congr(1),same_stats.incon(1),same_stats.nm_nm(1);...
    l2h_stats.congr(1),l2h_stats.incon(1),l2h_stats.nm_nm(1);...
    h2l_stats.congr(1),h2l_stats.incon(1),h2l_stats.nm_nm(1)].*100);
[ci1,ci2]=deal([]);
for f=["congr","incon","nm_nm"]
    ci1=[ci1,cellfun(@(x) x.(f)(2),{same_stats,l2h_stats,h2l_stats})];
    ci2=[ci2,cellfun(@(x) x.(f)(3),{same_stats,l2h_stats,h2l_stats})];
end
errorbar([bh.XEndPoints],[bh.YEndPoints],ci1.*100-[bh.YEndPoints],ci2.*100-[bh.YEndPoints],'k.');


bh(1).FaceColor='w';
bh(2).FaceColor=[0.5,0.5,0.5];
bh(3).FaceColor='k';

legend(bh,{'Same memory','Diff. memory','Non-memory'})
set(gca(),'XTick',1:3,'XTickLabel',{'Within reg.','Low->High','High->Low'},'XTickLabelRotation',30)
ylabel('Func. coupling probability (%)');
exportgraphics(fh,'conn_prob_bars_hier.pdf');

%% chisq test
chisq_3(hier_stats.same_stats.congr(4),hier_stats.same_stats.congr(5),...
    hier_stats.same_stats.incon(4),hier_stats.same_stats.incon(5),...
hier_stats.same_stats.nm_nm(4),hier_stats.same_stats.nm_nm(5))

chisq_3(hier_stats.l2h_stats.congr(4),hier_stats.l2h_stats.congr(5),...
    hier_stats.l2h_stats.incon(4),hier_stats.l2h_stats.incon(5),...
hier_stats.l2h_stats.nm_nm(4),hier_stats.l2h_stats.nm_nm(5))

chisq_3(hier_stats.h2l_stats.congr(4),hier_stats.h2l_stats.congr(5),...
    hier_stats.h2l_stats.incon(4),hier_stats.h2l_stats.incon(5),...
hier_stats.h2l_stats.nm_nm(4),hier_stats.h2l_stats.nm_nm(5))


chisq_3(hier_stats.same_stats.congr(4),hier_stats.same_stats.congr(5),...
    hier_stats.l2h_stats.congr(4),hier_stats.l2h_stats.congr(5),...
hier_stats.h2l_stats.congr(4),hier_stats.h2l_stats.congr(5))


chisq_3(hier_stats.same_stats.incon(4),hier_stats.same_stats.incon(5),...
    hier_stats.l2h_stats.incon(4),hier_stats.l2h_stats.incon(5),...
hier_stats.h2l_stats.incon(4),hier_stats.h2l_stats.incon(5))


chisq_3(hier_stats.same_stats.nm_nm(4),hier_stats.same_stats.nm_nm(5),...
    hier_stats.l2h_stats.nm_nm(4),hier_stats.l2h_stats.nm_nm(5),...
hier_stats.h2l_stats.nm_nm(4),hier_stats.h2l_stats.nm_nm(5))

end

