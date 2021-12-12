% input data from
% [sig,pair]=bz.load_sig_pair('pair',true)

function [hier_stats,fh,bh]=conn_prob_bars_hier(sig,pair,opt)
arguments
    sig (1,1) struct
    pair (1,1) struct
    opt.dist (1,1) double = 5
    opt.per_region_within (1,1) logical = false
    opt.save_int_data (1,1) logical = false
    opt.load_int_data (1,1) logical = true
    opt.bin_edge=0:8:24;
end

addpath('k:\code\align\')
sess_cnt=max(sig.sess);
same_stats=struct();
[same_stats.nm_nm,same_stats.congr,same_stats.incon,same_stats.mem_nm,same_stats.nm_mem]...
    =deal(nan(sess_cnt,1));
if opt.load_int_data
    load('conn_prob_bars_hier.mat')
else
    l2h_stats=same_stats;
    h2l_stats=same_stats;
    [~,sig_same,sig_h2l,sig_l2h]=bz.util.diff_at_level(sig.reg,'hierarchy',true);
    [~,pair_same,pair_h2l,pair_l2h]=bz.util.diff_at_level(pair.reg,'hierarchy',true);

    %save intermediate data
    if opt.save_int_data
        save('conn_prob_bars_hier.mat',...
            'h2l_stats',...
            'l2h_stats',...
            'pair_h2l',...
            'pair_l2h',...
            'pair_same',...
            'same_stats',...
            'sess_cnt',...
            'sig_h2l',...
            'sig_l2h',...
            'sig_same')
    end
end

sig_wave_id=[ephys.get_wave_id(sig.sess,sig.suid(:,1)),ephys.get_wave_id(sig.sess,sig.suid(:,2))];
pair_wave_id=[ephys.get_wave_id(pair.sess,pair.suid(:,1)),ephys.get_wave_id(pair.sess,pair.suid(:,2))];
% 
sig_wave_TCOM=[ephys.getTCOM(sig.sess,sig.suid(:,1)),ephys.getTCOM(sig.sess,sig.suid(:,2))];
pair_wave_TCOM=[ephys.getTCOM(pair.sess,pair.suid(:,1)),ephys.getTCOM(pair.sess,pair.suid(:,2))];

sess_sig_type=sig.mem_type(sig_same(:,opt.dist),:);
sess_pair_type=pair.mem_type(pair_same(:,opt.dist),:);
sig_wave=sig_wave_id(sig_same(:,opt.dist),:);
pair_wave=pair_wave_id(pair_same(:,opt.dist),:);

same_stats=bz.bars_util.get_ratio(sess_sig_type,sess_pair_type);
wavestats=bz.bars_util.get_wave_ratio(sig_wave,pair_wave);
% 

sig_congr=all(ismember(sig.mem_type,1:2),2) | all(ismember(sig.mem_type,3:4),2);
pair_congr=all(ismember(pair.mem_type,1:2),2) | all(ismember(pair.mem_type,3:4),2);

sig_TCOM=sig_wave_TCOM(sig_same(:,opt.dist) & sig_congr,:);
pair_TCOM=pair_wave_TCOM(pair_same(:,opt.dist) & pair_congr,:);
TCOM_profile=bz.bars_util.get_TCOM_profile(sig_TCOM,pair_TCOM);

% figure();imagesc(same_binmats.binmat(:,:,1).*100,[1,4])
% colormap('turbo');
% colorbar()

same_stats=cell2struct([struct2cell(wavestats);struct2cell(same_stats);{TCOM_profile}],...
    [fieldnames(wavestats);fieldnames(same_stats);'deltaTCOM_FC_rate']);

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

sig_wave=sig_wave_id(sig_h2l(:,opt.dist),:);
pair_wave=pair_wave_id(pair_h2l(:,opt.dist),:);

h2l_stats=bz.bars_util.get_ratio(sess_sig_type,sess_pair_type);
wavestats=bz.bars_util.get_wave_ratio(sig_wave,pair_wave);

sig_TCOM=sig_wave_TCOM(sig_h2l(:,opt.dist) & sig_congr,:);
pair_TCOM=pair_wave_TCOM(pair_h2l(:,opt.dist) & pair_congr,:);
TCOM_profile=bz.bars_util.get_TCOM_profile(sig_TCOM,pair_TCOM);

h2l_stats=cell2struct([struct2cell(wavestats);struct2cell(h2l_stats);{TCOM_profile}],...
    [fieldnames(wavestats);fieldnames(h2l_stats);'deltaTCOM_FC_rate']);

%l2h
sess_sig_type=sig.mem_type(sig_l2h(:,opt.dist),:);
sess_pair_type=pair.mem_type(pair_l2h(:,opt.dist),:);
sig_wave=sig_wave_id(sig_l2h(:,opt.dist),:);
pair_wave=pair_wave_id(pair_l2h(:,opt.dist),:);
l2h_stats=bz.bars_util.get_ratio(sess_sig_type,sess_pair_type);
wavestats=bz.bars_util.get_wave_ratio(sig_wave,pair_wave);
sig_TCOM=sig_wave_TCOM(sig_l2h(:,opt.dist) & sig_congr,:);
pair_TCOM=pair_wave_TCOM(pair_l2h(:,opt.dist) & pair_congr,:);
TCOM_profile=bz.bars_util.get_TCOM_profile(sig_TCOM,pair_TCOM);

l2h_stats=cell2struct([struct2cell(wavestats);struct2cell(l2h_stats);{TCOM_profile}],...
    [fieldnames(wavestats);fieldnames(l2h_stats);'deltaTCOM_FC_rate']);

hier_stats=struct('same_stats',same_stats,'l2h_stats',l2h_stats,'h2l_stats',h2l_stats);
assignin('base','hier_stats',hier_stats);


%%
fh=figure('Color','w','Position',[32,32,300,150]);
subplot(1,3,1)
hold on
mm=[same_stats.congr(1);same_stats.incon(1);same_stats.nm_nm(1)].*100;
for jj=1:3
    bh(jj)=bar(jj,mm(jj));
end
[ci1,ci2]=deal([]);
for f=["congr","incon","nm_nm"]
    ci1=[ci1,cellfun(@(x) x.(f)(2),{same_stats})];
    ci2=[ci2,cellfun(@(x) x.(f)(3),{same_stats})];
end
errorbar([bh.XEndPoints],[bh.YEndPoints],ci1.*100-[bh.YEndPoints],ci2.*100-[bh.YEndPoints],'k.');
set(gca(),'XTick',[])

bh(1).FaceColor='w';
bh(2).FaceColor=[0.5,0.5,0.5];
bh(3).FaceColor='k';
ylabel('FC rate (%)')
%%
subplot(1,3,2)
hold on
mm=[l2h_stats.congr(1);l2h_stats.incon(1);l2h_stats.nm_nm(1)].*100;
for jj=1:3
    bh(jj)=bar(jj,mm(jj));
end
[ci1,ci2]=deal([]);
for f=["congr","incon","nm_nm"]
    ci1=[ci1,cellfun(@(x) x.(f)(2),{l2h_stats})];
    ci2=[ci2,cellfun(@(x) x.(f)(3),{l2h_stats})];
end
errorbar([bh.XEndPoints],[bh.YEndPoints],ci1.*100-[bh.YEndPoints],ci2.*100-[bh.YEndPoints],'k.');
set(gca(),'XTick',[])
bh(1).FaceColor='w';
bh(2).FaceColor=[0.5,0.5,0.5];
bh(3).FaceColor='k';
ylabel('FC rate (%)')
%%
subplot(1,3,3)
hold on
mm=[h2l_stats.congr(1);h2l_stats.incon(1);h2l_stats.nm_nm(1)].*100;
for jj=1:3
    bh(jj)=bar(jj,mm(jj));
end
[ci1,ci2]=deal([]);
for f=["congr","incon","nm_nm"]
    ci1=[ci1,cellfun(@(x) x.(f)(2),{h2l_stats})];
    ci2=[ci2,cellfun(@(x) x.(f)(3),{h2l_stats})];
end
errorbar([bh.XEndPoints],[bh.YEndPoints],ci1.*100-[bh.YEndPoints],ci2.*100-[bh.YEndPoints],'k.');
set(gca(),'XTick',[])
bh(1).FaceColor='w';
bh(2).FaceColor=[0.5,0.5,0.5];
bh(3).FaceColor='k';
ylim([0,1])
ylabel('FC rate (%)')
% 
% legend(bh,{'Same memory','Diff. memory','Non-memory'})
% set(gca(),'XTick',1:3,'XTickLabel',{'Within reg.','Low->High','High->Low'},'XTickLabelRotation',30)
% ylabel('Func. coupling probability (%)');
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

