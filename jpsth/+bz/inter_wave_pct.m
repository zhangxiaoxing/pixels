function [fh,stats,t]=inter_wave_pct(pct_meta)
arguments
    pct_meta
end
% selstr=load('perm_sens.mat','sens_meta');

% persistent sig pair

[sig,pair]=bz.load_sig_sums_conn_file('pair',true);
sig=bz.join_fc_waveid(sig,pct_meta.wave_id);
pair=bz.join_fc_waveid(pair,pct_meta.wave_id);
%>>>>>>>>>>>>>>>>>>Skip hierarchy>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
[sig_diff,sig_same,~,~]=bz.util.diff_at_level(sig.reg,'hierarchy',false);
[pair_diff,pair_same,~,~]=bz.util.diff_at_level(pair.reg,'hierarchy',false);

[samestats,sameci]=statsOne(sig_same(:,5),pair_same(:,5),sig,pair);
[diffstats,diffci]=statsOne(sig_diff(:,5),pair_diff(:,5),sig,pair);
stats.samestats=samestats;
stats.sameci=sameci;
stats.diffstats=diffstats;
stats.diffci=diffci;


fh=figure('Color','w','Position',[32,32,1220,320]);
t=tiledlayout(1,3);
nexttile(2)
ax_same=plotOne(samestats,'scale',[min(samestats,[],"all"),max(samestats,[],"all")]);
nexttile(3)
ax_diff=plotOne(diffstats,'scale',[min(diffstats,[],"all"),max(diffstats,[],"all")]);
title(t,'Same-, cross-region FC rate')
% th=nexttile();
% ephys.util.figtable(fh,th,{'chisq-same';samep;'chisq-diff';diffp})

% exportgraphics(fh,'inter_wave_fc.pdf','ContentType','vector');
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end

function [cowavehat,cowaveci]=statsOne(sig_sel,pair_sel,sig,pair)
%% 4 Qdr FC rate loop
cowavehat=nan(4,4);
cowaveci=nan(4,4,2);
for fromqtr=1:4
    for toqtr=1:4
        sig_cowave=nnz(sig_sel & sig.waveid(:,1)==fromqtr & sig.waveid(:,2)==toqtr);
        pair_cowave=nnz(pair_sel & pair.waveid(:,1)==fromqtr & pair.waveid(:,2)==toqtr);
        [cowavehat(fromqtr,toqtr),cowaveci(fromqtr,toqtr,:)]=binofit(sig_cowave,pair_cowave);
    end
end
%% skipped statistics for now
% data=[cowavehat,cowaveci,crosswavehat,crosswaveci,nonmemhat,nonmemci];
% 
% [~,~,p]=crosstab([ones(pair_cowave,1); ...
%     4*ones(pair_cross_wave,1);...
%     5*ones(pair_nonmem,1)],...
%     [(1:(pair_cowave))>(sig_cowave),...
%     (1:pair_cross_wave)>sig_cross_wave, ...
%     (1:pair_nonmem)>sig_nonmem].');
% stats=p;
end

function ax=plotOne(stats,opt)
arguments
    stats
    opt.scale
end
ax=imagesc(stats,opt.scale);
cbh=colorbar();
cbh.Label.String='F.C. Rate (%)';
cbh.TickLabels=cbh.Ticks*100;
set(gca(),'YDir','normal')
colormap('turbo')
xlabel('From class # neuron');
ylabel('To class # neuron');
set(gca(),'XTick',1:4,'YTick',1:4)
end

function [data1,data2,stats]=statsTwo(sig_sel1,sig_sel2,pair_sel1,pair_sel2,sig,pair)
% >>>>>>>>>>>>>>>>>>>>>>>>both>>>>>>>>>>>>>>>>>>>>
sig_both1=nnz(sig_sel1 & (all(sig.waveid==5,2) | all(sig.waveid==6,2)));
pair_both1=nnz(pair_sel1 & (all(pair.waveid==5,2) | all(pair.waveid==6,2)));

sig_both2=nnz(sig_sel2 & (all(sig.waveid==5,2) | all(sig.waveid==6,2)));
pair_both2=nnz(pair_sel2 & (all(pair.waveid==5,2) | all(pair.waveid==6,2)));

[bothhat1,bothci1]=binofit(sig_both1,pair_both1);
[bothhat2,bothci2]=binofit(sig_both2,pair_both2);
% >>>>>>>>>>>>>>>>>>>>common wave>>>>>>>>>>>>>>>>
sig_sel_common=all(ismember(sig.waveid,[1 5]),2)...
    | all(ismember(sig.waveid,[3 5]),2)...
    | all(ismember(sig.waveid,[2 6]),2)...
    | all(ismember(sig.waveid,[4 6]),2);
sig_cowave1=nnz(sig_sel1 &sig_sel_common);
sig_cowave2=nnz(sig_sel2 &sig_sel_common);

pair_sel_common=all(ismember(pair.waveid,[1 5]),2)...
    | all(ismember(pair.waveid,[3 5]),2)...
    | all(ismember(pair.waveid,[2 6]),2)...
    | all(ismember(pair.waveid,[4 6]),2);
pair_cowave1=nnz(pair_sel1 & pair_sel_common);
pair_cowave2=nnz(pair_sel2 & pair_sel_common);

[cowavehat1,cowaveci1]=binofit(sig_cowave1,pair_cowave1);
[cowavehat2,cowaveci2]=binofit(sig_cowave2,pair_cowave2);

% >>>>>>>>>>>>>>>>>>>>>>> nonmem >>>>>>>>>>>>>>>>>>>>
sig_nonmem1=nnz(sig_sel1 & all(sig.waveid==0,2));
sig_nonmem2=nnz(sig_sel2 & all(sig.waveid==0,2));
pair_nonmem1=nnz(pair_sel1 & all(pair.waveid==0,2));
pair_nonmem2=nnz(pair_sel2 & all(pair.waveid==0,2));

[nonmemhat1,nonmemci1]=binofit(sig_nonmem1,pair_nonmem1);
[nonmemhat2,nonmemci2]=binofit(sig_nonmem2,pair_nonmem2);

data1=[bothhat1,bothci1,cowavehat1,cowaveci1,nonmemhat1,nonmemci1];
data2=[bothhat2,bothci2,cowavehat2,cowaveci2,nonmemhat2,nonmemci2];

%>>>>>>>>>>>>>>>>>>>>>> STATS >>>>>>>>>>>>>>>>>>>>
[~,~,p(1)]=crosstab([ones(pair_both1,1);2*ones(pair_both2,1)],...
    [(1:pair_both1)>sig_both1,(1:pair_both2)>sig_both2].');

[~,~,p(2)]=crosstab([ones(pair_cowave1,1);2*ones(pair_cowave2,1)],...
    [(1:pair_cowave1)>sig_cowave1,(1:pair_cowave2)>sig_cowave2].');

[~,~,p(3)]=crosstab([ones(pair_nonmem1,1);2*ones(pair_nonmem2,1)],...
    [(1:pair_nonmem1)>sig_nonmem1,(1:pair_nonmem2)>sig_nonmem2].');

stats=p;
end

function plotTwo(stats1,stats2,opt)
arguments
    stats1
    stats2
    opt.ylim = [0,0.02]
end

hold on
bh=bar([stats1((1:3)*3-2);stats2((1:3)*3-2)].','EdgeColor','k');
bh(1).FaceColor='w';
bh(2).FaceColor='k';

errorbar(bh(1).XEndPoints,bh(1).YEndPoints, ...
    diff([stats1(1:3:end);stats1(2:3:end)]), ...
    diff([stats1(1:3:end);stats1(3:3:end)]), ...
    'k.','CapSize',4)

errorbar(bh(2).XEndPoints,bh(2).YEndPoints, ...
    diff([stats2(1:3:end);stats2(2:3:end)]), ...
    diff([stats2(1:3:end);stats2(3:3:end)]), ...
    'k.','CapSize',4)

legend(bh,{'Upward','Downward'},'Location','northoutside','Orientation','horizontal');

if ~isempty(opt.ylim), ylim(opt.ylim);end
ylabel('F.C. Rate (%)')
set(gca(),'XTick',1:3,'XTickLabel',...
    {'both','same wave','nonmem' },...
    'XTickLabelRotation',90,...
    'YTickLabel',100.*get(gca(),'YTick'));
end

