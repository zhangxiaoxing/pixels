function fh=inter_wave(sens_meta)
split_diff_samp=false;
% selstr=load('perm_sens.mat','sens_meta');

% persistent sig pair

[sig,pair]=bz.load_sig_sums_conn_file('pair',true);
sig=bz.join_fc_waveid(sig,sens_meta.wave_id);
pair=bz.join_fc_waveid(pair,sens_meta.wave_id);
%>>>>>>>>>>>>>>>>>>Skip hierarchy>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
[sig_diff,sig_same,~,~]=bz.util.diff_at_level(sig.reg,'hierarchy',false);
[pair_diff,pair_same,~,~]=bz.util.diff_at_level(pair.reg,'hierarchy',false);

[samestats,samep]=statsOne(sig_same(:,5),pair_same(:,5),sig,pair,split_diff_samp);
[diffstats,diffp]=statsOne(sig_diff(:,5),pair_diff(:,5),sig,pair,split_diff_samp);

fh=figure('Color','w','Position',[32,32,600,225]);
t=tiledlayout(1,3);
same_ax=plotOne(samestats,'ylim',[0,0.055]);
diff_ax=plotOne(diffstats,'ylim',[0,0.02]);
title(t,'same-, diff- region FC rate')
th=nexttile();
ephys.util.figtable(fh,th,{'chisq-same';samep;'chisq-diff';diffp})

% exportgraphics(fh,'inter_wave_fc.pdf','ContentType','vector');
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%>>>>>>>>>>>>>>>>>>Hierarchy>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

[~,sig_same,sig_h2l,sig_l2h]=bz.util.diff_at_level(sig.reg,'hierarchy',true,'hiermap','AON','descend',true);
[~,pair_same,pair_h2l,pair_l2h]=bz.util.diff_at_level(pair.reg,'hierarchy',true,'hiermap','AON','descend',true);

[l2hstats,h2lstats,p]=statsTwo(sig_l2h(:,5),sig_h2l(:,5),pair_l2h(:,5),pair_h2l(:,5),sig,pair);
fh=figure('Color','w','Position',[32,32,440,250]);
nexttile()
plotTwo(l2hstats,h2lstats,'ylim',[0,0.02]);
th=nexttile();
ephys.util.figtable(fh,th,p)
title(t,'same-, diff- region FC rate')

end

function [data,stats]=statsOne(sig_sel,pair_sel,sig,pair,split)

% both
sig_both=nnz(sig_sel & (all(sig.waveid==5,2) | all(sig.waveid==6,2)));
pair_both=nnz(pair_sel & (all(pair.waveid==5,2) | all(pair.waveid==6,2)));
[bothhat,bothci]=binofit(sig_both,pair_both);
% common wave
sig_cowave=nnz(sig_sel & (all(ismember(sig.waveid,[1 5]),2)...
    | all(ismember(sig.waveid,[3 5]),2)...
    | all(ismember(sig.waveid,[2 6]),2)...
    | all(ismember(sig.waveid,[4 6]),2)));

pair_cowave=nnz(pair_sel & (all(ismember(pair.waveid,[1 5]),2)...
    | all(ismember(pair.waveid,[3 5]),2)...
    | all(ismember(pair.waveid,[2 6]),2)...
    | all(ismember(pair.waveid,[4 6]),2)));
[cowavehat,cowaveci]=binofit(sig_cowave,pair_cowave);
% same sample, diff wave
sig_cosample=nnz(sig_sel & ((any(sig.waveid==1,2) & any(sig.waveid==3,2))...
    | (any(sig.waveid==2,2) & any(sig.waveid==4,2))));
pair_cosample=nnz(pair_sel & ((any(pair.waveid==1,2) & any(pair.waveid==3,2))...
    | (any(pair.waveid==2,2) & any(pair.waveid==4,2))));
[cosamplehat,cosampleci]=binofit(sig_cosample,pair_cosample);
% nonmem
sig_nonmem=nnz(sig_sel & all(sig.waveid==0,2));
pair_nonmem=nnz(pair_sel & all(pair.waveid==0,2));
[nonmemhat,nonmemci]=binofit(sig_nonmem,pair_nonmem);

if split
    %diff sample, same duration
    sig_diffsample=nnz(sig_sel & ...
        ((any(ismember(sig.waveid,[1 5]),2) & any(ismember(sig.waveid,[2 6]),2)) ...
        |(any(ismember(sig.waveid,[3 5]),2) & any(ismember(sig.waveid,[4 6]),2))));

    pair_diffsample=nnz(pair_sel & ...
        ((any(ismember(pair.waveid,[1 5]),2) & any(ismember(pair.waveid,[2 6]),2)) ...
        |(any(ismember(pair.waveid,[3 5]),2) & any(ismember(pair.waveid,[4 6]),2))));

    [diffsample_samedur_hat,diffsample_samedur_ci]= ...
        binofit(sig_diffsample,pair_diffsample);

    % diff sample, diff duration
    sig_diffsample_diff_dur=nnz(sig_sel & (...
        (any(sig.waveid==1,2) & any(sig.waveid==4,2)) ...
        |(any(sig.waveid==3,2) & any(sig.waveid==2,2))));

    pair_diffsample_diff_dur=nnz(pair_sel & (...
        (any(pair.waveid==1,2) & any(pair.waveid==4,2)) ...
        |(any(pair.waveid==3,2) & any(pair.waveid==2,2))));

    [diffsample_diffdur_hat,diffsample_diffdur_ci]=binofit(sig_diffsample_diff_dur,pair_diffsample_diff_dur);

    data=[bothhat,bothci,cowavehat,cowaveci,cosamplehat,cosampleci,diffsample_samedur_hat,diffsample_samedur_ci,diffsample_diffdur_hat,diffsample_diffdur_ci,nonmemhat,nonmemci];

    [tbl,chi2,p]=crosstab([ones(pair_both,1);2*ones(pair_cowave,1);3*ones(pair_cosample,1);...
        4*ones(pair_diffsample,1);...
        5*ones(pair_diffsample_diff_dur,1);...
        6*ones(pair_nonmem,1)],...
        [(1:pair_both)>sig_both,(1:pair_cowave)>sig_cowave,(1:pair_cosample)>sig_cosample,...
        (1:pair_diffsample)>sig_diffsample, ...
        (1:pair_diffsample_diff_dur)>sig_diffsample_diff_dur, ...
        (1:pair_nonmem)>sig_nonmem].');

    stats=p;
else
    %diff sample
    sig_diffsample=nnz(sig_sel & ...
        (any(ismember(sig.waveid,[1 3 5]),2) & any(ismember(sig.waveid,[2 4 6]),2)));
    pair_diffsample=nnz(pair_sel & ...
        (any(ismember(pair.waveid,[1 3 5]),2) & any(ismember(pair.waveid,[2 4 6]),2))); ...
    [diffsample_samedur_hat,diffsample_samedur_ci]= ...
        binofit(sig_diffsample,pair_diffsample);

    data=[bothhat,bothci,cowavehat,cowaveci,cosamplehat,cosampleci,diffsample_samedur_hat,diffsample_samedur_ci,nonmemhat,nonmemci];

    [tbl,chi2,p]=crosstab([ones(pair_both,1);2*ones(pair_cowave,1);3*ones(pair_cosample,1);...
        4*ones(pair_diffsample,1);...
        5*ones(pair_nonmem,1)],...
        [(1:pair_both)>sig_both,(1:pair_cowave)>sig_cowave,(1:pair_cosample)>sig_cosample,...
        (1:pair_diffsample)>sig_diffsample, ...
        (1:pair_nonmem)>sig_nonmem].');
    stats=p;
end
end

function ax=plotOne(stats,opt)
arguments
    stats
    opt.ylim = []
end
ax=nexttile();
hold on
for ii=1:numel(stats)./3
    bh(ii)=bar(ii,stats(ii*3-2),'FaceColor','w','EdgeColor','k');
    errorbar(ii,stats(ii*3-2),diff(stats(ii*3+[-2,-1]),1,2),diff(stats(ii*3+[-2,0]),1,2),'k.')
end
bh(1).FaceColor=[1,0.5,0.5];
bh(2).FaceColor=[1,0.5,1];
bh(3).FaceColor=[0.5,0.5,1];
bh(4).FaceColor=[0.5,0.5,0.5];
if numel(bh)==6
    bh(5).FaceColor='w';
    bh(6).FaceColor='k';
else
    bh(5).FaceColor='k';
end
if ~isempty(opt.ylim), ylim(opt.ylim);end
ylabel('FC Rate (%)')
set(gca(),'XTick',[],'YTickLabel',100.*get(gca(),'YTick'));
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

