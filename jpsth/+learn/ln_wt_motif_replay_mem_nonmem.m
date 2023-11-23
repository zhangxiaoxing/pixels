[wt_c_fh,wt_c_yy,wt_c_gg]=wave.motif_freq_mem_vs_nonmem(type='chain',criteria='WT');
[wt_l_fh,wt_l_yy,wt_l_gg]=wave.motif_freq_mem_vs_nonmem(type='loop',criteria='WT');

[ln_c_fh,ln_c_yy,ln_c_gg]=wave.motif_freq_mem_vs_nonmem(type='chain',criteria='Learning');
[ln_l_fh,ln_l_yy,ln_l_gg]=wave.motif_freq_mem_vs_nonmem(type='loop',criteria='Learning');


bcix=@(x,y,z) bootci(1000, @median, x(isfinite(x) & y.'==z));

fini_splitapply=@(f,x,y) splitapply(f,x(isfinite(x)),y(isfinite(x)));

wtcmm=fini_splitapply(@median,wt_c_yy,wt_c_gg.');
lncmm=fini_splitapply(@median,ln_c_yy,ln_c_gg.');

wtlmm=fini_splitapply(@median,wt_l_yy,wt_l_gg.');
lnlmm=fini_splitapply(@median,ln_l_yy,ln_l_gg.');

lncci=cell2mat(arrayfun(@(x) bcix(ln_c_yy,ln_c_gg,x),1:10,'UniformOutput',false));
wtcci=cell2mat(arrayfun(@(x) bcix(wt_c_yy,wt_c_gg,x),1:10,'UniformOutput',false));

lnlci=cell2mat(arrayfun(@(x) bcix(ln_l_yy,ln_l_gg,x),1:10,'UniformOutput',false));
wtlci=cell2mat(arrayfun(@(x) bcix(wt_l_yy,wt_l_gg,x),1:10,'UniformOutput',false));

mem_c_p=arrayfun(@(x) ranksum(wt_c_yy(wt_c_gg==x),ln_c_yy(ln_c_gg==x)),1:2:10);
nonmem_c_p=arrayfun(@(x) ranksum(wt_c_yy(wt_c_gg==x),ln_c_yy(ln_c_gg==x)),2:2:10);
mem_l_p=arrayfun(@(x) ranksum(wt_l_yy(wt_l_gg==x),ln_l_yy(ln_l_gg==x)),1:2:10);
nonmem_l_p=arrayfun(@(x) ranksum(wt_l_yy(wt_l_gg==x),ln_l_yy(ln_l_gg==x)),2:2:10);

p2str=@(p) subsref(["***","**","*","ns"],substruct('()',{find([p<0.001;p<0.005;p<0.05;true],1,'first')}));

fh=figure();
tiledlayout(2,2)
nexttile()
hold on
bh=bar([lncmm(1:2:10);wtcmm(1:2:10)].');
errorbar(bh(1).XEndPoints,bh(1).YEndPoints,lncci(1,1:2:10)-bh(1).YEndPoints,lncci(2,1:2:10)-bh(1).YEndPoints,'k.')
errorbar(bh(2).XEndPoints,bh(2).YEndPoints,wtcci(1,1:2:10)-bh(2).YEndPoints,wtcci(2,1:2:10)-bh(2).YEndPoints,'k.')
set(gca(),'XTick',1:5,'XTickLabel',{'Delay','NPDelay','ITI','Before','After'})
ylabel('Unit motif freq (Hz)')
text(1:5,mean(ylim())*ones(1,5),arrayfun(@(x) p2str(x),mem_c_p),'HorizontalAlignment','center')
title('memory-chain')

nexttile()
hold on
bh=bar([lncmm(2:2:10);wtcmm(2:2:10)].');
errorbar(bh(1).XEndPoints,bh(1).YEndPoints,lncci(1,2:2:10)-bh(1).YEndPoints,lncci(2,2:2:10)-bh(1).YEndPoints,'k.')
errorbar(bh(2).XEndPoints,bh(2).YEndPoints,wtcci(1,2:2:10)-bh(2).YEndPoints,wtcci(2,2:2:10)-bh(2).YEndPoints,'k.')
set(gca(),'XTick',1:5,'XTickLabel',{'Delay','NPDelay','ITI','Before','After'})
ylabel('Unit motif freq (Hz)')
text(1:5,mean(ylim())*ones(1,5),arrayfun(@(x) p2str(x),nonmem_c_p),'HorizontalAlignment','center')
title('nonmemory-chain')

nexttile()
hold on
bh=bar([lnlmm(1:2:10);wtlmm(1:2:10)].');
errorbar(bh(1).XEndPoints,bh(1).YEndPoints,lnlci(1,1:2:10)-bh(1).YEndPoints,lnlci(2,1:2:10)-bh(1).YEndPoints,'k.')
errorbar(bh(2).XEndPoints,bh(2).YEndPoints,wtlci(1,1:2:10)-bh(2).YEndPoints,wtlci(2,1:2:10)-bh(2).YEndPoints,'k.')
set(gca(),'XTick',1:5,'XTickLabel',{'Delay','NPDelay','ITI','Before','After'})
ylabel('Unit motif freq (Hz)')
text(1:5,mean(ylim())*ones(1,5),arrayfun(@(x) p2str(x),mem_l_p),'HorizontalAlignment','center')
title('memory-loop')

nexttile()
hold on
bh=bar([lnlmm(2:2:10);wtlmm(2:2:10)].');
errorbar(bh(1).XEndPoints,bh(1).YEndPoints,lnlci(1,2:2:10)-bh(1).YEndPoints,lnlci(2,2:2:10)-bh(1).YEndPoints,'k.')
errorbar(bh(2).XEndPoints,bh(2).YEndPoints,wtlci(1,2:2:10)-bh(2).YEndPoints,wtlci(2,2:2:10)-bh(2).YEndPoints,'k.')
set(gca(),'XTick',1:5,'XTickLabel',{'Delay','NPDelay','ITI','Before','After'})
ylabel('Unit motif freq (Hz)')
text(1:5,mean(ylim())*ones(1,5),arrayfun(@(x) p2str(x),nonmem_l_p),'HorizontalAlignment','center')
title('nonmemory-loop')

savefig(fh,fullfile("binary","ln_wt_motif_replay_mem_nonmem.fig"))





