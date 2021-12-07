keyboard();
[hier_stats,fh,bh]=bz.conn_prob_bars_hier(sig,pair,"bin_edge",0:12:24);

fh=figure('Color','w');
hold on;
pidx=0;
for fn=["same_stats","l2h_stats","h2l_stats"]
    wavemat=hier_stats.(fn).congr_wave;
    [phat36,pci36]=binofit(sum(wavemat(5:6,5)),sum(wavemat(5:6,6)));
    [phat3,pci3]=binofit(sum(wavemat(1:2,5)),sum(wavemat(1:2,6)));
    [phat6,pci6]=binofit(sum(wavemat(3:4,5)),sum(wavemat(3:4,6)));
    mm=[phat36,phat3,phat6,hier_stats.(fn).congr_overlapped_wave(1),hier_stats.(fn).congr_inter_wave(1)].*100;
    lbound=[pci36(1),pci3(1),pci6(1),hier_stats.(fn).congr_overlapped_wave(2),hier_stats.(fn).congr_inter_wave(2)].*100-mm;
    ubound=[pci36(2),pci3(2),pci6(2),hier_stats.(fn).congr_overlapped_wave(3),hier_stats.(fn).congr_inter_wave(3)].*100-mm;
    bar((1:5)+pidx,mm,...
        'FaceColor','w','EdgeColor','k')
    errorbar((1:5)+pidx,mm,lbound,ubound,'k.')

    ph=plot([0.5,5.5]+pidx,[1,1].*hier_stats.(fn).congr(1).*100,'--r');
    pidx=pidx+6;
end

ylabel('FC Rate (%)')
set(gca,'XTick',[1:5,7:11,13:17],...
    'XTickLabel',cellstr(repmat(["3 and 6","3 only","6 only","partial coactive","indepedent"],1,3)),...
    'XTickLabelRotation',90)
legend([ph],{'Same-memory w/o classification'});
exportgraphics(gcf,'inter_wave_fc.pdf','ContentType','vector');


fh=figure('Color','w','Position',[32,32,900,200]);
plotOneBar(hier_stats.same_stats,1,'Same region',[0,3])
plotOneBar(hier_stats.l2h_stats,2,'Olfactory to Motor',[0,1.1])
plotOneBar(hier_stats.h2l_stats,3,'Motor to Olfactory',[0,1.1])
exportgraphics(fh,'per_bin_FC_bars.pdf','ContentType','vector')

function plotOne(stats,idx,t)
subplot(1,3,idx)
hold on
imagesc(stats.per_bin_FC(:,:,1).*100,[0,4])
colormap('turbo');
cbh=colorbar();
set(gca,'XTick',0.5:1:3.5,'XTickLabel',0:2:6,'YTick',0.5:1:3.5,'YTickLabel',0:2:6)
ylabel('Leading neuron FRTC (s)')
xlabel('Following neuron FRTC (s)')
xlim([0.5,3.5]);
ylim([0.5,3.5]);
cbh.Label.String='FC Rate (%)';
title(t);
end


function plotOneBar(stats,idx,t,yspan)
subplot(1,3,idx)
hold on
mm=[stats.per_bin_FC(1,1,1),stats.per_bin_FC(2,2,1),stats.per_bin_FC(1,2,1),stats.per_bin_FC(2,1,1)].*100;
uci=[stats.per_bin_FC(1,1,2),stats.per_bin_FC(2,2,2),stats.per_bin_FC(1,2,2),stats.per_bin_FC(2,1,2)].*100;
bci=[stats.per_bin_FC(1,1,3),stats.per_bin_FC(2,2,3),stats.per_bin_FC(1,2,3),stats.per_bin_FC(2,1,3)].*100;
n=sum(stats.per_bin_FC(:,:,4),'all');
bar(mm,'FaceColor','w','EdgeColor','k')
errorbar(1:4,mm,uci-mm,bci-mm,'k.')
set(gca(),'XTick',1:4,'XTickLabel',{'Early-early','Late-late','Early-late','Late-early'})
ylabel('FC Rate (%)');
ylim(yspan);
title(sprintf('%s,n=%d',t,n));
end