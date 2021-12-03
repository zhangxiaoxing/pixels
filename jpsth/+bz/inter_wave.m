[hier_stats,fh,bh]=conn_prob_bars_hier(sig,pair);

figure();
hold on;
pidx=0;
for fn=["same_stats","l2h_stats","h2l_stats"]
    wavemat=hier_stats.(fn).congr_wave;
    [phat36,pci36]=binofit(sum(wavemat(5:6,5)),sum(wavemat(5:6,6)));
    [phat3,pci3]=binofit(sum(wavemat(1:2,5)),sum(wavemat(1:2,6)));
    [phat6,pci6]=binofit(sum(wavemat(3:4,5)),sum(wavemat(3:4,6)));
    
    bar((1:5)+pidx,[phat36,phat3,phat6,hier_stats.(fn).congr_overlapped_wave(1),hier_stats.(fn).congr_inter_wave(1)].*100)
    plot([1,5]+pidx,[1,1].*hier_stats.(fn).congr(1).*100,'--r')
    pidx=pidx+6;
end