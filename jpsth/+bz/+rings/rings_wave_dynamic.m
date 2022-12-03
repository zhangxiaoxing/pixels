load('bzdata\sums_ring_stats_all.mat')
% ring in wave

% classify rings based on su-waveid combination

% tag ring spks with trial type, trial time

% histogram of loop activity duration
mm_dur=cell(1,3);

for rsize=1:3
    one_rsize=sums_all{rsize};
    for ridx=1:size(one_rsize,1)
        rfreq_str=one_rsize{ridx,5};
        mm_dur{rsize}=[mm_dur{rsize};mean(rfreq_str.durs)./30]; % spkTS unit, i.e. 1/30 ms
    end
end

xx=5:10:95;
edges=0:10:100;
figure()
hold on
plot(xx,histcounts(mm_dur{1},edges));
plot(xx,histcounts(mm_dur{2},edges))
plot(xx,histcounts(mm_dur{3},edges))

edges=0:10:200;
per_ring_hist=zeros(1,20);
for rsize=1:3
    one_rsize=sums_all{rsize};
    for ridx=1:size(one_rsize,1)
        rfreq_str=one_rsize{ridx,5};
        per_ring_hist=per_ring_hist+histcounts((rfreq_str.durs)./30,edges); % spkTS unit, i.e. 1/30 ms
    end
end

xx=5:10:195;
figure()
hold on
plot(xx,per_ring_hist./sum(per_ring_hist,'all'));
set(gca(),'YScale','log')
xlabel('Time (ms)')
ylabel('Probability')