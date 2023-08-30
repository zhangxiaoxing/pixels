function ssloop_time_const(ring_replay,opt)
arguments
    ring_replay= []
    opt.skip_save = true
end
if isempty(ring_replay)
    load(fullfile('binary','motif_replay.mat'),'ring_replay');
end
runlength=[];
for rii=1:size(ring_replay,1)
    trl_align=ring_replay.trl_align{rii};
    pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
    runlength=[runlength;cellfun(@(x) x(end)-x(1), ring_replay.ts{rii}(pref_delay))./30];
end

singlehist=histcounts(runlength,[0:9,10:10:100],'Normalization','pdf');

fh=figure();
hold on
sh=plot([0.5:9.5,15:10:95],singlehist,'-k');
set(gca(),'XScale','log','YScale','log');
ylim([1e-6,1]);
xlim([1,500])
xlabel('Time (ms)');
ylabel('Probability density');
qtrs=prctile(runlength,[2.5,50,97.5]);
xline(qtrs,'--k',["2.5pct ","50pct ","97.5pct "]+string(num2cell(qtrs)));
title('single spike loops')
if ~opt.skip_save
    savefig(fh,fullfile('binary','ssloop_time_constant.fig'));
end
end