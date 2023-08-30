function sschain_time_const(chain_replay,opt)
arguments
    chain_replay= []
    opt.skip_save = true
end
if isempty(chain_replay)
    load(fullfile('binary','motif_replay.mat'),'chain_replay');
end
runlength=[];
for cii=1:size(chain_replay,1)
    trl_align=chain_replay.trl_align{cii};
    pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
    runlength=[runlength;arrayfun(@(x) diff(chain_replay.ts{cii}(x,[1 end]),1,2), find(pref_delay))./30];
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
title('single spike chains')
if ~opt.skip_save
    savefig(fh,fullfile('binary','sschain_time_constant.fig'));
end
end