function runlength=sschain_time_const(chain_replay,opt)
arguments
    chain_replay= []
    opt.skip_save = true
    opt.skip_plot = false
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
end

% Learning implementation WIP
if isempty(chain_replay)
    load(fullfile('binary','motif_replay.mat'),'chain_replay');
end
runlength=[];
for cii=1:size(chain_replay,1)
    trl_align=chain_replay.trl_align{cii};
    pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
    runlength=[runlength;arrayfun(@(x) diff(chain_replay.ts{cii}(x,[1 end]),1,2), find(pref_delay))./30];
end

if false
    chain.per_sequence_run_length_in_ms=runlength;
    fid=fopen(fullfile('binary','upload','F2J_chain_activity_runlength.json'),'w');
    fprintf(fid,jsonencode(chain));
    fclose(fid);
end

singlehist=histcounts(runlength,[0:9,10:10:100],'Normalization','pdf');
if ~opt.skip_plot
    fh=figure();
    hold on
    sh=plot([0.5:9.5,15:10:95],singlehist,'-k');
    set(gca(),'XScale','log','YScale','log');
    ylim([1e-6,1]);
    xlim([1,500]);
    xlabel('Time (ms)');
    ylabel('Probability density');
    qtrs=prctile(runlength,[1,25,50,75,99]);
    xline(qtrs,'--k',["1% ","25% ","50% ","75% ","99% "]+string(num2cell(qtrs)));
    title('single spike chains')
    if ~opt.skip_save
        savefig(fh,fullfile('binary','sschain_time_constant.fig'));
    end
end
end