function runlength=ssloop_time_const(ring_replay,opt)
arguments
    ring_replay= []
    opt.skip_save = true
    opt.skip_plot = false
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
end
if isempty(ring_replay)
    switch opt.criteria
        case 'WT'
            load(fullfile('binary','motif_replay.mat'),'ring_replay');
        case 'Learning'
            load(fullfile('binary','LN_motif_replay.mat'),'ring_replay');
        otherwise
            keyboard()
    end
end

runlength=[];
for rii=1:size(ring_replay,1)
    trl_align=ring_replay.trl_align{rii};
    pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
    runlength=[runlength;cellfun(@(x) x(end)-x(1), ring_replay.ts{rii}(pref_delay))./30];
end
if false
    loop.per_sequence_run_length_in_ms=runlength;
    fid=fopen(fullfile('binary','upload','F2M_loop_activity_runlength.json'),'w');
    fprintf(fid,jsonencode(loop));
    fclose(fid);
end
singlehist=histcounts(runlength,[0:9,10:10:100],'Normalization','pdf');
if ~opt.skip_plot
    fh=figure();
    hold on
    sh=plot([0.5:9.5,15:10:95],singlehist,'-k');
    set(gca(),'XScale','log','YScale','log');
    ylim([1e-6,1]);
    xlim([1,500])
    xlabel('Time (ms)');
    ylabel('Probability density');
    qtrs=prctile(runlength,[1,25,50,75,99]);
    xline(qtrs,'--k',["1% ","25% ","50% ","75% ","99% "]+string(num2cell(qtrs)));
    title('single spike loops')
    if ~opt.skip_save
        savefig(fh,fullfile('binary','ssloop_time_constant.fig'));
    end
end
end