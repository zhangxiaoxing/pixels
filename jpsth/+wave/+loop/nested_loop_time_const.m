function plot_nested_loop_runlength()
load(fullfile('binary','delay_iti_runlength_covered.mat'),'run_length')
fh=figure();
chained_loops_pdf=histcounts(run_length.delay(:,2),[0:19,20:20:500],'Normalization','pdf');
plot([0.5:19.5,30:20:490],chained_loops_pdf,'-k');
xlim([1,500])
qtrs=prctile(run_length.delay(:,2),[1,25,50,75,99]);
xline(qtrs,'--k',["1% ","25% ","50% ","75% ","99% "]+qtrs);
ylim([1e-6,1])
set(gca(),'XScale','log','YScale','log')
xlabel('Time (ms)')
ylabel('Probability density')
title('nest loops')
savefig(fh,fullfile('binary','nested_loop_time_constant.fig'));

if false
    fid=fopen(fullfile('binary','upload','F4G_Nested_loops_runlength.json'),'w');
    dout.nested_loops_runlength=run_length.delay(:,2);
    fprintf(fid,jsonencode(dout))
    fclose(fid)
end