function fh=plot_delay_vs_iti(run_length)
delay_pdf=histcounts(run_length.delay(:,2),[0:19,20:20:300],'Normalization','pdf');
iti_pdf=histcounts(run_length.iti(:,2),[0:19,20:20:300],'Normalization','pdf');
pre_post_pdf=histcounts(run_length.pre_post(:,2),[0:19,20:20:300],'Normalization','pdf');

fh=figure('Position',[500,100,400,300]);
hold on
dh=plot([0.5:19.5,30:20:290],delay_pdf,'-r');
ih=plot([0.5:19.5,30:20:290],iti_pdf,'-k');
oh=plot([0.5:19.5,30:20:290],pre_post_pdf,'-b');

xlim([1,500])

% qtrs=prctile(run_length,[10,50,90]);
% xline(qtrs,'--k',string(qtrs)) % 17 24 35

ylim([1e-7,0.1])
set(gca(),'XScale','log','YScale','log')
xlabel('Time (ms)')
ylabel('Probability density')
legend([dh,ih,oh],{'Delay','ITI','Surround'})
end

