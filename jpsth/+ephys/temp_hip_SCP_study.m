% Temporary study on HIP SCP replay dynamics. Unlikely to be truly useful
% 2023.09.24

load(fullfile('binary','per_reg_SC_replay.mat'));

HIP_SCPF=scregs.from.HIP.DelaySP./scregs.from.HIP.DelaySec;
HIP_NPratio=SPratio.lead.HIP.NPDelaySP./SPratio.lead.HIP.NPDelaySPK;
HIP_ITIratio=SPratio.lead.HIP.ITISP./SPratio.lead.HIP.ITISPK;
NP_ITI_index=(HIP_NPratio-HIP_ITIratio)./(HIP_NPratio+HIP_ITIratio);

figure()
tiledlayout(1,2)
nexttile
scatter(HIP_SCPF,NP_ITI_index,16,'black','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
set(gca,'XScale','log')
xline(nanmedian(HIP_SCPF),'r--')
yline(0,'k:')
yline(nanmedian(NP_ITI_index),'r--')
title('From / Lead')
xlabel('SCP freq (Hz)')
ylabel('NP-ITI/NP+ITI')
ylim([-1,1]);

HIP_SCPF=scregs.to.HIP.DelaySP./scregs.to.HIP.DelaySec;
HIP_NPratio=SPratio.follow.HIP.NPDelaySP./SPratio.follow.HIP.NPDelaySPK;
HIP_ITIratio=SPratio.follow.HIP.ITISP./SPratio.follow.HIP.ITISPK;
NP_ITI_index=(HIP_NPratio-HIP_ITIratio)./(HIP_NPratio+HIP_ITIratio);

nexttile
scatter(HIP_SCPF,NP_ITI_index,16,'black','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
set(gca,'XScale','log')
xline(nanmedian(HIP_SCPF),'r--')
yline(0,'k:')
yline(nanmedian(NP_ITI_index),'r--')
title('To / Follow')
xlabel('SCP freq (Hz)')
ylabel('NP-ITI/NP+ITI')
ylim([-1,1]);