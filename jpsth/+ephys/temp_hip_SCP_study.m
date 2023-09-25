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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HIP_SCPF=scregs.from.HIP.DelaySP./scregs.from.HIP.DelaySec;
HIP_NPratio=SPratio.lead.HIP.NPDelaySP./SPratio.lead.HIP.NPDelaySPK;
HIP_ITIratio=SPratio.lead.HIP.ITISP./SPratio.lead.HIP.ITISPK;
NP_ITI_index=HIP_NPratio./HIP_ITIratio;

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
ylabel('NP frac/ITI frac')
ylim([0,4.1])

HIP_SCPF=scregs.to.HIP.DelaySP./scregs.to.HIP.DelaySec;
HIP_NPratio=SPratio.follow.HIP.NPDelaySP./SPratio.follow.HIP.NPDelaySPK;
HIP_ITIratio=SPratio.follow.HIP.ITISP./SPratio.follow.HIP.ITISPK;
NP_ITI_index=HIP_NPratio./HIP_ITIratio;

nexttile
scatter(HIP_SCPF,NP_ITI_index,16,'black','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
set(gca,'XScale','log')
xline(nanmedian(HIP_SCPF),'r--')
yline(0,'k:')
yline(nanmedian(NP_ITI_index),'r--')
title('To / Follow')
xlabel('SCP freq (Hz)')
ylabel('NP frac/ITI frac')
ylim([0,4.1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% HIP_SCPF=scregs.from.HIP.DelaySP./scregs.from.HIP.DelaySec;
HIP_NPratio=SPratio.lead.HIP.NPDelaySP./SPratio.lead.HIP.NPDelaySPK;
HIP_ITIratio=SPratio.lead.HIP.ITISP./SPratio.lead.HIP.ITISPK;
NP_ITI_index=(HIP_NPratio-HIP_ITIratio)./(HIP_NPratio+HIP_ITIratio);

from_hip_sel=strcmp(reg_pairs(:,1),'HIP');
toreg=reg_pairs(from_hip_sel,2);

[G,ureg]=findgroups(toreg);
mdg=splitapply(@nanmedian,NP_ITI_index,G);

[~,sidx]=sort(mdg);

figure()
hold on
boxplot(NP_ITI_index,toreg,'Colors','k','Whisker',inf,'GroupOrder',ureg(sidx));
yline(0,'r:')
ylabel('Coupled fraction index (NP-ITI/NP+ITI)')
title('SCP from HIP to various regions')


HIP_SCPF=scregs.from.HIP.DelaySP./scregs.from.HIP.DelaySec;
HIP_NPratio=SPratio.lead.HIP.NPDelaySP./SPratio.lead.HIP.NPDelaySPK;
HIP_ITIratio=SPratio.lead.HIP.ITISP./SPratio.lead.HIP.ITISPK;
NP_ITI_index=(HIP_NPratio-HIP_ITIratio)./(HIP_NPratio+HIP_ITIratio);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HIP_SCPF=scregs.to.HIP.DelaySP./scregs.to.HIP.DelaySec;
HIP_NPratio=SPratio.follow.HIP.NPDelaySP./SPratio.follow.HIP.NPDelaySPK;
HIP_ITIratio=SPratio.follow.HIP.ITISP./SPratio.follow.HIP.ITISPK;
NP_ITI_index=(HIP_NPratio-HIP_ITIratio)./(HIP_NPratio+HIP_ITIratio);

from_hip_sel=strcmp(reg_pairs(:,2),'HIP');
toreg=reg_pairs(from_hip_sel,1);

[G,ureg]=findgroups(toreg);
mdg=splitapply(@nanmedian,NP_ITI_index,G);

[~,sidx]=sort(mdg);

figure()
hold on
boxplot(NP_ITI_index,toreg,'Colors','k','Whisker',inf,'GroupOrder',ureg(sidx));
yline(0,'r:')
ylabel('Coupled fraction index (NP-ITI/NP+ITI)')
title('SCP from various regions to HIP')




