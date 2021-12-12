keyboard();
[hier_stats,fh,bh]=bz.conn_prob_bars_hier(sig,pair,"bin_edge",0:12:24);


for fn=["same_stats","l2h_stats","h2l_stats"]
    fh=figure('Color','w','Position',[32,32,100,150]);
    hold on;
    pidx=0;
    wavemat=hier_stats.(fn).congr_wave;
    [phat36,pci36]=binofit(sum(wavemat(5:6,5)),sum(wavemat(5:6,6)));
    [phat3,pci3]=binofit(sum(wavemat(1:2,5)),sum(wavemat(1:2,6)));
    [phat6,pci6]=binofit(sum(wavemat(3:4,5)),sum(wavemat(3:4,6)));
    mm=[phat36,phat3,phat6,hier_stats.(fn).congr_overlapped_wave(1),hier_stats.(fn).congr_inter_wave(1)].*100;
    lbound=[pci36(1),pci3(1),pci6(1),hier_stats.(fn).congr_overlapped_wave(2),hier_stats.(fn).congr_inter_wave(2)].*100-mm;
    ubound=[pci36(2),pci3(2),pci6(2),hier_stats.(fn).congr_overlapped_wave(3),hier_stats.(fn).congr_inter_wave(3)].*100-mm;
    for jj=1:5
        bh(jj)=bar(jj,mm(jj),...
            'FaceColor','w','EdgeColor','k');
    end
    bh(2).FaceColor=[1,0.5,0.5];
    bh(3).FaceColor=[0.5,0.5,1];
    bh(4).FaceColor=[0.5,0.5,0.5];
    bh(5).FaceColor='k';

    errorbar((1:5)+pidx,mm,lbound,ubound,'k.')
    

    ph=plot([0.5,5.5]+pidx,[1,1].*hier_stats.(fn).congr(1).*100,'--r');
    pidx=pidx+6;
    ylabel('FC Rate (%)')
    if strcmp(fn,"h2l_stats"), ylim([0,2]);end
    p=chisqInter(hier_stats.(fn));
    text(mean(xlim()),max(ylim()),num2str(p),'HorizontalAlignment','center','VerticalAlignment','top')

%     set(gca,'XTick',[1:5,7:11,13:17],...
% %         'XTickLabel',cellstr(repmat(["3 and 6","3 only","6 only","partial coactive","indepedent"],1,3)),...
%     'XTickLabel',{"Wave 1","Wave 2","Wave 3","W1 and Wave2/3","indepedent"},...
%         'XTickLabelRotation',90)
%     legend([ph],{'Same-memory w/o classification'});
    exportgraphics(fh,sprintf('inter_wave_fc_%s.pdf',fn),'ContentType','vector');
end

fh=figure('Color','w','Position',[32,32,900,200]);
plotOneHist(hier_stats.same_stats,1,'Same region',[0,4])
plotOneHist(hier_stats.l2h_stats,2,'Olfactory to motor',[0,1.5])
plotOneHist(hier_stats.h2l_stats,3,'Motor to olfactory',[0,1.5])
exportgraphics(fh,'dTCOM_FC_rate.pdf','ContentType','vector')

% fh=figure('Color','w','Position',[32,32,900,200]);
% plotOneBar(hier_stats.same_stats,1,'Same region',[0,3])
% p=chisqOne(hier_stats.same_stats);
% text(mean(xlim()),max(ylim()),num2str(p),'HorizontalAlignment','center','VerticalAlignment','top')
% plotOneBar(hier_stats.l2h_stats,2,'Olfactory to Motor',[0,1.1])
% p=chisqOne(hier_stats.l2h_stats);
% text(mean(xlim()),max(ylim()),num2str(p),'HorizontalAlignment','center','VerticalAlignment','top')
% plotOneBar(hier_stats.h2l_stats,3,'Motor to Olfactory',[0,1.1])
% p=chisqOne(hier_stats.h2l_stats);
% text(mean(xlim()),max(ylim()),num2str(p),'HorizontalAlignment','center','VerticalAlignment','top')
% exportgraphics(fh,'per_bin_FC_bars.pdf','ContentType','vector')

% pHier=[chisqHier(hier_stats.same_stats),chisqHier(hier_stats.l2h_stats),chisqHier(hier_stats.h2l_stats)]

function plotOne(stats,idx,t)
subplot(1,3,idx)
hold on
imagesc(stats.per_bin_FC(:,:,1).*100,[0,4])
colormap('turbo');
cbh=colorbar();
set(gca,'XTick',0.5:1:3.5,'XTickLabel',0:2:6,'YTick',0.5:1:3.5,'YTickLabel',0:2:6)
ylabel('Leading neuron FRTC (s)')
xlabel('Following neuron FRTC (s)')
xlim([0.5,3.5]);
ylim([0.5,3.5]);
cbh.Label.String='FC Rate (%)';
title(t);
end


function plotOneBar(stats,idx,t,yspan)
subplot(1,3,idx)
hold on
mm=[stats.per_bin_FC(1,1,1),stats.per_bin_FC(2,2,1),stats.per_bin_FC(1,2,1),stats.per_bin_FC(2,1,1)].*100;
uci=[stats.per_bin_FC(1,1,2),stats.per_bin_FC(2,2,2),stats.per_bin_FC(1,2,2),stats.per_bin_FC(2,1,2)].*100;
bci=[stats.per_bin_FC(1,1,3),stats.per_bin_FC(2,2,3),stats.per_bin_FC(1,2,3),stats.per_bin_FC(2,1,3)].*100;
n=sum(stats.per_bin_FC(:,:,4),'all');
for jj=1:4
    bh(jj)=bar(jj,mm(jj),'FaceColor','w','EdgeColor','k');
end
bh(2).FaceColor='r';
bh(3).FaceColor='b';
bh(4).FaceColor='k';
errorbar(1:4,mm,uci-mm,bci-mm,'k.')
% set(gca(),'XTick',1:4,'XTickLabel',{'Early-early','Late-late','Early-late','Late-early'})
ylabel('FC Rate (%)');
ylim(yspan);
title(sprintf('%s,n=%d',t,n));
end


function plotOneHist(stats,idx,t,yspan)
subplot(1,3,idx)
hold on
for jj=1:size(stats.deltaTCOM_FC_rate,2)
    [phat(jj),pci(jj,:)]=binofit(stats.deltaTCOM_FC_rate(1,jj),stats.deltaTCOM_FC_rate(2,jj));
end
fill([stats.deltaTCOM_FC_rate(3,:),stats.deltaTCOM_FC_rate(3,end:-1:1)],[pci(:,1);flip(pci(:,2))].*100,'r','FaceAlpha',0.1,'EdgeColor','none')
plot(stats.deltaTCOM_FC_rate(3,:),phat.*100,'-r')
ns=sum(stats.deltaTCOM_FC_rate(1,:));
np=sum(stats.deltaTCOM_FC_rate(2,:));
ylabel('FC Rate (%)');
ylim(yspan);
title(sprintf('%s,%d of %d',t,ns,np));
set(gca,'XTick',-20:20:20,'XTickLabel',-5:5:5)
xline(0,'--k')
xlabel('follow-lead dTCOM')
end

function p=chisqOne(stats)
    vec1=[ones(stats.per_bin_FC(1,1,5),1);...
        2*ones(stats.per_bin_FC(1,2,5),1);...
        3*ones(stats.per_bin_FC(2,1,5),1);...
        4*ones(stats.per_bin_FC(2,2,5),1);...
        ];

    vec2=[(1:stats.per_bin_FC(1,1,5)).'>stats.per_bin_FC(1,1,4);...
        (1:stats.per_bin_FC(1,2,5)).'>stats.per_bin_FC(1,2,4);...
        (1:stats.per_bin_FC(2,1,5)).'>stats.per_bin_FC(2,1,4);...
        (1:stats.per_bin_FC(2,2,5)).'>stats.per_bin_FC(2,2,4);...
        ];
    [~,~,p]=crosstab(vec1,vec2);
end


function p=chisqInter(stats)
    vec1=[ones(sum(stats.congr_wave(5:6,6)),1);...
        2*ones(sum(stats.congr_wave(1:2,6)),1);...
        3*ones(sum(stats.congr_wave(3:4,6)),1);...
        4*ones(stats.congr_inter_wave(5),1);...
        5*ones(stats.congr_overlapped_wave(5),1);...
        ];

    vec2=[(1:sum(stats.congr_wave(5:6,6))).'>sum(stats.congr_wave(5:6,5));...
        (1:sum(stats.congr_wave(1:2,6))).'>sum(stats.congr_wave(1:2,5));...
        (1:sum(stats.congr_wave(3:4,6))).'>sum(stats.congr_wave(3:4,5));...
        (1:stats.congr_inter_wave(5)).'>stats.congr_inter_wave(4);...
        (1:stats.congr_overlapped_wave(5)).'>stats.congr_overlapped_wave(4)...
        ];
    [~,~,p]=crosstab(vec1,vec2);
end

function p=chisqHier(stats)

vec1=[2*ones(stats.per_bin_FC(1,2,5),1);...
    3*ones(stats.per_bin_FC(2,1,5),1);...
    ];

vec2=[...
    (1:stats.per_bin_FC(1,2,5)).'>stats.per_bin_FC(1,2,4);...
    (1:stats.per_bin_FC(2,1,5)).'>stats.per_bin_FC(2,1,4);...
    ];
[~,~,p]=crosstab(vec1,vec2);
end
