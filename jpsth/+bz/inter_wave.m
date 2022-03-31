keyboard();
[sig,pair]=bz.load_sig_pair('pair',true);
[sig_same, sig_h2l, sig_l2h,pair_same, pair_h2l, pair_l2h]=hier.get_reg_hier_relation();

[samestats,samep]=statsOne(sig_same,pair_same,sig,pair);
[l2hstats,l2hp]=statsOne(sig_l2h,pair_l2h,sig,pair);
[h2lstats,h2lp]=statsOne(sig_h2l,pair_h2l,sig,pair);


fh=figure('Color','w','Position',[32,32,400,225]);
plotOne(samestats,1)
plotOne(l2hstats,2,'ylim',[0,0.02])
plotOne(h2lstats,3,'ylim',[0,0.02])
exportgraphics(fh,'inter_wave_fc.pdf','ContentType','vector');

return

function [data,stats]=statsOne(sig_sel,pair_sel,sig,pair)
% both
sig_both=nnz(sig_sel & (all(sig.waveid==5,2) | all(sig.waveid==6,2)));
pair_both=nnz(pair_sel & (all(pair.waveid==5,2) | all(pair.waveid==6,2)));
[bothhat,bothci]=binofit(sig_both,pair_both);
% common wave
sig_cowave=nnz(sig_sel & (all(ismember(sig.waveid,[1 5]),2)...
    | all(ismember(sig.waveid,[3 5]),2)...
    | all(ismember(sig.waveid,[2 6]),2)...
    | all(ismember(sig.waveid,[4 6]),2)));

pair_cowave=nnz(pair_sel & (all(ismember(pair.waveid,[1 5]),2)...
    | all(ismember(pair.waveid,[3 5]),2)...
    | all(ismember(pair.waveid,[2 6]),2)...
    | all(ismember(pair.waveid,[4 6]),2)));
[cowavehat,cowaveci]=binofit(sig_cowave,pair_cowave);
% same sample, diff wave
sig_cosample=nnz(sig_sel & ((any(sig.waveid==1,2) & any(sig.waveid==1,2))...
    | (any(sig.waveid==2,2) & any(sig.waveid==4,2))));
pair_cosample=nnz(pair_sel & ((any(pair.waveid==1,2) & any(pair.waveid==1,2))...
    | (any(pair.waveid==2,2) & any(pair.waveid==4,2))));
[cosamplehat,cosampleci]=binofit(sig_cosample,pair_cosample);

% diff sample
sig_diffsample=nnz(sig_sel & any(ismember(sig.waveid,[1 3 5]),2)...
    & any(ismember(sig.waveid,[2 4 6]),2));
pair_diffsample=nnz(pair_sel & any(ismember(pair.waveid,[1 3 5]),2)...
    & any(ismember(pair.waveid,[2 4 6]),2));
[diffsamplehat,diffsampleci]=binofit(sig_diffsample,pair_diffsample);

% nonmem
sig_nonmem=nnz(sig_sel & all(sig.waveid==0,2));
pair_nonmem=nnz(pair_sel & all(pair.waveid==0,2));
[nonmemhat,nonmemci]=binofit(sig_nonmem,pair_nonmem);

data=[bothhat,bothci,cowavehat,cowaveci,cosamplehat,cosampleci,diffsamplehat,diffsampleci,nonmemhat,nonmemci];

[tbl,chi2,p]=crosstab([ones(pair_both,1);2*ones(pair_cowave,1);3*ones(pair_cosample,1);4*ones(pair_diffsample,1);5*ones(pair_nonmem,1)],...
    [(1:pair_both)>sig_both,(1:pair_cowave)>sig_cowave,(1:pair_cosample)>sig_cosample,(1:pair_diffsample)>sig_diffsample,(1:pair_nonmem)>sig_nonmem].');
stats=p;
end


% 




function fh=plotOne(stats,sidx,opt)
arguments
    stats
    sidx
    opt.ylim = []
end
    subplot(1,3,sidx)
    hold on
    for ii=1:5
        bh(ii)=bar(ii,stats(ii*3-2),'FaceColor','w','EdgeColor','k');
        errorbar(ii,stats(ii*3-2),diff(stats(ii*3+[-2,-1]),1,2),diff(stats(ii*3+[-2,0]),1,2),'k.')
    end
    bh(2).FaceColor=[1,0.5,0.5];
    bh(3).FaceColor=[0.5,0.5,1];
    bh(4).FaceColor=[0.5,0.5,0.5];
    bh(5).FaceColor='k';
    if ~isempty(opt.ylim), ylim(opt.ylim);end
    ylabel('FC Rate (%)')
    set(gca(),'XTick',[],'YTickLabel',100.*get(gca(),'YTick'));
%     if strcmp(fn,"h2l_stats"), ylim([0,2]);end
%     p=chisqInter(hier_stats.(fn));
%     text(mean(xlim()),max(ylim()),num2str(p),'HorizontalAlignment','center','VerticalAlignment','top')
% 
% %     set(gca,'XTick',[1:5,7:11,13:17],...
% % %         'XTickLabel',cellstr(repmat(["3 and 6","3 only","6 only","partial coactive","indepedent"],1,3)),...
% %     'XTickLabel',{"Wave 1","Wave 2","Wave 3","W1 and Wave2/3","indepedent"},...
% %         'XTickLabelRotation',90)
% %     legend([ph],{'Same-memory w/o classification'});
%     exportgraphics(fh,sprintf('inter_wave_fc_%s.pdf',fn),'ContentType','vector');
end

% 
% fh=figure('Color','w','Position',[32,32,900,200]);
% plotOneHist(hier_stats.same_stats,1,'Same region',[0,4])
% plotOneHist(hier_stats.l2h_stats,2,'Olfactory to motor',[0,1.5])
% plotOneHist(hier_stats.h2l_stats,3,'Motor to olfactory',[0,1.5])
% exportgraphics(fh,'dTCOM_FC_rate.pdf','ContentType','vector')
% 
% % fh=figure('Color','w','Position',[32,32,900,200]);
% % plotOneBar(hier_stats.same_stats,1,'Same region',[0,3])
% % p=chisqOne(hier_stats.same_stats);
% % text(mean(xlim()),max(ylim()),num2str(p),'HorizontalAlignment','center','VerticalAlignment','top')
% % plotOneBar(hier_stats.l2h_stats,2,'Olfactory to Motor',[0,1.1])
% % p=chisqOne(hier_stats.l2h_stats);
% % text(mean(xlim()),max(ylim()),num2str(p),'HorizontalAlignment','center','VerticalAlignment','top')
% % plotOneBar(hier_stats.h2l_stats,3,'Motor to Olfactory',[0,1.1])
% % p=chisqOne(hier_stats.h2l_stats);
% % text(mean(xlim()),max(ylim()),num2str(p),'HorizontalAlignment','center','VerticalAlignment','top')
% % exportgraphics(fh,'per_bin_FC_bars.pdf','ContentType','vector')
% 
% % pHier=[chisqHier(hier_stats.same_stats),chisqHier(hier_stats.l2h_stats),chisqHier(hier_stats.h2l_stats)]
% 
% function plotOne(stats,idx,t)
% subplot(1,3,idx)
% hold on
% imagesc(stats.per_bin_FC(:,:,1).*100,[0,4])
% colormap('turbo');
% cbh=colorbar();
% set(gca,'XTick',0.5:1:3.5,'XTickLabel',0:2:6,'YTick',0.5:1:3.5,'YTickLabel',0:2:6)
% ylabel('Leading neuron FRTC (s)')
% xlabel('Following neuron FRTC (s)')
% xlim([0.5,3.5]);
% ylim([0.5,3.5]);
% cbh.Label.String='FC Rate (%)';
% title(t);
% end
% 
% 
% function plotOneBar(stats,idx,t,yspan)
% subplot(1,3,idx)
% hold on
% mm=[stats.per_bin_FC(1,1,1),stats.per_bin_FC(2,2,1),stats.per_bin_FC(1,2,1),stats.per_bin_FC(2,1,1)].*100;
% uci=[stats.per_bin_FC(1,1,2),stats.per_bin_FC(2,2,2),stats.per_bin_FC(1,2,2),stats.per_bin_FC(2,1,2)].*100;
% bci=[stats.per_bin_FC(1,1,3),stats.per_bin_FC(2,2,3),stats.per_bin_FC(1,2,3),stats.per_bin_FC(2,1,3)].*100;
% n=sum(stats.per_bin_FC(:,:,4),'all');
% for jj=1:4
%     bh(jj)=bar(jj,mm(jj),'FaceColor','w','EdgeColor','k');
% end
% bh(2).FaceColor='r';
% bh(3).FaceColor='b';
% bh(4).FaceColor='k';
% errorbar(1:4,mm,uci-mm,bci-mm,'k.')
% % set(gca(),'XTick',1:4,'XTickLabel',{'Early-early','Late-late','Early-late','Late-early'})
% ylabel('FC Rate (%)');
% ylim(yspan);
% title(sprintf('%s,n=%d',t,n));
% end
% 
% 
% function plotOneHist(stats,idx,t,yspan)
% subplot(1,3,idx)
% hold on
% for jj=1:size(stats.deltaTCOM_FC_rate,2)
%     [phat(jj),pci(jj,:)]=binofit(stats.deltaTCOM_FC_rate(1,jj),stats.deltaTCOM_FC_rate(2,jj));
% end
% fill([stats.deltaTCOM_FC_rate(3,:),stats.deltaTCOM_FC_rate(3,end:-1:1)],[pci(:,1);flip(pci(:,2))].*100,'r','FaceAlpha',0.1,'EdgeColor','none')
% plot(stats.deltaTCOM_FC_rate(3,:),phat.*100,'-r')
% ns=sum(stats.deltaTCOM_FC_rate(1,:));
% np=sum(stats.deltaTCOM_FC_rate(2,:));
% ylabel('FC Rate (%)');
% ylim(yspan);
% title(sprintf('%s,%d of %d',t,ns,np));
% set(gca,'XTick',-20:20:20,'XTickLabel',-5:5:5)
% xline(0,'--k')
% xlabel('follow-lead dTCOM')
% end
% 
% function p=chisqOne(stats)
%     vec1=[ones(stats.per_bin_FC(1,1,5),1);...
%         2*ones(stats.per_bin_FC(1,2,5),1);...
%         3*ones(stats.per_bin_FC(2,1,5),1);...
%         4*ones(stats.per_bin_FC(2,2,5),1);...
%         ];
% 
%     vec2=[(1:stats.per_bin_FC(1,1,5)).'>stats.per_bin_FC(1,1,4);...
%         (1:stats.per_bin_FC(1,2,5)).'>stats.per_bin_FC(1,2,4);...
%         (1:stats.per_bin_FC(2,1,5)).'>stats.per_bin_FC(2,1,4);...
%         (1:stats.per_bin_FC(2,2,5)).'>stats.per_bin_FC(2,2,4);...
%         ];
%     [~,~,p]=crosstab(vec1,vec2);
% end
% 
% 
% function p=chisqInter(stats)
%     vec1=[ones(sum(stats.congr_wave(5:6,6)),1);...
%         2*ones(sum(stats.congr_wave(1:2,6)),1);...
%         3*ones(sum(stats.congr_wave(3:4,6)),1);...
%         4*ones(stats.congr_inter_wave(5),1);...
%         5*ones(stats.congr_overlapped_wave(5),1);...
%         ];
% 
%     vec2=[(1:sum(stats.congr_wave(5:6,6))).'>sum(stats.congr_wave(5:6,5));...
%         (1:sum(stats.congr_wave(1:2,6))).'>sum(stats.congr_wave(1:2,5));...
%         (1:sum(stats.congr_wave(3:4,6))).'>sum(stats.congr_wave(3:4,5));...
%         (1:stats.congr_inter_wave(5)).'>stats.congr_inter_wave(4);...
%         (1:stats.congr_overlapped_wave(5)).'>stats.congr_overlapped_wave(4)...
%         ];
%     [~,~,p]=crosstab(vec1,vec2);
% end
% 
% function p=chisqHier(stats)
% 
% vec1=[2*ones(stats.per_bin_FC(1,2,5),1);...
%     3*ones(stats.per_bin_FC(2,1,5),1);...
%     ];
% 
% vec2=[...
%     (1:stats.per_bin_FC(1,2,5)).'>stats.per_bin_FC(1,2,4);...
%     (1:stats.per_bin_FC(2,1,5)).'>stats.per_bin_FC(2,1,4);...
%     ];
% [~,~,p]=crosstab(vec1,vec2);
% end
