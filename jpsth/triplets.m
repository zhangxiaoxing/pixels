load('motif_count.mat')
% load('0831_selec_conn_chain_duo_6s_1_2.mat')

% su_set=unique(pair_chain(:));
% tripCandidate=tripCount(su_set);
for msize=3:5
    plotOne(msize,motif_count,motif_count_shuf,motif_baseline_count,motif_baseline_count_shuf)
end


function plotOne(msize,countm, count_shuf,countb,countb_shuf)
midx=msize-2;
mmshuf=mean(count_shuf(midx,:,:),3);
% stdshuf=std(count_shuf(midx,:,:),0,3);
% (triplet_count'-mmshuf)./stdshuf
mmbase=mean(countb_shuf(midx,:));
fh=figure('Color','w','Position',[100,100,215,235]);
hold on;

yyaxis left
sbh=bar(-1-0.1,mmbase,0.6,'FaceColor','k','FaceAlpha',1);
rbh=bar(-1+0.1,countb(midx),0.6,'FaceColor','w','FaceAlpha',0.5);
ylabel('Baseline inactive congruent rings');

yyaxis right
sh=bar((1:6)-0.1,mmshuf,0.6,'FaceColor','k','FaceAlpha',1);
rh=bar((1:6)+0.1,countm(midx,:),0.6,'FaceColor','w','FaceAlpha',0.5);
% ylim([0,0.0004])
xlim([-1.75,6.75])
set(gca,'XTick',[1 5]);
xlabel('Time bin (sec)')
ylabel('Delay active congruent rings')
% legend([ih,oh,nh],{'recorded','shuffled','non-sel.'});
legend([rh,sh],{'recorded','shuffled'});
keyboard
exportgraphics(fh,sprintf('motif_of_%d.pdf',msize),'ContentType','vector');
end




function out=tripCount(su_set)
out=0;
for i=1:114
    subsetCount=nnz(su_set>i*100000 & su_set<(i+1)*100000);
    if subsetCount>=3
        out=out+nchoosek(subsetCount,3);
    end
end
end