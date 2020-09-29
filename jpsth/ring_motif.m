load ringsX.mat
load('candidate_count.mat')
shufrpt=size(rings_shuf,1);
for midx=1:3
    for bin=1:6
        motif_count(midx,bin)=size(cell2mat([(rings(midx,:,bin,1)),(rings(midx,:,bin,2))]'),1);
    end
end
motif_count_shuf=nan(3,6,shufrpt);
for midx=1:3
    for bin=1:6
        for rpt=1:shufrpt
            motif_count_shuf(midx,bin,rpt)=size(cell2mat([squeeze(rings_shuf(rpt,midx,:,bin,1));squeeze(rings_shuf(rpt,midx,:,bin,2))]),1);
        end
    end
end
%% baseline
for midx=1:3
    base_motif_count(midx)=size(cell2mat([(base_rings(midx,:,1)),(base_rings(midx,:,2))]'),1);
end

base_motif_count_shuf=nan(3,shufrpt);
for midx=1:3
    for rpt=1:shufrpt
        base_motif_count_shuf(midx,rpt)=size(cell2mat([squeeze(rings_shuf(rpt,midx,:,1));squeeze(rings_shuf(rpt,midx,:,2))]),1);
    end
end
% shuf_mm=mean(shuf_count,3);
%% get figure
close all
for msize=3:5
    plotOne(msize,motif_count,motif_count_shuf,base_motif_count,base_motif_count_shuf,candidate_delay,candidate_count_base);
end


return


function plotOne(msize,motif_delay, motif_shuf,motif_base,motif_base_shuf,candidate_delay,candidate_count_base)
midx=msize-2;
mmshufC=mean(motif_shuf(midx,:,:),3);
mmbase_shuf_C=mean(motif_base_shuf(midx,:));
% mmbase_shuf=mmbase_shuf_C./sum(candidate_count_base(:,midx));
% mmshuf=mmshufC./sum(candidate_delay(:,:,midx));
fh=figure('Color','w','Position',[100,100,215,235]);
yyaxis left
hold on;
rbh=bar(-1,motif_base(midx)./mmbase_shuf_C,0.6,'FaceColor','w','FaceAlpha',0.5);
rh=bar((1:6),motif_delay(midx,:)./mmshufC,0.6,'FaceColor','w','FaceAlpha',0.5);
xlim([-1.75,6.75])
yline(1,'--','Color',[0.5,0.5,0.5])
set(gca,'XTick',[-1,1 5],'XTickLabel',{'ITI','1','5'},'YScale','log','YColor','k','FontSize',10);
ylim([0.5,500])
xlabel('Time bin (sec)')
ylabel('Observation / expectation')


yyaxis right
hold on
plot([-1,1:6],[motif_base(midx)./sum(candidate_count_base(:,midx)),...
    motif_delay(midx,:)./sum(candidate_delay(:,:,midx))],':ro','MarkerSize',3)
set(gca,'YScale','log','YColor','r','FontSize',10);
ylim([0.0000099,0.01]);
ylabel('Normalized density')
% legend([rh,sh],{'recorded','shuffled'});
% keyboard
% exportgraphics(fh,sprintf('motif_of_%d.pdf',msize),'ContentType','vector');
end
