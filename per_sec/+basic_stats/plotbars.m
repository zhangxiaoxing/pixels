function fh=plotbars(reg,count,depth)
arguments
    reg (:,1) cell
    count (:,3) double
    depth (1,1) double {mustBeMember(depth,1:6)}
end

sust=count(:,2)./count(:,1);
trans=count(:,3)./count(:,1);
stats=[trans,sust]*100;
[~,idces]=sort(trans,'descend');
fh=figure('Color','w','Position',[100,100,numel(reg)*20+50,300]);
bh=bar(stats(idces,:));
set(gca(),'XTick',1:numel(reg),'XTickLabel',reg(idces),'XTickLabelRotation',90)
ylabel('Fraction of selective neuron (%)')
legend(bh,{'Transient','Sustained'});
title(sprintf('Branch level %d',depth+2));
end