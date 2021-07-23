frac.sust=ephys.per_region_fraction('memtype','sust');
frac.trans=ephys.per_region_fraction('memtype','trans');
coord=[];
fh=figure('Color','w');
hold on;
for ridx=1:size(frac.trans,1)
    if frac.trans{ridx,3}==5 && frac.trans{ridx,4}>40 && ~isempty(frac.trans{ridx,2})
        ssel=strcmp(frac.sust(:,2),frac.trans(ridx,2));
        xx=frac.trans{ridx,1}.*100;
        yy=frac.sust{ssel,1}.*100;
        coord=[coord;xx,yy];
        plot(xx,yy,'o','MarkerFaceColor','r','MarkerEdgeColor','none');
        text(xx,yy,frac.trans{ridx,2},'HorizontalAlignment','center','VerticalAlignment','top');
    end
end
[r,p]=corr(coord(:,1),coord(:,2));
text(min(xlim()),max(ylim()),sprintf('r=%.3f,p=%.3f',r,p),'HorizontalAlignment','left','VerticalAlignment','top')
ylabel('Sustained fraction (%)');
xlabel('Transient fraction (%)');
% ylim([0,4])