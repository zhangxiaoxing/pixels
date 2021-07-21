fcom=load('per_region_com_collection.mat','collection');
ffrac=load('per_region_fraction_collection.mat','collection');

ureg=fcom.collection(cell2mat(fcom.collection(:,4))>20 & cell2mat(fcom.collection(:,3))==5,2);

fh=figure('Color','w');
hold on;
coord=[];
regs=[];

for ri=1:numel(ureg)
    if isempty(deblank(ureg{ri})), continue;end
    fridx=find(strcmp(ffrac.collection(:,2),ureg(ri)));
    if ~isempty(fridx) || ffrac.collection{fridx,4}>40
        comidx=find(strcmp(fcom.collection(:,2),ureg(ri)));
        xx=fcom.collection{comidx,1}./4;
        yy=ffrac.collection{fridx,1};
        coord=[coord;xx,yy];
        regs=[regs,fcom.collection{comidx,2}];
        plot(xx,yy,'wo','MarkerFaceColor',[255,127,127]./255);
        text(xx,yy,fcom.collection{comidx,2},'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end
ylim([0.1,0.6])
xlim([2.4,3.5]);
set(gca(),'XTick',2.5:0.5:3.5);
ylabel('Fraction of delay selective neuron')
xlabel('F.R. center of mass (s)')
[r,p]=corr(coord(:,1),coord(:,2));
text(max(xlim()),max(ylim()),sprintf('r = %.3f, p = %.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top');