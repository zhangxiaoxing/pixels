%%
figure('Color','w','Position',[32,32,315,300]);
hold on;
corrmat=[dursel(dursel(:,6)<0.05,9),dursel(dursel(:,6)<0.05,10)];
corrmat(:,3)=1;
scatter(corrmat(:,1),corrmat(:,2),9,'k','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.1);
[r,p]=corr(corrmat(:,1),corrmat(:,2));
b=corrmat(:,2)\corrmat(:,1);
plot(xlim(),xlim().*b,'--r')
xlim([-0.5,0.5])
ylim([-0.5,0.5])
xlabel('Selectivity index in correct trial');
ylabel('Selectivity index in error trial');
text(max(xlim()),max(ylim()),sprintf('r = %.2f, p < 0.001',r),'HorizontalAlignment','right','VerticalAlignment','top')
%%