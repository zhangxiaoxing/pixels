for ppow=0:5
    node(ppow+1)=10^ppow;
    edge(ppow+1)=(10^ppow)^2*0.01;
end
figure('Color','w')
hold on
fill([384,1536,1536,384],[1,1,10^8,10^8],'r','EdgeColor','none','FaceAlpha',0.1)
fill([16,128,128,16],[1,1,10^8,10^8],'k','EdgeColor','none','FaceAlpha',0.1)
xline(10^4,'--r')
% xline(1.5*10^5,'--k')
pnh=plot(node*2,node,'-k');
pfh=plot(node*2,edge,'-r');
set(gca,'XScale','log','YScale','log')
ylabel('Measurable units')
xlabel('Recording-sites/session')
ylim([1,10^8])
xlim([1,1.5*10^5])
legend([pnh,pfh],{'Neurons','Functional couplings'},'Location','northoutside','Orientation','horizontal')