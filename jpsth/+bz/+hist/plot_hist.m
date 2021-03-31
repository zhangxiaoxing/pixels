function fh=plot_hist(congru,nonmem)
congci=fliplr(bootci(1000,@(x) mean(x),congru(:,2:end)))*100;
nmci=fliplr(bootci(1000,@(x) mean(x),nonmem(:,2:end)))*100;
congmm=flip(mean(congru(:,2:end)))*100;
nmmm=flip(mean(nonmem(:,2:end)))*100;
fh=figure('Color','w','Position',[100,100,400,300]);
hold on;
fill([10:20:200,fliplr(10:20:200)],[congci(1,:),fliplr(congci(2,:))],'r','EdgeColor','none','FaceAlpha',0.2);
fill([10:20:200,fliplr(10:20:200)],[nmci(1,:),fliplr(nmci(2,:))],'k','EdgeColor','none','FaceAlpha',0.2);
hcong=plot(10:20:200,congmm,'-r');
hnm=plot(10:20:200,nmmm,'-k');
xlabel('T(obervation) - T(pre spike) (ms)');
ylabel('Post neuron firing increase (%)');
grid on
legend([hcong,hnm],{'Congruent FC','NonMem FC'});
% exportgraphics(fh,'prehistory_post_firing.pdf');
end