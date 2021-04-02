function fh=plot_hist(congru,incongru,nonmem,opt)
arguments
    congru (:,11) double
    incongru (:,11) double
    nonmem (:,11) double
    opt.type (1,:) char {mustBeMember(opt.type,{'spk','fc_eff','fc_prob'})} = 'spk'
end
congci=fliplr(bootci(1000,@(x) mean(x),congru(:,2:end)))*100;
incongci=fliplr(bootci(1000,@(x) mean(x),incongru(:,2:end)))*100;
nmci=fliplr(bootci(1000,@(x) mean(x),nonmem(:,2:end)))*100;
congmm=flip(mean(congru(:,2:end)))*100;
incongmm=flip(mean(incongru(:,2:end)))*100;
nmmm=flip(mean(nonmem(:,2:end)))*100;
fh=figure('Color','w','Position',[100,100,400,300]);
hold on;
fill([10:20:200,fliplr(10:20:200)],[congci(1,:),fliplr(congci(2,:))],'r','EdgeColor','none','FaceAlpha',0.2);
fill([10:20:200,fliplr(10:20:200)],[incongci(1,:),fliplr(incongci(2,:))],'b','EdgeColor','none','FaceAlpha',0.2);
fill([10:20:200,fliplr(10:20:200)],[nmci(1,:),fliplr(nmci(2,:))],'k','EdgeColor','none','FaceAlpha',0.2);
hcong=plot(10:20:200,congmm,'-r');
hincong=plot(10:20:200,incongmm,'-b');
hnm=plot(10:20:200,nmmm,'-k');
switch opt.type
    case 'spk'
        xlabel('T(obervation) - T(pre spike) (ms)');
        ylabel('Post neuron firing increase (%)');
    case 'fc_eff'
        xlabel('Pre-spike latency (ms)');
        ylabel('Functional coupling efficacy increase (%)');
    case 'fc_prob'
        xlabel('Pre-spike latency (ms)');
        ylabel('Functional coupling event increase (%)');
        
end
grid on
legend([hcong,hincong,hnm],{'Congruent FC','Incongruent FC','NonMem FC'});
% exportgraphics(fh,'prehistory_post_firing.pdf');
end