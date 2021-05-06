function fh=plot_hist(congru,incongru,nonmem,opt)
arguments
    congru (:,11) double
    incongru (:,11) double
    nonmem (:,11) double
    opt.type (1,:) char {mustBeMember(opt.type,{'spk','fc_eff','fc_prob'})} = 'spk'
    opt.title (1,:) char = '' %figure title
    opt.binw (1,1) double = 20
    opt.ylim (1,2) double {mustBeNonnegative}= [0,0]
end
congci=fliplr(bootci(1000,@(x) mean(x),congru(:,2:end)))*100;
incongci=fliplr(bootci(1000,@(x) mean(x),incongru(:,2:end)))*100;
nmci=fliplr(bootci(1000,@(x) mean(x),nonmem(:,2:end)))*100;
congmm=flip(mean(congru(:,2:end)))*100;
incongmm=flip(mean(incongru(:,2:end)))*100;
nmmm=flip(mean(nonmem(:,2:end)))*100;
minbin=opt.binw./2;
mxbin=(size(congru,2)-1)*opt.binw;

fh=figure('Color','w','Position',[100,100,400,300]);
hold on;
fill([minbin:opt.binw:mxbin,fliplr(minbin:opt.binw:mxbin)],[congci(1,:),fliplr(congci(2,:))],'r','EdgeColor','none','FaceAlpha',0.2);
fill([minbin:opt.binw:mxbin,fliplr(minbin:opt.binw:mxbin)],[incongci(1,:),fliplr(incongci(2,:))],'b','EdgeColor','none','FaceAlpha',0.2);
fill([minbin:opt.binw:mxbin,fliplr(minbin:opt.binw:mxbin)],[nmci(1,:),fliplr(nmci(2,:))],'k','EdgeColor','none','FaceAlpha',0.2);
hcong=plot(minbin:opt.binw:mxbin,congmm,'-r');
hincong=plot(minbin:opt.binw:mxbin,incongmm,'-b');
hnm=plot(minbin:opt.binw:mxbin,nmmm,'-k');
switch opt.type
    case 'spk'
        xlabel('Pre - post latency (ms)');
        ylabel('Post spike increase (%)');
    case 'fc_eff'
        xlabel('Pre-spike latency (ms)');
        ylabel('Functional coupling efficacy increase (%)');
    case 'fc_prob'
        xlabel('Pre-spike latency (ms)');
        ylabel('Functional coupling event increase (%)');
        
end
grid on
legend([hcong,hincong,hnm],{'Congruent FC','Incongruent FC','Non-Memory FC'});
if ~isempty(opt.title), title(opt.title);end
if max(opt.ylim)>0, ylim(opt.ylim);end
% exportgraphics(fh,'prehistory_post_firing.pdf');
end