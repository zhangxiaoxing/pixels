% likely obsolete as of 22-06-15

function plot_fc_decoding(sel_meta,delay,opt)
arguments
    sel_meta
    delay (1,1) double {mustBeMember(delay,[3 6])} = 6
    opt.trials (1,1) double {mustBeMember(opt.trials, [20 30 40])} = 20
    opt.debug (1,1) logical = false
    opt.rpts (1,1) double =10
    opt.decoder (1,:) char {mustBeMember(opt.decoder, {'LDA','SVM','NB'})}='SVM'
end

out=cellfun(@(x) ...
    fc.dec.dec(...
    fc.get_fc_decoding(sel_meta,'rpts',opt.rpts,'type',x,'trials',opt.trials,'delay',delay,'debug',opt.debug),...
    'decoder',opt.decoder),...
    {'none','both','pre','post'});

out(end+1).cvcorr=mat2cell(cell2mat(cat(1,out(:).shufcorr)),... % mat % merge shuffle
    size(out(1).shufcorr{1},1)*size(out,2),... % n row
    ones(1,size(out(1).shufcorr,2))); % n col

out=cell2struct([{out(:).cvcorr};... % cvcorr % none,both,shuffle
    cellfun(@(dd) bootci(100,@(x) mean(x),cell2mat(dd.cvcorr)),num2cell(out),'UniformOutput',false);... % ci
    {'b','r','c','m','k'}],... % color
    {'cvcorr';'ci';'color'},1);

fh=figure('Position',[100,100,400,400],'Color','w');
hold on;
xspan=(-3:10)+0.5;
cellfun(@(dd)...
    fill([xspan,fliplr(xspan)],...
    [dd.ci(1,:),fliplr(dd.ci(2,:))]*100,dd.color,'FaceAlpha',0.2,'EdgeColor','none'),...
    num2cell(out));

fhs=cellfun(@(dd)...
    plot(xspan,mean(cell2mat(dd.cvcorr))*100,'-','Color',dd.color),...
    num2cell(out),'UniformOutput',false);

arrayfun(@(x) xline(x,':k'),[0,1,delay+1,delay+2]);

ylabel('Classification accuracy (%)');
xlabel('Time (s)');
xlim([-2,10])
legend([fhs{[5,1,2,3,4]}],{'Shuffle','Neither rate-select','Both rate-select','pre rate-select','post rate-select'},'Location','northoutside')
title(sprintf('%ds delay',delay));
end