function ax=plotOne(data,subidx,ttext,opt)
arguments
    data double
    subidx (1,1) double {mustBePositive,mustBeInteger}
    ttext (1,:) char
    opt.com double = []
end
subplot(2,2,subidx)
hold on
persistent gk
if isempty(gk)
    gk = fspecial('gaussian', [3 3], 1);
end
imagesc(conv2(data,gk,'same'),[-1 1])
colormap('jet')
ax=gca();
ax.YDir='normal';
ylabel('SU #')
ax.XTick=[12.5,32.5,52.5];
ax.XTickLabel={'0','5','10'};
xlabel('Time (s)');
arrayfun(@(x) xline(x,'--w'),[12.5,16.5,40.5,44.5]);
ylim([0.5,size(data,1)+0.5]);
xlim([9,40]);
title(ttext);
if ~isempty(opt.com)
    plot(opt.com,1:size(data,1),'--k');
end