function ax=plotOne(data,subidx,ttext)
subplot(1,4,subidx)
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
title(ttext);