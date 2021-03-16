function plotCCG(range)
mono=evalin('base','mono');
if ~exist('range','var')
    range=1:size(mono.sig_con,1);
elseif numel(range)==1
    range=range:size(mono.sig_con,1);
elseif numel(range)~=2
    disp('Error parsing range')
    return
end
for j=range
    fh=figure('Position',[32,32,320,240]);
    plot(mono.ccgR(:,mono.sig_con(j,1),mono.sig_con(j,2)),'r-');
    xlim([1,501])
    arrayfun(@(x) xline(x,'b:'),[251-25,251,251+25])
    set(gca(),'XTick',[1,251-25,251,251+25,501],'XTickLabel',[-100,-10,0,10,100])
    xlim([151,351])
    title(j)
    keyboard();
    close(fh)
end
end