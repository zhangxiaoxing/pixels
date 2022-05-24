function plotseq(in,msize,tidx,ids,sessIdx,ssidx,cocount)
fh=figure('Color','w','Position',[100,100,400,400]);
hold on
color={'r','b','c','m','g'};
ph=matlab.graphics.chart.primitive.Line.empty(0,5);
for bin=1:6
    for i=1:msize
        j=i;
        selp=in(:,2)==i & in(:,3)==1 & in(:,1)>=bin & in(:,1)<bin+1;
        seln=in(:,2)==i & in(:,3)==0 & in(:,1)>=bin & in(:,1)<bin+1;
        tph=plot(repmat(in(selp,1)',2,1)-bin,(bin-1)*(msize+1)+repmat([j-0.6;j+0.6],1,nnz(selp)),'-','Color',color{i});
%         ph(i)=tph(1);
        plot(repmat(in(seln,1)',2,1)-bin,(bin-1)*(msize+1)+repmat([j-0.6;j+0.6],1,nnz(seln)),'-','Color',[0.8,0.8,0.8]);

    end
end
xlabel('Time within time-bin(s)')
ylabel('Time-bin')
set(gca,'XTick',0:0.5:1,'YTick',(msize/2):(msize+1):((msize+1)*6),'YTickLabel',1:6)
if msize==3
    title(sprintf('Sess#%d, Trial#%d, SU%d, %d, %d',sessIdx,tidx,ids(1),ids(2),ids(3)));
    %     print(fh,sprintf('spike_seq_%d_%d_%d.png',ids(1),ids(2),ids(3)),'-dpng');
elseif msize==4
    title(sprintf('Sess#%d, Trial#%d, SU%d, %d, %d, %d',sessIdx,tidx,ids(1),ids(2),ids(3),ids(4)));
    %     print(fh,sprintf('spike_seq_%d_%d_%d_%d.png',ids(1),ids(2),ids(3),ids(4)),'-dpng');
elseif msize==4
    title(sprintf('Sess#%d, Trial#%d',sessIdx,tidx));    
end
print(fh,sprintf('r%d_spike_seq_%03d_%d.png',msize,cocount,ssidx),'-dpng');
if false
%     legend(ph,arrayfun(@(x) ['cell ',num2str(x)],1:5,'UniformOutput',false),'NumColumns',3,'Location','northoutside')
    ylim([0,6*(msize+1)]);
    fh.Position(3:4)=[225,225];
    set(gca,'FontSize',10);
    exportgraphics(fh,sprintf('r%d_spike_seq_%d_%d.pdf',msize,cocount,ssidx),'ContentType','vector');
end

close(fh);
end