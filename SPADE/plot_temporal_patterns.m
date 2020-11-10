close all
fh=figure('Color','w','Position',[100,100,1000,200]);
fidx=1;
keysseq=[];
for lag2=0:3
    for lag3=lag2:3
        for lag4=lag3:3
            subplot(2,10,fidx)
            hold on
            plot([0,0],1+[-0.4,0.4],'r-')
            plot([lag2,lag2],2+[-0.4,0.4],'r-')
            plot([lag3,lag3],3+[-0.4,0.4],'r-')
            plot([lag4,lag4],4+[-0.4,0.4],'r-')
            text(3,1,sprintf('#%d',fidx),'HorizontalAlignment','right');
            xlim([-0.5,3.5])
            set(gca,'XTick',0:3,'XTickLabel',0:4:12,'YTick',[])%,1:4,'YTickLabel',{'Neu1','Neu2','Neu3','Neu4'})
            fidx=fidx+1;
            keysseq(end+1)=lag2*100+lag3*10+lag4;
        end
    end
end
exportgraphics(fh,'4_su_temporal_pattern.pdf')
