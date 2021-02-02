plotSpikeRaster()
keyboard()

if false
    binCount=[];
    for bin=1:6
        load(sprintf('nonsel_silent_fc_%d.mat',bin));
        %%sums
        binmat=cell2mat(sums(:,[2,5]));
        binmat(:,12)=bin;
        binCount=[binCount;binmat];
    end
    
    binCntIdx=find(sum(binCount(:,3:11)<0.01,2)>=4);
end


for bin=0%1:6
    load(sprintf('nonsel_silent_fc_%d.mat',bin));
    idces=find(sum(cell2mat(sums(:,5))<0.01,2)>=4);
    for idx=9%idces'
        keyboard();
        stats=sums{idx,4};
        trials=sums{idx,3};
        sel_6s_S1=trials(:,5)==4 & trials(:,8)== 6;
        sel_6s_S2=trials(:,5)==8 & trials(:,8)== 6;
        
        sel_3s_S1=trials(:,5)==4 & trials(:,8)== 3;
        sel_3s_S2=trials(:,5)==8 & trials(:,8)== 3;
        
        fh=figure('Color','w','Position',[32,32,900,800]);
        sgtitle(sprintf('bin %d, Idx %d',bin,idx));
        
        
        subplot(3,2,1);
        plotOne(3)
        ylabel('Events Freq (Hz)')
        title ('1 -> 2, 10ms window');
        
        subplot(3,2,4);
        plotOne(1)
        ylabel('Pre cell FR (Hz)')
        
        subplot(3,2,6)
        plotOne(2)
        ylabel('Post cell FR (Hz)')
        % legend([h61,h31,h62,h32],{'S1 6s','S1 3s','S2 6s','S2 3s'},'Location','eastoutside')
        
        
        subplot(3,2,3)
        plotOne(4)
        ylabel('Events Freq (Hz)')
        title ('2 -> 1, 10ms window');
        
        
        subplot(3,2,5)
        plotOne(5)
        ylabel('Events Freq (Hz)')
        title ('1 -> 2, 50ms-60ms window');
        
    end
end


function ax=plotOne(statIdx)
stats=evalin('base','stats');
sel_6s_S1=evalin('base','sel_6s_S1');
sel_3s_S1=evalin('base','sel_3s_S1');
sel_6s_S2=evalin('base','sel_6s_S2');
sel_3s_S2=evalin('base','sel_3s_S2');


hold on;

ci61=bootci(1000,@(x) mean(x), squeeze(stats(statIdx,sel_6s_S1,:)));
ci31=bootci(1000,@(x) mean(x), squeeze(stats(statIdx,sel_3s_S1,:)));
ci62=bootci(1000,@(x) mean(x), squeeze(stats(statIdx,sel_6s_S2,:)));
ci32=bootci(1000,@(x) mean(x), squeeze(stats(statIdx,sel_3s_S2,:)));

fill([-3:10,10:-1:-3]+0.5,[ci61(1,:),fliplr(ci61(2,:))],'r','EdgeColor','none','FaceAlpha',0.1)
fill([-3:10,10:-1:-3]+0.5,[ci31(1,:),fliplr(ci31(2,:))],'m','EdgeColor','none','FaceAlpha',0.1)
fill([-3:10,10:-1:-3]+0.5,[ci62(1,:),fliplr(ci62(2,:))],'b','EdgeColor','none','FaceAlpha',0.1)
fill([-3:10,10:-1:-3]+0.5,[ci32(1,:),fliplr(ci32(2,:))],'c','EdgeColor','none','FaceAlpha',0.1)

h61=plot((-3:10)+0.5,squeeze(mean(stats(statIdx,sel_6s_S1,:),2)),'-r');
h31=plot((-3:10)+0.5,squeeze(mean(stats(statIdx,sel_3s_S1,:),2)),'-m');
h62=plot((-3:10)+0.5,squeeze(mean(stats(statIdx,sel_6s_S2,:),2)),'-b');
h32=plot((-3:10)+0.5,squeeze(mean(stats(statIdx,sel_3s_S2,:),2)),'-c');
arrayfun(@(x) xline(x,'k:'),[0 1 4 5 7 8]);
xlabel('Time (s)')
xlim([-2,10]);
end

function plotPertrial()

figure('Position',[32,32,600,400])
hold on
s1t=sums{idx,4}(3,sel_3s_S1,7)
s2t=sums{idx,4}(3,sel_3s_S2,7)
plot(find(sel_3s_S1),s1t,'ro')
plot(find(sel_3s_S2),s2t,'bo')


figure('Position',[32,32,600,400])
hold on
s1t=sums{idx,4}(3,sel_3s_S1,7)./sqrt(sums{idx,4}(1,sel_3s_S1,7).*sums{idx,4}(2,sel_3s_S1,7))
s2t=sums{idx,4}(3,sel_3s_S2,7)./sqrt(sums{idx,4}(1,sel_3s_S2,7).*sums{idx,4}(2,sel_3s_S2,7))
plot(find(sel_3s_S1),s1t,'ro')
plot(find(sel_3s_S2),s2t,'bo')


s1t=sums{idx,4}(3,sel_3s_S1,7);
s1s=sqrt(sums{idx,4}(1,sel_3s_S1,7).*sums{idx,4}(2,sel_3s_S1,7));
s1i=find(sel_3s_S1);
s2t=sums{idx,4}(3,sel_3s_S2,7);
s2s=sqrt(sums{idx,4}(1,sel_3s_S2,7).*sums{idx,4}(2,sel_3s_S2,7));
s2i=find(sel_3s_S2);

figure('Position',[32,32,600,400])
hold on
for i=1:nnz(sel_3s_S1)
    text(s1t(i),s1s(i),num2str(s1i(i)));
end
grid on
plot([0,12],[0,45],'r:')
xlim([0,12])
ylim([0,45])

figure('Position',[32,32,600,400])
hold on
for i=1:nnz(sel_3s_S2)
    text(s2t(i),s2s(i),num2str(s2i(i)));
end
grid on
plot([0,12],[0,45],'r:')
xlim([0,12])
ylim([0,45])



plot(find(sel_3s_S1),s1t,'ro')
plot(find(sel_3s_S2),s2t,'bo')

figure('Position',[32,32,600,400])
hold on
histogram(s1t,0:12,'Normalization','probability')
histogram(s2t,0:12,'Normalization','probability')

end
