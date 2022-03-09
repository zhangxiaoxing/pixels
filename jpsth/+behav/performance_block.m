%TODO performance in block of 4+4 trials
[~,~,sessmap]=ephys.sessid2path(0);
sesses=cell2mat(sessmap.keys());
all_perf=cell(0,10);
for sess=reshape(sesses,1,[])
    path=ephys.sessid2path(sess);
    MID=regexp(path,'(?<=M)\d{2,3}(?=_)','match','once');
    NID=regexp(path,'(?<=[_\\])\d{2,3}(?=_)','match','once');
    if xor(isempty(MID),isempty(NID))
        mmid=str2double([MID,NID]);
    elseif isequal(MID,NID)
        mmid=str2double(MID);
    else
        error('Failed to id mouse');
    end

    sess_perf=cell(1,8);
    [~,~,trials]=ephys.getSPKID_TS(sess,'skip_spike',true);
    
    dur_resp=behav.tag_block(trials,'wt',true);
    
    for bi=1:4
        sess_perf{bi}=histcounts(dur_resp(dur_resp(:,1)==3 & dur_resp(:,3)==bi,2),-0.5:1:1.5);
        sess_perf{bi+4}=histcounts(dur_resp(dur_resp(:,1)==6 & dur_resp(:,3)==bi,2),-0.5:1:1.5);
    end
    all_perf(end+1,:)=[mmid,sess,sess_perf];
end
%% per mouse
per_mouse=[];
for umid=reshape(unique(cell2mat(all_perf(:,1))),1,[])
    perfmat=sum(cell2mat(all_perf(cell2mat(all_perf(:,1))==umid,3:end)),1);
    per_mouse=[per_mouse;umid,perfmat(2:2:16)./(perfmat(2:2:16)+perfmat(1:2:15))];
end

fh=figure('Color','w');
subplot(1,3,1)
hold on
scatter(per_mouse(:,5),per_mouse(:,6))
plot([0.4,1],[0.4,1],'--r')
xlim([0.4,1])
ylim([0.4,1])
xlabel('Last 3s-delay trials')
ylabel('First 6s-delay trials')
[r,p]=corr(per_mouse(:,5),per_mouse(:,6));
text(min(xlim()),min(ylim()),sprintf('r=%.2f,p=%.2f',r,p),'HorizontalAlignment','left','VerticalAlignment','bottom')

subplot(1,3,2)
hold on
scatter(per_mouse(:,9),per_mouse(:,2))
plot([0.4,1],[0.4,1],'--r')
xlim([0.4,1])
ylim([0.4,1])
xlabel('Last 6s-delay trials')
ylabel('First 3s-delay trials')
[r,p]=corr(per_mouse(:,9),per_mouse(:,2));
text(min(xlim()),min(ylim()),sprintf('r=%.2f,p=%.2f',r,p),'HorizontalAlignment','left','VerticalAlignment','bottom')

subplot(1,3,3)
hold on
scatter(per_mouse(:,7),per_mouse(:,8))
plot([0.4,1],[0.4,1],'--r')
xlim([0.4,1])
ylim([0.4,1])
xlabel('2nd 6s-delay trials')
ylabel('3rd 6s-delay trials')
[r,p]=corr(per_mouse(:,7),per_mouse(:,8));
text(min(xlim()),min(ylim()),sprintf('r=%.2f,p=%.2f',r,p),'HorizontalAlignment','left','VerticalAlignment','bottom')


return

%% session scatter
fh=figure('Color','w');
subplot(1,3,1)
hold on
scatter(all_perf(:,4),all_perf(:,5))
plot([0.4,1],[0.4,1],'--r')
xlim([0.4,1])
ylim([0.4,1])
xlabel('Last 3s-delay trials')
ylabel('First 6s-delay trials')
[r,p]=corr(all_perf(:,4),all_perf(:,5));
text(min(xlim()),min(ylim()),sprintf('r=%.2f,p=%.2f',r,p),'HorizontalAlignment','left','VerticalAlignment','bottom')

subplot(1,3,2)
hold on
scatter(all_perf(:,8),all_perf(:,1))
plot([0.4,1],[0.4,1],'--r')
xlim([0.4,1])
ylim([0.4,1])
xlabel('Last 6s-delay trials')
ylabel('First 3s-delay trials')
[r,p]=corr(all_perf(:,8),all_perf(:,1));
text(min(xlim()),min(ylim()),sprintf('r=%.2f,p=%.2f',r,p),'HorizontalAlignment','left','VerticalAlignment','bottom')

subplot(1,3,3)
hold on
scatter(all_perf(:,6),all_perf(:,7))
plot([0.4,1],[0.4,1],'--r')
xlim([0.4,1])
ylim([0.4,1])
xlabel('2nd 6s-delay trials')
ylabel('3rd 6s-delay trials')
[r,p]=corr(all_perf(:,6),all_perf(:,7));
text(min(xlim()),min(ylim()),sprintf('r=%.2f,p=%.2f',r,p),'HorizontalAlignment','left','VerticalAlignment','bottom')
%%




if false %bar graph
    fh=figure('Color','w','Position',[32,32,175,150]);
    hold on
    mm=mean(all_perf).*100;
    sem=std(all_perf)./sqrt(size(all_perf,1)).*100;
    bar(mm,'FaceColor','w','EdgeColor','k');
    errorbar(1:8,mm,sem,'k.')
    ylim([50,100])
    set(gca(),'XTick',1:8,'XTickLabel',[1:4,1:4])
    ylabel('Correct rate (%)')
    anova1(all_perf)
    exportgraphics(fh,'perf_block.pdf','ContentType','vector')
end