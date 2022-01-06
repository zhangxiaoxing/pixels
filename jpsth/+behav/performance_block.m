%TODO performance in block of 4+4 trials
[~,~,sessmap]=ephys.sessid2path(0);
sesses=cell2mat(sessmap.keys());
all_perf=[];
for sess=reshape(sesses,1,[])
    sess_perf=nan(1,8);
    [~,~,trials]=ephys.getSPKID_TS(sess,'skip_spike',true);
    
    dur_resp=behav.tag_block(trials);
    
    for bi=1:4
        sess_perf(bi)=mean(dur_resp(dur_resp(:,1)==3 & dur_resp(:,3)==bi,2));
        sess_perf(bi+4)=mean(dur_resp(dur_resp(:,1)==6 & dur_resp(:,3)==bi,2));
    end
    all_perf=[all_perf;sess_perf];
end

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
