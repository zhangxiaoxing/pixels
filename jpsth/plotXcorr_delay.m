stats=cell(0);
for sidx=1:size(sums,1)
    xc=sums{sidx,5};
    sustCount=numel(sums{sidx,3});
    for si=1:size(xc.xcorr,1)
        if si<=sustCount
            su1='sust';
        else
            su1='transient';
        end
        for sj=1:size(xc.xcorr,2)
            if sj==si
                su2='auto-corr';
            elseif sj<=sustCount
                su2='sust';
            else
                su2='transient';
            end
            totalCount=nansum(squeeze(xc.xcorr(si,sj,:)));
            if totalCount<50
                stats(end+1,:)={sidx,si,sj,su1,su2,totalCount,NaN};
                continue
            end
            leadingProp=nansum(squeeze(xc.xcorr(si,sj,1:40)))/totalCount;
            stats(end+1,:)={sidx,si,sj,su1,su2,totalCount,leadingProp};
            
            if leadingProp>0.75 || leadingProp<0.25
                fh=figure('Color','w','Position',[100,100,400,300]);
                stem(xc.time*1000,squeeze(xc.xcorr(si,sj,:)),'Marker','none','LineWidth',10)
                title(sprintf('%s-%s,%d,%d,%d',su1,su2,sidx,si,sj))
                xlabel('time (ms)')
                ylabel('spike-pair count')
                print(fh,sprintf('%s_%s_%d_%d_%d.png',su1,su2,sidx,si,sj),'-dpng')

%                 pause
                close(fh)
            end
        end
    end
end

autoCount=sum(strcmp(stats(:,5),'auto-corr'));
stCount=sum(strcmp(stats(:,4),'sust') & strcmp(stats(:,5),'transient'));
stCount50=sum(strcmp(stats(:,4),'sust') & strcmp(stats(:,5),'transient') & cell2mat(stats(:,6))>=50);
% ssCount=sum(strcmp(stats(:,4),'sust') & strcmp(stats(:,5),'sust'));
ttCount=sum(strcmp(stats(:,4),'transient') & strcmp(stats(:,5),'transient'))/2;
ttCount50=sum(strcmp(stats(:,4),'transient') & strcmp(stats(:,5),'transient')& cell2mat(stats(:,6))>=50)/2;

