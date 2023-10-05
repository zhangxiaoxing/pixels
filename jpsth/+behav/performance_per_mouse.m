[~,~,sessmap]=ephys.sessid2path(0);
sesses=cell2mat(sessmap.keys());
all_perf=[];
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

    all_perf=[all_perf;repmat(mmid,size(trials(:,1))),repmat(sess,size(trials(:,1))),trials];
end

[perfsum,hrsum,crsum]=deal([]);
for mm=reshape(double(unique(all_perf(:,1))),1,[])
    perfsum=[perfsum;mm,mean(all_perf(all_perf(:,1)==mm & all_perf(:,11)==1,12))];
    hrsum=[hrsum;mm,mean(all_perf(all_perf(:,1)==mm & all_perf(:,11)==1 & all_perf(:,7)~=all_perf(:,8),12))];
    crsum=[crsum;mm,mean(all_perf(all_perf(:,1)==mm & all_perf(:,11)==1 & all_perf(:,7)==all_perf(:,8),12))];
end

figure
tiledlayout(1,3)
nexttile
hold on
scatter(ones(size(perfsum(:,2))),perfsum(:,2),'k','filled','MarkerFaceAlpha',0.2)
errorbar(1,mean(perfsum(:,2)),std(perfsum(:,2))./sqrt(size(perfsum,1)),'ro','LineWidth',1,'CapSize',12)
ylim([0.4,1])
set(gca(),'XTick',[],'YTick',0.5:0.25:1,'YTickLabel',50:25:100)
ylabel('Performace (%)')

nexttile
hold on
scatter(ones(size(hrsum(:,2))),hrsum(:,2),'k','filled','MarkerFaceAlpha',0.2)
errorbar(1,mean(hrsum(:,2)),std(hrsum(:,2))./sqrt(size(hrsum,1)),'ro','LineWidth',1,'CapSize',12)
ylim([0.4,1])
set(gca(),'XTick',[],'YTick',0.5:0.25:1,'YTickLabel',50:25:100)
ylabel('Hit rate (%)')

nexttile
hold on
scatter(ones(size(crsum(:,2))),crsum(:,2),'k','filled','MarkerFaceAlpha',0.2)
errorbar(1,mean(crsum(:,2)),std(crsum(:,2))./sqrt(size(crsum,1)),'ro','LineWidth',1,'CapSize',12)
ylim([0.4,1])
set(gca(),'XTick',[],'YTick',0.5:0.25:1,'YTickLabel',50:25:100)
ylabel('Correct rejection (%)')

fid=fopen(fullfile('binary','upload','SF1D_behavi.json'),'w');
dout.performace=array2table(perfsum,'VariableNames',{'MouseID','Performace_correct_rate'})
dout.hit_rate=array2table(hrsum,'VariableNames',{'MouseID','Hit_rate'})
dout.correct_rejection_rate=array2table(crsum,'VariableNames',{'MouseID','Correct_rejection_rate'})
fprintf(fid,jsonencode(dout))
fclose(fid)
