function previous_trial_glm()
load(fullfile('binary','trials_dict.mat'),'trials_dict');
fittbl=[];
for sessid=reshape(trials_dict.keys,1,[])
    trials=cell2mat(trials_dict(sessid));
    for wtt=reshape(setdiff(find(trials(:,9)),1:3),1,[])
        fittbl=[fittbl;...
            trials(wtt,5)~=trials(wtt,6),...
            trials(wtt-1,5)~=trials(wtt,6),...
            trials(wtt-2,5)~=trials(wtt,6),...
            trials(wtt-3,5)~=trials(wtt,6),...
            trials(wtt,7)>0];
    end
end
if false
glmstr=cell2struct(mat2cell(fittbl,size(fittbl,1),ones(size(fittbl,2),1)).',...
    {'Current_trial_pair_match','Previous_trial_pair_match','Second_to_current_trial_pair_match','Third_to_current_trial_pair_match','Lick_response'});
fid=fopen(fullfile('binary','upload','SF7B_Previous_trial_effects_on_behavior.json'),'w');
fprintf(fid,jsonencode(glmstr));
fclose(fid)
end
[b,dev,stats]=glmfit(fittbl(:,1:4),fittbl(:,5),'binomial','link','identity');

fh=figure();
hold on
bar(flip(b),'FaceColor',[0.5,0.5,0.5])
errorbar(1:5,flip(b),flip(stats.se),'k.','CapSize',12)
ylabel('Weigh (b) in predicting lick responses')
set(gca,'XTick',1:5,'XTickLabel',{'3rd to','2nd to','previous','current','constant'})
savefig(fh,fullfile('binary','behav_previous_trls_glm.fig'))