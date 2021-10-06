%bool_congru+,bool_corss+,sess,rsize,loops_id,cids (+up to 2 fill-in), cid,delay correct,iti,delay error

if ~exist('loops_sums','var')
    load('loops_coding.mat','loops_sums')  %from \+bz\+rings\rings_coding.m
end
if ~exist('meta','var')
    meta=ephys.util.load_meta();
end

mtypes=["congru","nonmem"];
regtypes=["cross","within"];
rsizes=3:5;
rings_conding_sums=[];
for rsize=rsizes
    for mtype=mtypes
        for rtype=regtypes
            if ~isfield(loops_sums,mtype), continue;end
            if ~isfield(loops_sums.(mtype),sprintf('%s_%d',rtype,rsize)), continue;end
            loops=cell2mat(loops_sums.(mtype).(sprintf('%s_%d',rtype,rsize)).');
            for li=1:numel(loops)
%                 if loops(li).sessid==55
%                     keyboard();
%                 end
                trials=loops(li).trialinfo;
                pref_trials=find(trials(:,5)==4 & trials(:,8)==6 & trials(:,9)==1 & trials(:,10)==1);
                nonp_trials=find(trials(:,5)==8 & trials(:,8)==6 & trials(:,9)==1 & trials(:,10)==1);
                pref_err_trls=find(trials(:,5)==4 & trials(:,8)==6 & trials(:,10)==0);
                nonp_err_trls=find(trials(:,5)==8 & trials(:,8)==6 & trials(:,10)==0);
                memtype=meta.mem_type(meta.sess==loops(li).sessid & meta.allcid==str2double(loops(li).label{1}));
                if strcmp(mtype,'congru') && memtype>2
                    [pref_trials,nonp_trials]=deal(nonp_trials,pref_trials);
                    [pref_err_trls,nonp_err_trls]=deal(nonp_err_trls,pref_err_trls);
                end
                %bool_congru+,bool_corss+,sess,rsize,loops_id,cids (+up to 2 fill-in), cid,iti sel, delay correct, delay error
                cidfill=nan(1,5);
                cidfill(1:rsize)=str2double(loops(li).label);
                for ci=1:numel(loops(li).label)
                    pref_delay=nnz(ismember(loops(li).trial{ci},pref_trials) & loops(li).time{ci} >= 1 & loops(li).time{ci} < 7);
                    nonp_delay=nnz(ismember(loops(li).trial{ci},nonp_trials) & loops(li).time{ci} >= 1 & loops(li).time{ci} < 7);
                    pref_iti=nnz(ismember(loops(li).trial{ci},pref_trials) & loops(li).time{ci} >= -3 & loops(li).time{ci} < -1);
                    nonp_iti=nnz(ismember(loops(li).trial{ci},nonp_trials) & loops(li).time{ci} >= -3 & loops(li).time{ci} < -1);
                    pref_err=nnz(ismember(loops(li).trial{ci},pref_err_trls) & loops(li).time{ci} >= 1 & loops(li).time{ci} < 7);
                    nonp_err=nnz(ismember(loops(li).trial{ci},nonp_err_trls) & loops(li).time{ci} >= 1 & loops(li).time{ci} < 7);
                    rings_conding_sums=[rings_conding_sums;...
                        strcmp(mtype,'congru'),strcmp(rtype,'cross'),...
                        loops(li).sessid,rsize,li,cidfill,cidfill(ci),...
                        pref_delay,nnz(pref_trials),...
                        nonp_delay,nnz(nonp_trials),...
                        pref_iti,nnz(pref_trials),...
                        nonp_iti,nnz(nonp_trials),...
                        pref_err,nnz(pref_err_trls),...
                        nonp_err,nnz(nonp_err_trls)];
                end
            end
        end
    end
end

per_trl=[rings_conding_sums(:,12)./rings_conding_sums(:,13),rings_conding_sums(:,14)./rings_conding_sums(:,15),...
    rings_conding_sums(:,16)./rings_conding_sums(:,17),rings_conding_sums(:,18)./rings_conding_sums(:,19),...
    rings_conding_sums(:,20)./rings_conding_sums(:,21),rings_conding_sums(:,22)./rings_conding_sums(:,23)];

sel_idx=[(per_trl(:,1)-per_trl(:,2))./(per_trl(:,1)+per_trl(:,2)),...
    (per_trl(:,3)-per_trl(:,4))./(per_trl(:,3)+per_trl(:,4)),...
    (per_trl(:,5)-per_trl(:,6))./(per_trl(:,5)+per_trl(:,6))];

freq_sel=max(per_trl(:,1:2),[],2)>=0.5 & min(rings_conding_sums(:,[13,15,21,23]),[],2)>=5;

congru_sel=rings_conding_sums(:,1)==1;
cong_c=histcounts(sel_idx(congru_sel & freq_sel,1),-1:0.1:1,'Normalization','probability');
non_c=histcounts(abs(sel_idx(~congru_sel & freq_sel,1)),-1:0.1:1,'Normalization','probability');
cong_e=histcounts(sel_idx(congru_sel & freq_sel,3),-1:0.1:1,'Normalization','probability');

fh=figure('Color','w','Position',[32,32,225,225]);
hold on;
% bh=bar([cong_c;cong_e;non_c].',1);
bh=bar([cong_c;cong_e].',1);
bh(1).FaceColor='r';
bh(2).FaceColor='k';
% bh(3).FaceColor='c';
for ii=1:2, bh(ii).EdgeColor='none';end
set(gca,'XTick',0.5:5:20.5,'XTickLabel',-1:0.5:1)
xlabel('Selectivity index')
ylabel('Probability');

exportgraphics(fh,'loops_coding.pdf');

mm=[mean(sel_idx(congru_sel & freq_sel,1)),...
    nanmean(sel_idx(congru_sel & freq_sel,3))];%,...
    %mean(abs(sel_idx(~congru_sel & freq_sel,1)))];
ci=[bootci(500,@(x) mean(x),sel_idx(congru_sel & freq_sel,1)),...
    bootci(500,@(x) nanmean(x),sel_idx(congru_sel & freq_sel,3))];%,...
    %bootci(500,@(x) mean(x),abs(sel_idx(~congru_sel & freq_sel,1)))];

fh=figure('Color','w','Position',[32,32,225,225]);
hold on;
bar(1,mm(1),0.8,'FaceColor','r','EdgeColor','none');
bar(2,mm(2),0.8,'FaceColor','k','EdgeColor','none');
% bar(3,mm(3),0.8,'FaceColor','c','EdgeColor','none');
errorbar(1:2,mm,ci(1:2:3)-mm,ci(2:2:4)-mm,'k.','LineWidth',1,'CapSize',20)
xlim([0.5,2.5])
set(gca(),'XTick',[])
ylabel('Mean selectivity index')
exportgraphics(fh,'loops_coding_mean.pdf');

p=anova1([sel_idx(congru_sel & freq_sel,1);...
    sel_idx(congru_sel & freq_sel & ~isnan(sel_idx(:,3)),3);...
    abs(sel_idx(~congru_sel & freq_sel,1))],...
    [ones(nnz(congru_sel & freq_sel),1);...
    2*ones(nnz(congru_sel & freq_sel & ~isnan(sel_idx(:,3))),1);...
    3*ones(nnz(~congru_sel & freq_sel),1)])
    
