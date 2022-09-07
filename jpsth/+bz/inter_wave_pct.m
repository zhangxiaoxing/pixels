function fh=inter_wave_pct(pct_meta,opt)
arguments
    pct_meta
    opt.min_pair_per_session (1,1) double = 20
    opt.per_sess (1,1) logical = false
end
% selstr=load('perm_sens.mat','sens_meta');

% persistent sig pair

[sig,pair]=bz.load_sig_sums_conn_file('pair',true);
sig=bz.join_fc_waveid(sig,pct_meta.wave_id);
sig=bz.join_fc_waveid(sig,pct_meta.mat_coord,'pct_mat',true); %lead_olf, lead_dur,folo_olf,folo_dur
pair=bz.join_fc_waveid(pair,pct_meta.wave_id);
pair=bz.join_fc_waveid(pair,pct_meta.mat_coord,"pct_mat",true);
%>>>>>>>>>>>>>>>>>>Skip hierarchy>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
[sig_diff,sig_same,~,~]=bz.util.diff_at_level(sig.reg,'hierarchy',false);
[pair_diff,pair_same,~,~]=bz.util.diff_at_level(pair.reg,'hierarchy',false);
% 
% [samestats,sameci]=statsMixed(sig_same(:,5),pair_same(:,5),sig,pair,opt.min_pair_per_session);
% [diffstats,diffci]=statsMixed(sig_diff(:,5),pair_diff(:,5),sig,pair,opt.min_pair_per_session );

fh.samereg=stats_congru(sig_same(:,5),pair_same(:,5),sig,pair,opt.min_pair_per_session,opt.per_sess);
fh.crossreg=stats_congru(sig_diff(:,5),pair_diff(:,5),sig,pair,opt.min_pair_per_session,opt.per_sess);

%% single-mixed
statsMat(sig,pair)
statsMixed(sig,pair)

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end


function stats=statsMat(sig,pair)

stats=struct();
% from 
struc_from_tag={'from_mix','from_olf','from_dur','from_nonmem'};
lead_member={1:4,5:6,7:8,0};
for ii=1:4
    pair_cnt=nan(5,5);
    sig_cnt=nan(5,5);

    pair_lead_sel=ismember(pair.waveid(:,1),lead_member{ii});
    sig_lead_sel=ismember(sig.waveid(:,1),lead_member{ii});
    for oidx=1:5
        for didx=1:5
            pair_coord_sel=pair.pct_coord(:,3)==oidx ...
                & pair.pct_coord(:,4)==didx;
            sig_coord_sel=sig.pct_coord(:,3)==oidx ...
                & sig.pct_coord(:,4)==didx;
            pair_cnt(oidx,didx)=nnz(pair_lead_sel & pair_coord_sel);
            sig_cnt(oidx,didx)=nnz(sig_lead_sel & sig_coord_sel);
        end
    end
    stats.(struc_from_tag{ii}).sig_cnt=sig_cnt;
    stats.(struc_from_tag{ii}).pair_cnt=pair_cnt;
end
%to
struc_to_tag={'to_mix','to_olf','to_dur','to_nonmem'};
lead_member={1:4,5:6,7:8,0};
for ii=1:4
    pair_cnt=nan(5,5);
    sig_cnt=nan(5,5);

    pair_lead_sel=ismember(pair.waveid(:,2),lead_member{ii});
    sig_lead_sel=ismember(sig.waveid(:,2),lead_member{ii});
    for oidx=1:5
        for didx=1:5
            pair_coord_sel=pair.pct_coord(:,1)==oidx ...
                & pair.pct_coord(:,2)==didx;
            sig_coord_sel=sig.pct_coord(:,1)==oidx ...
                & sig.pct_coord(:,2)==didx;
            pair_cnt(oidx,didx)=nnz(pair_lead_sel & pair_coord_sel);
            sig_cnt(oidx,didx)=nnz(sig_lead_sel & sig_coord_sel);
        end
    end
    stats.(struc_to_tag{ii}).sig_cnt=sig_cnt;
    stats.(struc_to_tag{ii}).pair_cnt=pair_cnt;
end

%% matrix
figure('Color','w')
tiledlayout(2,4);

for ii=1:4
    nexttile(ii)
    imagesc((stats.(struc_from_tag{ii}).sig_cnt./stats.(struc_from_tag{ii}).pair_cnt).*100,[0.4,1]);
    title(struc_from_tag{ii},'Interpreter','none')
    colormap(flip(colormap('gray')));
    set(gca(),'YDir','normal','XTick',0.5:1:10.5,'XTickLabel',0:20:100,...
    'YTick',0.5:1:10.5,'YTickLabel',0:20:100);
    xlabel('Duration rank (%)')
    ylabel('Olfactory rank (%)')
    colorbar()
end

for ii=1:4
    nexttile(4+ii)
    imagesc((stats.(struc_to_tag{ii}).sig_cnt./stats.(struc_to_tag{ii}).pair_cnt).*100,[0.4,1]);
    title(struc_to_tag{ii},'Interpreter','none')
    colormap(flip(colormap('gray')));
    set(gca(),'YDir','normal','XTick',0.5:1:5.5,'XTickLabel',0:20:100,...
    'YTick',0.5:1:5.5,'YTickLabel',0:20:100);
    xlabel('Duration rank (%)')
    ylabel('Olfactory rank (%)')
    colorbar()
end

%% margin
figure('Color','w')
tiledlayout(4,4);

for ii=1:4
    nexttile(ii)
    bar(sum(stats.(struc_from_tag{ii}).sig_cnt)./sum(stats.(struc_from_tag{ii}).pair_cnt).*100,'grouped','black');
    title('collapse olf')
    xlabel('to dur rank')
    ylabel('fc rate')
    ylim([0,1.05])
    nexttile(4+ii)
    bar(sum(stats.(struc_from_tag{ii}).sig_cnt,2)./sum(stats.(struc_from_tag{ii}).pair_cnt,2).*100,'grouped','black','Horizontal','on');
    title('collapse dur')
    xlim([0,1.05])
    ylabel('to olf rank')
    xlabel('fc rate')
end

for ii=1:4
    nexttile(8+ii)
    bar(sum(stats.(struc_to_tag{ii}).sig_cnt)./sum(stats.(struc_to_tag{ii}).pair_cnt).*100,'grouped','black');
    title('collapse olf')
    xlabel('from dur rank')
    ylabel('fc rate')
    ylim([0,1.05])
    nexttile(12+ii)
    bar(sum(stats.(struc_to_tag{ii}).sig_cnt,2)./sum(stats.(struc_to_tag{ii}).pair_cnt,2).*100,'grouped','black','Horizontal','on');
    title('collapse dur')
    ylabel('from olf rank')
    xlabel('fc rate')
    xlim([0,1.05])
end

end


function statsMixed(sig,pair)
%% 4 Qdr FC rate loop

pair_sel_olf=((pair.waveid(:,1)==1 | pair.waveid(:,1)==2) & pair.waveid(:,2)==5) ...
    |((pair.waveid(:,1)==3 | pair.waveid(:,1)==4) & pair.waveid(:,2)==6); 
% pair info was flipped and duplicated, so mix-> single always equals single->mix
% in number.

pair_sel_dur=((pair.waveid(:,1)==1 | pair.waveid(:,1)==3) & pair.waveid(:,2)==7) ...
    |((pair.waveid(:,1)==2 | pair.waveid(:,1)==4) & pair.waveid(:,2)==8);

ms_sig_sel_olf=((sig.waveid(:,1)==1 | sig.waveid(:,1)==2) & sig.waveid(:,2)==5) ...
    |((sig.waveid(:,1)==3 | sig.waveid(:,1)==4) & sig.waveid(:,2)==6);

ms_sig_sel_dur=((sig.waveid(:,1)==1 | sig.waveid(:,1)==3) & sig.waveid(:,2)==7) ...
    |((sig.waveid(:,1)==2 | sig.waveid(:,1)==4) & sig.waveid(:,2)==8);

sm_sig_sel_olf=((sig.waveid(:,2)==1 | sig.waveid(:,2)==2) & sig.waveid(:,1)==5) ...
    |((sig.waveid(:,2)==3 | sig.waveid(:,2)==4) & sig.waveid(:,1)==6);

sm_sig_sel_dur=((sig.waveid(:,2)==1 | sig.waveid(:,2)==3) & sig.waveid(:,1)==7) ...
    |((sig.waveid(:,2)==2 | sig.waveid(:,2)==4) & sig.waveid(:,1)==8);


olfstats=[nnz(sm_sig_sel_olf),nnz(ms_sig_sel_olf),nnz(pair_sel_olf)];
durstats=[nnz(sm_sig_sel_dur),nnz(ms_sig_sel_dur),nnz(pair_sel_dur)];

%% mix-mix, single-single
mm_pair=pair.waveid(:,1)==pair.waveid(:,2) & ismember(pair.waveid(:,1),1:4);
mm_sig=sig.waveid(:,1)==sig.waveid(:,2) & ismember(sig.waveid(:,1),1:4);

olf_pair=pair.waveid(:,1)==pair.waveid(:,2) & ismember(pair.waveid(:,1),5:6);
olf_sig=sig.waveid(:,1)==sig.waveid(:,2) & ismember(sig.waveid(:,1),5:6);

dur_pair=pair.waveid(:,1)==pair.waveid(:,2) & ismember(pair.waveid(:,1),7:8);
dur_sig=sig.waveid(:,1)==sig.waveid(:,2) & ismember(sig.waveid(:,1),7:8);

[mm_hat,mm_ci]=binofit(nnz(mm_sig),nnz(mm_pair));
[olf_hat,olf_ci]=binofit(nnz(olf_sig),nnz(olf_pair));
[dur_hat,dur_ci]=binofit(nnz(dur_sig),nnz(dur_pair));

mmm=[mm_hat,olf_hat,dur_hat];

[~,smci_olf]=binofit(olfstats(1),olfstats(3));
[~,msci_olf]=binofit(olfstats(2),olfstats(3));

[~,smci_dur]=binofit(durstats(1),durstats(3));
[~,msci_dur]=binofit(durstats(2),durstats(3));

[~,~,polf]=crosstab([zeros(1,nnz(olf_pair)),ones(1,olfstats(3)),ones(1,olfstats(3)).*2],...
    [1:nnz(olf_pair)>nnz(olf_sig),1:olfstats(3)>olfstats(1),1:olfstats(3)>olfstats(2)]);

[~,~,pdur]=crosstab([zeros(1,nnz(dur_pair)),ones(1,durstats(3)),ones(1,durstats(3)).*2],...
    [1:nnz(dur_pair)>nnz(dur_sig),1:durstats(3)>durstats(1),1:durstats(3)>durstats(2)]);

[~,~,povall]=crosstab([zeros(1,nnz(olf_pair)),ones(1,olfstats(3)),ones(1,olfstats(3)).*2,ones(1,nnz(dur_pair)).*3,ones(1,durstats(3)).*4,ones(1,durstats(3)).*5,ones(1,nnz(mm_pair)).*6],...
    [1:nnz(olf_pair)>nnz(olf_sig),1:olfstats(3)>olfstats(1),1:olfstats(3)>olfstats(2),1:nnz(dur_pair)>nnz(dur_sig),1:durstats(3)>durstats(1),1:durstats(3)>durstats(2),1:(nnz(mm_pair))>nnz(mm_sig)]);


figure('Color','w','Position',[100,100,275,235])
hold on
bh1=bar(1,mmm(1),0.25);
bh2=bar(2,[mmm(2),olfstats(1:2)./olfstats(3)]);
bh3=bar(3,[mmm(3),durstats(1:2)./durstats(3)]);
errorbar([bh1.XEndPoints,bh2.XEndPoints,bh3.XEndPoints],...
    [bh1.YEndPoints,bh2.YEndPoints,bh3.YEndPoints],...
    [mm_ci(1),olf_ci(1),smci_olf(1),msci_olf(1),dur_ci(1),smci_dur(1),msci_dur(1)]-[bh1.YEndPoints,bh2.YEndPoints,bh3.YEndPoints],...
    [mm_ci(2),olf_ci(2),smci_olf(2),msci_olf(2),dur_ci(2),smci_dur(2),msci_dur(2)]-[bh1.YEndPoints,bh2.YEndPoints,bh3.YEndPoints],...
    'k.');
    
[bh2(1).FaceColor,bh3(1).FaceColor]=deal('k');
[bh2(2).FaceColor,bh3(2).FaceColor]=deal([0.5,0.5,0.5]);
[bh2(3).FaceColor,bh3(3).FaceColor]=deal('w');

% legend(bh,{'Single->Multiplexed','Multiplexed->Single'},'Location','northoutside')
set(gca(),'XTick',1:3,'XTickLabel',{'m2m','Olfactory','Duration'},'YTickLabel',get(gca(),'YTick').*100)
ylabel('F.C. rate (%)')
title (sprintf('p=%.3f,%.3f,%.3f',povall,polf,pdur));










%%
% cowavehat=nan(5,5);
% cowaveci=nan(5,5,2);
% 
% for fromqtr=0:4
%     for toqtr=0:4
%         sig_cowave=sig.sess(sig_sel & sig.waveid(:,1)==fromqtr & sig.waveid(:,2)==toqtr);
%         pair_cowave=pair.sess(pair_sel & pair.waveid(:,1)==fromqtr & pair.waveid(:,2)==toqtr);
%         %         [cowavehat(fromqtr,toqtr),cowaveci(fromqtr,toqtr,:)]=binofit(sig_cowave,pair_cowave);
%         sess_vec=[histcounts(sig_cowave,(0:116)+0.5);histcounts(pair_cowave,(0:116)+0.5)];
%         cowavehat(fromqtr+1,toqtr+1)=mean(sess_vec(1,sess_vec(2,:)>min_pair_per_session)./sess_vec(2,sess_vec(2,:)>min_pair_per_session));
%     end
% end
end



function fh=stats_congru(sig_sel,pair_sel,sig,pair,min_pair_per_session,per_sess)
% [fromhat,tohat,fromsem,tosem]=deal(nan(5,1));

nonmem_sig=pct.su_pairs.get_nonmem(sig.waveid);
nonmem_pair=pct.su_pairs.get_nonmem(pair.waveid);

congru_sig=pct.su_pairs.get_congru(sig.waveid);
congru_pair=pct.su_pairs.get_congru(pair.waveid);

incong_sig=pct.su_pairs.get_incongru(sig.waveid);
incong_pair=pct.su_pairs.get_incongru(pair.waveid);

if per_sess
    sig_nonmem=sig.sess(sig_sel & nonmem_sig);
    pair_nonmem=pair.sess(pair_sel & nonmem_pair);
    sess_vec=[histcounts(sig_nonmem,(0:116)+0.5);histcounts(pair_nonmem,(0:116)+0.5)];
    inc_vec=sess_vec(1,sess_vec(2,:)>min_pair_per_session)./sess_vec(2,sess_vec(2,:)>min_pair_per_session);
    nonmem_hat=mean(inc_vec);
    nonmem_sem=std(inc_vec)./sqrt(numel(inc_vec));
    finivec=[reshape(inc_vec,[],1),repmat(0,numel(inc_vec),1)];

    sig_congru=sig.sess(sig_sel & congru_sig);
    pair_congru=pair.sess(pair_sel & congru_pair);
    sess_vec=[histcounts(sig_congru,(0:116)+0.5);histcounts(pair_congru,(0:116)+0.5)];
    inc_vec=sess_vec(1,sess_vec(2,:)>min_pair_per_session)./sess_vec(2,sess_vec(2,:)>min_pair_per_session);
    congru_hat=mean(inc_vec);
    congru_sem=std(inc_vec)./sqrt(numel(inc_vec));
    finivec=[finivec;reshape(inc_vec,[],1),repmat(1,numel(inc_vec),1)];

    sig_incong=sig.sess(sig_sel & incong_sig);
    pair_incong=pair.sess(pair_sel & incong_pair);
    sess_vec=[histcounts(sig_incong,(0:116)+0.5);histcounts(pair_incong,(0:116)+0.5)];
    inc_vec=sess_vec(1,sess_vec(2,:)>min_pair_per_session)./sess_vec(2,sess_vec(2,:)>min_pair_per_session);
    incong_hat=mean(inc_vec);
    incong_sem=std(inc_vec)./sqrt(numel(inc_vec));
    finivec=[finivec;reshape(inc_vec,[],1),repmat(2,numel(inc_vec),1)];

%     [h,p]=ttest2(finivec(finivec(:,2)==0,1),finivec(finivec(:,2)==2,1))

    p=anovan(finivec(:,1),finivec(:,2),'display','off');
    disp("From bars anova p="+num2str(p));

else
    sig_nonmem=nnz(sig_sel & nonmem_sig);
    pair_nonmem=nnz(pair_sel & nonmem_pair);
    [nonmem_hat,nonmem_sem]=binofit(sig_nonmem,pair_nonmem);

    sig_congru=nnz(sig_sel & congru_sig);
    pair_congru=nnz(pair_sel & congru_pair);
    [congru_hat,congru_sem]=binofit(sig_congru,pair_congru);

    sig_incong=nnz(sig_sel & incong_sig);
    pair_incong=nnz(pair_sel & incong_pair);
    [incong_hat,incong_sem]=binofit(sig_incong,pair_incong);

    [~,~,p]=crosstab([zeros(pair_nonmem,1);ones(pair_congru,1);2.*ones(pair_incong,1)],...
    [(1:pair_nonmem)>sig_nonmem,(1:pair_congru)>sig_congru,(1:pair_incong)>sig_incong].');
end


fh=figure('Color','w','Position',[100,100,235,235]);
hold on
bh=bar(1:3,diag([nonmem_hat,incong_hat,congru_hat]),'stacked');
bh(1).FaceColor='k';
bh(2).FaceColor='b';
bh(3).FaceColor='r';
if per_sess
    errorbar(1:3,[nonmem_hat,incong_hat,congru_hat],...
        [nonmem_sem,incong_sem,congru_sem],...
        'k.');
else
errorbar(1:3,[nonmem_hat,incong_hat,congru_hat],...
    [nonmem_sem(1),incong_sem(1),congru_sem(1)]-[nonmem_hat,incong_hat,congru_hat],...
    [nonmem_sem(2),incong_sem(2),congru_sem(2)]-[nonmem_hat,incong_hat,congru_hat],...
    'k.');
end

set(gca(),'XTick',1:3,'XTickLabel',{'NONMEM','INCONG','CONGRU'},'XTickLabelRotation',90,...
    'YTickLabel',100.*get(gca(),'YTick'));

title("chisq-p "+num2str(p))
xlim([0.4,3.6]);
ylabel('F.C. rate (%)')

end

function [fromhat,fromsem,tohat,tosem]=statsX(sig_sel,pair_sel,sig,pair,min_pair_per_session)
[fromhat,tohat,fromsem,tosem]=deal(nan(5,1));
anovastats=[];
for fromqtr=0:8
    %% UNFINISHED
    sig_cowave=sig.sess(sig_sel & sig.waveid(:,1)==fromqtr & ismember(sig.waveid(:,2),0:8));
    pair_cowave=pair.sess(pair_sel & pair.waveid(:,1)==fromqtr & ismember(pair.waveid(:,2),0:8));
    %     [fromhat(fromqtr),fromci(fromqtr,:)]=binofit(sig_cowave,pair_cowave);
    sess_vec=[histcounts(sig_cowave,(0:116)+0.5);histcounts(pair_cowave,(0:116)+0.5)];
    fromhat(fromqtr+1)=mean(sess_vec(1,sess_vec(2,:)>min_pair_per_session)./sess_vec(2,sess_vec(2,:)>min_pair_per_session));
    fromsem(fromqtr+1)=std(sess_vec(1,sess_vec(2,:)>min_pair_per_session)./sess_vec(2,sess_vec(2,:)>min_pair_per_session))./sqrt(nnz(sess_vec(2,:)>min_pair_per_session));
    finivec=reshape(sess_vec(1,sess_vec(2,:)>min_pair_per_session),[],1);
    anovastats=cat(1,anovastats,[finivec,fromqtr.*ones(numel(finivec),1)]);
end
p=anovan(anovastats(:,1),anovastats(:,2),'display','off');
disp("From bars anova p="+num2str(p))

anovastats=[];
for toqtr=0:8
    sig_cowave=sig.sess(sig_sel  & ismember(sig.waveid(:,1),0:8) & sig.waveid(:,2)==toqtr);
    pair_cowave=pair.sess(pair_sel  & ismember(pair.waveid(:,1),0:8) & pair.waveid(:,2)==toqtr);
    sess_vec=[histcounts(sig_cowave,(0:116)+0.5);histcounts(pair_cowave,(0:116)+0.5)];
    tohat(toqtr+1)=mean(sess_vec(1,sess_vec(2,:)>min_pair_per_session)./sess_vec(2,sess_vec(2,:)>min_pair_per_session));
    tosem(toqtr+1)=std(sess_vec(1,sess_vec(2,:)>min_pair_per_session)./sess_vec(2,sess_vec(2,:)>min_pair_per_session))./sqrt(nnz(sess_vec(2,:)>min_pair_per_session));
    finivec=reshape(sess_vec(1,sess_vec(2,:)>min_pair_per_session),[],1);
    anovastats=cat(1,anovastats,[finivec,toqtr.*ones(numel(finivec),1)]);
end

p=anovan(anovastats(:,1),anovastats(:,2),'display','off');
disp("To bars anova p="+num2str(p))


%% skipped statistics for now
% data=[cowavehat,cowaveci,crosswavehat,crosswaveci,nonmemhat,nonmemci];
%
% [~,~,p]=crosstab([ones(pair_cowave,1); ...
%     4*ones(pair_cross_wave,1);...
%     5*ones(pair_nonmem,1)],...
%     [(1:(pair_cowave))>(sig_cowave),...
%     (1:pair_cross_wave)>sig_cross_wave, ...
%     (1:pair_nonmem)>sig_nonmem].');
% stats=p;
end

function ax=plotOne(stats,opt)
arguments
    stats
    opt.scale
end
ax=imagesc(stats,opt.scale);
cbh=colorbar();
cbh.Label.String='F.C. Rate (%)';
cbh.TickLabels=cbh.Ticks*100;
set(gca(),'YDir','normal')
colormap(flip(colormap('gray')))
xlabel('From class # neuron');
ylabel('To class # neuron');
% set(gca(),'XTick',1:4,'YTick',1:4)
end

function plotOneBar(t,stats,ci)
nexttile(t)
hold on
bh=bar(stats,'FaceColor','w');
errorbar(1:9,stats,ci,'k.');
set(gca(),'YTickLabel',get(gca(),'YTick').*100);
set(gca(),'XTick',1:4,'XTickLabel',{'Non','Olf','Dur','Mix'});
end
