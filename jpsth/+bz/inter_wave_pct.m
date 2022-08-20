function fh=inter_wave_pct(pct_meta,opt)
arguments
    pct_meta
    opt.min_pair_per_session (1,1) double = 10
end
% selstr=load('perm_sens.mat','sens_meta');

% persistent sig pair

[sig,pair]=bz.load_sig_sums_conn_file('pair',true);
sig=bz.join_fc_waveid(sig,pct_meta.wave_id);
pair=bz.join_fc_waveid(pair,pct_meta.wave_id);
%>>>>>>>>>>>>>>>>>>Skip hierarchy>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
[sig_diff,sig_same,~,~]=bz.util.diff_at_level(sig.reg,'hierarchy',false);
[pair_diff,pair_same,~,~]=bz.util.diff_at_level(pair.reg,'hierarchy',false);
% 
% [samestats,sameci]=statsMixed(sig_same(:,5),pair_same(:,5),sig,pair,opt.min_pair_per_session);
% [diffstats,diffci]=statsMixed(sig_diff(:,5),pair_diff(:,5),sig,pair,opt.min_pair_per_session );

fh.samereg=stats_congru(sig_same(:,5),pair_same(:,5),sig,pair,opt.min_pair_per_session);
fh.crossreg=stats_congru(sig_diff(:,5),pair_diff(:,5),sig,pair,opt.min_pair_per_session );

% 
% stats.samestats=samestats;
% stats.sameci=sameci;
% stats.diffstats=diffstats;
% stats.diffci=diffci;
% 
% fh=struct();
% 
% fh.mat=figure('Color','w','Position',[32,32,1220,320]);
% t=tiledlayout(1,3);
% nexttile(2)
% ax_same=plotOne(samestats.','scale',[min(samestats,[],"all"),max(samestats,[],"all")]);
% nexttile(3)
% ax_diff=plotOne(diffstats.','scale',[min(diffstats,[],"all"),max(diffstats,[],"all")]);
% title(t,'Same-, cross-region FC rate')
% 
% fh.bar=figure('Color','w','Position',[100,100,1600,330]);
% t=tiledlayout(1,4);
% plotOneBar(t,samefromhat,samefromci);
% plotOneBar(t,sametohat,sametoci);
% plotOneBar(t,difffromhat,difffromci);
% plotOneBar(t,difftohat,difftoci);


% th=nexttile();
% ephys.util.figtable(fh,th,{'chisq-same';samep;'chisq-diff';diffp})

% exportgraphics(fh,'inter_wave_fc.pdf','ContentType','vector');
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end

function [cowavehat,cowaveci]=statsMixed(sig_sel,pair_sel,sig,pair,min_pair_per_session)
%% 4 Qdr FC rate loop
cowavehat=nan(5,5);
cowaveci=nan(5,5,2);

for fromqtr=0:4
    for toqtr=0:4
        sig_cowave=sig.sess(sig_sel & sig.waveid(:,1)==fromqtr & sig.waveid(:,2)==toqtr);
        pair_cowave=pair.sess(pair_sel & pair.waveid(:,1)==fromqtr & pair.waveid(:,2)==toqtr);
        %         [cowavehat(fromqtr,toqtr),cowaveci(fromqtr,toqtr,:)]=binofit(sig_cowave,pair_cowave);
        sess_vec=[histcounts(sig_cowave,(0:116)+0.5);histcounts(pair_cowave,(0:116)+0.5)];
        cowavehat(fromqtr+1,toqtr+1)=mean(sess_vec(1,sess_vec(2,:)>min_pair_per_session)./sess_vec(2,sess_vec(2,:)>min_pair_per_session));
    end
end
end



function fh=stats_congru(sig_sel,pair_sel,sig,pair,min_pair_per_session)
% [fromhat,tohat,fromsem,tosem]=deal(nan(5,1));

nonmem_sig=pct.su_pairs.get_nonmem(sig.waveid);
nonmem_pair=pct.su_pairs.get_nonmem(pair.waveid);

congru_sig=pct.su_pairs.get_congru(sig.waveid);
congru_pair=pct.su_pairs.get_congru(pair.waveid);

incong_sig=pct.su_pairs.get_incongru(sig.waveid);
incong_pair=pct.su_pairs.get_incongru(pair.waveid);

sig_nonmem=nnz(sig_sel & nonmem_sig);
pair_nonmem=nnz(pair_sel & nonmem_pair);
[nonmem_hat,nonmem_sem]=binofit(sig_nonmem,pair_nonmem);


% 
% sig_nonmem=sig.sess(sig_sel & nonmem_sig);
% pair_nonmem=pair.sess(pair_sel & nonmem_pair);
% sess_vec=[histcounts(sig_nonmem,(0:116)+0.5);histcounts(pair_nonmem,(0:116)+0.5)];
% inc_vec=sess_vec(1,sess_vec(2,:)>min_pair_per_session)./sess_vec(2,sess_vec(2,:)>min_pair_per_session);
% nonmem_hat=mean(inc_vec);
% nonmem_sem=std(inc_vec)./sqrt(numel(inc_vec));
% finivec=[reshape(inc_vec,[],1),repmat(0,numel(inc_vec),1)];

sig_congru=nnz(sig_sel & congru_sig);
pair_congru=nnz(pair_sel & congru_pair);
[congru_hat,congru_sem]=binofit(sig_congru,pair_congru);

% sig_congru=sig.sess(sig_sel & congru_sig);
% pair_congru=pair.sess(pair_sel & congru_pair);
% sess_vec=[histcounts(sig_congru,(0:116)+0.5);histcounts(pair_congru,(0:116)+0.5)];
% inc_vec=sess_vec(1,sess_vec(2,:)>min_pair_per_session)./sess_vec(2,sess_vec(2,:)>min_pair_per_session);
% congru_hat=mean(inc_vec);
% congru_sem=std(inc_vec)./sqrt(numel(inc_vec));
% finivec=[finivec;reshape(inc_vec,[],1),repmat(1,numel(inc_vec),1)];

sig_incong=nnz(sig_sel & incong_sig);
pair_incong=nnz(pair_sel & incong_pair);
[incong_hat,incong_sem]=binofit(sig_incong,pair_incong);
% sig_incong=sig.sess(sig_sel & incong_sig);
% pair_incong=pair.sess(pair_sel & incong_pair);
% sess_vec=[histcounts(sig_incong,(0:116)+0.5);histcounts(pair_incong,(0:116)+0.5)];
% inc_vec=sess_vec(1,sess_vec(2,:)>min_pair_per_session)./sess_vec(2,sess_vec(2,:)>min_pair_per_session);
% incong_hat=mean(inc_vec);
% incong_sem=std(inc_vec)./sqrt(numel(inc_vec));
% finivec=[finivec;reshape(inc_vec,[],1),repmat(2,numel(inc_vec),1)];

% [h,p]=ttest2(finivec(finivec(:,2)==0,1),finivec(finivec(:,2)==2,1))
% 
% p=anovan(finivec(:,1),finivec(:,2),'display','off');
% disp("From bars anova p="+num2str(p));

fh=figure('Color','w','Position',[100,100,235,235]);
hold on
bh=bar(1:3,diag([nonmem_hat,incong_hat,congru_hat]),'stacked');
bh(1).FaceColor='k';
bh(2).FaceColor='b';
bh(3).FaceColor='r';
errorbar(1:3,[nonmem_hat,incong_hat,congru_hat],...
    [nonmem_sem(1),incong_sem(1),congru_sem(1)]-[nonmem_hat,incong_hat,congru_hat],...
    [nonmem_sem(2),incong_sem(2),congru_sem(2)]-[nonmem_hat,incong_hat,congru_hat],...
    'k.');

[~,~,p]=crosstab([zeros(pair_nonmem,1);ones(pair_congru,1);2.*ones(pair_incong,1)],...
    [(1:pair_nonmem)>sig_nonmem,(1:pair_congru)>sig_congru,(1:pair_incong)>sig_incong].');

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
