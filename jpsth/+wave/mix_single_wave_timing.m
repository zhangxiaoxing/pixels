global_init;
sumeta=ephys.util.load_meta('skip_stats',true);
wrs_mux_meta=ephys.get_wrs_mux_meta();
com_map=wave.get_pct_com_map(wrs_mux_meta,'curve',true);

% function com_map2meta_like()
pref2waveid=containers.Map({'s1d3','s1d6','s2d3','s2d6','olf_s1','olf_s2','dur_d3','dur_d6'},...
    num2cell(1:8));
com_meta=struct();
[com_meta.sess,com_meta.allcid,com_meta.wave_id,com_meta.com]=deal([]);
com_meta.reg_tree=cell(0);
sess=fieldnames(com_map);
for si=1:numel(sess)
    sessid=str2double(replace(sess{si},'s',''));
    sesscid=sumeta.allcid(sumeta.sess==sessid);
    sessreg=sumeta.reg_tree(:,sumeta.sess==sessid);
    prefs=fieldnames(com_map.(sess{si}));
    for pi=1:numel(prefs)
        cmap=com_map.(sess{si}).(prefs{pi}).com;
        keys=cell2mat(cmap.keys());
        values=cell2mat(cmap.values());
        com_meta.sess=[com_meta.sess;...
            repmat(sessid,numel(keys),1)];
        com_meta.allcid=[com_meta.allcid;...
            reshape(keys,[],1)];
        com_meta.com=[com_meta.com;...
            reshape(values,[],1)];
        com_meta.wave_id=[com_meta.wave_id;...
            repmat(pref2waveid(prefs{pi}),numel(keys),1)];

        [~,sidx]=ismember(uint16(keys),sesscid);
        com_meta.reg_tree=[com_meta.reg_tree,sessreg(:,sidx)];
    end
end

% end

%% mean
mm=[mean(com_meta.com(ismember(com_meta.wave_id,5:6))),...
mean(com_meta.com(ismember(com_meta.wave_id,7:8))),...
mean(com_meta.com(ismember(com_meta.wave_id,1:4)))];

sem=[std(com_meta.com(ismember(com_meta.wave_id,5:6)))./sqrt(nnz(ismember(com_meta.wave_id,5:6))),...
mean(com_meta.com(ismember(com_meta.wave_id,7:8)))./sqrt(nnz(ismember(com_meta.wave_id,7:8))),...
mean(com_meta.com(ismember(com_meta.wave_id,1:4)))./sqrt(nnz(ismember(com_meta.wave_id,1:4)))];
figure()
hold on
bh=bar(diag(mm),'stacked');
errorbar(bh(3).XEndPoints,bh(3).YEndPoints,sem,'k.')
xlim([0.5,3.5])
set(gca(),'YTick',0:2:12,'YTickLabel',0:0.5:3,...
    'XTick',1:3,'XTickLabel',{'Olf','Dur','Mix'})
ylabel('Wave-specific mean TCOM (s)')
[bh.FaceColor]=deal('r','b','w');
anovap=anovan([com_meta.com(ismember(com_meta.wave_id,5:6));...
    com_meta.com(ismember(com_meta.wave_id,7:8));...
    com_meta.com(ismember(com_meta.wave_id,1:4))], ...
    [zeros(nnz(ismember(com_meta.wave_id,5:6)),1);...
    ones(nnz(ismember(com_meta.wave_id,7:8)),1);...
    2*ones(nnz(ismember(com_meta.wave_id,1:4)),1)],...
    'display','off');

omp=ranksum(com_meta.com(ismember(com_meta.wave_id,1:4)),...
    com_meta.com(ismember(com_meta.wave_id,5:6)));
odp=ranksum(com_meta.com(ismember(com_meta.wave_id,7:8)),...
    com_meta.com(ismember(com_meta.wave_id,5:6)));
mdp=ranksum(com_meta.com(ismember(com_meta.wave_id,7:8)),...
    com_meta.com(ismember(com_meta.wave_id,1:4)));
title(sprintf('anova %.3f, wrs %.3f,%.3f,%.3f',anovap,omp,odp,mdp))
ylim([4,8])


%% cdf in specific selective subpopulations
mixcdf=histcounts(com_meta.com(ismember(com_meta.wave_id,1:4)),0:0.25:12,'Normalization','cdf');
olfcdf=histcounts(com_meta.com(ismember(com_meta.wave_id,5:6)),0:0.25:12,'Normalization','cdf');
durcdf=histcounts(com_meta.com(ismember(com_meta.wave_id,7:8)),0:0.25:12,'Normalization','cdf');

figure()
hold on
mh=plot(mixcdf,'-k');
oh=plot(olfcdf,'-r');
dh=plot(durcdf,'-b');
set(gca(),'XTick',0:16:48,'XTickLabel',0:1:3,'YLim',[0,1]);

ylabel('Cumulated density (%)')
xlabel('Time (s)')
legend([oh,dh,mh],{'Olfactory','Duration','Mixed'},'Location','northoutside','Orientation','horizontal')

%% percentage of total population

mixcuc=histcounts(com_meta.com(ismember(com_meta.wave_id,1:4)),0:0.25:12,'Normalization','cumcount');
olfcuc=histcounts(com_meta.com(ismember(com_meta.wave_id,5:6)),0:0.25:12,'Normalization','cumcount');
durcuc=histcounts(com_meta.com(ismember(com_meta.wave_id,7:8)),0:0.25:12,'Normalization','cumcount');

figure()
hold on
mh=plot(mixcuc,'-k');
oh=plot(olfcuc,'-r');
dh=plot(durcuc,'-b');
ylim([0,3303]);
set(gca(),'XTick',0:16:48,'XTickLabel',0:1:3,...
    'YTick',0:(330.3*5):3303,'YTickLabel',0:5:10);
xlabel('Time (s)');
ylabel('Percentage of total neuron (%)')
xlabel('Time (s)')


ylabel('Cumulated density (%)')
xlabel('Time (s)')
legend([oh,dh,mh],{'Olfactory','Duration','Mixed'},'Location','northoutside','Orientation','horizontal')




%%
% olf_com_per_reg=struct();
% [olf_com_per_reg.collection,olf_com_per_reg.com_meta]=wave.per_region_COM(...
%     com_map,'sel_type','olf');
% % [mean/media,reg,depth,count,std/iqr]
% 
% depth_sel=cell2mat(olf_com_per_reg.collection(:,3))==5 ...
%     & ~ismissing(olf_com_per_reg.collection(:,2));
% sortrows(olf_com_per_reg.collection(depth_sel,:),1)
%%



