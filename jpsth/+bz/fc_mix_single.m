keyboard();
global_init;
meta=ephys.util.load_meta('skip_stats',true);
wrs_mux_meta=ephys.get_wrs_mux_meta();
com_map=wave.get_pct_com_map(wrs_mux_meta,'curve',true);

%%
fc_rate_mix_simple(wrs_mux_meta,'condense_plot',true)

%%
FC_TCOM(wrs_mux_meta,com_map,tcom_maps)

if false
    olfgrp=struct('tag','OLF','lst',{{'AON','TT','DP','PIR'}});
    pfcgrp=struct('tag','PFC','lst',{{'PL','ILA','ACA','ORB','AI'}});
    thsgrp=struct('tag','THS','lst',{{'VENT','ILM','LAT','MED'}});

    fc_rate_mix_simple(wrs_mux_meta,olfgrp,olfgrp)
    fc_rate_mix_simple(wrs_mux_meta,pfcgrp,pfcgrp)
    fc_rate_mix_simple(wrs_mux_meta,thsgrp,thsgrp)

    fc_rate_mix_simple(wrs_mux_meta,olfgrp,pfcgrp)
    fc_rate_mix_simple(wrs_mux_meta,pfcgrp,thsgrp)
    fc_rate_mix_simple(wrs_mux_meta,olfgrp,thsgrp)
end

if false
    [olf_com_per_reg.collection,olf_com_per_reg.com_meta]=wave.per_region_COM(...
        com_map,'sel_type','olf');
    depth_sel=cell2mat(olf_com_per_reg.collection(:,3))==5 ...
        & ~ismissing(olf_com_per_reg.collection(:,2));
    reg_tcom=sortrows(olf_com_per_reg.collection(depth_sel,:),1);
    subtotal=sum([reg_tcom{:,4}]);
    cumu_cnt=arrayfun(@(x) sum([reg_tcom{1:x,4}]),1:size(reg_tcom,1));

    %% 3 zones
    [min1_3,min1_3id]=min(abs(cumu_cnt-subtotal./3));
    [min2_3,min2_3id]=min(abs(cumu_cnt-2.*subtotal./3));

    early3grp=struct('tag','early 3 groups','lst',{reg_tcom(1:min1_3id,2)});
    mid3grp=struct('tag','middle 3 groups','lst',{reg_tcom((min1_3id+1):min2_3id,2)});
    late3grp=struct('tag','late 3 groups','lst',{reg_tcom((min2_3id+1):end,2)});

    fc_rate_mix_simple(wrs_mux_meta,early3grp,early3grp)
    fc_rate_mix_simple(wrs_mux_meta,mid3grp,mid3grp)
    fc_rate_mix_simple(wrs_mux_meta,late3grp,late3grp)

    fc_rate_mix_simple(wrs_mux_meta,early3grp,mid3grp)
    fc_rate_mix_simple(wrs_mux_meta,mid3grp,late3grp)
    fc_rate_mix_simple(wrs_mux_meta,early3grp,late3grp)
    %% 2 zones
    [~,half_idx]=min(abs(cumu_cnt-subtotal./2));

    early2grp=struct('tag','early 2 groups','lst',{reg_tcom(1:half_idx,2)});
    late2grp=struct('tag','late 2 groups','lst',{reg_tcom((half_idx+1):end,2)});

    fc_rate_mix_simple(wrs_mux_meta,early2grp,early2grp)
    fc_rate_mix_simple(wrs_mux_meta,late2grp,late2grp)
    fc_rate_mix_simple(wrs_mux_meta,early2grp,late2grp)
end


%%

function FC_TCOM(wrs_mux_meta,com_map,tcom_maps)
%% timing
fc_stats=fc.fc_com_reg_wave(wrs_mux_meta,com_map,tcom_maps); % produce two figures
congrusel=pct.su_pairs.get_congru(fc_stats(:,10:11));

mix_olf_sel=ismember(fc_stats(:,10),1:4) & ismember(fc_stats(:,11),5:6) & congrusel;
olf_mix_sel=ismember(fc_stats(:,11),1:4) & ismember(fc_stats(:,10),5:6) & congrusel;
mix_dur_sel=ismember(fc_stats(:,10),1:4) & ismember(fc_stats(:,11),7:8) & congrusel;
dur_mix_sel=ismember(fc_stats(:,11),1:4) & ismember(fc_stats(:,10),7:8) & congrusel;

olf_dur_sel=ismember(fc_stats(:,10),5:6) & ismember(fc_stats(:,11),7:8); % undefined congruent-incongruent relationship
dur_olf_sel=ismember(fc_stats(:,11),5:6) & ismember(fc_stats(:,10),7:8); % likewise ^^

mix_dur_tcom=fc_stats(mix_dur_sel,4:5);
mddf=diff(mix_dur_tcom,1,2);

mix_olf_tcom=fc_stats(mix_olf_sel,4:5);
modf=diff(mix_olf_tcom,1,2);

dur_mix_tcom=fc_stats(dur_mix_sel,4:5);
dmdf=diff(dur_mix_tcom,1,2);

olf_mix_tcom=fc_stats(olf_mix_sel,4:5);
omdf=diff(olf_mix_tcom,1,2);

olf_dur_tcom=fc_stats(olf_dur_sel,4:5);
oddf=diff(olf_dur_tcom,1,2);

dur_olf_tcom=fc_stats(dur_olf_sel,4:5);
dodf=diff(dur_olf_tcom,1,2);

mon=nnz(mix_olf_sel);
omn=nnz(olf_mix_sel);
pmo=binocdf(min([mon;omn]),sum([mon;omn]),0.5)*2;

mdn=nnz(mix_dur_sel);
dmn=nnz(dur_mix_sel);
pmd=binocdf(min([mdn;dmn]),sum([mdn;dmn]),0.5)*2;

odn=nnz(olf_dur_sel);
don=nnz(dur_olf_sel);
pod=binocdf(min([odn;don]),sum([odn;don]),0.5)*2;

if false %latency histogram
    figure()
    tiledlayout(1,2)
    nexttile()
    histogram(mddf)
    nexttile()
    histogram(modf)
end
if false %TCOM pdf cdf curve
    figure()
    hold on
    edges=-10:1:10;
    mdcdf=histcounts(mddf,edges,'Normalization','cdf');
    mocdf=histcounts(modf,edges,'Normalization','cdf');
    dmcdf=histcounts(dmdf,edges,'Normalization','cdf');
    omcdf=histcounts(omdf,edges,'Normalization','cdf');


    plot(edges,[0,mdcdf],'-bo');
    plot(edges,[0,mocdf],'-ro');
    plot(edges,[0,dmcdf],'--bo');
    plot(edges,[0,omcdf],'--ro');
    yline(0.5,'k:');
    xline(0,'k:');
end
% rratio=@(x) x./(sum(x,'all'));
ebarsort=@(x) subsref(x,struct(type={'()'},subs={{[1,3,2,4]}}));
%%
figure() % mix vs olfactory, mo,om
tiledlayout(1,3)
nexttile()
mocount=[nnz(modf>0),nnz(modf<=0)];
omcount=[nnz(omdf>0),nnz(omdf<=0)];
pmo=binocdf(min(mocount),numel(modf),0.5)*2;
pom=binocdf(min(omcount),numel(omdf),0.5)*2;
[mohat,moci]=binofit(mocount,sum(mocount));
[omhat,omci]=binofit(omcount,sum(omcount));
hold on
bh=bar([mohat;omhat],'grouped');
errorbar([bh.XEndPoints],[bh.YEndPoints],...
    ebarsort([moci(:,1),omci(:,1)])-[bh.YEndPoints],...
    ebarsort([moci(:,2),omci(:,2)])-[bh.YEndPoints],'k.');
% title('mix-olf, olf-mix')
set(gca(),'YTick',0.3:0.1:0.7,'YTickLabel',30:10:70,'XTick',1:2,...
    'XTickLabel',{'Mixed>Olf.','Olf.>Mixed'},'YLim',[0.3,0.7])
yline(0.5,'k:')
ylabel('Consistent (%)')

nexttile()
mdcount=[nnz(mddf>0),nnz(mddf<=0)];
dmcount=[nnz(dmdf>0),nnz(dmdf<=0)];
pmd=binocdf(min(mdcount),numel(mddf),0.5)*2;
pdm=binocdf(min(dmcount),numel(dmdf),0.5)*2;
[mdhat,mdci]=binofit(mdcount,sum(mdcount));
[dmhat,dmci]=binofit(dmcount,sum(dmcount));
hold on
bh=bar([mdhat;dmhat],'grouped');
errorbar([bh.XEndPoints],[bh.YEndPoints],...
    ebarsort([mdci(:,1),dmci(:,1)])-[bh.YEndPoints],...
    ebarsort([mdci(:,2),dmci(:,2)])-[bh.YEndPoints],'k.');
set(gca(),'YTick',0.3:0.1:0.7,'YTickLabel',30:10:70,'XTick',1:2,...
    'XTickLabel',{'Mixed>Dur.','Dur.>Mixed'},'YLim',[0.3,0.7])
yline(0.5,'k:')
ylabel('Consistent (%)')

nexttile()
odcount=[nnz(oddf>0),nnz(oddf<=0)];
docount=[nnz(dodf>0),nnz(dodf<=0)];
pod=binocdf(min(odcount),numel(oddf),0.5)*2;
pdo=binocdf(min(docount),numel(dodf),0.5)*2;
[odhat,odci]=binofit(odcount,sum(odcount));
[dohat,doci]=binofit(docount,sum(docount));
hold on
bh=bar([odhat;dohat],'grouped');
errorbar([bh.XEndPoints],[bh.YEndPoints],...
    ebarsort([odci(:,1),doci(:,1)])-[bh.YEndPoints],...
    ebarsort([odci(:,2),doci(:,2)])-[bh.YEndPoints],'k.');
set(gca(),'YTick',0.3:0.1:0.7,'YTickLabel',30:10:70,'XTick',1:2,...
    'XTickLabel',{'Olf.>Dur.','Dur.>Olf.'},'YLim',[0.3,0.7])
yline(0.5,'k:')
ylabel('Consistent (%)')
sgtitle(sprintf('%.3f, ',pmo,pom,pmd,pdm,pod,pdo))

%%
end

function split_class_stats(sel_meta)

%s1d3-s1d3, s1d3-s1d6,s1-s1d3 s1-s1d6
rrate=@(x,y) [x,y,x./y];
[sig,pair]=bz.load_sig_sums_conn_file('pair',true);
sig=bz.join_fc_waveid(sig,sel_meta.wave_id);
pair=bz.join_fc_waveid(pair,sel_meta.wave_id);
s1d3_s1d3_rate=rrate(nnz(sig.waveid(:,1)==1 & sig.waveid(:,2)==1),nnz(pair.waveid(:,1)==1 & pair.waveid(:,2)==1));
s1d6_s1d6_rate=rrate(nnz(sig.waveid(:,1)==2 & sig.waveid(:,2)==2),nnz(pair.waveid(:,1)==2 & pair.waveid(:,2)==2));
s1sg_s1sg_rate=rrate(nnz(sig.waveid(:,1)==5 & sig.waveid(:,2)==5),nnz(pair.waveid(:,1)==5 & pair.waveid(:,2)==5));

s1d3_s1d6_rate=rrate(nnz(sig.waveid(:,1)==1 & sig.waveid(:,2)==2),nnz(pair.waveid(:,1)==1 & pair.waveid(:,2)==2));
s1d6_s1d3_rate=rrate(nnz(sig.waveid(:,1)==2 & sig.waveid(:,2)==1),nnz(pair.waveid(:,1)==2 & pair.waveid(:,2)==1));

s1d3_s1sg_rate=rrate(nnz(sig.waveid(:,1)==1 & sig.waveid(:,2)==5),nnz(pair.waveid(:,1)==1 & pair.waveid(:,2)==5));
s1sg_s1d3_rate=rrate(nnz(sig.waveid(:,1)==5 & sig.waveid(:,2)==1),nnz(pair.waveid(:,1)==5 & pair.waveid(:,2)==1));

s1d6_s1sg_rate=rrate(nnz(sig.waveid(:,1)==2 & sig.waveid(:,2)==5),nnz(pair.waveid(:,1)==2 & pair.waveid(:,2)==5));
s1sg_s1d6_rate=rrate(nnz(sig.waveid(:,1)==5 & sig.waveid(:,2)==2),nnz(pair.waveid(:,1)==5 & pair.waveid(:,2)==2));

% TODO total input, i.e. rate * lead freq

end

%% activity_stats, fcsp, fc_coding, etc
function activity_stats()

mix_sel=ismember(wrs_mux_meta.wave_id,1:4);
olf_sel=ismember(wrs_mux_meta.wave_id,5:6);
dur_sel=ismember(wrs_mux_meta.wave_id,7:8);

fl=dir(fullfile('fccoding','*.mat'));

[mix_olf_sum,mix_dur_sum,olf_mix_sum,dur_mix_sum]=deal([]);

for fi=1:numel(fl)
    fstr=load(fullfile(fl(fi).folder,fl(fi).name));

    fcsuids=uint16(cell2mat(fstr.onesess.fc(:,1)));
    sesssel=meta.sess==fstr.onesess.fidx;
    sesssuids=meta.allcid(sesssel);
    [~,sessidx]=ismember(fcsuids,sesssuids);

    sess_wave_id=wrs_mux_meta.wave_id(sesssel);
    fc_wave_id=sess_wave_id(sessidx);

    congrusel=pct.su_pairs.get_congru(fc_wave_id);

    mix_su=meta.allcid(sesssel & mix_sel);
    olf_su=meta.allcid(sesssel & olf_sel);
    dur_su=meta.allcid(sesssel & dur_sel);

    trials=fstr.onesess.trials;

    %% mix->olf vs mix->dur
    mix_olf_sel=ismember(fcsuids(:,1),mix_su) & ismember(fcsuids(:,2),olf_su) & congrusel;
    mix_dur_sel=ismember(fcsuids(:,1),mix_su) & ismember(fcsuids(:,2),dur_su) & congrusel;

    for fci=reshape(find(mix_olf_sel),1,[])
        curr_suid=fcsuids(fci,:);
        curr_wave=fc_wave_id(fci,:);
        % consistent trial inconsistent trial
        [consis,incons]=consist_incons_trial(trials,curr_wave);
        fc_data=cell2mat(fstr.onesess.fc(fci,[2 8 9])); %fcsp, lead_spks, follow_spks
        consis_data=fc_data(consis,:);
        incons_data=fc_data(incons,:);

        mix_olf_sum=[mix_olf_sum;...
            mean(consis_data),...
            mean(consis_data(consis_data(:,2)>0,1)./consis_data(consis_data(:,2)>0,2)),...
            mean(consis_data(consis_data(:,3)>0,1)./consis_data(consis_data(:,3)>0,3)),...
            mean(incons_data),...
            mean(incons_data(incons_data(:,2)>0,1)./incons_data(incons_data(:,2)>0,2)),...
            mean(incons_data(incons_data(:,3)>0,1)./incons_data(incons_data(:,3)>0,3))...
            ];
    end


    for fci=reshape(find(mix_dur_sel),1,[])
        curr_suid=fcsuids(fci,:);
        curr_wave=fc_wave_id(fci,:);
        % consistent trial inconsistent trial
        [consis,incons]=consist_incons_trial(trials,curr_wave);
        fc_data=cell2mat(fstr.onesess.fc(fci,[2 8 9]));
        consis_data=fc_data(consis,:);
        incons_data=fc_data(incons,:);

        mix_dur_sum=[mix_dur_sum;...
            mean(consis_data),...
            mean(consis_data(consis_data(:,2)>0,1)./consis_data(consis_data(:,2)>0,2)),...
            mean(consis_data(consis_data(:,3)>0,1)./consis_data(consis_data(:,3)>0,3)),...
            mean(incons_data),...
            mean(incons_data(incons_data(:,2)>0,1)./incons_data(incons_data(:,2)>0,2)),...
            mean(incons_data(incons_data(:,3)>0,1)./incons_data(incons_data(:,3)>0,3))...
            ];
    end


    %% olf-mix vs dur->mix
    olf_mix_sel=ismember(fcsuids(:,2),mix_su) & ismember(fcsuids(:,1),olf_su) & congrusel;
    dur_mix_sel=ismember(fcsuids(:,2),mix_su) & ismember(fcsuids(:,1),dur_su) & congrusel;

    for fci=reshape(find(olf_mix_sel),1,[])
        curr_suid=fcsuids(fci,:);
        curr_wave=fc_wave_id(fci,:);
        % consistent trial inconsistent trial
        [consis,incons]=consist_incons_trial(trials,curr_wave);
        fc_data=cell2mat(fstr.onesess.fc(fci,[2 8 9]));
        consis_data=fc_data(consis,:);
        incons_data=fc_data(incons,:);

        olf_mix_sum=[olf_mix_sum;...
            mean(consis_data),...
            mean(consis_data(consis_data(:,2)>0,1)./consis_data(consis_data(:,2)>0,2)),...
            mean(consis_data(consis_data(:,3)>0,1)./consis_data(consis_data(:,3)>0,3)),...
            mean(incons_data),...
            mean(incons_data(incons_data(:,2)>0,1)./incons_data(incons_data(:,2)>0,2)),...
            mean(incons_data(incons_data(:,3)>0,1)./incons_data(incons_data(:,3)>0,3))...
            ];
    end


    for fci=reshape(find(dur_mix_sel),1,[])
        curr_suid=fcsuids(fci,:);
        curr_wave=fc_wave_id(fci,:);
        % consistent trial inconsistent trial
        [consis,incons]=consist_incons_trial(trials,curr_wave);
        fc_data=cell2mat(fstr.onesess.fc(fci,[2 8 9]));
        consis_data=fc_data(consis,:);
        incons_data=fc_data(incons,:);

        dur_mix_sum=[dur_mix_sum;...
            mean(consis_data),...
            mean(consis_data(consis_data(:,2)>0,1)./consis_data(consis_data(:,2)>0,2)),...
            mean(consis_data(consis_data(:,3)>0,1)./consis_data(consis_data(:,3)>0,3)),...
            mean(incons_data),...
            mean(incons_data(incons_data(:,2)>0,1)./incons_data(incons_data(:,2)>0,2)),...
            mean(incons_data(incons_data(:,3)>0,1)./incons_data(incons_data(:,3)>0,3))...
            ];
    end


end

mix_olf_spk_mm=mean(mix_olf_sum(:,[1:3,6:8])./3);
mix_dur_spk_mm=mean(mix_dur_sum(:,[1:3,6:8])./3);
olf_mix_spk_mm=mean(olf_mix_sum(:,[1:3,6:8])./3);
dur_mix_spk_mm=mean(dur_mix_sum(:,[1:3,6:8])./3);

mix_olf_spk_sem=std(mix_olf_sum(:,[1:3,6:8])./3)./sqrt(size(mix_olf_sum,1));
mix_dur_spk_sem=std(mix_dur_sum(:,[1:3,6:8])./3)./sqrt(size(mix_dur_sum,1));
olf_mix_spk_sem=std(olf_mix_sum(:,[1:3,6:8])./3)./sqrt(size(olf_mix_sum,1));
dur_mix_spk_sem=std(dur_mix_sum(:,[1:3,6:8])./3)./sqrt(size(dur_mix_sum,1));


[mosp,mdsp,omsp,dmsp]=deal(nan(1,3));
for ii=1:3
    mosp(ii)=ranksum(mix_olf_sum(:,ii),mix_olf_sum(:,ii+5));
    mdsp(ii)=ranksum(mix_dur_sum(:,ii),mix_dur_sum(:,ii+5));
    omsp(ii)=ranksum(olf_mix_sum(:,ii),olf_mix_sum(:,ii+5));
    dmsp(ii)=ranksum(dur_mix_sum(:,ii),dur_mix_sum(:,ii+5));
end

mix_olf_rate_mm=mean(mix_olf_sum(:,[4,5,9,10]));
mix_dur_rate_mm=mean(mix_dur_sum(:,[4,5,9,10]));
olf_mix_rate_mm=nanmean(olf_mix_sum(:,[4,5,9,10])); % zero follow spk in inconsistent trials
dur_mix_rate_mm=mean(dur_mix_sum(:,[4,5,9,10]));


mix_olf_rate_sem=std(mix_olf_sum(:,[4,5,9,10]))./sqrt(size(mix_olf_sum,1));
mix_dur_rate_sem=std(mix_dur_sum(:,[4,5,9,10]))./sqrt(size(mix_dur_sum,1));
olf_mix_rate_sem=nanstd(olf_mix_sum(:,[4,5,9,10]))./sqrt(sum(isfinite(olf_mix_sum(:,[4,5,9,10])))); % zero follow spk in inconsistent trials
dur_mix_rate_sem=std(dur_mix_sum(:,[4,5,9,10]))./sqrt(size(dur_mix_sum,1));

[morp,mdrp,omrp,dmrp]=deal(nan(1,2));
for ii=1:2
    morp(ii)=ranksum(mix_olf_sum(:,ii+3),mix_olf_sum(:,ii+8));
    mdrp(ii)=ranksum(mix_dur_sum(:,ii+3),mix_dur_sum(:,ii+8));
    omrp(ii)=ranksum(olf_mix_sum(:,ii+3),olf_mix_sum(:,ii+8));
    dmrp(ii)=ranksum(dur_mix_sum(:,ii+3),dur_mix_sum(:,ii+8));
end

pv2str=@(x) subsref({'N.S.','*','**','***'},...
    struct(type='{}',subs={{[x>=0.05,x>=0.01 & x<0.05,x>=0.001 & x<0.01, x<0.001]}}));

%% spike count plot
figure('Color','w')
tiledlayout(2,2)
nexttile()
hold on
bh=bar(mix_olf_spk_mm([1:3;4:6].'));
errorbar([bh.XEndPoints],[bh.YEndPoints],mix_olf_spk_sem,'k.')
title('mix->olf')
set(gca,'XTick',1:3,'XTickLabel',{'#FC','#Lead','#Follow'});
legend(bh,{'both preferred','mix nonpref, simple preferred'},'Location','northoutside','Orientation','horizontal')
ymax=max(ylim());
for ii=1:3, text(ii,ymax,pv2str(mosp(ii)),'HorizontalAlignment','center','VerticalAlignment','middle');end
nexttile()
hold on
bh=bar(mix_dur_spk_mm([1:3;4:6].'));
errorbar([bh.XEndPoints],[bh.YEndPoints],mix_dur_spk_sem,'k.');
ymax=max(ylim());
for ii=1:3, text(ii,ymax,pv2str(mdsp(ii)),'HorizontalAlignment','center','VerticalAlignment','middle');end
title('mix->dur')
set(gca,'XTick',1:3,'XTickLabel',{'#FC','#Lead','#Follow'});
nexttile()
hold on
bh=bar(olf_mix_spk_mm([1:3;4:6].'));
errorbar([bh.XEndPoints],[bh.YEndPoints],olf_mix_spk_sem,'k.');
ymax=max(ylim());
for ii=1:3, text(ii,ymax,pv2str(omsp(ii)),'HorizontalAlignment','center','VerticalAlignment','middle');end
title('olf->mix')
set(gca,'XTick',1:3,'XTickLabel',{'#FC','#Lead','#Follow'});
nexttile()
hold on
bh=bar(dur_mix_spk_mm([1:3;4:6].'));
errorbar([bh.XEndPoints],[bh.YEndPoints],dur_mix_spk_sem,'k.');
ymax=max(ylim());
for ii=1:3, text(ii,ymax,pv2str(dmsp(ii)),'HorizontalAlignment','center','VerticalAlignment','middle');end
title('dur->mix')
set(gca,'XTick',1:3,'XTickLabel',{'#FC','#Lead','#Follow'});
sgtitle('Activity (SPK)')

%% fcsp/spike rate plot
figure('Color','w')
tiledlayout(2,2)
nexttile()
hold on
bh=bar(mix_olf_rate_mm([1,3;2,4]));
errorbar([bh.XEndPoints],[bh.YEndPoints],mix_olf_rate_sem,'k.');
ymax=max(ylim());
for ii=1:2, text(ii,ymax,pv2str(morp(ii)),'HorizontalAlignment','center','VerticalAlignment','middle');end
title('mix->olf')
set(gca,'XTick',1:2,'XTickLabel',{'Lead','Follow'});
legend(bh,{'both preferred','mix nonpref, simple preferred'},'Location','northoutside','Orientation','horizontal')
nexttile()
hold on
bh=bar(mix_dur_rate_mm([1,3;2,4]));
errorbar([bh.XEndPoints],[bh.YEndPoints],mix_dur_rate_sem,'k.');
ymax=max(ylim());
for ii=1:2, text(ii,ymax,pv2str(mdrp(ii)),'HorizontalAlignment','center','VerticalAlignment','middle');end
title('mix->dur')
set(gca,'XTick',1:2,'XTickLabel',{'Lead','Follow'});
nexttile()
hold on
bh=bar(olf_mix_rate_mm([1,3;2,4]));
errorbar([bh.XEndPoints],[bh.YEndPoints],olf_mix_rate_sem,'k.');
ymax=max(ylim());
for ii=1:2, text(ii,ymax,pv2str(omrp(ii)),'HorizontalAlignment','center','VerticalAlignment','middle');end
title('olf->mix')
set(gca,'XTick',1:2,'XTickLabel',{'Lead','Follow'});
nexttile()
hold on
bh=bar(dur_mix_rate_mm([1,3;2,4]));
errorbar([bh.XEndPoints],[bh.YEndPoints],dur_mix_rate_sem,'k.');
ymax=max(ylim());
for ii=1:2, text(ii,ymax,pv2str(dmrp(ii)),'HorizontalAlignment','center','VerticalAlignment','middle');end
title('dur->mix')
set(gca,'XTick',1:2,'XTickLabel',{'Lead','Follow'});
sgtitle('Efficiency (FCSP/SPK)')
end

function [consis,incons]=consist_incons_trial(trials,wave_id)
arguments
    trials (:,10) double
    wave_id (1,2) double
end
valid_sel=all(trials(:,9:10)~=0,2); % & ismember(trials(:,8),[3,6]) & ismember(trials(:,5),[4,8]);
if all(ismember(wave_id,[1 5]),'all')
    consis=trials(:,5)==4 & trials(:,8)==3 & valid_sel;
    incons=trials(:,5)==4 & trials(:,8)==6 & valid_sel;
elseif all(ismember(wave_id,[2 5]),'all')
    consis=trials(:,5)==4 & trials(:,8)==6 & valid_sel;
    incons=trials(:,5)==4 & trials(:,8)==3 & valid_sel;
elseif all(ismember(wave_id,[3 6]),'all')
    consis=trials(:,5)==8 & trials(:,8)==3 & valid_sel;
    incons=trials(:,5)==8 & trials(:,8)==6 & valid_sel;
elseif all(ismember(wave_id,[4 6]),'all')
    consis=trials(:,5)==8 & trials(:,8)==6 & valid_sel;
    incons=trials(:,5)==8 & trials(:,8)==3 & valid_sel;
elseif all(ismember(wave_id,[1 7]),'all')
    consis=trials(:,5)==4 & trials(:,8)==3 & valid_sel;
    incons=trials(:,5)==8 & trials(:,8)==3 & valid_sel;
elseif all(ismember(wave_id,[2 8]),'all')
    consis=trials(:,5)==4 & trials(:,8)==6 & valid_sel;
    incons=trials(:,5)==8 & trials(:,8)==6 & valid_sel;
elseif all(ismember(wave_id,[3 7]),'all')
    consis=trials(:,5)==8 & trials(:,8)==3 & valid_sel;
    incons=trials(:,5)==4 & trials(:,8)==3 & valid_sel;
elseif all(ismember(wave_id,[4 8]),'all')
    consis=trials(:,5)==8 & trials(:,8)==6 & valid_sel;
    incons=trials(:,5)==4 & trials(:,8)==6 & valid_sel;
end
end

function fc_rate_mix_simple(sel_meta,leadgrp,followgrp,opt)
arguments
    sel_meta
    leadgrp = []
    followgrp = []
    opt.condense_plot (1,1) logical = false
end
[sig,pair]=bz.load_sig_sums_conn_file('pair',true);
sig=bz.join_fc_waveid(sig,sel_meta.wave_id);
pair=bz.join_fc_waveid(pair,sel_meta.wave_id);

[inhibit_sig,inhibit_pair]=bz.load_sig_sums_conn_file('pair',true,'inhibit',true);
inhibit_sig=bz.join_fc_waveid(inhibit_sig,sel_meta.wave_id);
inhibit_pair=bz.join_fc_waveid(inhibit_pair,sel_meta.wave_id);

if exist("leadgrp","var") && exist("followgrp","var") && ~isempty(leadgrp) && ~isempty(followgrp)
    sig=reg_group(sig,leadgrp.lst,followgrp.lst);
    pair=reg_group(pair,leadgrp.lst,followgrp.lst);
    inhibit_sig=reg_group(inhibit_sig,leadgrp.lst,followgrp.lst);
    inhibit_pair=reg_group(inhibit_pair,leadgrp.lst,followgrp.lst);
end %if non-exist, skip selection


%% olfactory >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% s1d3,s1d6,s1
s1_out=triplet_stats(sig,pair,1,2,5);
inhibit_s1_out=triplet_stats(inhibit_sig,inhibit_pair,1,2,5);

s2_out=triplet_stats(sig,pair,3,4,6);
inhibit_s2_out=triplet_stats(inhibit_sig,inhibit_pair,3,4,6);

d3_out=triplet_stats(sig,pair,1,3,7);
inhibit_d3_out=triplet_stats(inhibit_sig,inhibit_pair,1,3,7);

d6_out=triplet_stats(sig,pair,2,4,8);
inhibit_d6_out=triplet_stats(inhibit_sig,inhibit_pair,2,4,8);

m2o_count=[s1_out.a2c_c2a(1,1:2);...
    s1_out.b2c_c2b(1,1:2);...
    s2_out.a2c_c2a(1,1:2);...
    s2_out.b2c_c2b(1,1:2)];
inhibit_m2o_count=[inhibit_s1_out.a2c_c2a(1,1:2);...
    inhibit_s1_out.b2c_c2b(1,1:2);...
    inhibit_s2_out.a2c_c2a(1,1:2);...
    inhibit_s2_out.b2c_c2b(1,1:2)];

o2m_count=[s1_out.a2c_c2a(2,1:2);...
    s1_out.b2c_c2b(2,1:2);...
    s2_out.a2c_c2a(2,1:2);...
    s2_out.b2c_c2b(2,1:2)];
inhibit_o2m_count=[inhibit_s1_out.a2c_c2a(2,1:2);...
    inhibit_s1_out.b2c_c2b(2,1:2);...
    inhibit_s2_out.a2c_c2a(2,1:2);...
    inhibit_s2_out.b2c_c2b(2,1:2)];

% m2m_count=[s1_out.a2b_b2a(:,1:2);...
%     s2_out.a2b_b2a(:,1:2)];
% inhibit_m2m_count=[inhibit_s1_out.a2b_b2a(:,1:2);...
%     inhibit_s2_out.a2b_b2a(:,1:2)];

if false % global cross mix sig, pair, currently not in use
    sig_mix_N=nnz(all(ismember(sig.waveid,1:4),2) & sig.waveid(:,1)~=sig.waveid(:,2));
    pair_mix_N=nnz(all(ismember(pair.waveid,1:4),2) & pair.waveid(:,1)~=pair.waveid(:,2));

    inhibit_sig_mix_N=nnz(all(ismember(inhibit_sig.waveid,1:4),2) & inhibit_sig.waveid(:,1)~=inhibit_sig.waveid(:,2));
    inhibit_pair_mix_N=nnz(all(ismember(inhibit_pair.waveid,1:4),2) & inhibit_pair.waveid(:,1)~=inhibit_pair.waveid(:,2));
end

count_sum=[sum(m2o_count);...
    sum(o2m_count);...
%     sum(m2m_count);...
%     sig_mix_N,pair_mix_N;...
    sum(inhibit_m2o_count);...
    sum(inhibit_o2m_count)];%;...
%     inhibit_sig_mix_N,inhibit_pair_mix_N];
%     sum(inhibit_m2m_count)];

p_m_o=[];
for bias=[0,2]
    [~,~,p]=crosstab([zeros(count_sum(1+bias,2),1);...
        ones(count_sum(2+bias,2),1)],...
        [1:count_sum(1+bias,2)>count_sum(1+bias,1),...
    1:count_sum(2+bias,2)>count_sum(2+bias,1)]);
    p_m_o=[p_m_o;p];
end

[mohat,moci]=binofit(count_sum(:,1),count_sum(:,2));

%% duration >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
m2d_count=[d3_out.a2c_c2a(1,1:2);...
    d3_out.b2c_c2b(1,1:2);...
    d6_out.a2c_c2a(1,1:2);...
    d6_out.b2c_c2b(1,1:2)];
inhibit_m2d_count=[inhibit_d3_out.a2c_c2a(1,1:2);...
    inhibit_d3_out.b2c_c2b(1,1:2);...
    inhibit_d6_out.a2c_c2a(1,1:2);...
    inhibit_d6_out.b2c_c2b(1,1:2)];

d2m_count=[d3_out.a2c_c2a(2,1:2);...
    d3_out.b2c_c2b(2,1:2);...
    d6_out.a2c_c2a(2,1:2);...
    d6_out.b2c_c2b(2,1:2)];
inhibit_d2m_count=[inhibit_d3_out.a2c_c2a(2,1:2);...
    inhibit_d3_out.b2c_c2b(2,1:2);...
    inhibit_d6_out.a2c_c2a(2,1:2);...
    inhibit_d6_out.b2c_c2b(2,1:2)];

% m2m_d_count=[d3_out.a2b_b2a(:,1:2);d6_out.a2b_b2a(:,1:2)];
% inhibit_m2m_d_count=[inhibit_d3_out.a2b_b2a(:,1:2);...
%     inhibit_d6_out.a2b_b2a(:,1:2)];

d_count_sum=[sum(m2d_count);...
    sum(d2m_count);...
    sum(inhibit_m2d_count);...
    sum(inhibit_d2m_count)];

[mdhat,mdci]=binofit(d_count_sum(:,1),d_count_sum(:,2));

p_m_d=[];
for bias=[0,2]
    [~,~,p]=crosstab([zeros(d_count_sum(1+bias,2),1);...
        ones(d_count_sum(2+bias,2),1)],...
        [1:d_count_sum(1+bias,2)>d_count_sum(1+bias,1),... % isequal(1:count_sum(1,2)>count_sum(1,1),(1:count_sum(1,2))>count_sum(1,1))
        1:d_count_sum(2+bias,2)>d_count_sum(2+bias,1)]);
    p_m_d=[p_m_d;p];
end
%% olf <-> dur >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

sig_o2d=nnz(ismember(sig.waveid(:,1),5:6) & ismember(sig.waveid(:,2),7:8));
sig_d2o=nnz(ismember(sig.waveid(:,1),7:8) & ismember(sig.waveid(:,2),5:6));
pair_o2d=nnz(ismember(pair.waveid(:,1),5:6) & ismember(pair.waveid(:,2),7:8)); % pair_d2o is the same by design

inhibit_sig_o2d=nnz(ismember(inhibit_sig.waveid(:,1),5:6) & ismember(inhibit_sig.waveid(:,2),7:8));
inhibit_sig_d2o=nnz(ismember(inhibit_sig.waveid(:,1),7:8) & ismember(inhibit_sig.waveid(:,2),5:6));
inhibit_pair_o2d=nnz(ismember(inhibit_pair.waveid(:,1),5:6) & ismember(inhibit_pair.waveid(:,2),7:8)); % pair_d2o is the same by design

[odhat,odci]=binofit([sig_o2d;sig_d2o;inhibit_sig_o2d;inhibit_sig_d2o],...
    [pair_o2d;pair_o2d;inhibit_pair_o2d;inhibit_pair_o2d]);

[~,~,p_o_d(1)]=crosstab([zeros(pair_o2d,1);...
    ones(pair_o2d,1)],...
    [1:pair_o2d>sig_o2d,... % isequal(1:count_sum(1,2)>count_sum(1,1),(1:count_sum(1,2))>count_sum(1,1))
    1:pair_o2d>sig_d2o]);
[~,~,p_o_d(2)]=crosstab([zeros(inhibit_pair_o2d,1);...
    ones(inhibit_pair_o2d,1)],...
    [1:inhibit_pair_o2d>inhibit_sig_o2d,... % isequal(1:count_sum(1,2)>count_sum(1,1),(1:count_sum(1,2))>count_sum(1,1))
    1:inhibit_pair_o2d>inhibit_sig_d2o]);


if opt.condense_plot
    figure('Color','w')
    hold on
    bh=bar(1:6,[mohat(1:2);mdhat(1:2);odhat(1:2)]);
    errorbar([bh.XEndPoints],[bh.YEndPoints],...
        [moci(1:2,1);mdci(1:2,1);odci(1:2,1)]-bh.YEndPoints.',...
        [moci(1:2,2);mdci(1:2,2);odci(1:2,2)]-bh.YEndPoints.',...
        'k.');
    set(gca(),'XTick',1:6,'XTickLabel',{'M2O','O2M','M2D','D2M','O2D','D2O'})
    if exist("leadgrp","var") && exist("followgrp","var") && ~isempty(leadgrp) && ~isempty(followgrp)
        set(gca(),'YLim',[0,0.025],'YTick',0:0.005:0.025,'YTickLabel',0:0.5:2.5);
    else
        set(gca(),'YLim',[0,0.015],'YTick',0:0.005:0.015,'YTickLabel',0:0.5:1.5);
    end
    ylabel('F.C. rate (%)')
    keyboard()
%     legend(bh,{'Mixed to olf','Olf. to mixed','Mixed to mixed'},'Location','northoutside','Orientation','horizontal')
    title(sprintf('%.3f, ',p_m_o(1),p_m_d(1),p_o_d(1)))

else
    figure('Color','w')
    tiledlayout(1,3)
    nexttile()
    hold on
    bh=bar([mohat(1:2).';mohat(3:4).']);
    errorbar([bh.XEndPoints],[bh.YEndPoints],moci([1 3 2 4],1),moci([1 3 2 4],2),'k.')
    set(gca(),'XTick',1:2,'XTickLabel',{'Excitatory','Inhibitory'})
    if exist("leadgrp","var") && exist("followgrp","var") && ~isempty(leadgrp) && ~isempty(followgrp)
        set(gca(),'YLim',[0,0.025],'YTick',0:0.005:0.025,'YTickLabel',0:0.5:2.5);
    else
        set(gca(),'YLim',[0,0.015],'YTick',0:0.005:0.015,'YTickLabel',0:0.5:1.5);
    end
    ylabel('F.C. rate (%)')
    legend(bh,{'Mixed to olf','Olf. to mixed','Mixed to mixed'},'Location','northoutside','Orientation','horizontal')
    title(sprintf('%.3f, ',p_m_o))


    % figure('Color','w')
    nexttile()
    hold on
    bh=bar([mdhat(1:2).';mdhat(3:4).']);
    errorbar([bh.XEndPoints],[bh.YEndPoints],mdci([1 3 2 4],1),mdci([1 3 2 4],2),'k.')
    set(gca(),'XTick',1:2,'XTickLabel',{'Excitatory','Inhibitory'})
    if exist("leadgrp","var") && exist("followgrp","var") && ~isempty(leadgrp) && ~isempty(followgrp)
        set(gca(),'YLim',[0,0.025],'YTick',0:0.005:0.025,'YTickLabel',0:0.5:2.5);
    else
        set(gca(),'YLim',[0,0.015],'YTick',0:0.005:0.015,'YTickLabel',0:0.5:1.5);
    end
    ylabel('F.C. rate (%)')
    legend(bh,{'Mixed to Dur.','Dur. to mixed','Mixed to mixed'},'Location','northoutside','Orientation','horizontal')
    title(sprintf('%.3f, ',p_m_d))


    if false % global cross mix sig, pair, currently not in use
        sig_mix_N=nnz(all(ismember(sig.waveid,1:4),2) & sig.waveid(:,1)~=sig.waveid(:,2));
        pair_mix_N=nnz(all(ismember(pair.waveid,1:4),2) & pair.waveid(:,1)~=pair.waveid(:,2));

        inhibit_sig_mix_N=nnz(all(ismember(inhibit_sig.waveid,1:4),2) & inhibit_sig.waveid(:,1)~=inhibit_sig.waveid(:,2));
        inhibit_pair_mix_N=nnz(all(ismember(inhibit_pair.waveid,1:4),2) & inhibit_pair.waveid(:,1)~=inhibit_pair.waveid(:,2));
    end

    % figure('Color','w')
    nexttile()
    hold on
    bh=bar(odhat([1,2;3,4]));
    errorbar([bh.XEndPoints],[bh.YEndPoints],...
        odci([1 3 2 4],1).'-[bh.YEndPoints],...
        odci([1 3 2 4],2).'-[bh.YEndPoints],...
        'k.');
    set(gca(),'XTick',1:2,'XTickLabel',{'Excitatory','Inhibitory'})
    if exist("leadgrp","var") && exist("followgrp","var") && ~isempty(leadgrp) && ~isempty(followgrp)
        set(gca(),'YLim',[0,0.025],'YTick',0:0.005:0.025,'YTickLabel',0:0.5:2.5);
    else
        set(gca(),'YLim',[0,0.015],'YTick',0:0.005:0.015,'YTickLabel',0:0.5:1.5);
    end
    ylabel('F.C. rate (%)')
    legend(bh,{'Olf. to Dur.','Dur. to Olf.'},'Location','northoutside','Orientation','horizontal')
    title(sprintf('%.3f, ',p_o_d))

    if exist("leadgrp","var") && exist("followgrp","var") && ~isempty(leadgrp) && ~isempty(followgrp)
        sgtitle([leadgrp.tag,' ',followgrp.tag])
    end
end
end

function out=triplet_stats(sig,pair,a,b,c)
rratio=@(x,y) [x,y,x./y];
out.a2b_b2a=[rratio(nnz(sig.waveid(:,1)==a & sig.waveid(:,2)==b),nnz(pair.waveid(:,1)==a & pair.waveid(:,2)==b));...
    rratio(nnz(sig.waveid(:,1)==b & sig.waveid(:,2)==a),nnz(pair.waveid(:,1)==b & pair.waveid(:,2)==a))];

out.a2c_c2a=[rratio(nnz(sig.waveid(:,1)==a & sig.waveid(:,2)==c),nnz(pair.waveid(:,1)==a & pair.waveid(:,2)==c));...
    rratio(nnz(sig.waveid(:,1)==c & sig.waveid(:,2)==a),nnz(pair.waveid(:,1)==c & pair.waveid(:,2)==a))];

out.b2c_c2b=[rratio(nnz(sig.waveid(:,1)==b & sig.waveid(:,2)==c),nnz(pair.waveid(:,1)==b & pair.waveid(:,2)==c));...
    rratio(nnz(sig.waveid(:,1)==c & sig.waveid(:,2)==b),nnz(pair.waveid(:,1)==c & pair.waveid(:,2)==b))];

out.a2a=rratio(nnz(sig.waveid(:,1)==a & sig.waveid(:,2)==a),nnz(pair.waveid(:,1)==a & pair.waveid(:,2)==a));
out.b2b=rratio(nnz(sig.waveid(:,1)==b & sig.waveid(:,2)==b),nnz(pair.waveid(:,1)==b & pair.waveid(:,2)==b));
out.c2c=rratio(nnz(sig.waveid(:,1)==c & sig.waveid(:,2)==c),nnz(pair.waveid(:,1)==c & pair.waveid(:,2)==c));
end


% function plot_triplet(excite,inhibit)
% figure()
% hold on
% annotation('arrow',[2,4]./5,[4.5,4.5]./5)
% annotation('arrow',[2,4]./5,[4,4]./5,'HeadStyle','rectangle')
% 
% annotation('arrow',[4,2]./5,[3.5,3.5]./5)
% annotation('arrow',[4,2]./5,[3,3]./5,'HeadStyle','rectangle')
% 
% 
% annotation('arrow',[1,1]./5,[2,3]./5)
% annotation('arrow',[1.5,1.5]./5,[2,3]./5,'HeadStyle','rectangle')
% 
% annotation('arrow',[4,4]./5,[2,3]./5)
% annotation('arrow',[4.5,4.5]./5,[2,3]./5,'HeadStyle','rectangle')
% end

% leadgrp={'AON','TT','DP','PIR'}
function out=reg_group(dinput,leadgrp,followgrp)
persistent idmap
if isempty(idmap)
    idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
end

leadgrp_ccfid=cell2mat(idmap.reg2ccfid.values(leadgrp));
followgrp_ccfid=cell2mat(idmap.reg2ccfid.values(followgrp));


% currgrp=zeros(size(dinput.suid));
lead_sel=ismember(squeeze(dinput.reg(:,5,1)),leadgrp_ccfid);
follow_sel=ismember(squeeze(dinput.reg(:,5,2)),followgrp_ccfid);

out=struct();
out.suid=dinput.suid(lead_sel & follow_sel,:);
out.sess=dinput.sess(lead_sel & follow_sel);
out.reg=dinput.reg(lead_sel & follow_sel,:,:);
out.waveid=dinput.waveid(lead_sel & follow_sel,:);

end

