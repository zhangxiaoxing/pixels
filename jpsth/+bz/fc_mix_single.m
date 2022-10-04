global_init;
meta=ephys.util.load_meta('skip_stats',true);
wrs_mux_meta=ephys.get_wrs_mux_meta();
com_map=wave.get_pct_com_map(wrs_mux_meta,'curve',true);

function conn_stat(wrs_mux_meta,com_map,tcom_maps)

fc_stats=fc.fc_com_reg_wave(wrs_mux_meta,com_map,tcom_maps);
congrusel=pct.su_pairs.get_congru(fc_stats(:,10:11));

mix_olf_sel=ismember(fc_stats(:,10),1:4) & ismember(fc_stats(:,11),5:6) & congrusel;
olf_mix_sel=ismember(fc_stats(:,11),1:4) & ismember(fc_stats(:,10),5:6) & congrusel;
mix_dur_sel=ismember(fc_stats(:,10),1:4) & ismember(fc_stats(:,11),7:8) & congrusel;
dur_mix_sel=ismember(fc_stats(:,11),1:4) & ismember(fc_stats(:,10),7:8) & congrusel;

mix_dur_tcom=fc_stats(mix_dur_sel,4:5);
mddf=diff(mix_dur_tcom,1,2);

mix_olf_tcom=fc_stats(mix_olf_sel,4:5);
modf=diff(mix_olf_tcom,1,2);

dur_mix_tcom=fc_stats(dur_mix_sel,4:5);
dmdf=diff(dur_mix_tcom,1,2);

olf_mix_tcom=fc_stats(olf_mix_sel,4:5);
omdf=diff(olf_mix_tcom,1,2);

mon=nnz(mix_olf_sel);
omn=nnz(olf_mix_sel);
pmo=binocdf(min([mon;omn]),sum([omn;mon]),0.5)*2;

mdn=nnz(mix_dur_sel);
dmn=nnz(dur_mix_sel);
pmd=binocdf(min([mdn;dmn]),sum([dmn;mdn]),0.5)*2;


figure()
tiledlayout(1,2)
nexttile()
histogram(mddf)
nexttile()
histogram(modf)

%TCOM pdf cdf curve
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
s1d6_s1sg_rate=rrate(nnz(sig.waveid(:,1)==5 & sig.waveid(:,2)==2),nnz(pair.waveid(:,1)==5 & pair.waveid(:,2)==2));

% TODO total input, i.e. rate * lead freq

end



















function activity_stats

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
        fc_data=cell2mat(fstr.onesess.fc(fci,[2 8 9]));
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
legend(bh,{'both preferred','mix nonpref, single preferred'},'Location','northoutside','Orientation','horizontal')
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
legend(bh,{'both preferred','mix nonpref, single preferred'},'Location','northoutside','Orientation','horizontal')
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

