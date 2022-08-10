warning('gather all plots from top')
keyboard()
clear all
close all
%% >>>>>>>>>>>>>>> statistic options >>>>>>>>>>>>>>>
global gather_config
gather_config=struct();
gather_config.fc_win=10;
gather_config.adjust_white_matter=true;
gather_config.corr_type='Pearson';
gather_config.fnsuffix='_10ms_adj_pearson';

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if strcmp(gather_config.corr_type,'Pearson')
    corr_ln_log='PearsonLinearLog';
    corr_log_log='PearsonLogLog';
else
    corr_ln_log='Spearman';
    corr_log_log='Spearman';
end

%% RANKSUM1 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

th=figure('Color','w','Position',[32,32,1600,900]);
blame=vcs.blame();

str=strjoin({'Stats:RANKSUM1',blame.datetime,blame.git_hash,blame.git_status},'=====\n');
annotation('textbox',[0.05,0.05,0.9,0.9],...
    'String',str,'FitBoxToText', ...
    'on','Interpreter','none');
% exportgraphics(th,'collections.pdf','ContentType','vector','Append',true)

% sens_meta=ephys.get_sens_meta('load_file',false,'permutation',true,'perm_repeat',1000,'save_file',true,'uneven_duration',true);
sens_meta=ephys.get_sens_meta('load_file',true);
stats_type='RANKSUM_per_bin';


% regulated by regional SU number > 100 in getGreRegs
[sens_map_cells,sens_reg_bar_fh]=ephys.Both_either_reg_bars( ...
    'stats_model','RANKSUM', ...
    'skip_plot',false, ...
    'waveid',sens_meta.wave_id, ...
    'range','grey', ...
    'data_type','sensory', ...
    'stats_type',stats_type);
sens_reg_bar_fh.Children.Children(5).YLim=[1e-3,0.5];
sens_GLM_fh=wave.connectivity_proportion_GLM( ...
    sens_map_cells, ...
    corr_log_log, ...
    'range','grey', ...
    'data_type','sensory', ...
    'stats_type','Ranksum',...
    'feat_tag',{'Context-independent','Context-dependent','All selective neurons'});
% sens_reg_bar_fh.reg_bar.Children(2).YLim=[0.001,0.5];



% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> sens dur FC >>>>>>>>>>>>>>>

fh_interwave_sens=bz.inter_wave(sens_meta);
fh_interwave_dur=bz.inter_wave(dur_meta,true);
fh_interwave_dur.Children.Children(3).YLim=[0,0.03];
fh_interwave_dur.Children.Children(2).YLim=[0,0.01];
fh_interwave_dur.Children.Children(2).YTickLabel=0:0.5:1;
fh_interwave_dur.Children.Children(2).YTick=0:0.005:0.01;

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Sensory wave part 1 >>>>>>>>>>>>>>>>>>>>>>>>
sust_trans_fh=ephys.sust_trans_bar();
wave_half_half_fh=wave.plot_wave_half_half(sens_meta);
stats_half_half_fh=wave.COM_half_half(sens_meta);
[wave_fh,stats]=wave.plot_wave_3_6(sens_meta);

%>>>>>>>>>>>>>>>>>>>>>>>>>>>> Duration distribution >>>>>>>>>>>>>>>>>>>>>>>>>>>>
fh=behav.per_sess_duration_coding();
% dur_meta=ephys.get_dur_meta('load_file',false,'merge_bin',false,'save_file',true,'perm_repeat',1000,'permutation',true);
dur_meta=ephys.get_dur_meta('load_file',true);
[dur_dec_fh,~]=wave.get_duration_decoding(dur_meta);

[dur_map_cells,dur_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','RANKSUM','skip_plot',false,'waveid',dur_meta.wave_id,'range','grey','data_type','duration','stats_type',stats_type);
dur_reg_bar_fh.Children.Children(5).YLim=[5e-5,0.5];
dur_reg_bar_fh.Children.Children(5).YTick=10.^[-4:0];
dur_GLM_fhgrey=wave.connectivity_proportion_GLM(dur_map_cells,corr_ln_log, ...
    'range','grey','data_type','duration','stats_type',stats_type,'feat_tag',{'Context-independent','Context-dependent','All selective neurons'});
dur_GLM_fhch=wave.connectivity_proportion_GLM(dur_map_cells,corr_ln_log, ...
    'range','CH','data_type','duration','stats_type',    stats_type,'feat_tag',{'Context-independent','Context-dependent','All selective neurons'});
dur_GLM_fhctx=wave.connectivity_proportion_GLM(dur_map_cells,corr_ln_log, ...
    'range','CTX','data_type','duration','stats_type',  stats_type,'feat_tag',{'Context-independent','Context-dependent','All selective neurons'});


% plot Neuron SC
if false
    %     dur_wo_sens=find(ismember(dur_meta.wave_id,5:6));
    dur_intact=find(anova2meta.dur & ~anova2meta.sens & anova2meta.interact & dur_meta.wave_id>0 &dur_meta.wave_id<5);
    meta=ephys.util.load_meta('skip_stats',true);
    for ii=[9071,21530,3387,2617]
        %         fh=ephys.sens_dur_SC(dur_intact(ii),meta,sens_meta,dur_meta,'skip_raster',false);
        scfh=ephys.sens_dur_SC(ii,meta,sens_meta,dur_meta,'skip_raster',false);%
        %             11241,20621,9071
    end
end
% dur_reg_bar_fh.reg_bar.Children(2).YLim=[0.0001,0.5];
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%>>>>>>>>>>>>>>>>>>> Single mix modality >>>>>>>>>>>>>>>>>>>>>>>>>
dur_sens_corr_fh=hier.sens_dur_corr(dur_map_cells{3},sens_map_cells{3});

single_mix_meta.single1=sens_meta.wave_id>0;
single_mix_meta.single2=dur_meta.wave_id>0;
[singlemix_map_cells,singlemix_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','SINGLE_MIX','skip_plot',false,'meta',single_mix_meta,'range','grey','data_type','single_mix','stats_type',stats_type);
% {indep, dep, I or D}
singlemix_GLM_fh=wave.connectivity_proportion_GLM(singlemix_map_cells,corr_log_log, ...
    'range','grey','data_type','single_mix','stats_type',stats_type,'feat_tag',{'Single modality','Mixed-modality','All selective neurons'});
% TODO mix density vs. sens density | dur density


perm_meta.sens_only=sens_meta.wave_id>0;
perm_meta.dur_only=dur_meta.wave_id>0;
perm_meta.mixed=sens_meta.wave_id>0 & dur_meta.wave_id>0;
inter_wave_fh=bz.inter_wave_ext_bars(perm_meta);  % dur, olf vs Isocortex, Straitum and Midbrain


%>>>>>>>>>>>>>>>>>>>>>> Sensory wave Part2 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% sense, 6s
sens_tcom_fh=struct();
tcom_maps=cell(1,2);
% sens_meta_even=ephys.sens_uneven2even(sens_meta);
for currdelay=[6,3]
    sens_com.(['d',num2str(currdelay)])=wave.get_com_map(sens_meta, ...
        'wave',['anyContext',num2str(currdelay/3)], ... %indep+dep
        'delay',currdelay);
    [fcom.(['d',num2str(currdelay)]).collection,...
        fcom.(['d',num2str(currdelay)]).com_meta]=wave.per_region_COM(...
        sens_com.(['d',num2str(currdelay)]),'stats_method','mean');

    ureg=intersect(ephys.getGreyRegs('range','grey'),...
        fcom.(['d',num2str(currdelay)]).collection(:,2));
    ffrac.(['d',num2str(currdelay)]).collection=...
        [num2cell(cellfun(@(x) x(1),sens_map_cells{1}.values(ureg))),...
        ureg,...
        num2cell(ones(numel(ureg),1)*5),...
        num2cell(cellfun(@(x) x(3),sens_map_cells{1}.values(ureg)))];
    
    [sens_tcom_fh.(['d',num2str(currdelay)]),t]=wave.per_region_COM_frac(fcom.(['d',num2str(currdelay)]),ffrac.(['d',num2str(currdelay)]),'hier_reg','AON','corr',gather_config.corr_type);
    t.Title.String=['Sense wave, delay=',num2str(currdelay)];
    [~,tcidx]=ismember(ureg,fcom.(['d',num2str(currdelay)]).collection(:,2));
    tcom_maps{currdelay/3}=containers.Map(ureg,num2cell(cellfun(@(x) x/4, fcom.(['d',num2str(currdelay)]).collection(tcidx,1))));
end

[comdiff_stats,com_pair]=wave.fc_com_pvsst(sens_com.d3,sens_com.d6,sens_meta,'hiermap','AON','tbl_title','Olfactory-AON-indep','mem_type','indep');
[comdiff_stats,com_pair]=wave.fc_com_pvsst(sens_com.d3,sens_com.d6,sens_meta,'hiermap','AON','tbl_title','Olfactory-AON-mixed','mem_type','mixed');
[comdiff_stats,com_pair]=wave.fc_com_pvsst(sens_com.d3,sens_com.d6,sens_meta,'hiermap','AON','tbl_title','Olfactory-AON-congru');
[comdiff_stats,com_pair]=wave.fc_com_pvsst(sens_com.d3,sens_com.d6,sens_meta,'hiermap','AON','tbl_title','Olfactory-AON-incong','mem_type','incong');


sensory_TCOM_GLM_fh=wave.connectivity_proportion_GLM(tcom_maps,corr_ln_log, ...
    'range','grey','data_type','3s_6s_sensory_TCOM','stats_type',stats_type,'feat_tag',{'Sens TCOM 3s','Sens TCOM 6s'});

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Dur TCOM>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
for currcue=[4,8]
    dur_com.(['s',num2str(currcue)])=wave.get_com_map(dur_meta, ...
        'wave',['anyContext',num2str(currcue/4)], ... %indep+dep
        'sens_cue',currcue,'curve',true);
end

%>>>>>>>>>>>> Create merged-odor dataset>>>>>
padOdor=@(y,z) num2cell(cellfun(@(x) x+z*100000,y));
usess=unique([fieldnames(dur_com.s4);fieldnames(dur_com.s8)]);
feat_fields=fieldnames(dur_com.s4.s1);
for sess=reshape(usess,1,[])
    for fn=reshape(feat_fields,1,[])
        fnkeys=[padOdor(dur_com.s4.(char(sess)).(fn{1}).keys,4),padOdor(dur_com.s8.(char(sess)).(fn{1}).keys,8)];
        if ~isempty(fnkeys)
            dur_com.summed.(char(sess)).(fn{1})=...
                containers.Map(fnkeys, ...
                [dur_com.s4.(char(sess)).(fn{1}).values,dur_com.s8.(char(sess)).(fn{1}).values]);
        end
    end
    end

dur_wave_fh=wave.plot_dur_wave(dur_com);
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
[fcom.dur.collection,fcom.dur.com_meta]=wave.per_region_COM(...
    dur_com.summed,'stats_method','mean','min_su',20);
for range=["grey","CH","CTX"]
    ureg=intersect(ephys.getGreyRegs('range',range),...
        fcom.dur.collection(:,2));
    ffrac.dur.collection=...
        [num2cell(cellfun(@(x) x(1),dur_map_cells{3}.values(ureg))),...
        ureg,...
        num2cell(ones(numel(ureg),1)*5),...
        num2cell(cellfun(@(x) x(3),dur_map_cells{3}.values(ureg)))];

    [dur_tcom_fh.(range),t]=wave.per_region_COM_frac(fcom.dur,ffrac.dur,'hier_reg','COA','corr',gather_config.corr_type,'density_scale','linear');
    t.Title.String=sprintf('Duration wave-%s',range);

    [~,tcidx]=ismember(ureg,fcom.dur.collection(:,2));
    dur_tcom_maps={containers.Map(ureg,num2cell(cellfun(@(x) x/4, fcom.dur.collection(tcidx,1))))};
    dur_TCOM_GLM_grey_fh=wave.connectivity_proportion_GLM(dur_tcom_maps,corr_ln_log, ...
        'range',range,'data_type','s4_s8_duration_TCOM','stats_type',stats_type,'feat_tag',"Duration TCOM-"+range);
end

[comdiff_stats,com_pair]=wave.fc_com_pvsst(dur_com.s4,dur_com.s8,dur_meta,'hiermap','ATN','tbl_title','Duration-ATN');
%========================================================================

sens_dur_TCOM_corr_fh=wave.sens_dur_TCOM_corr(fcom);

tcom_corr_bar_fh=wave.sens_dur_wave_bar(sens_map_cells,dur_map_cells,fcom);
% >>>>>>>>>>>>>>>>>>>>>>> FC DECODING >>>>>>>>>>>>>>>>>>>>>>>
bz.fccoding.plot_coding(sens_meta)
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% >>>>>>>>>>>>>>>>>>>>> Movie >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if false
    [~,out]=wave.reg_wave_movie(fcom,'spatial_coord',false,'write_video',true,'color_type','wavefront','export_data',true,'waves','both');
    [~,~]=wave.reg_wave_movie(fcom,'spatial_coord',false,'write_video',true,'color_type','wavefront','export_data',false,'waves','sens');
    [~,~]=wave.reg_wave_movie(fcom,'spatial_coord',false,'write_video',true,'color_type','wavefront','export_data',false,'waves','dur');
end
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




if exist(sprintf('collections%s.pdf',gather_config.fnsuffix),'file')
    delete(sprintf('collections%s.pdf',gather_config.fnsuffix))
end

fhandles=get(groot(),'Children');
for hc=reshape(fhandles,1,[])
    %exportgraphics(gcf(),sprintf('collections%s.pdf',gather_config.fnsuffix),'ContentType','vector','Append',true)
    exportgraphics(hc,sprintf('collections%s.pdf',gather_config.fnsuffix),'ContentType','vector','Append',true)
end
savefig(fhandles,sprintf('Ranksum1%s.fig',gather_config.fnsuffix));

return
