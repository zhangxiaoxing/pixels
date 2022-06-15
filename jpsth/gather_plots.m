warning('gather all plots from top')
keyboard()
clear all
close all

%% RANKSUM1
fnsuffix='_10ms_adj_pearson';
corr_type='Pearson';

if strcmp(corr_type,'Pearson')
    corr_ln_log='PearsonLinearLog';
    corr_log_log='PearsonLogLog';
else
    corr_ln_log='Spearman';
    corr_log_log='Spearman';
end

th=figure('Color','w','Position',[32,32,1600,900]);
blame=vcs.blame();
str=strjoin({'Stats:RANKSUM1',blame.datetime,blame.git_hash,blame.git_status},'=====\n');
annotation('textbox',[0.05,0.05,0.9,0.9],...
    'String',str,'FitBoxToText', ...
    'on','Interpreter','none');
% exportgraphics(th,'collections.pdf','ContentType','vector','Append',true)

% sens_meta=ephys.get_sens_meta();
% sens_waveid=ephys.get_wave_id(meta);
stats_type='RANKSUM_per_bin';
load perm_sens.mat

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

% dur_meta=ephys.get_dur_meta();
load perm_dur.mat
[dur_map_cells,dur_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','RANKSUM','skip_plot',false,'waveid',dur_meta.wave_id,'range','grey','data_type','duration','stats_type',stats_type);
dur_reg_bar_fh.Children.Children(5).YLim=[5e-5,0.5];
dur_reg_bar_fh.Children.Children(5).YTick=10.^[-4:0];
dur_GLM_fhgrey=wave.connectivity_proportion_GLM(dur_map_cells,corr_ln_log, ...
    'range','grey','data_type','duration','stats_type',stats_type,'feat_tag',{'Context-independent','Context-dependent','All selective neurons'});
dur_GLM_fhch=wave.connectivity_proportion_GLM(dur_map_cells,corr_ln_log, ...
    'range','CH','data_type','duration','stats_type',    stats_type,'feat_tag',{'Context-independent','Context-dependent','All selective neurons'});
dur_GLM_fhctx=wave.connectivity_proportion_GLM(dur_map_cells,corr_ln_log, ...
    'range','CTX','data_type','duration','stats_type',  stats_type,'feat_tag',{'Context-independent','Context-dependent','All selective neurons'});


% plot duration SC
if false
%     dur_wo_sens=find(ismember(dur_meta.wave_id,5:6));
    dur_intact=find(anova2meta.dur & ~anova2meta.sens & anova2meta.interact & dur_meta.wave_id>0 &dur_meta.wave_id<5);
    meta=ephys.util.load_meta('skip_stats',true);
    for ii=[11241]
%         fh=ephys.sens_dur_SC(dur_intact(ii),meta,sens_meta,dur_meta,'skip_raster',false);
            fh=ephys.sens_dur_SC(ii,meta,sens_meta,dur_meta,'skip_raster',false);%
%             11241,20621,9071
        waitfor(fh);
    end
end


% dur_reg_bar_fh.reg_bar.Children(2).YLim=[0.0001,0.5];

dur_sens_corr_fh=hier.sens_dur_corr(dur_map_cells{3},sens_map_cells{3});

single_mix_meta.single1=sens_meta.wave_id>0;
single_mix_meta.single2=dur_meta.wave_id>0;
[singlemix_map_cells,singlemix_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','SINGLE_MIX','skip_plot',false,'meta',single_mix_meta,'range','grey','data_type','single_mix','stats_type',stats_type);
% {indep, dep, I or D}
singlemix_GLM_fh=wave.connectivity_proportion_GLM(singlemix_map_cells,corr_log_log, ...
    'range','grey','data_type','single_mix','stats_type',stats_type,'feat_tag',{'Single modality','Mixed-modality','All selective neurons'});
% TODO mix density vs. sens density | dur density

fh_interwave=bz.inter_wave(sens_meta);

perm_meta.sens_only=sens_meta.wave_id>0;
perm_meta.dur_only=dur_meta.wave_id>0;
perm_meta.mixed=sens_meta.wave_id>0 & dur_meta.wave_id>0;
inter_wave_fh=bz.inter_wave_ext_bars(perm_meta);



% TCOM
% sense, 6s
sens_tcom_fh=struct();
tcom_maps=cell(1,2);
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
    
    [sens_tcom_fh.(['d',num2str(currdelay)]),t]=wave.per_region_COM_frac(fcom.(['d',num2str(currdelay)]),ffrac.(['d',num2str(currdelay)]),'hier_reg','AON','corr',corr_type);
    t.Title.String=['Sense wave, delay=',num2str(currdelay)];
    [~,tcidx]=ismember(ureg,fcom.(['d',num2str(currdelay)]).collection(:,2));
    tcom_maps{currdelay/3}=containers.Map(ureg,num2cell(cellfun(@(x) x/4, fcom.(['d',num2str(currdelay)]).collection(tcidx,1))));
end
[comdiff_stats,com_pair]=wave.fc_com_pvsst(sens_com.d3,sens_com.d6,sens_meta,'hiermap','AON','tbl_title','Olfactory-AON-congru');
[comdiff_stats,com_pair]=wave.fc_com_pvsst(sens_com.d3,sens_com.d6,sens_meta,'hiermap','AON','tbl_title','Olfactory-AON-incong','mem_type','incong');


sensory_TCOM_GLM_fh=wave.connectivity_proportion_GLM(tcom_maps,corr_ln_log, ...
    'range','grey','data_type','3s_6s_sensory_TCOM','stats_type',stats_type,'feat_tag',{'Sens TCOM 3s','Sens TCOM 6s'});

%=============================Dur TCOM===========================
for currcue=[4,8]
    dur_com.(['s',num2str(currcue)])=wave.get_com_map(dur_meta, ...
        'wave',['anyContext',num2str(currcue/4)], ... %indep+dep
        'sens_cue',currcue);
end

padOdor=@(y,z) num2cell(cellfun(@(x) x+z*100000,y));
usess=unique([fieldnames(dur_com.s4);fieldnames(dur_com.s8)]);
for sess=reshape(usess,1,[])
    c1keys=[padOdor(dur_com.s4.(char(sess)).c1.keys,4),padOdor(dur_com.s8.(char(sess)).c1.keys,8)];
    if ~isempty(c1keys)
        dur_com.summed.(char(sess)).c1=...
            containers.Map(c1keys, ...
            [dur_com.s4.(char(sess)).c1.values,dur_com.s8.(char(sess)).c1.values]);
    end
    c2keys=[padOdor(dur_com.s4.(char(sess)).c2.keys,4),padOdor(dur_com.s8.(char(sess)).c2.keys,8)];
    if ~isempty(c2keys)
        dur_com.summed.(char(sess)).c2=...
            containers.Map(c2keys, ...
            [dur_com.s4.(char(sess)).c2.values,dur_com.s8.(char(sess)).c2.values]);
    end
end

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

    [dur_tcom_fh.(range),t]=wave.per_region_COM_frac(fcom.dur,ffrac.dur,'hier_reg','COA','corr',corr_type,'density_scale','linear');
    t.Title.String=sprintf('Duration wave-%s',range);
end
 [~,tcidx]=ismember(ureg,fcom.dur.collection(:,2));
 dur_tcom_maps={containers.Map(ureg,num2cell(cellfun(@(x) x/4, fcom.dur.collection(tcidx,1))))};

dur_TCOM_GLM_grey_fh=wave.connectivity_proportion_GLM(dur_tcom_maps,corr_ln_log, ...
    'range','grey','data_type','s4_s8_duration_TCOM','stats_type',stats_type,'feat_tag',{'Duration TCOM-grey'});
dur_TCOM_GLM_ch_fh=wave.connectivity_proportion_GLM(dur_tcom_maps,corr_ln_log, ...
    'range','CH','data_type','s4_s8_duration_TCOM','stats_type',stats_type,'feat_tag',{'Duration TCOM-CH'});
dur_TCOM_GLM_ctx_fh=wave.connectivity_proportion_GLM(dur_tcom_maps,corr_ln_log, ...
    'range','CTX','data_type','s4_s8_duration_TCOM','stats_type',stats_type,'feat_tag',{'Duration TCOM-CTX'});

[comdiff_stats,com_pair]=wave.fc_com_pvsst(dur_com.s4,dur_com.s8,dur_meta,'hiermap','RSP','tbl_title','Duration-RSP');
%========================================================================

sens_dur_TCOM_corr_fh=wave.sens_dur_TCOM_corr(fcom);

tcom_corr_bar_fh=wave.sens_dur_wave_bar(sens_map_cells,dur_map_cells,fcom);
% >>>>>>>>>>>>>>>>>>>>>>> FC DECODING >>>>>>>>>>>>>>>>>>>>>>>
bz.fccoding.plot_coding(sens_meta)
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if exist(sprintf('collections%s.pdf',fnsuffix),'file')
    delete(sprintf('collections%s.pdf',fnsuffix))
end

fhandles=get(groot(),'Children');
for hc=reshape(fhandles,1,[])
    %exportgraphics(gcf(),sprintf('collections%s.pdf',fnsuffix),'ContentType','vector','Append',true)
    exportgraphics(hc,sprintf('collections%s.pdf',fnsuffix),'ContentType','vector','Append',true)
end
savefig(fhandles,sprintf('Ranksum1%s.fig',fnsuffix));

return

%% >>>>>>>>>>>>>>>  ANOVA2   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

anova2meta=ephys.selectivity_anova('merge_time_bin',true,'anova_model','full');

% %% ANOVA2ALT
% anova2_alt_meta=ephys.selectivity_anova('largest_varied_bin',true);

anovameta=anova2meta;
stats_type='ANOVA2';
[sens_map_cells,sens_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','ANOVA2','skip_plot',false,'meta',anovameta,'single_field','sens','range','grey','data_type','sensory','stats_type',stats_type);
sens_GLM_fh=wave.connectivity_proportion_GLM(sens_map_cells,'range','grey','data_type','sensory','stats_type',stats_type);

[dur_map_cells,dur_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','ANOVA2','skip_plot',false,'meta',anovameta,'single_field','dur','range','grey','data_type','duration','stats_type',stats_type);
dur_GLM_fh_grey=wave.connectivity_proportion_GLM(dur_map_cells,'range','grey','data_type','duration','stats_type',stats_type);
dur_GLM_fh_CH=wave.connectivity_proportion_GLM(dur_map_cells,'range','CH','data_type','duration','stats_type',stats_type);
dur_GLM_fh_CTX=wave.connectivity_proportion_GLM(dur_map_cells,'range','CTX','data_type','duration','stats_type',stats_type);

dur_sens_corr_fh=hier.sens_dur_corr(dur_map_cells{3},sens_map_cells{3});

single_mix_meta.single1=anovameta.sens;
single_mix_meta.single2=anovameta.dur;
[singlemix_map_cells,singlemix_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','SINGLE_MIX','skip_plot',false,'meta',single_mix_meta,'data_type','single_mix','stats_type',stats_type);
singlemix_GLM_fh_grey=wave.connectivity_proportion_GLM(singlemix_map_cells,'data_type','single_mix','stats_type',stats_type);
singlemix_GLM_fh_CH=wave.connectivity_proportion_GLM(singlemix_map_cells,'range','CH','data_type','single_mix','stats_type',stats_type);
singlemix_GLM_fh_CTX=wave.connectivity_proportion_GLM(singlemix_map_cells,'range','CTX','data_type','single_mix','stats_type',stats_type);

anova2meta.sens_only=anova2meta.sens;
anova2meta.dur_only=anova2meta.dur;
anova2meta.mixed=anova2meta.sens & anova2meta.dur;
inter_wave_fh=bz.inter_wave_ext_bars(anova2meta);

th=figure('Color','w','Position',[32,32,1600,900]);
blame=vcs.blame();
str=strjoin({'Stats:ANOVA2',blame.datetime,blame.git_hash,blame.git_status},'=====\n');
annotation('textbox',[0.05,0.05,0.9,0.9],'String',str,'FitBoxToText','on','Interpreter','none');
exportgraphics(th,'collections.pdf','ContentType','vector','Append',true)



fhandles={sens_reg_bar_fh,sens_GLM_fh,dur_reg_bar_fh,dur_GLM_fh_grey,dur_GLM_fh_CH,dur_GLM_fh_CTX,dur_sens_corr_fh,singlemix_reg_bar_fh,singlemix_GLM_fh_grey,singlemix_GLM_fh_CH,singlemix_GLM_fh_CTX,inter_wave_fh};

for hc=fhandles
    hh=hc{1};
    if isa(hh,'matlab.ui.Figure')
        disp('hh')
        exportgraphics(hh,'collections.pdf','ContentType','vector','Append',true)
    elseif isa(hh,"struct")
        fn=fieldnames(hh);
        disp(fn)
        for fh=reshape(fn,1,[])
            exportgraphics(hh.(fh{1}),'collections.pdf','ContentType','vector','Append',true)
        end
    end

end


% for fn=fieldnames(sens_reg_bar_fh)

%feat1: sens_only, feat2: mixed
%%F2
%TODO Heatmap session, corr stats

%TODO Heatmap population



%% ANOVA per bin + FDR
if false
    anova3meta=ephys.selectivity_anova('per_bin',true);
    stats_type='ANOVA_FDR';

    [sens_map_cells,sens_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','ANOVA2','skip_plot',false,'meta',anova3meta,'single_field','sens','range','grey','data_type','sensory','stats_type',stats_type);
    sens_GLM_fh=wave.connectivity_proportion_GLM(sens_map_cells,'range','grey','data_type','sensory','stats_type',stats_type);

    [dur_map_cells,dur_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','ANOVA2','skip_plot',false,'meta',anova3meta,'single_field','dur','range','grey','data_type','duration','stats_type',stats_type);
    dur_GLM_fh_grey=wave.connectivity_proportion_GLM(dur_map_cells,'range','grey','data_type','duration','stats_type',stats_type);
    dur_GLM_fh_CH=wave.connectivity_proportion_GLM(dur_map_cells,'range','CH','data_type','duration','stats_type',  stats_type);
    dur_GLM_fh_CTX=wave.connectivity_proportion_GLM(dur_map_cells,'range','CTX','data_type','duration','stats_type', stats_type);

    dur_sens_corr_fh=hier.sens_dur_corr(dur_map_cells{1},sens_map_cells{1});

    single_mix_meta.single1=anova3meta.sens;
    single_mix_meta.single2=anova3meta.dur;
    [singlemix_map_cells,singlemix_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','SINGLE_MIX','skip_plot',false,'meta',single_mix_meta);
    singlemix_GLM_fh=wave.connectivity_proportion_GLM(singlemix_map_cells);


    inter_wave_fh=bz.inter_wave_ext_bars(anovameta);

    th=figure('Color','w','Position',[32,32,1600,900]);
    blame=vcs.blame();
    str=strjoin({['St' ...
        'ats:ANOVA2'],blame.datetime,blame.git_hash,blame.git_status},'=====\n');
    annotation('textbox',[0.05,0.05,0.9,0.9],'String',str,'FitBoxToText','on','Interpreter','none');
    exportgraphics(th,'collections.pdf','ContentType','vector','Append',true)


    
    fhandles={sens_reg_bar_fh,sens_GLM_fh,dur_reg_bar_fh,dur_GLM_fh_grey,dur_GLM_fh_CH,dur_GLM_fh_CTX,dur_sens_corr_fh,singlemix_reg_bar_fh,singlemix_GLM_fh_grey,singlemix_GLM_fh_CH,singlemix_GLM_fh_CTX,inter_wave_fh};

    for hc=fhandles
        hh=hc{1};
        if isa(hh,'matlab.ui.Figure')
            disp('hh')
            exportgraphics(hh,'collections.pdf','ContentType','vector','Append',true)
        elseif isa(hh,"struct")
            fn=fieldnames(hh);
            disp(fn)
            for fh=reshape(fn,1,[])
                exportgraphics(hh.(fh{1}),'collections.pdf','ContentType','vector','Append',true)
            end
        end

    end




    %feat1: sens_only, feat2: mixed
    %%F2
    %TODO Heatmap session, corr stats

    %TODO Heatmap population

end

%% ANOVA3
if false
    anova3meta=ephys.selectivity_anova('anova_model','linear');

    [sens_map_cells,sens_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','ANOVA2','skip_plot',false,'meta',anova3meta,'single_field','sens');
    sens_GLM_fh=wave.connectivity_proportion_GLM(sens_map_cells);

    [dur_map_cells,dur_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','ANOVA2','skip_plot',false,'meta',anova3meta,'single_field','time_bin');
    dur_GLM_fh=wave.connectivity_proportion_GLM(dur_map_cells);

    dur_sens_corr_fh=hier.sens_dur_corr(dur_map_cells{1},sens_map_cells{1});

    single_mix_meta.single1=anova3meta.sens;
    single_mix_meta.single2=anova3meta.dur;
    [singlemix_map_cells,singlemix_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','SINGLE_MIX','skip_plot',false,'meta',single_mix_meta);
    singlemix_GLM_fh=wave.connectivity_proportion_GLM(singlemix_map_cells);

    inter_wave_fh=bz.inter_wave_ext_bars(anovameta);


    th=figure('Color','w','Position',[32,32,1600,900]);
    blame=vcs.blame();
    str=strjoin({'Stats:ANOVA2',blame.datetime,blame.git_hash,blame.git_status},'=====\n');
    annotation('textbox',[0.05,0.05,0.9,0.9],'String',str,'FitBoxToText','on','Interpreter','none');
    exportgraphics(th,'collections.pdf','ContentType','vector','Append',true)

    fhandles={sens_reg_bar_fh,sens_GLM_fh,dur_reg_bar_fh,dur_GLM_fh,dur_sens_corr_fh,inter_wave_fh};

    for hc=fhandles
        hh=hc{1};
        if isa(hh,'matlab.ui.Figure')
            disp('hh')
            exportgraphics(hh,'collections.pdf','ContentType','vector','Append',true)
        elseif isa(hh,"struct")
            fn=fieldnames(hh);
            disp(fn)
            for fh=reshape(fn,1,[])
                exportgraphics(hh.(fh{1}),'collections.pdf','ContentType','vector','Append',true)
            end
        end

    end


    % for fn=fieldnames(sens_reg_bar_fh)

    %feat1: sens_only, feat2: mixed
    %%F2
    %TODO Heatmap session, corr stats

    %TODO Heatmap population

    % TODO TCOM-density
    % TODO:get_com_map
    % TODO:get_region_COM
    % wave.per_region_COM_frac



end


%% RANKSUM with merged bins
if false
    sens_meta=ephys.get_sens_meta('merge_bin',true,'permutation',false);
    stats_type='RANKSUM_merge_bin';
    % load perm_sens.mat
    [sens_map_cells,sens_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','RANKSUM','skip_plot',false,'waveid',sens_meta.wave_id,'range','grey','data_type','sensory','stats_type',stats_type);
    sens_GLM_fh=wave.connectivity_proportion_GLM(sens_map_cells,'range','grey','data_type','sensory','stats_type','Ranksum');
    sens_reg_bar_fh.reg_bar.Children(2).YLim=[0.001,0.5];

    % dur_meta=ephys.get_dur_meta();
    dur_meta=ephys.get_dur_meta('merge_bin',true,'permutation',false);
    [dur_map_cells,dur_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','RANKSUM','skip_plot',false,'waveid',dur_meta.wave_id,'range','grey','data_type','duration','stats_type',stats_type);
    dur_GLM_fhgrey=wave.connectivity_proportion_GLM(dur_map_cells,'range','grey','data_type','duration','stats_type',stats_type);
    dur_GLM_fhch=wave.connectivity_proportion_GLM(dur_map_cells,'range','CH','data_type','duration','stats_type',    stats_type);
    dur_GLM_fhctx=wave.connectivity_proportion_GLM(dur_map_cells,'range','CTX','data_type','duration','stats_type',  stats_type);
end

%% list all regions regardless used or not

if false
    tv=idmap.reg2tree.values();
    for ii=1:numel(tv)
        if numel(tv{ii})>=7
            ureg=[ureg,tv{ii}{7}];
        end
    end
    ureg=unique(ureg);
    % grey_regs=ephys.getGreyRegs();
    for reg=reshape(ureg,1,[])
        fprintf("%s = %s\n",reg{1},char(idmap.reg2full(reg{1})))
    end
end


[tbl,chi,p]=crosstab((1:133)>62,[(1:62)>34,(1:71)>28])