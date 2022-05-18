warning('gather all plots from top')
keyboard()
%% ANOVA2
if exist('collections.pdf','file')
    delete('collections.pdf')
end
close all
% %% ANOVA2
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

% TODO TCOM-density
% TODO:get_com_map
% TODO:get_region_COM
% wave.per_region_COM_frac


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


    % TODO TCOM-density
    % TODO:get_com_map
    % TODO:get_region_COM
    wave.per_region_COM_frac




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


%% RANKSUM1
close all

th=figure('Color','w','Position',[32,32,1600,900]);
blame=vcs.blame();
str=strjoin({'Stats:RANKSUM1',blame.datetime,blame.git_hash,blame.git_status},'=====\n');
annotation('textbox',[0.05,0.05,0.9,0.9],...
    'String',str,'FitBoxToText', ...
    'on','Interpreter','none');
exportgraphics(th,'collections.pdf','ContentType','vector','Append',true)

% sens_meta=ephys.get_sens_meta();
% sens_waveid=ephys.get_wave_id(meta);
stats_type='RANKSUM_per_bin';
load perm_sens.mat
[sens_map_cells,sens_reg_bar_fh]=ephys.Both_either_reg_bars( ...
    'stats_model','RANKSUM', ...
    'skip_plot',false, ...
    'waveid',sens_meta.wave_id, ...
    'range','grey', ...
    'data_type','sensory', ...
    'stats_type',stats_type);
sens_GLM_fh=wave.connectivity_proportion_GLM( ...
    sens_map_cells, ...
    'range','grey', ...
    'data_type','sensory', ...
    'stats_type','Ranksum');
% sens_reg_bar_fh.reg_bar.Children(2).YLim=[0.001,0.5];

% dur_meta=ephys.get_dur_meta();
load perm_dur.mat
[dur_map_cells,dur_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','RANKSUM','skip_plot',false,'waveid',dur_meta.wave_id,'range','grey','data_type','duration','stats_type',stats_type);
dur_GLM_fhgrey=wave.connectivity_proportion_GLM(dur_map_cells,'range','grey','data_type','duration','stats_type',stats_type);
dur_GLM_fhch=wave.connectivity_proportion_GLM(dur_map_cells,'range','CH','data_type','duration','stats_type',    stats_type);
dur_GLM_fhctx=wave.connectivity_proportion_GLM(dur_map_cells,'range','CTX','data_type','duration','stats_type',  stats_type);


% plot duration SC
if false
    dur_wo_sens=find(ismember(dur_meta.wave_id,5:6));
    meta=ephys.util.load_meta('skip_stats',true);
    for ii=1:numel(dur_wo_sens)
        fh=ephys.sens_dur_SC(dur_wo_sens(ii),meta,sens_meta,dur_meta,'skip_raster',true);
        %     fh=ephys.sens_dur_SC(9071,meta,sens_meta,dur_meta,'skip_raster',false);%
        %     881,9071
        waitfor(fh);
    end
end


% dur_reg_bar_fh.reg_bar.Children(2).YLim=[0.0001,0.5];

dur_sens_corr_fh=hier.sens_dur_corr(dur_map_cells{3},sens_map_cells{3});

single_mix_meta.single1=sens_meta.wave_id>0;
single_mix_meta.single2=dur_meta.wave_id>0;
[singlemix_map_cells,singlemix_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','SINGLE_MIX','skip_plot',false,'meta',single_mix_meta,'range','grey','data_type','single_mix','stats_type',stats_type);
singlemix_GLM_fh=wave.connectivity_proportion_GLM(singlemix_map_cells,'range','grey','data_type','single_mix','stats_type',stats_type);

perm_meta.sens_only=sens_meta.wave_id>0;
perm_meta.dur_only=dur_meta.wave_id>0;
perm_meta.mixed=sens_meta.wave_id>0 & dur_meta.wave_id>0;
inter_wave_fh=bz.inter_wave_ext_bars(perm_meta);

% TCOM
% sense, 6s
sens_tcom_fh=struct();
tcom_maps=cell(1,2);
for currdelay=[6,3]
    sens_com=wave.get_com_map('sel_meta',sens_meta, ...
        'wave',['any',num2str(currdelay)], ... %indep+dep
        'delay',currdelay);
    [fcom.(['d',num2str(currdelay)]).collection,fcom.(['d',num2str(currdelay)]).com_meta]=wave.per_region_COM(sens_com,...
        'stats_method','mean');
    ureg=intersect(ephys.getGreyRegs('range','grey'),...
        fcom.(['d',num2str(currdelay)]).collection(:,2));
    ffrac.(['d',num2str(currdelay)]).collection=[num2cell(cellfun(@(x) x(1),sens_map_cells{1}.values(ureg))),...
        ureg,...
        num2cell(ones(numel(ureg),1)*5),...
        num2cell(cellfun(@(x) x(3),sens_map_cells{1}.values(ureg)))];
    
    [sens_tcom_fh.(['d',num2str(currdelay)]),t]=wave.per_region_COM_frac(fcom.(['d',num2str(currdelay)]),ffrac.(['d',num2str(currdelay)]),'hier_reg','AON');
    t.Title.String=['Sense wave, delay=',num2str(currdelay)];
    [~,tcidx]=ismember(ureg,fcom.(['d',num2str(currdelay)]).collection(:,2));
    tcom_maps{currdelay/3}=containers.Map(ureg,num2cell(cellfun(@(x) x/4, fcom.(['d',num2str(currdelay)]).collection(tcidx,1))));
end
sensory_TCOM_GLM_fh=wave.connectivity_proportion_GLM(tcom_maps,'range','grey','data_type','3s_6s_sensory_TCOM','stats_type',stats_type);


fh=cell(0);
dur_com=wave.get_dur_com_map(dur_meta);

[fcom.dur.collection,fcom.dur.com_meta]=wave.per_region_COM(dur_com,...
    'stats_method','mean');
ureg=intersect(ephys.getGreyRegs('range','grey'),...
    fcom.dur.collection(:,2));

%total number of duration neuron (indep+dep)
ffrac.dur.collection=[num2cell(cellfun(@(x) x(1),dur_map_cells{3}.values(ureg))),...
    ureg,...
    num2cell(ones(numel(ureg),1)*5),...
    num2cell(cellfun(@(x) x(3),dur_map_cells{3}.values(ureg)))];

[dur_tcom_fh,t]=wave.per_region_COM_frac(fcom.dur,ffrac.dur,'hier_reg','COA','range','grey');
t.Title.String='Duration wave';


[~,tcidx]=ismember(ureg,fcom.dur.collection(:,2));
dur_com_maps{1}=containers.Map(ureg,num2cell(cellfun(@(x) x/4, fcom.dur.collection(tcidx,1))));
duration_TCOM_GLM_fh=wave.connectivity_proportion_GLM(dur_com_maps,'range','grey','data_type','duration_TCOM','stats_type',stats_type);
duration_TCOM_GLM_fh=wave.connectivity_proportion_GLM(dur_com_maps,'range','CH','data_type','duration_TCOM','stats_type',stats_type);
duration_TCOM_GLM_fh=wave.connectivity_proportion_GLM(dur_com_maps,'range','CTX','data_type','duration_TCOM','stats_type',stats_type);

sens_dur_TCOM_corr_fh=wave.sens_dur_TCOM_corr(fcom);

fhandles=get(groot(),'Children');

for hc=reshape(fhandles,1,[])
    exportgraphics(hc,'collections.pdf','ContentType','vector','Append',true)
end
savefig(fhandles,'Ranksum1.fig');


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
