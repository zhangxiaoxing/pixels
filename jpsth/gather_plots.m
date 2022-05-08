%% ANOVA2
if exist('collections.pdf','file')
    delete('collections.pdf')
end
close all
% %% ANOVA2
anova2meta=ephys.selectivity_anova('merge_time_bin',true,'anova_model','full');
% 
% %% ANOVA2ALT
% anova2_alt_meta=ephys.selectivity_anova('largest_varied_bin',true);


% anovameta=ephys.selectivity_anova('merge_time_bin',true);

anovameta=anova2meta;

[sens_map_cells,sens_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','ANOVA2','skip_plot',false,'meta',anovameta,'single_field','sens');
sens_GLM_fh=wave.connectivity_proportion_GLM(sens_map_cells);

[dur_map_cells,dur_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','ANOVA2','skip_plot',false,'meta',anovameta,'single_field','dur');
dur_GLM_fh=wave.connectivity_proportion_GLM(dur_map_cells,'range','CTX');

dur_sens_corr_fh=hier.sens_dur_corr(dur_map_cells{1},sens_map_cells{1});

single_mix_meta.single1=anovameta.sens;
single_mix_meta.single2=anovameta.dur;
[singlemix_map_cells,singlemix_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','SINGLE_MIX','skip_plot',false,'meta',single_mix_meta);
singlemix_GLM_fh=wave.connectivity_proportion_GLM(singlemix_map_cells);

anova2meta.sens_only=anova2meta.sens & ~anova2meta.dur;
anova2meta.dur_only=~anova2meta.sens & anova2meta.dur;
anova2meta.mixed=anova2meta.sens & anova2meta.dur;
inter_wave_fh=bz.inter_wave_ext_bars(anova2meta);

th=figure('Color','w','Position',[32,32,1600,900]);
blame=vcs.blame();
str=strjoin({'Stats:ANOVA2',blame.datetime,blame.git_hash,blame.git_status},'=====\n');
annotation('textbox',[0.05,0.05,0.9,0.9],'String',str,'FitBoxToText','on');
exportgraphics(th,'collections.pdf','ContentType','vector','Append',true)

handles={sens_reg_bar_fh,sens_GLM_fh,dur_reg_bar_fh,dur_GLM_fh,dur_sens_corr_fh,inter_wave_fh};

for hc=handles
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


%% ANOVA per bin

anova3meta=ephys.selectivity_anova('per_bin',true);

[sens_map_cells,sens_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','ANOVA2','skip_plot',false,'meta',anova3meta,'single_field','sens');
sens_GLM_fh=wave.connectivity_proportion_GLM(sens_map_cells);

[dur_map_cells,dur_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','ANOVA2','skip_plot',false,'meta',anova3meta,'single_field','dur');
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
annotation('textbox',[0.05,0.05,0.9,0.9],'String',str,'FitBoxToText','on');
exportgraphics(th,'collections.pdf','ContentType','vector','Append',true)

handles={sens_reg_bar_fh,sens_GLM_fh,dur_reg_bar_fh,dur_GLM_fh,dur_sens_corr_fh,inter_wave_fh};

for hc=handles
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


% ANOVA3
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
annotation('textbox',[0.05,0.05,0.9,0.9],'String',str,'FitBoxToText','on');
exportgraphics(th,'collections.pdf','ContentType','vector','Append',true)

handles={sens_reg_bar_fh,sens_GLM_fh,dur_reg_bar_fh,dur_GLM_fh,dur_sens_corr_fh,inter_wave_fh};

for hc=handles
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






%% RANKSUM1 %TODO pivot to permutation
close all

th=figure('Color','w','Position',[32,32,1600,900]);
blame=vcs.blame();
str=strjoin({'Stats:RANKSUM1',blame.datetime,blame.git_hash,blame.git_status},'=====\n');
annotation('textbox',[0.05,0.05,0.9,0.9],'String',str,'FitBoxToText','on');
exportgraphics(th,'collections.pdf','ContentType','vector','Append',true)

% sens_meta=ephys.get_sens_meta();
% sens_waveid=ephys.get_wave_id(meta);
load perm_sens.mat
[sens_map_cells,sens_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','RANKSUM','skip_plot',false,'waveid',sens_meta.wave_id);
sens_GLM_fh=wave.connectivity_proportion_GLM(sens_map_cells);
sens_reg_bar_fh.reg_bar.Children(2).YLim=[0.001,0.5];

% dur_meta=ephys.get_dur_meta();
load perm_dur.mat
[dur_map_cells,dur_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','RANKSUM','skip_plot',false,'waveid',dur_meta.wave_id);
dur_GLM_fhgrey=wave.connectivity_proportion_GLM(dur_map_cells,'range','grey');
dur_GLM_fhch=wave.connectivity_proportion_GLM(dur_map_cells,'range','CH');
dur_GLM_fhctx=wave.connectivity_proportion_GLM(dur_map_cells,'range','CTX');

dur_reg_bar_fh.reg_bar.Children(2).YLim=[0.0001,0.5];

dur_sens_corr_fh=hier.sens_dur_corr(dur_map_cells{2},sens_map_cells{1});

single_mix_meta.single1=sens_meta.wave_id>0;
single_mix_meta.single2=dur_meta.wave_id>0;
[singlemix_map_cells,singlemix_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','SINGLE_MIX','skip_plot',false,'meta',single_mix_meta);
singlemix_GLM_fh=wave.connectivity_proportion_GLM(singlemix_map_cells);

perm_meta.sens_only=sens_meta.wave_id>0 & dur_meta.wave_id==0;
perm_meta.dur_only=sens_meta.wave_id==0 & dur_meta.wave_id>0;
perm_meta.mixed=sens_meta.wave_id>0 & dur_meta.wave_id>0;
inter_wave_fh=bz.inter_wave_ext_bars(perm_meta);

handles={sens_reg_bar_fh,sens_GLM_fh,dur_reg_bar_fh,dur_GLM_fh,dur_sens_corr_fh};

for hc=handles
    hh=hc{1};
    if isa(hh,'matlab.ui.Figure')
%         disp('hh')
        exportgraphics(hh,'collections.pdf','ContentType','vector','Append',true)
    elseif isa(hh,"struct")
        fn=fieldnames(hh);
%         disp(fn)
        for fh=reshape(fn,1,[])
            exportgraphics(hh.(fh{1}),'collections.pdf','ContentType','vector','Append',true)
        end
    end
    
end




dur_waveid=dur_meta.wave_id;
nnz(ismember(sens_waveid,1:4))
nnz(ismember(dur_waveid,1:4))
nnz(ismember(sens_waveid,1:4) & ismember(dur_waveid,1:4))

nnz(ismember(sens_waveid,5:6) & ismember(dur_waveid,5:6))

nnz(ismember(sens_waveid,1:4) & ismember(dur_waveid,5:6))

nnz(ismember(sens_waveid,5:6) & ismember(dur_waveid,1:4))

histogram(dur_waveid(ismember(sens_waveid,1:4)),-0.5:6.5)
histogram(sens_waveid(ismember(dur_waveid,1:4)),-0.5:6.5)

%% RANKSUM2 TODO:merge into ranksum1
close all

th=figure('Color','w','Position',[32,32,1600,900]);
blame=vcs.blame();
str=strjoin({'Stats:RANKSUM2',blame.datetime,blame.git_hash,blame.git_status},'=====\n');
annotation('textbox',[0.05,0.05,0.9,0.9],'String',str,'FitBoxToText','on');
exportgraphics(th,'collections.pdf','ContentType','vector','Append',true)

meta=ephys.util.load_meta();
dur_meta=ephys.get_dur_meta();

sens_waveid=ephys.get_wave_id(meta.sess,meta.allcid);
dur_waveid=dur_meta.wave_id;
waveid=zeros(size(sens_waveid));
waveid(sens_waveid>0 & dur_waveid==0)=1;
waveid(sens_waveid>0 & dur_waveid>0)=2;

[sens_map_cells,sens_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','RANKSUM2','skip_plot',false,'waveid',waveid);
sens_GLM_fh=wave.connectivity_proportion_GLM(sens_map_cells);

waveid=zeros(size(sens_waveid));
waveid(sens_waveid==0 & dur_waveid>0)=1;
waveid(sens_waveid>0 & dur_waveid>0)=2;

[dur_map_cells,dur_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','RANKSUM2','skip_plot',false,'waveid',waveid);
dur_GLM_fh=wave.connectivity_proportion_GLM(dur_map_cells);

dur_sens_corr_fh=hier.sens_dur_corr(dur_map_cells{1},sens_map_cells{1});

wrsmeta.sens_only=sens_waveid>0 & dur_waveid==0;
wrsmeta.dur_only=sens_waveid==0 & dur_waveid>0;
wrsmeta.mixed=sens_waveid>0 & dur_waveid>0;
inter_wave_fh=bz.inter_wave_ext_bars(wrsmeta);

handles={sens_reg_bar_fh,sens_GLM_fh,dur_reg_bar_fh,dur_GLM_fh,dur_sens_corr_fh,inter_wave_fh};

for hc=handles
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

exportgraphics(gcf(),'collections.pdf','ContentType','vector','Append',true)


%%


ureg=cell(0);
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
