%% ANOVA2
if exist('collections.pdf','file')
    delete('collections.pdf')
end
[sens_map_cells,sens_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','ANOVA2','skip_plot',false);
sens_GLM_fh=wave.connectivity_proportion_GLM(sens_map_cells);

[dur_map_cells,dur_reg_bar_fh]=ephys.duration_reg_bars('stats_model','ANOVA2','skip_plot',false);
dur_GLM_fh=wave.connectivity_proportion_GLM(dur_map_cells);

dur_sens_corr_fh=hier.sens_dur_corr(dur_map_cells{1},sens_map_cells{1});

inter_wave_fh=bz.inter_wave_ext_bars(anovameta);


figure('Color','w')
str='Stats:ANOVA2'
annotation('textbox',[0.2,0.2,0.6,0.6],'String',str,'FitBoxToText','on');

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

%% RANKSUM1
meta=ephys.util.load_meta();
sens_waveid=ephys.get_wave_id(meta.sess,meta.allcid);
[sens_map_cells,sens_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','RANKSUM','skip_plot',false,'waveid',sens_waveid);
sens_GLM_fh=wave.connectivity_proportion_GLM(sens_map_cells);

dur_meta=ephys.get_dur_meta();
[dur_map_cells,dur_reg_bar_fh]=ephys.Both_either_reg_bars('stats_model','RANKSUM','skip_plot',false,'waveid',dur_meta.wave_id);
dur_GLM_fh=wave.connectivity_proportion_GLM(dur_map_cells);

dur_sens_corr_fh=hier.sens_dur_corr(dur_map_cells{1},sens_map_cells{1});

%% RANKSUM2
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

