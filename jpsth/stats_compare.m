%% ANOVA2
anova2meta=ephys.selectivity_anova('merge_time_bin',true);

%% ANOVA3
anova3meta=ephys.selectivity_anova('merge_time_bin',false);
%% RANKSUM1
meta=ephys.util.load_meta();
sens_waveid1=ephys.get_wave_id(meta.sess,meta.allcid);
dur_meta=ephys.get_dur_meta();
dur_waveid1=dur_meta.wave_id;


%% RANKSUM2

% meta=ephys.util.load_meta();
% dur_meta=ephys.get_dur_meta();
% sens_waveid=ephys.get_wave_id(meta.sess,meta.allcid);
% dur_waveid=dur_meta.wave_id;

waveid2=zeros(size(sens_waveid1));
waveid2(sens_waveid1>0 & dur_waveid1==0)=1;
waveid2(sens_waveid1==0 & dur_waveid1>0)=2;
waveid2(sens_waveid1>0 & dur_waveid1>0)=4;


%% plot SC
sens_wo_dur=find(ismember(sens_waveid1,1:4) & dur_waveid1==0);
nnz(ismember(sens_waveid,5:6) & ismember(dur_waveid,5:6))
nnz(ismember(sens_waveid,1:4) & ismember(dur_waveid,5:6))
nnz(ismember(sens_waveid,5:6) & ismember(dur_waveid,1:4))