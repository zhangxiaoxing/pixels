%% ANOVA2
anova2meta=ephys.selectivity_anova('merge_time_bin',true);

%% ANOVA2ALT
% anova2_alt_meta=ephys.selectivity_anova('largest_varied_bin',true);

%% RANKSUM1
meta=ephys.util.load_meta();
sens_waveid1=ephys.get_wave_id(meta);
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
sens_wo_dur=find(ismember(sens_waveid1,1:4) & ismember(dur_waveid1,1:4));

for ii=1:numel(sens_wo_dur)
    fh=ephys.sens_dur_SC(sens_wo_dur(ii));
    waitfor(fh);
end
%% mixed but independent
sens_meta.allcid=meta.allcid;
sens_meta.allpath=meta.allpath;
sens_meta.sess=meta.sess;

dur_waveid1=dur_meta.wave_id;
dur_dep=find(sens_meta.wave_id>4 & dur_meta.wave_id>4);
% dur_dep=find(ismember(dur_waveid1,1:4) & sens_waveid1==0);
for ii=1:numel(dur_dep)
    fh=ephys.sens_dur_SC(dur_dep(ii),sens_meta,dur_meta,'anova_meta',anova2meta);
    waitfor(fh);
end


%% dependent but single-modality
dur_dep=find(ismember(sens_meta.wave_id,1:4) & dur_meta.wave_id==0);
% dur_dep=find(ismember(dur_waveid1,1:4) & sens_waveid1==0);
for ii=1:numel(dur_dep)
    fh=ephys.sens_dur_SC(dur_dep(ii),sens_meta,dur_meta,'anova_meta',anova2meta);
    waitfor(fh);
end




%%
sranksel={ismember(sens_meta.wave_id,5:6),ismember(sens_meta.wave_id,1:4),sens_meta.wave_id==0};
anovasel={anova2meta.anovap(:,1)<0.05,anova2meta.anovap(:,2)<0.05,anova2meta.anovap(:,3)<0.05,anova2meta.non_mem};

for ridx=1:3
disp(arrayfun(@(x) nnz(sranksel{ridx} & anovasel{x}),1:4))
end


%%
dranksel={ismember(dur_meta.wave_id,5:6),ismember(dur_meta.wave_id,1:4),dur_meta.wave_id==0};
for ridx=1:3
disp(arrayfun(@(x) nnz(dranksel{ridx} & anovasel{x}),1:4))
end


%%
for ridx=1:3
disp(arrayfun(@(x) nnz(sranksel{ridx} & dranksel{x}),1:3))
end


