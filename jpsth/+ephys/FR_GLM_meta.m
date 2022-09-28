%TODO brain region filter, olfaction filter.
function out=FR_GLM_meta(opt)
arguments
    opt.su_list % su of interest
    opt.sess
    opt.save_file (1,1) logical = false
    opt.fn (1,:) char = 'GLM_meta.mat'
    opt.collect_file (1,1) logical = false
end
if isunix && opt.collect_file
    fs=split(ls('GLM_meta*.mat'));
    fs=fs(~ismissing(fs));
    fstr=load(fs{1});
    out=fstr.out;
    fldn=fieldnames(fstr.out);
    for fid=2:numel(fs)
        fstr=load(fs{fid});
        for fldidx=1:numel(fldn)
            out.(fldn{fldidx})=[out.(fldn{fldidx});...
                fstr.out.(fldn{fldidx})];
        end
    end
    blame=vcs.blame();
    save('GLM_summed_meta.mat','out','blame');
    return
end
zthres=norminv(0.95);

%% gen data
curr_sess=-1;
meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
homedir=ephys.util.getHomedir('type','raw');

if isfield(opt,'su_list') && ~isempty(opt.su_list)
    globalidx=reshape(opt.su_list,1,[]);
elseif isfield(opt,'sess') && ~isempty(opt.sess)
    globalidx=reshape(find(meta.sess==opt.sess),1,[]);
    opt.su_list=globalidx;
    opt.save_file=true;
    opt.fn=sprintf('GLM_meta_sess%03d.mat',opt.sess);
else
    globalidx=1:numel(meta.allcid);
end

out=struct();
out.global_idx=nan(numel(globalidx),1);
out.samp_shuf_mm=nan(numel(globalidx),3);
out.samp_shuf_std=nan(numel(globalidx),3);
out.dur_shuf_mm =nan(numel(globalidx),3);
out.dur_shuf_std =nan(numel(globalidx),3);
out.full_mm =nan(numel(globalidx),3);

out.samp_sel=nan(numel(globalidx),3);
out.dur_sel =nan(numel(globalidx),3);

for gidx=1:numel(globalidx)
    ii=opt.su_list(gidx);
    if ispc
        disp(gidx)
    end
    if meta.sess(ii)~=curr_sess
        disp(meta.sess(ii))
        fpath=fullfile(homedir,ephys.sessid2path(meta.sess(ii)),"FR_All_1000.hdf5");
        fpath=replace(fpath,'\',filesep());
        fr=h5read(fpath,'/FR_All'); %Trial x SU x time-bin
        trials=h5read(fpath,'/Trials');
        suid=h5read(fpath,'/SU_id');
        curr_sess=meta.sess(ii);          
        csel=trials(:,9)~=0 & trials(:,10)~=0 & ismember(trials(:,8),[3 6]) & ismember(trials(:,5),[4 8]);
        % keep only correct well trained trials for stratified CV
        trials=trials(csel,:);
        fr=fr(csel,:,:);
        
        cvgrp=nan(size(fr,1),1);
        cvgrp(trials(:,5)==4 & trials(:,8)==3)=1;
        cvgrp(trials(:,5)==8 & trials(:,8)==3)=2;
        cvgrp(trials(:,5)==4 & trials(:,8)==6)=3;
        cvgrp(trials(:,5)==8 & trials(:,8)==6)=4;
    end
    suidx=(suid==meta.allcid(ii));
    
    pcc_full=nan(20,3); % rpt*kfold, bins
    cv=cvpartition(cvgrp,'KFold',2);
    
    for rpt=1:10
        for kf=1:2 % 2-fold cv real data 10 repeat
            x_samp=trials(training(cv,kf),5);
            x_dur=trials(training(cv,kf),8);
            for bin=5:7 % 3-bin (optional TODO: 6-bin)
                y=fr(training(cv,kf),suidx,bin);
                % 4-class FR-full regressor GLM model
                mdl=fitglm([x_samp,x_dur],y,'linear','CategoricalVars',[1,2]);
                pcc_full(rpt*2-2+kf,bin-4)=corr(y,mdl.Fitted.Response);
            end
        end
        cv=repartition(cv);
    end
    % 400 shuffle regressor
    pcc_wo_samp=nan(400,3); % rpt*kfold, bins
    pcc_wo_dur=nan(400,3); % rpt*kfold, bins
    for rpt=1:200
        for kf=1:2 % 2-fold cv real data 10 repeat
            x_samp=trials(training(cv,kf),5);
            x_dur=trials(training(cv,kf),8);
            for bin=5:7 % 3-bin (optional TODO: 6-bin)
                y=fr(training(cv,kf),suidx,bin);
                % 4-class FR-full regressor GLM model
                mdl_wo_samp=fitglm([datasample(x_samp,numel(x_samp),'Replace',false),x_dur],...
                    y,'linear','CategoricalVars',[1,2]);
                mdl_wo_dur=fitglm([x_samp,datasample(x_dur,numel(x_dur),'Replace',false)],...
                    y,'linear','CategoricalVars',[1,2]);
                pcc_wo_samp(rpt*2-2+kf,bin-4)=corr(y,mdl_wo_samp.Fitted.Response);
                pcc_wo_dur(rpt*2-2+kf,bin-4)=corr(y,mdl_wo_dur.Fitted.Response);
            end
        end
        cv=repartition(cv);
    end

    wosmm=mean(pcc_wo_samp);
    wosstd=std(pcc_wo_samp);

    wodmm=mean(pcc_wo_dur);
    wodstd=std(pcc_wo_dur);

    pccmm=mean(pcc_full);

    shuf_wo_samp=wosmm+zthres.*wosstd;
    shuf_wo_dur=wodmm+zthres.*wodstd;

    out.global_idx(gidx)=opt.su_list(gidx);
    out.samp_shuf_mm(gidx,:)=wosmm;
    out.samp_shuf_std(gidx,:)=wosstd;
    out.dur_shuf_mm(gidx,:)=wodmm;
    out.dur_shuf_std(gidx,:)=wodstd;
    out.full_mm(gidx,:)=pccmm;
    
    out.samp_sel(gidx,:)=pccmm>shuf_wo_samp;
    out.dur_sel(gidx,:)=pccmm>shuf_wo_dur;
end

if opt.save_file
    blame=vcs.blame();
    save(opt.fn,'out','blame')
end

end
