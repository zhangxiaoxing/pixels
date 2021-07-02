function [spk_out,skip,pp,rsq]=history_coeff_one(sessid,suid,opt)
arguments
    sessid (1,1) int32
    suid (1,2) int32
    opt.tsbin_size (1,1) double = 600
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO','MY'})}='neupix'
    opt.laser (1,:) char {mustBeMember(opt.laser,{'on','off','any'})} = 'any'
    opt.epoch (1,:) char {mustBeMember(opt.epoch,{'delay','ITI','any'})} = 'any'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
    opt.correct_error (1,:) char {mustBeMember(opt.correct_error,{'correct','error','any'})} = 'any'
end
persistent bitmask X
if isempty(bitmask) || isempty(X)
    bitmask=2.^(0:9)';
    X=buildX();
end

%whole session
skip=true;
if strcmp(opt.correct_error,'any')
    if strcmp(opt.type,'neupix') || strcmp(opt.type,'MY')
        [spkID,spkTS,~,~,~]=ephys.getSPKID_TS(sessid,...
            'criteria',opt.criteria);
        tspre=spkTS(spkID==suid(1));
        tspost=spkTS(spkID==suid(2));
    else
        [tspre,tspost]=ephys.getSPKID_TS_HEM(sessid,suid(1),suid(2),'laser',opt.laser);
    end
    
    if (~isempty(tspre)) && ~isempty(tspost)
        tmax=max([tspre;tspost]);
        histpre=histcounts(tspre,1:opt.tsbin_size:tmax)>0;
        histpost=histcounts(tspost,1:opt.tsbin_size:tmax)>0;
        post_spike_prob=zeros(1024,2);
        %TODO: per trial
        for i=1:(length(histpre)-10)
            if any(histpre(i:i+9))
                hist_type=histpre(i:i+9)*bitmask; %will supply last bin later
            else
                hist_type=0;
            end
            
            post_spike_prob(hist_type+1,2)=post_spike_prob(hist_type+1,2)+histpost(i+10); % blind detect post spike
            post_spike_prob(hist_type+1,1)=post_spike_prob(hist_type+1,1)+1;
            
        end
    end
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% correct error trials %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    post_spike_prob=zeros(1024,2);
    [~,~,trials,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,...
        'criteria',opt.criteria,'keep_trial',true);
    if strcmp(opt.correct_error,'correct')
        target_trials=find(all(trials(:,9:10)>0,2));
    elseif strcmp(opt.correct_error,'error')
        target_trials=find(~(trials(:,10)));
    end
    presel=str2double(FT_SPIKE.label)==suid(1);
    postsel=str2double(FT_SPIKE.label)==suid(2);
    for tt=reshape(target_trials,1,[])
        tspre=FT_SPIKE.timestamp{presel}(FT_SPIKE.trial{presel}==tt);
        tspost=FT_SPIKE.timestamp{postsel}(FT_SPIKE.trial{postsel}==tt);
        if isempty(tspre) || isempty(tspost)
            continue
        end
        skip=false;
        histpre=histcounts(tspre,trials(tt,1)-30000*3:opt.tsbin_size:trials(tt,1)+30000*11)>0;
        histpost=histcounts(tspost,trials(tt,1)-30000*3:opt.tsbin_size:trials(tt,1)+30000*11)>0;
        %TODO: per trial
        for i=1:(length(histpre)-10)
            if any(histpre(i:i+9))
                hist_type=histpre(i:i+9)*bitmask; %will supply last bin later
            else
                hist_type=0;
            end
            post_spike_prob(hist_type+1,2)=post_spike_prob(hist_type+1,2)+histpost(i+10); % blind detect post spike
            post_spike_prob(hist_type+1,1)=post_spike_prob(hist_type+1,1)+1;
        end
    end
end
%% 

spksel=post_spike_prob(:,1)>0;
if (~skip) && nnz(spksel)>1
    glmopt=statset('fitglm');
    glmopt.MaxIter=1000;
    spk_mdl=fitglm(X(spksel,:),post_spike_prob(spksel,[2,1]),'Distribution','binomial','Link','identity','Options',glmopt);
    spk_out=spk_mdl.Coefficients.Estimate;
    skip=checkWarning();
    pp=spk_mdl.coefTest();
    rsq=spk_mdl.Rsquared.Ordinary;
else
    spk_out=zeros(1,11);
    skip=true;
    pp=1;
    rsq=0;
end
end

function X=buildX()
X=zeros(1024,10);
for row=1:1024
    for col=1:10
        if bitand(row-1,2.^(col-1))
            X(row,col)=1;
        end
    end
end
end

function out=checkWarning()
[~,wid]=lastwarn('','zx:unset');
out=strcmp(wid,'stats:glmfit:IterationLimit');
end
