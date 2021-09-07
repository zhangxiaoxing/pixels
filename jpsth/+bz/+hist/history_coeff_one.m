function [postspk_,fc_]=history_coeff_one(sessid,suid,opt)
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
if strcmp(opt.correct_error,'any')
    if strcmp(opt.type,'neupix') || strcmp(opt.type,'MY')
        [spkID,spkTS,~,~,~]=ephys.getSPKID_TS(sessid,...
            'criteria',opt.criteria);
        tspre=spkTS(spkID==suid(1));
        tspost=spkTS(spkID==suid(2));
    else
        [tspre,tspost]=ephys.getSPKID_TS_HEM(sessid,suid(1),suid(2),'laser',opt.laser);
    end
    
    tmax=max([tspre;tspost]);
    histpre=histcounts(tspre,1:opt.tsbin_size:tmax)>0;
    histpost=histcounts(tspost,1:opt.tsbin_size:tmax)>0;
    post_spike_prob=zeros(1024,2);
    %per trial
    for i=1:(length(histpre)-10)
        if any(histpre(i:i+9))
            hist_type=histpre(i:i+9)*bitmask; %will supply last bin later
        else
            hist_type=0;
        end
        
        post_spike_prob(hist_type+1,2)=post_spike_prob(hist_type+1,2)+histpost(i+10); % blind detect post spike
        post_spike_prob(hist_type+1,1)=post_spike_prob(hist_type+1,1)+1;
    end
    
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% correct error trials %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    post_spike_prob=zeros(1024,2);
    homo_plastic=zeros(1024,2);
    
    [~,~,trials,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,...
        'criteria',opt.criteria,'keep_trial',true);
    if strcmp(opt.correct_error,'correct')
        target_trials=find(all(trials(:,9:10)>0,2));
    elseif strcmp(opt.correct_error,'error')
        target_trials=find(~(trials(:,10)));
    end
    presel=str2double(FT_SPIKE.label)==suid(1);
    postsel=str2double(FT_SPIKE.label)==suid(2);
    % per trial
    for tt=reshape(target_trials,1,[])
        %         disp(tt)
        %         if tt==28, keyboard;end
        tspre=FT_SPIKE.timestamp{presel}(FT_SPIKE.trial{presel}==tt);
        tspost=FT_SPIKE.timestamp{postsel}(FT_SPIKE.trial{postsel}==tt);
        tonset=trials(tt,1)-30000*3;
        toffset=trials(tt,1)+30000*11;
        tedges=tonset:opt.tsbin_size:toffset;
        histpre=histcounts(tspre,tedges)>0;
        histpost=histcounts(tspost,tedges)>0;
        for i=1:(length(histpre)-10)
            if any(histpre(i:i+9))
                hist_type=histpre(i:i+9)*bitmask; %skipped last bin for balanced design
            else
                hist_type=0;
            end
            post_spike_prob(hist_type+1,2)=post_spike_prob(hist_type+1,2)+histpost(i+10); % fast-forward detect post spike
            post_spike_prob(hist_type+1,1)=post_spike_prob(hist_type+1,1)+1;
            
            %% fc effi
            % TODO post sel
            tpresel=tspre> tedges(i+10) & tspre<=tedges(i+11);
            tpostsel=tspost> tedges(i+10) & tspost<=tedges(i+11);
            if nnz(tpostsel)>0
                [fchit,fcmiss]=getfc(tspre(tpresel),tspost(tpostsel));
            else
                fchit=0;fcmiss=numel(tpresel);
            end
            homo_plastic(hist_type+1,1)=homo_plastic(hist_type+1,1)+fchit+fcmiss;
            homo_plastic(hist_type+1,2)=homo_plastic(hist_type+1,2)+fchit;
        end
        %% end of fc effi
        
    end
end
%% postspk
spksel=post_spike_prob(:,1)>0;
postspk_=struct();
if nnz(spksel)>20
    try
        glmopt=statset('fitglm');
        glmopt.MaxIter=1000;
        glmopt.Display='iter';
        spk_mdl=fitglm(X(spksel,:),post_spike_prob(spksel,[2,1]),'linear','Distribution','binomial','Link','identity','Options',glmopt);
        postspk_.coeff=spk_mdl.Coefficients.Estimate;
        postspk_.skip=checkWarning();
        postspk_.pp=spk_mdl.coefTest();
        postspk_.rsq=spk_mdl.Rsquared.Ordinary;
        %         if postspk_.rsq<0,keyboard();end
    catch ME
        postspk_.coeff=zeros(1,11);
        postspk_.skip=true;
        postspk_.pp=-1;
        postspk_.rsq=0;
    end
else
    postspk_.coeff=zeros(1,11);
    postspk_.skip=true;
    postspk_.pp=1;
    postspk_.rsq=0;
end
%% fc effi
if ~strcmp(opt.correct_error,'any')
    fcsel=homo_plastic(:,1)>0;
    if nnz(fcsel)>20
        try
            glmopt=statset('fitglm');
            glmopt.MaxIter=1000;
            fc_eff_mdl=fitglm(X(fcsel,:),homo_plastic(fcsel,[2,1]),'Distribution','binomial','Link','identity','Options',glmopt);
            fc_.coeff=fc_eff_mdl.Coefficients.Estimate;
            fc_.skip=checkWarning();
            fc_.pp=fc_eff_mdl.coefTest();
            fc_.rsq=fc_eff_mdl.Rsquared.Ordinary;
        catch ME
            fc_.coeff=zeros(1,11);
            fc_.skip=true;
            fc_.pp=-1;
            fc_.rsq=0;
            %         if fc_.rsq<0,keyboard();end
        end
    else
        fc_.coeff=zeros(1,11);
        fc_.skip=true;
        fc_.pp=1;
        fc_.rsq=0;
    end
end
%% fc effi end

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

function [hit,miss]=getfc(pre,post)
%TODO proper handling pre spikes
lag=arrayfun(@(x) post-(x),pre,'UniformOutput',false);
hit=nnz(cellfun(@(x) any(x>=24 & x<300,'all'),lag)); %assuming 30K sps,0.8ms min latency according to English 2017 code
miss=numel(pre)-hit;
end