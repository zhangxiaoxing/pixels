function [spk_out,fc_eff_out,fc_prob_out,maxiter]=history_coeff(sessid,suid,opt)
arguments
    sessid (1,1) int32
    suid (1,2) int32
    opt.tsbin_size (1,1) double = 600
    opt.postspike (1,1) logical = true
    opt.fc_effi (1,1) logical = false
    opt.fc_prob (1,1) logical = false
end
persistent bitmask X
if isempty(bitmask) || isempty(X)
    bitmask=2.^(0:9)';
    X=buildX();
end
[spkID,spkTS,~,~,~]=ephys.getSPKID_TS(sessid);
tspre=spkTS(spkID==suid(1));
tspost=spkTS(spkID==suid(2));
tmax=max([tspre;tspost]);
histpre=histcounts(tspre,0:opt.tsbin_size:tmax)>0;
histpost=histcounts(tspost,0:opt.tsbin_size:tmax)>0;
binpre=discretize(tspre,0:opt.tsbin_size:tmax);
binpost=discretize(tspost,0:opt.tsbin_size:tmax);
post_spike_prob=zeros(1024,2);
fc_effi=zeros(1024,2);
fc_prob=zeros(1024,1);
pre_ptr=1;precnt=numel(tspre);%for performance optimizaiton
post_ptr=1;postcnt=numel(tspost);
for i=1:(length(histpre)-10)
    if any(histpre(i:i+9))
        hist_type=histpre(i:i+9)*bitmask;
    else
        hist_type=0;
    end
    post_spike_prob(hist_type+1,1)=post_spike_prob(hist_type+1,1)+1;
    post_spike_prob(hist_type+1,2)=post_spike_prob(hist_type+1,2)+histpost(i+9); % blind detect post spike
    
    while pre_ptr<=precnt && binpre(pre_ptr)<i+10, pre_ptr=pre_ptr+1; end
    while post_ptr<=postcnt && binpost(post_ptr)<i+10, post_ptr=post_ptr+1; end
    
    if pre_ptr>precnt || binpre(pre_ptr)>i+10 %not threre yet
        continue
    elseif binpre(pre_ptr)==i+10
        presel=[];postsel=[];
        while pre_ptr<=precnt && binpre(pre_ptr)==i+10
            presel=[presel,pre_ptr];
            pre_ptr=pre_ptr+1;
        end
        while post_ptr<=postcnt && binpost(post_ptr)==i+10
            postsel=[postsel,post_ptr];
            post_ptr=post_ptr+1;
        end
        if nnz(postsel)>0
            [fchit,fcmiss]=getfc(tspre(presel),tspost(postsel));
        else
            fchit=0;fcmiss=numel(presel);
        end
        fc_effi(hist_type+1,1)=fc_effi(hist_type+1,1)+fchit+fcmiss;
        fc_effi(hist_type+1,2)=fc_effi(hist_type+1,2)+fchit;
        if fchit>0
            fc_prob(hist_type+1,1)=fc_prob(hist_type+1,1)+1;
        end
    end
end
maxiter=false(1,3);

if opt.postspike
    spksel=post_spike_prob(:,1)>0;
    spk_mdl=fitglm(X(spksel,:),post_spike_prob(spksel,[2,1]),'Distribution','binomial','Link','identity');
    spk_out=spk_mdl.Coefficients.Estimate;
    maxiter(1)=checkWarning();
else
    spk_out=zeros(1,11);
end
if opt.fc_effi
    fceffsel=fc_effi(:,1)>0;
    fc_eff_mdl=fitglm(X(fceffsel,:),fc_effi(fceffsel,[2,1]),'Distribution','binomial','Link','identity');
    fc_eff_out=fc_eff_mdl.Coefficients.Estimate;
    maxiter(2)=checkWarning();
else
    fc_eff_out=zeros(1,11);
end

if opt.fc_prob
    fc_mdl=fitglm(X(spksel,:),[fc_prob(spksel),post_spike_prob(spksel,1)],'Distribution','binomial','Link','identity');
    fc_prob_out=fc_mdl.Coefficients.Estimate;
    maxiter(3)=checkWarning();
else
    fc_prob_out=zeros(1,11);
end


end

function [hit,miss]=getfc(pre,post)
lag=post-(pre.');
hit=nnz(any(lag>15 & lag<300,1)); %assuming 30K sps.
miss=numel(pre)-hit;
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