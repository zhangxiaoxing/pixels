function [spk_out,skip]=history_coeff(sessid,suid,opt)
arguments
    sessid (1,1) int32
    suid (1,2) int32
    opt.tsbin_size (1,1) double = 600
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO'})}='neupix'
    opt.laser (1,:) char {mustBeMember(opt.laser,{'on','off','any'})} = 'any'
    opt.epoch (1,:) char {mustBeMember(opt.epoch,{'delay','ITI','any'})} = 'any'
end
persistent bitmask X
if isempty(bitmask) || isempty(X)
    bitmask=2.^(0:9)';
    X=buildX();
end
if strcmp(opt.type,'neupix')
    [spkID,spkTS,~,~,~]=ephys.getSPKID_TS(sessid,'epoch',opt.epoch);
    tspre=spkTS(spkID==suid(1));
    tspost=spkTS(spkID==suid(2));
else
    [tspre,tspost]=ephys.getSPKID_TS_HEM(sessid,suid(1),suid(2),'laser',opt.laser);
end

if isempty(tspre) || isempty(tspost)
    spk_out=zeros(1,11);
    skip=true;
    return
end

tmax=max([tspre;tspost]);
histpre=histcounts(tspre,1:opt.tsbin_size:tmax)>0;
histpost=histcounts(tspost,1:opt.tsbin_size:tmax)>0;
post_spike_prob=zeros(1024,2);
for i=1:(length(histpre)-10)
    if any(histpre(i:i+9))
        hist_type=histpre(i:i+9)*bitmask; %will supply last bin later
    else
        hist_type=0;
    end

    post_spike_prob(hist_type+1,2)=post_spike_prob(hist_type+1,2)+histpost(i+10); % blind detect post spike
    post_spike_prob(hist_type+1,1)=post_spike_prob(hist_type+1,1)+1;
   
end
skip=true;

glmopt=statset('fitglm');
glmopt.MaxIter=1000;

spksel=post_spike_prob(:,1)>0;
if nnz(spksel)>1
    spk_mdl=fitglm(X(spksel,:),post_spike_prob(spksel,[2,1]),'Distribution','binomial','Link','identity','Options',glmopt);
    spk_out=spk_mdl.Coefficients.Estimate;
    skip=checkWarning();
else
    spk_out=zeros(1,11);
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
