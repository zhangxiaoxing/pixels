function [spk_out,fc_out]=history_coeff(sessid,suid)
arguments
    sessid (1,1) int32
    suid (1,2) int32
end
persistent bitmask X
if isempty(bitmask) || isempty(X)
    bitmask=2.^(0:9)'; %TODO constant to persistent?
    %TODO spinoff model matrix
    X=zeros(1024,10);
    for row=1:1024
        for col=1:10
            if bitand(row-1,2.^(col-1))
                X(row,col)=1;
            end
        end
    end
end

[spkID,spkTS,~,~,~]=ephys.getSPKID_TS(sessid);
tspre=spkTS(spkID==suid(1));
tspost=spkTS(spkID==suid(2));
tmax=max([tspre;tspost]);
histpre=histcounts(tspre,0:600:tmax)>0;
histpost=histcounts(tspost,0:600:tmax)>0;
binpre=discretize(tspre,0:600:tmax);
binpost=discretize(tspost,0:600:tmax);
post_spike_prob=zeros(1024,2);
fc_prob=zeros(1024,2);
for i=1:(length(histpre)-10)
    hist_type=histpre(i:i+9)*bitmask;
    post_spike_prob(hist_type+1,1)=post_spike_prob(hist_type+1,1)+1;
    post_spike_prob(hist_type+1,2)=post_spike_prob(hist_type+1,2)+histpost(i+9); % blind detect post spike
    presel=binpre==i+10;
    postsel=binpost==i+10;
    if nnz(presel)>0
        if nnz(postsel)>0
            [fchit,fcmiss]=getfc(tspre(presel),tspost(postsel));
        else
            fchit=0;fcmiss=nnz(presel);
        end
        fc_prob(hist_type+1,1)=fc_prob(hist_type+1,1)+fchit+fcmiss;
        fc_prob(hist_type+1,2)=fc_prob(hist_type+1,2)+fchit;
    end
end

spksel=post_spike_prob(:,1)>0;
fcsel=fc_prob(:,1)>0;
spk_mdl=fitglm(X(spksel,:),post_spike_prob(spksel,[2,1]),'Distribution','binomial','Link','identity');
fc_mdl=fitglm(X(fcsel,:),fc_prob(fcsel,[2,1]),'Distribution','binomial','Link','identity');
spk_out=spk_mdl.Coefficients.Estimate;
fc_out=fc_mdl.Coefficients.Estimate;
end

function [hit,miss]=getfc(pre,post)
lag=post-(pre.');
hit=nnz(any(lag>15 & lag<300,1)); %assuming 30K sps.
miss=numel(pre)-hit;
end