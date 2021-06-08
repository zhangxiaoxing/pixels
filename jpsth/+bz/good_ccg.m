function out=good_ccg(one_ccg)
%adapted from jpsth\+ephys\+waveform\goodWaveform.m
arguments
    one_ccg (501,1) double {mustBeNonempty}
end
%criteria 1
out=[-1,-1,-1,-1];
one_ccg=one_ccg-mean(one_ccg(1:50));

[mmin,~]=min(one_ccg);
[mmax,xi]=max(one_ccg);
% [mmax,xi]=max(one_ccg((mi+1):end));
if mmax>-mmin
    out(1)=1;
    out(2)=xi;
else
    return
end
%criteria 2
[lc_pk,~]=findpeaks(one_ccg,'MinPeakProminence',0.2*mmax);
out(3)=numel(lc_pk);
scale=max(abs(one_ccg));
one_ccg=one_ccg./scale;

%fwhm
lcross=find(one_ccg>0.5,1);
rcross=find(one_ccg(xi:end)<0.5,1)+xi;
if numel(lcross)==1 && numel(rcross)==1
    out(4)=rcross-lcross;
end
end
