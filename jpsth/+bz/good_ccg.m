function out=good_ccg(one_ccg)
%adapted from jpsth\+ephys\+waveform\goodWaveform.m
arguments
    one_ccg (501,1) double {mustBeNonempty}
end
%criteria 1
out=ones(1,6)*-1;
one_ccg=one_ccg-mean(one_ccg(1:50));

[mmin,~]=min(one_ccg);
[mmax,xi]=max(one_ccg);
% [mmax,xi]=max(one_ccg((mi+1):end));
if mmax>-mmin
    out(1)=1;
end
out(2)=xi;
%criteria 2
if max(one_ccg)>0
    one_ccg=one_ccg./max(one_ccg);
else
    one_ccg=one_ccg-min(one_ccg);
    one_ccg=one_ccg./max(one_ccg);
end
[lc_pk,~]=findpeaks(one_ccg,'MinPeakHeight',0.25,'MinPeakDistance',12);
out(3)=numel(lc_pk);

%fwhm
lcross=find(one_ccg(226:xi)>0.5,1)+225;
rcross=find(one_ccg(xi:end)<0.5,1)+xi;
if numel(lcross)==1 && numel(rcross)==1
    out(4)=rcross-lcross;
    out(5)=lcross;
    out(6)=rcross;
end

end
