function [sust,trans]=get_selective(opt)
%TODO optional delay duration
arguments
    opt.homedir (1,:) char = fullfile('K:','code','per_sec');
end
persistent sust_ trans_
if isempty(sust_)
    wrs_p=h5read(fullfile(opt.homedir,'transient_6.hdf5'),'/wrs_p');
    sust_=all(wrs_p(5:10,:)<0.05,1);
    trans_=any(wrs_p(5:10,:)<0.05,1) & ~sust_;
end
sust=sust_;
trans=trans_;
end