function out=get_TCOM_profile(sig_bin,pair_bin,opt)

arguments
    sig_bin (:,2) double
    pair_bin (:,2) double
    opt.edge (1,:) double = [-24,-12:2:12,24]
end

pair_diff=diff(pair_bin,1,2);
sig_diff=diff(sig_bin,1,2);
pair_hist=histcounts(pair_diff,opt.edge);
sig_hist=histcounts(sig_diff,opt.edge);
out=[sig_hist;pair_hist;opt.edge(1:end-1)+diff(opt.edge,1,2)./2];

