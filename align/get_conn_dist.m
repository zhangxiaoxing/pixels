function [s_dist,p_dist]=get_conn_dist(sig,pair)
arguments
    sig (1,1) struct
    pair (1,1) struct
end

persistent s_dist_ p_dist_ sigsize pairsize

if isempty(s_dist_) || sigsize~=numel(sig.sess) || pairsize~=numel(pair.sess)
    disp('Calculating distance');
    addpath('k:\code\align\')
    s_dist_=reg_tree_dist(sig.reg);
    p_dist_=reg_tree_dist(pair.reg);
    sigsize=numel(sig.sess);
    pairsize=numel(pair.sess);
end

s_dist=s_dist_;
p_dist=p_dist_;
end
