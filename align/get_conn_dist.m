function [s_dist,p_dist]=get_conn_dist(sig,pair)
arguments
    sig (1,1) struct
    pair (1,1) struct
end

persistent s_dist_ p_dist_

if isempty(s_dist_)
    disp('Calculating distance');
    addpath('k:\code\align\')
    s_dist_=reg_tree_dist(sig.reg);
    p_dist_=reg_tree_dist(pair.reg);
end

s_dist=s_dist_;
p_dist=p_dist_;
end