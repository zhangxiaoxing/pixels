function connmat
homedir = fullfile('K:','code','per_sec');
reg_tree=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/reg_tree'));
trial_counts=h5read(fullfile(homedir,'transient_6.hdf5'),'/trial_counts');
trial_sel=all(trial_counts>=20,1);
% wf_sel=(wf_good>0)';
reg_sel=strcmp(reg_tree(1,:),'BS') | strcmp(reg_tree(1,:),'CH');
subsel=trial_sel & reg_sel;

% for depth=3
%     per_dep=bz.net_at_depth(sig,pair,depth,'type','all','subsel',subsel);
%     bz.plot_conn_mat(per_dep,depth)
% end

for depth=5
    per_dep_m=bz.net_at_depth(sig,pair,depth,'type','memory','subsel',subsel);
    bz.plot_conn_mat(per_dep_m,depth,'Memory')
end

for depth=5
    per_dep_n=bz.net_at_depth(sig,pair,depth,'type','nonmem','subsel',subsel);
    bz.plot_conn_mat(per_dep_n,depth,'Nonmemory')
end

end