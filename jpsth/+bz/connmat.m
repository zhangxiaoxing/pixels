% function connmat


homedir = fullfile('K:','code','per_sec');
reg_tree=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/reg_tree'));
trial_counts=h5read(fullfile(homedir,'transient_6.hdf5'),'/trial_counts');
trial_sel=all(trial_counts>=20,1);
% wf_sel=(wf_good>0)';
reg_sel=strcmp(reg_tree(1,:),'BS') | strcmp(reg_tree(1,:),'CH');
addpath('K:\code\per_sec\');
[sust_sel,trans_sel]=basic_stats.get_selective();
subsel=trial_sel & reg_sel & (sust_sel | trans_sel);

for depth=3:5
    per_dep_a=bz.net_at_depth(sig,pair,depth,'type','all','subsel',subsel,'overwrite',true);
end

for depth=3:5
    per_dep_m=bz.net_at_depth(sig,pair,depth,'type','memory','subsel',subsel);
    fh=bz.plot_conn_mat(per_dep_m,depth,'Memory');
    
%     exportgraphics(fh,sprintf('Brainwide_memory_%d.pdf',depth))
end

for depth=4
    per_dep_n=bz.net_at_depth(sig,pair,depth,'type','nonmem','subsel',subsel);
    fh=bz.plot_conn_mat(per_dep_n,depth,'Nonmemory');
%     exportgraphics(fh,sprintf('Brainwide_nonmem_%d.pdf',depth))
end
for depth=3:5
    per_dep_m=bz.net_at_depth(sig,pair,depth,'type','memory','subsel',subsel);
    per_dep_n=bz.net_at_depth(sig,pair,depth,'type','nonmem','subsel',subsel);
    per_dep_a=bz.net_at_depth(sig,pair,depth,'type','all','subsel',subsel);
    bz.gehpi.bz2gephi(per_dep_m,per_dep_n,per_dep_a,depth);
end

% end