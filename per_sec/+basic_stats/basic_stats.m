homedir= fullfile('K:','code','per_sec');
trial_counts=h5read(fullfile(homedir,'transient_6.hdf5'),'/trial_counts');
% allpath=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/path'));
% allcid=h5read(fullfile(homedir,'transient_6.hdf5'),'/cluster_id');
reg_tree=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/reg_tree'));
% reg_tree_depth=h5read(fullfile(homedir,'transient_6.hdf5'),'/reg_tree_depth');
wf_good=h5read(fullfile(homedir,'transient_6.hdf5'),'/wf_good');

wf_sel=(wf_good>0)';
trial_sel=all(trial_counts>=20,1);
reg_sel=strcmp(reg_tree(1,:),'BS') | strcmp(reg_tree(1,:),'CH');
for dd=3
    [reg,count]=basic_stats.ratio_tree_depth(reg_tree,dd,'subsel',reg_sel & trial_sel);
    basic_stats.plotbars(reg,count,dd)
end