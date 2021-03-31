homedir=ephys.util.getHomedir();
mem_type=h5read(fullfile(homedir,'transient_6.hdf5'),'/mem_type');
trial_counts=h5read(fullfile(homedir,'transient_6.hdf5'),'/trial_counts');
wf_good=h5read(fullfile(homedir,'transient_6.hdf5'),'/wf_good')';
trial_sel=all(trial_counts>=20,1);
reg_tree=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/reg_tree'));
reg_sel=strcmp(reg_tree(1,:),'BS') | strcmp(reg_tree(1,:),'CH');
path=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/path'));

wrs_p=h5read(fullfile(homedir,'transient_6.hdf5'),'/wrs_p')';
selectivity=h5read(fullfile(homedir,'transient_6.hdf5'),'/selectivity')';
cid=h5read(fullfile(homedir,'transient_6.hdf5'),'/cluster_id')';
typesel=(mem_type==2 | mem_type==4) & all(trial_counts>=20,1) & wf_good & reg_sel;
[GC,GR]=groupcounts(path(typesel));

[~,sessidx]=sort(GC,'descend');
sesspath=GR(sessidx(1));

sesssel=strcmp(path,sesspath) & typesel';

sess_wrs=wrs_p(sesssel,:);
sess_cid=cid(sesssel);
sess_selec=selectivity(sesssel,:);
for bin=5:10
    [~,binidx]=max(abs(sess_selec(:,bin)));
    keyboard
end