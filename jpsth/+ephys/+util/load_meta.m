function out=load_meta()
persistent meta_str
if isempty(meta_str)
    homedir=ephys.util.getHomedir();
    meta_str.trial_counts=h5read(fullfile(homedir,'transient_6.hdf5'),'/trial_counts');
    meta_str.wrs_p=h5read(fullfile(homedir,'transient_6.hdf5'),'/wrs_p');
    meta_str.selec=h5read(fullfile(homedir,'transient_6.hdf5'),'/selectivity');
    meta_str.allpath=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/path'));
    meta_str.allcid=h5read(fullfile(homedir,'transient_6.hdf5'),'/cluster_id');
    meta_str.reg_tree=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/reg_tree'));
    meta_str.mem_type=h5read(fullfile(homedir,'transient_6.hdf5'),'/mem_type');
end
out=meta_str;
end