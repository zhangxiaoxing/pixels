function out=ucid2memtype(id,opt)
arguments
    id (1,1) int32
    opt.subsel (1,1) logical = false
end

persistent ucid2memtype
if isempty(ucid2memtype)
    homedir = fullfile('K:','code','per_sec');
    reg_tree=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/reg_tree'));
    trial_counts=h5read(fullfile(homedir,'transient_6.hdf5'),'/trial_counts');
    cid=h5read(fullfile(homedir,'transient_6.hdf5'),'/cluster_id');
    allpath=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/path'));
    mem_type=h5read(fullfile(homedir,'transient_6.hdf5'),'/mem_type');

    trial_sel=all(trial_counts>=20,1);
    % wf_sel=(wf_good>0)';
    reg_sel=strcmp(reg_tree(1,:),'BS') | strcmp(reg_tree(1,:),'CH');
    subsel=trial_sel & reg_sel;
    if opt.subsel
        cid=int32(cid(subsel));
        allpath=allpath(subsel);
        mem_type=mem_type(subsel)';
    end
    sess=arrayfun(@(x) ephys.path2sessid(allpath{x}),1:numel(allpath));
    ucid=sess'.*100000+int32(cid);

    ucid2memtype=containers.Map(ucid,mem_type);
end
out=ucid2memtype(id);
end