homedir=fullfile('K:','code','per_sec');
trial_counts=h5read(fullfile(homedir,'transient_6.hdf5'),'/trial_counts');
wrs_p=h5read(fullfile(homedir,'transient_6.hdf5'),'/wrs_p');

allpath=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/path'));
allcid=h5read(fullfile(homedir,'transient_6.hdf5'),'/cluster_id');
reg_tree=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/reg_tree'));
reg_tree_depth=h5read(fullfile(homedir,'transient_6.hdf5'),'/reg_tree_depth');
wf_good=h5read(fullfile(homedir,'transient_6.hdf5'),'/wf_good');

[B,BG]=groupcounts(allpath);

wf_sel=(wf_good~=0)';
reg_sel=strcmp(reg_tree(1,:),'BS') | strcmp(reg_tree(1,:),'CH');
trial_sel=all(trial_counts>20,1);

goodsu=(wf_sel & reg_sel & trial_sel);

sust_sel=all(wrs_p(5:10,:)<0.05,1);
trans_sel=any(wrs_p(5:10,:)<0.05,1);


disp([nnz(goodsu),nnz(sust_sel & goodsu),nnz(trans_sel & goodsu)]);

ureg_3=unique(reg_tree(3,reg_sel & trial_sel));
dep_stat=[];
for i=1:length(ureg_3)
    if isempty(ureg_3{i})
        continue
    end
    dep_sel=strcmp(reg_tree(3,:),ureg_3{i});
    dep_stat=cat(1,...
        dep_stat,...
        [nnz(trial_sel & dep_sel),nnz(sust_sel & trial_sel & dep_sel),nnz(trans_sel & trial_sel & dep_sel)]);
end

dep_stat(:,4)=dep_stat(:,2)./dep_stat(:,1);
dep_stat(:,5)=dep_stat(:,3)./dep_stat(:,1);
[ureg_3(2:end)',  num2cell(dep_stat)]

ureg_4=unique(reg_tree(4,reg_sel & trial_sel));
dep_stat=[];
for i=1:length(ureg_4)
    if isempty(ureg_4{i})
        continue
    end
    dep_sel=strcmp(reg_tree(4,:),ureg_4{i});
    dep_stat=cat(1,...
        dep_stat,...
        [nnz(trial_sel & dep_sel),nnz(sust_sel & trial_sel & dep_sel),nnz(trans_sel & trial_sel & dep_sel)]);
end

dep_stat(:,4)=dep_stat(:,2)./dep_stat(:,1);
dep_stat(:,5)=dep_stat(:,3)./dep_stat(:,1);
[ureg_4(2:end)',  num2cell(dep_stat)]