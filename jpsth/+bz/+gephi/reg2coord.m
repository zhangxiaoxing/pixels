function out=reg2coord(opt)
arguments
    opt.write_file (1,1) logical = false
end
persistent reg_coord_map opt_
if isempty(reg_coord_map) || ~isequaln(opt,opt_)
    homedir=fullfile('..','per_sec');
    reg_tree=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/reg_tree'));
    load('K:\neupix\track_meta\sucoords318.mat');
    reg_coord_map=containers.Map('KeyType','char','ValueType','any');
    idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));
    csvcell={'Id','Label','AP','DV','ML','depth','L1','L2','L3','L4','L5','L6'};
    for dep=1:6
        ureg=unique(reg_tree(dep,:));
        ureg=ureg(~cellfun('isempty',ureg));
        for i=1:numel(ureg)
            rsel=strcmp(reg_tree(dep,:),ureg{i});
            tree=cell(1,6);
            branch=idmap.reg2tree(ureg{i});
            branch=branch(3:end);
            tree(1:dep)=branch;
            coordmm=nanmean(coord(rsel,:),1);
            reg_coord_map(ureg{i})=coordmm;
            csvcell(end+1,:)=[{char(ureg{i}),...
                char(ureg{i}),...
                coordmm(1),...
                1000-coordmm(2),...
                coordmm(3),...
                dep},...
                tree];
        end
    end

    if opt.write_file
        save(fullfile('K:','neupix','track_meta','reg2coord.mat'),'reg_coord_map');
        writecell(csvcell,fullfile('..','..','neupix','track_meta','reg2coord.csv'));
        writecell(csvcell,fullfile('bzdata','reg2coord.csv'));
    end
end
out=reg_coord_map;