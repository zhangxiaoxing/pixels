function [out,map_]=path2sessid(path,opt)
arguments
    path (1,:) char
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO','MY'})}='neupix'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
end
persistent map type criteria

if isempty(map) || ~strcmp(opt.type,type) || ~strcmp(opt.criteria,criteria)
    homedir=ephys.util.getHomedir('dtype',opt.type);
    if strcmp(opt.type,'neupix') || strcmp(opt.type,'MY')
        if ~strcmp(opt.criteria,'Learning')
            allpath=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/path'));
        else
            allpath=deblank(h5read(fullfile(homedir,'transient_6_complete.hdf5'),'/path'));
        end
        else
        fullpath=deblank(h5read(fullfile(homedir,'Selectivity_AIopto_0419.hdf5'),'/path'));
        allpath=regexp(fullpath,'.*(?=\\.*)','match','once');
    end
    upath=unique(allpath);
    [upath,uidx]=sort(upath);
    map=containers.Map('KeyType','char','ValueType','int32');
    for i=1:numel(uidx)
        map(upath{i})=uidx(i);
    end
    type=opt.type;
    criteria=opt.criteria;
    
end
p=replace(path,'/home/zx/neupix/SPKINFO/','');
p=replace(p,'/','\');
out=map(p);
map_=map;
end
