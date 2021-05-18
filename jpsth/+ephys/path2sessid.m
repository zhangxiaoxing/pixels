function out=path2sessid(path,opt)
arguments
    path (1,:) char
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO'})}='neupix'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
end
persistent map

if isempty(map)
    homedir=ephys.util.getHomedir('dtype',opt.type);
    if strcmp(opt.type,'neupix')
        if strcmp(opt.criteria,'Learning')
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
end
out=map(path);

end
