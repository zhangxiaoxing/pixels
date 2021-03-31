function out=path2sessid(path)
arguments
    path (1,:) char
end
persistent map

if isempty(map)
    homedir=fullfile('K:','code','per_sec');
    allpath=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/path'));
    upath=unique(allpath);
    [upath,uidx]=sort(upath);
    map=containers.Map('KeyType','char','ValueType','int32');
    for i=1:numel(uidx)
        map(upath{i})=uidx(i);
    end
end

out=map(path);

end