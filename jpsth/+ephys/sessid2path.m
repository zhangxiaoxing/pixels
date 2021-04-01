function out=sessid2path(sessid)
arguments
    sessid (1,1) double {mustBeInteger,mustBePositive}
end
persistent map

if isempty(map)
    homedir=ephys.util.getHomedir();
    allpath=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/path'));
    upath=unique(allpath);
    [upath,uidx]=sort(upath);
    map=containers.Map('KeyType','int32','ValueType','char');
    for i=1:numel(uidx)
        map(uidx(i))=upath{i};
    end
end

out=map(sessid);

end
