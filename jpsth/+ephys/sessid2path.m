function [out,homedir]=sessid2path(sessid,opt)
arguments
    sessid (1,1) double {mustBeInteger,mustBePositive}
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO'})}='neupix'
end
persistent map
homedir=ephys.util.getHomedir('dtype',opt.type);
if isempty(map)
    
    if strcmp(opt.type,'neupix')
        allpath=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/path'));
    else
        fullpath=deblank(h5read(fullfile(homedir,'Selectivity_AIopto_0419.hdf5'),'/path'));
        allpath=regexp(fullpath,'.*(?=\\.*)','match','once');
    end
    upath=unique(allpath);
    [upath,uidx]=sort(upath);
    map=containers.Map('KeyType','int32','ValueType','char');
    for i=1:numel(uidx)
        map(uidx(i))=upath{i};
    end
end

out=map(sessid);

end
