function [out,homedir,map_]=sessid2path(sessid,opt)
arguments
    sessid (1,1) double {mustBeInteger}
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO'})}='neupix'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
end
persistent map type_ criteria_
homedir=ephys.util.getHomedir('dtype',opt.type);
if isempty(map) || ~strcmp(opt.type,type_) || ~strcmp(opt.criteria,criteria_)
    
    if strcmp(opt.type,'neupix')
        if strcmp(opt.criteria,'WT')
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
    map=containers.Map('KeyType','int32','ValueType','char');
    for i=1:numel(uidx)
        map(uidx(i))=upath{i};
    end
    type_=opt.type;
    criteria_=opt.criteria;
end
if map.isKey(sessid)
    out=map(sessid);
else
    out=[];
end
map_=map;
end
