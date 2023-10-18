function [out,map_]=path2sessid(path,opt)
arguments
    path (1,:) char
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO','MY'})}='neupix'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
end
persistent opt_ map

if isempty(map) || ~isequaln(opt,opt_)
    if strcmp(opt.type,'neupix')
        if strcmp(opt.criteria,'WT')
            su_meta=ephys.util.load_meta("save_file",false,"adjust_white_matter",true);
        elseif 
        else
            error("Unfinished");
        end
        map=containers.Map(su_meta.allpath,su_meta.sess);
    else
        error("Unfinished");
    end
    opt_=opt;
end
if map.isKey(path)
    out=map(path);
else
    out=[];
end
map_=map;
end
