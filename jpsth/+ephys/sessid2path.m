function [out,homedir,map_]=sessid2path(sessid,opt)
arguments
    sessid (1,1) double {mustBeInteger}
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO'})}='neupix'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
end
homedir='';
persistent opt_ map
if isempty(map) || ~isequaln(opt,opt_)
    if strcmp(opt.type,'neupix')
        su_meta=ephys.util.load_meta("save_file",false,"adjust_white_matter",true,'criteria',opt.criteria);
        map=containers.Map(num2cell(su_meta.sess),su_meta.allpath);
    else
        error("Unfinished");
    end
    opt_=opt;
end
if map.isKey(sessid)
    out=map(sessid);
else
    out=[];
end
map_=map;
end
