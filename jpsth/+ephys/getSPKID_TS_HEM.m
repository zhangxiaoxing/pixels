function [tspre,tspost]=getSPKID_TS_HEM(sessid,preid,postid)

persistent fstr sessid_
if isempty(fstr) || sessid_~=sessid
    homedir=ephys.util.getHomedir('type','raw','dtype','AIOPTO');
    fpath=ephys.sessid2path(sessid,'type','AIOPTO');
    fstr=load(fullfile(homedir,fpath,'FT_SPIKE.mat'));
    fstr.FT_SPIKE.labelid=cellfun(@(x) str2double(x),fstr.FT_SPIKE.label);
    sessid_=sessid;
end
tspre=fstr.FT_SPIKE.timestamp{find(fstr.FT_SPIKE.labelid==preid,1)}.';
tspost=fstr.FT_SPIKE.timestamp{find(fstr.FT_SPIKE.labelid==postid,1)}.';
end