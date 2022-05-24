function [tspre,tspost]=getSPKID_TS_HEM(sessid,preid,postid,opt)
arguments
    sessid (1,1) double {mustBeInteger,mustBePositive}
    preid (1,1) double {mustBeInteger,mustBePositive}
    postid (1,1) double {mustBeInteger,mustBePositive}
    opt.laser (1,:) char {mustBeMember(opt.laser,{'on','off','any'})} = 'any'
end

persistent fstr sessid_
if isempty(fstr) || sessid_~=sessid
    homedir=ephys.util.getHomedir('type','raw','dtype','AIOPTO');
    fpath=ephys.sessid2path(sessid,'type','AIOPTO');
    fstr=load(fullfile(homedir,fpath,'FT_SPIKE.mat'));
    fstr.FT_SPIKE.labelid=cellfun(@(x) str2double(x),fstr.FT_SPIKE.label);
    sessid_=sessid;
end
switch opt.laser
    case 'any'
        tspre=fstr.FT_SPIKE.timestamp{find(fstr.FT_SPIKE.labelid==preid,1)}.';
        tspost=fstr.FT_SPIKE.timestamp{find(fstr.FT_SPIKE.labelid==postid,1)}.';
    case 'on'
        tspre=laserProcess(fstr.FT_SPIKE.timestamp{find(fstr.FT_SPIKE.labelid==preid,1)}.',fstr.FT_SPIKE.trialinfo,'on');
        tspost=laserProcess(fstr.FT_SPIKE.timestamp{find(fstr.FT_SPIKE.labelid==postid,1)}.',fstr.FT_SPIKE.trialinfo,'on');
    case 'off'
        tspre=laserProcess(fstr.FT_SPIKE.timestamp{find(fstr.FT_SPIKE.labelid==preid,1)}.',fstr.FT_SPIKE.trialinfo,'off');
        tspost=laserProcess(fstr.FT_SPIKE.timestamp{find(fstr.FT_SPIKE.labelid==postid,1)}.',fstr.FT_SPIKE.trialinfo,'off');
end
end

function out=laserProcess(spkTS,trialInfo,laserType)
tssel=false(size(spkTS));
if strcmp(laserType,'off')
    t=trialInfo(trialInfo(:,9)<0,:);
else
    t=trialInfo(trialInfo(:,9)>0,:);
end
for i=1:size(t,1)
    tssel(spkTS>=(trialInfo(i,2)-45000) &...
        spkTS<trialInfo(i,2))=true;
end
out=spkTS(tssel);
end