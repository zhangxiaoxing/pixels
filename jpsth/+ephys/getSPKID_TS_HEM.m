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

if strcmp(laserType,'off')
    tssel=true(size(spkTS));
    for i=1:size(trialInfo,1)
        if trialInfo(i,9)>0
            tssel(spkTS>=(trialInfo(i,10)-30000) &...
                spkTS<(trialInfo(i,10)+300000))=false;
        end
    end
elseif strcmp(laserType,'on')
    tssel=false(size(spkTS));
    for i=1:size(trialInfo,1)
        if trialInfo(i,9)>0
            tssel(spkTS>=(trialInfo(i,10)+15000) &...
                spkTS<(trialInfo(i,10)+300000))=true;
        end
    end
end
out=spkTS(tssel);
end