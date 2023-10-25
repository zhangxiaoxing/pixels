function rec_dur=sessid2length(sessidx,opt)
arguments
    sessidx
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
end
persistent sessidx_ rec_dur_
if isempty(sessidx_) || sessidx~=sessidx_
    metaf=dir(fullfile(...
        ephys.util.getHomedir(),...
        replace(ephys.sessid2path(sessidx,'criteria',opt.criteria),'\',filesep()),...
        "*.ap.meta"));
    ts=textscan(...
        fopen(fullfile(metaf.folder,metaf.name)),...
        '%s','Delimiter',{'\n'});
    rec_dur=str2double(replace(ts{1}{startsWith(ts{1},'fileSizeBytes')},'fileSizeBytes=',''))...
        ./385/2;
    rec_dur_=rec_dur;
    sessidx_=sessidx;
    fclose('all');
else
    rec_dur=rec_dur_;
end
end

