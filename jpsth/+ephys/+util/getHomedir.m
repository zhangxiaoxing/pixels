function homedir=getHomedir(opt)
arguments
    opt.type (1,:) char {mustBeMember(opt.type,{'sums','raw'})} = 'sums'
end

if ispc
    if strcmp(opt.type,'sums')
        homedir = fullfile('K:','code','per_sec');
    elseif strcmp(opt.type,'raw')
        homedir = fullfile('K:','neupix','SPKINFO');
    end
elseif isunix
    if strcmp(opt.type,'sums')
        homedir = fullfile('~','pixels','per_sec');
    elseif strcmp(opt.type,'raw')
        homedir = fullfile('~','neupix','SPKINFO');
    end
end
