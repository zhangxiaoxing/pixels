function homedir=getHomedir(opt)
arguments
    opt.type (1,:) char {mustBeMember(opt.type,{'sums','raw'})} = 'sums'
    opt.dtype (1,:) char {mustBeMember(opt.dtype,{'neupix','AIOPTO'})}='neupix'
    
end
if strcmp(opt.dtype,'neupix')
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
else
    if ispc
        if strcmp(opt.type,'sums')
            homedir = fullfile('K:','neupix','AIOPTO','META');
        elseif strcmp(opt.type,'raw')
            homedir = fullfile('K:','neupix','AIOPTO','RECDATA');
        end
    elseif isunix
        if strcmp(opt.type,'sums')
            homedir = fullfile('~','neupix','AIOPTO','META');
        elseif strcmp(opt.type,'raw')
            homedir = fullfile('~','neupix','AIOPTO','RECDATA');
        end
    end
    
end
