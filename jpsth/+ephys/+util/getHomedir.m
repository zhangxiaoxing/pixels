function homedir=getHomedir(opt)
arguments
    opt.dtype (1,:) char {mustBeMember(opt.dtype,{'neupix','AIOPTO','MY'})}='neupix'
end
if strcmp(opt.dtype,'neupix') ||  strcmp(opt.dtype,'MY')
    if ispc
        if endsWith(pwd(),[filesep(),'jpsth'])
            homedir = fullfile('..','..','neupix','SPKINFO');
            [avail,attstr]=fileattrib(homedir);
            if avail
                homedir=attstr.Name;
            else
                keyboard();
            end
        else
            keyboard();
        end
    elseif isunix
        homedir = fullfile('..','..','neupix','SPKINFO');
    end
else
    if ispc
        homedir = fullfile('K:','neupix','AIOPTO','RECDATA');
    elseif isunix
        homedir = fullfile('~','neupix','AIOPTO','RECDATA');
    end
end
