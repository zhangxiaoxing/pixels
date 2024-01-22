function [homedir,homedir_all]=getHomedir(opt)
arguments
    opt.dtype (1,:) char {mustBeMember(opt.dtype,{'neupix','AIOPTO','MY'})}='neupix'
end
if strcmp(opt.dtype,'neupix') ||  strcmp(opt.dtype,'MY')
    if endsWith(pwd(),[filesep(),'jpsth'])
        homedir = fullfile('..','..','neupix','npdata_out');
        dirall = struct2table(dir(fullfile('..','..','neupix','npdata_out*')));
        homedir_all = table2cell(rowfun(@(x,y) fullfile(x,y),dirall(:,["folder","name"])));
        
        [avail,attstr]=fileattrib(homedir);
        if avail
            homedir=attstr.Name;
        else
            error("Error listing data files")
        end
    else
        error("Unexpected working dir")
    end
else
    if ispc
        homedir = fullfile('K:','neupix','AIOPTO','RECDATA');
    elseif isunix
        homedir = fullfile('~','neupix','AIOPTO','RECDATA');
    end
end
