function [avail,folder]=loaded(mono,folder,opt)
arguments
    mono struct = []
    folder char = []
    opt.mono (1,1) logical = true
    opt.folder (1,1) logical = true
end
avail=true;
if opt.mono && isempty(mono)
    avail=false;
end

if opt.folder && isempty(folder)
    avail=false;
end
if avail && startsWith(folder,'/home') && ispc
    folder=replace(replace(folder,'/','\'),'\home\zx','K:');
end
end