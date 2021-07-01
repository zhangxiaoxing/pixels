function out=blame()
out=struct();
out.datetime=char(datetime('now','TimeZone','Asia/Shanghai'));
if ispc
[~, out.git_hash] = system('"C:/Program Files/Git/mingw64/bin/git"  rev-parse --verify HEAD'); 
[~, out.git_status] = system('"C:/Program Files/Git/mingw64/bin/git" --no-pager diff --no-color'); 
elseif isunix
    [~, out.git_hash] = system('git rev-parse --verify HEAD')
    [~, out.git_status] = system('git --no-pager diff --no-color')
end
