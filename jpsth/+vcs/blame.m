function out=blame()
out=struct();
out.datetime=char(datetime('now','TimeZone','Asia/Shanghai'));
if ispc
    if exist("C:/Program Files/Git/mingw64/bin/git.exe","file")
        [~, out.git_hash] = system('"C:/Program Files/Git/mingw64/bin/git"  rev-parse --verify HEAD');
        [~, out.git_status] = system('"C:/Program Files/Git/mingw64/bin/git" --no-pager diff --no-color');
    elseif exist("C:/msys64/usr/bin/git.exe","file")
        [~, out.git_hash] = system('"C:/msys64/usr/bin/git.exe"  rev-parse --verify HEAD');
        [~, out.git_status] = system('"C:/msys64/usr/bin/git.exe" --no-pager diff --no-color');    
    end
elseif isunix
    [~, out.git_hash] = system('git rev-parse --verify HEAD');
    [~, out.git_status] = system('git --no-pager diff --no-color');
end
