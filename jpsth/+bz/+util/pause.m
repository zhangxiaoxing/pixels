function pause(i,fpath)
if exist(fpath,'file')
    disp('paused by file')
    disp(i)
    keyboard()
end