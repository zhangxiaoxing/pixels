if ~exist('mono','var')
    disp('Missing mono variable')
    return
end
% cids=mono.sig_con(sigidx,:);
if ~exist('folder','var')
    disp('Missing folder variable')
    return
end
if startsWith(folder,'/home') && ispc
    folder=replace(replace(folder,'/','\'),'\home\zx','K:');
end