if isunix
    addpath(fullfile('/','home','zx','Lib','npy-matlab-master','npy-matlab'))
    addpath(fullfile('/','home','zx','Lib','fieldtrip-20200320'))
    homedir=fullfile('/','home','zx','neupix','wyt','DataSum');
else
    addpath(fullfile('K:','Lib','npy-matlab-master','npy-matlab'))
    addpath(fullfile('K:','Lib','fieldtrip-20200320'))
    homedir=fullfile('K:','neupix','wyt','DataSum');
end
ft_defaults
