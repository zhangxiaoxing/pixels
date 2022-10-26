function dependency(opt)
arguments
    opt.buz (1,1) logical = true;
    opt.ft (1,1) logical = true;
end
if ispc
    % dPCA
    %addpath('..\..\Lib\dPCA\matlab\')


    if opt.buz
        addpath('K:\Lib\buzcode\buzcode\io')
        addpath('K:\Lib\buzcode\buzcode\utilities\')
        addpath('K:\Lib\buzcode\buzcode\analysis\spikes\correlation\')
        addpath('K:\Lib\buzcode\buzcode\analysis\spikes\functionalConnectionIdentification\')
        addpath('K:\Lib\buzcode\buzcode\visualization\')
        addpath('K:\Lib\buzcode\buzcode\externalPackages\FMAToolbox\General\');
        addpath('K:\Lib\buzcode\buzcode\externalPackages\FMAToolbox\Helpers\');
    end
    if opt.ft
        if exist(fullfile('K:','Lib','fieldtrip-20200320'),"dir")
%         addpath(fullfile('K:','Lib','npy-matlab-master','npy-matlab'))
            addpath(fullfile('K:','Lib','fieldtrip-20200320'))
        elseif exist(fullfile('..','Lib','fieldtrip-20221022'),"dir")
            addpath(fullfile('..','Lib','fieldtrip-20221022'))
        end
        ft_defaults
    end
else
    if opt.buz
        addpath('~/Lib/buzcode/buzcode/io')
        addpath('~/Lib/buzcode/buzcode/utilities/')
        addpath('~/Lib/buzcode/buzcode/analysis/spikes/correlation/')
        addpath('~/Lib/buzcode/buzcode/analysis/spikes/functionalConnectionIdentification/')
        addpath('~/Lib/buzcode/buzcode/visualization/')
        addpath('~/Lib/buzcode/buzcode/externalPackages/FMAToolbox/General/');
        addpath('~/Lib/buzcode/buzcode/externalPackages/FMAToolbox/Helpers/');
    end
    if opt.ft
        addpath(fullfile('~','Lib','npy-matlab-master','npy-matlab'))
        addpath(fullfile('~','Lib','fieldtrip-20200320'))
        ft_defaults
    end
end
