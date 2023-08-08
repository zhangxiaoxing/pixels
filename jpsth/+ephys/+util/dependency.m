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
end

if opt.ft
    if exist(fullfile('..','..','Lib','fieldtrip'),"dir")
        addpath(fullfile('..','..','Lib','fieldtrip'))
    elseif exist(fullfile('..','..','..','Lib','fieldtrip'),"dir")
        addpath(fullfile('..','..','..','Lib','fieldtrip'))
    end
    ft_defaults
end
