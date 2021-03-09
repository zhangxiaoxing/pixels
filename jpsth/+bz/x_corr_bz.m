if ~exist('i','var')
    disp('missing i value')
    if isunix
        quit(0)
    else
        return
    end
end

if ~exist('debug','var')
    debug=false;
end

prefix='0203';
if isfile(sprintf('%s_BZ_XCORR_duo_f%d.mat',prefix,i))
    quit(0)
end
disp(i);

bz.util.dependency
dualprobe=dir(fullfile(homedir,'**','spike_info.mat'));

error_list=cell(0);
% for i=102%11:numel(dualprobe)
    bz.util.pause(i,'xcorrpause');
    folder=dualprobe(i).folder;
    if isempty(behav.procPerf(h5read(fullfile(folder,'events.hdf5'),'/trials')',''))
        if isunix
            quit(0)
        else
            return
        end
    end
%     [avail,spktrial]=bz.pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,currmodel); % posix
    fstr=load(fullfile(dualprobe(i).folder,dualprobe(i).name));
    spkID=[fstr.spike_info{1}{1};fstr.spike_info{1}{2}];
    spkTS=[fstr.spike_info{2}{1};fstr.spike_info{2}{2}];
    
    suids=ephys.goodCid(dualprobe(i).folder);
    if debug
        suids=suids(1:20);
        keyboard()
    end
    susel=ismember(spkID,suids);
    spkID=double(spkID(susel));
    spkTS=double(spkTS(susel));
    
    
    mono=bz.sortSpikeIDz(spkTS,spkID);
    
    if debug && false
        bz.util.plotCCG
    end
    save(sprintf('%s_BZ_XCORR_duo_f%d.mat',prefix,i),'mono','-v7.3','folder')
% end
if isunix
    quit(0)
end

