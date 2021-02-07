debug=false;
bz.util.dependency
prefix='0203';
dualprobe=dir(fullfile(homedir,'**','spike_info.mat'));
singleprobe=dir(fullfile(homedir,'DataSum','singleProbe','**','spike_times.npy'));

error_list=cell(0);
% for i=102%11:numel(dualprobe)
    bz.util.pause(i,'xcorrpause');
    folder=dualprobe(i).folder;
    if isempty(bz.util.procPerf(h5read(fullfile(folder,'events.hdf5'),'/trials')',''))
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
    
    suids=bz.util.goodCid(dualprobe(i).folder);
    if debug
        suids=suids(1:20);
    end
    susel=ismember(spkID,suids);
    spkID=double(spkID(susel));
    spkTS=double(spkTS(susel));
    
    
    mono=bz.sortSpikeIDz(spkTS,spkID);
    
    if debug && false
        keyboard();
        for j=1:size(mono.sig_con,1)
            fh=figure('Position',[32,32,320,240]);
            plot(mono.ccgR(:,mono.sig_con(j,1),mono.sig_con(j,2)),'r-');
            xlim([1,501])
            arrayfun(@(x) xline(x,'b:'),[251-25,251,251+25])
            set(gca(),'XTick',[1,251-25,251,251+25,501],'XTickLabel',[-100,-10,0,10,100])
            title(j)
            close(fh)
        end
        
    end
    save(sprintf('%s_BZ_XCORR_duo_f%d.mat',prefix,i),'mono','-v7.3','folder')
% end
if isunix
    quit(0)
end

