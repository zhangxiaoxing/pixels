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

prefix='0315';
if isfile(sprintf('%s_BZ_XCORR_duo_f%d.mat',prefix,i))
    quit(0)
end
disp(i);

bz.util.dependency
sess=dir(fullfile(homedir,'**','spike_info.hdf5'));

error_list=cell(0);
% for i=102%11:numel(sess)
    bz.util.pause(i,'xcorrpause');
    folder=sess(i).folder;
    trials=h5read(fullfile(folder,'FR_All_1000.hdf5'),'/Trials');
    SU_id=h5read(fullfile(folder,'FR_All_1000.hdf5'),'/SU_id');
    if sum(trials(:,9))<40 || numel(SU_id)<2
        if isunix
            quit(0)
        else
            return
        end
    end
    
    cstr=h5info(fullfile(sess(i).folder,sess(i).name));
    spkID=[];spkTS=[];
    for prb=1:size(cstr.Groups,1)
        prbName=cstr.Groups(prb).Name;
        spkID=cat(1,spkID,h5read(fullfile(sess(i).folder,sess(i).name),[prbName,'/clusters']));
        spkTS=cat(1,spkTS,h5read(fullfile(sess(i).folder,sess(i).name),[prbName,'/times']));
    end
    
    susel=ismember(spkID,SU_id);
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

