fpath='K:\LiukaiData\20230904\01-wam-8-2seq-ag-20230904_124235-1.mountainsort';
fl=dir(fullfile(fpath,'*_offset.mat'));
spkts=[];
spkid=[];
for fi=1:numel(fl)
    tfstr=load(fullfile(fl(fi).folder,fl(fi).name));
    susel=arrayfun(@(x) ~isempty(tfstr.firings_offset(x).unit_indx),1:numel(tfstr.firings_offset));
    if ~any(susel)
        continue
    end
    for leadii=reshape(find(susel),1,[])
        nsu=numel(tfstr.firings_offset(leadii).unit_indx);
        for uii=1:nsu
            ts=tfstr.firings_offset(leadii).ts{uii}.';
            spkts=[spkts;double(ts)];
            id=tfstr.firings_offset(leadii).unit_indx(uii)+1000*leadii+1000000*fi;
            spkid=[spkid;double(repmat(id,numel(ts),1))];
        end
    end
end
ephys.util.dependency()
rez=bz.sortSpikeIDz(spkts,spkid);

figure()
tiledlayout('flow')
for ii=1:size(rez.sig_con,1)
    nexttile
    plot(rez.ccgR(:,rez.sig_con(ii,1),rez.sig_con(ii,2)))
end
