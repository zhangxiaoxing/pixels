function mono_rez=sortSpikeIDz(spkTS,spkID)
    spkts=spkTS./30000;
    spkid(:,2)=spkID;
    spkid(:,1)=0;
    [~,~,UC]=unique(spkid(:,2));
    spkid(:,3)=UC;
    mono_rez = bz.zx_MonoSynConvClick (spkid,spkts,'alpha',0.05);
end