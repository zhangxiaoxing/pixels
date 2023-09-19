function mono_rez=sortSpikeIDz(spkTS,spkID,opt)
arguments
    spkTS (:,1) double {mustBeNonempty,mustBeNonmissing}
    spkID (:,1) double {mustBeNonempty,mustBeNonmissing}
    opt.negccg (1,1) logical = false
    opt.wallclocktime (1,1) logical = false
end

spkid(:,2)=spkID;
spkid(:,1)=0;
[~,~,UC]=unique(spkid(:,2));
spkid(:,3)=UC;

if opt.wallclocktime
    mono_rez = bz.zx_MonoSynConvClick (spkid,spkTS,'alpha',0.05,'negccg',opt.negccg);
else
    spkts=spkTS./30000;
    mono_rez = bz.zx_MonoSynConvClick (spkid,spkts,'alpha',0.05,'negccg',opt.negccg);
end
end