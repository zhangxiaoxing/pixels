oo=[];
for nn=reshape(find(wrs_mux_meta.wave_id==5),1,[])
    bins=wrs_mux_meta.p_olf(nn,:)<0.05;
    oo=[oo,reshape(wrs_mux_meta.class_fr(nn,1:2,bins),1,[])];
end
for nn=reshape(find(wrs_mux_meta.wave_id==6),1,[])
    bins=wrs_mux_meta.p_olf(nn,:)<0.05;
    oo=[oo,reshape(wrs_mux_meta.class_fr(nn,3:4,bins),1,[])];
end

dd=[];
for nn=reshape(find(wrs_mux_meta.wave_id==7),1,[])
    bins=wrs_mux_meta.p_dur(nn,:)<0.05;
    dd=[dd,reshape(wrs_mux_meta.class_fr(nn,[1 3],bins),1,[])];
end
for nn=reshape(find(wrs_mux_meta.wave_id==8),1,[])
    bins=wrs_mux_meta.p_dur(nn,:)<0.05;
    dd=[dd,reshape(wrs_mux_meta.class_fr(nn,[2 4],bins),1,[])];
end

mm=[];
for classs=1:4
for nn=reshape(find(wrs_mux_meta.wave_id==classs),1,[])
    bins=wrs_mux_meta.p_mux(nn,:)<0.05 | wrs_mux_meta.p_olf(nn,:)<0.05 | wrs_mux_meta.p_dur(nn,:)<0.05;
    mm=[mm,reshape(wrs_mux_meta.class_fr(nn,classs,bins),1,[])];
end
end

persess=[];
for ii=1:116
    disp(ii);
    [spkID,spkTS,~,~,~,FT_SPK]=ephys.getSPKID_TS(ii,'keep_trial',true);
    persess=[persess,cellfun(@(x) nnz(x<-0.5 & x>=-1.5),FT_SPK.time)./size(FT_SPK.trialinfo,1)];
end

mean([oo,dd,mm])
std([oo,dd,mm])./sqrt(numel([oo,dd,mm]))

mean(persess)
std(persess)./sqrt(numel(persess))