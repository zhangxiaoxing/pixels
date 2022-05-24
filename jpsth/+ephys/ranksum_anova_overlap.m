anovameta=wave.get_dur_waveid();
asel=any(anovameta.anovap(:,[1,4,5,7])<0.05,2);
meta=ephys.util.load_meta();
greysel=ismember(meta.reg_tree(1,:),{'CH','BS'}).';

d6sel=any(meta.fdr_6(1:4,:)<0.05).';
d3sel=any(meta.fdr_3(1:4,:)<0.05).';

% d6sel=any(meta.wrs_p_6(4:7,:)<0.05).';
% d3sel=any(meta.wrs_p_3(4:7,:)<0.05).';

overlap_rate=nnz(asel & (d3sel | d6sel))/nnz((d3sel|d6sel))
grey_overlap_rate=nnz(asel & (d3sel | d6sel) & greysel)./nnz((d3sel|d6sel) & greysel)


waveid=ephys.get_wave_id(meta.sess,meta.allcid);
anovameta=ephys.selectivity_anova();
% dur_waveid_=zeros(size(anovameta.sess));
dur_all=any(anovameta.anovap(:,[2 4 6 7])<0.05,2);

nnz(waveid>0 & waveid<5)
nnz(dur_all)

nnz(waveid>0 & waveid<5 & dur_all)
