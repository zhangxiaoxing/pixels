asel=any(anovameta.anovap(:,[1,4,5,7])<0.05,2);
meta=ephys.util.load_meta();
greysel=ismember(meta.reg_tree(1,:),{'CH','BS'}).';

d6sel=any(meta.fdr_6(1:4,:)<0.05).';
d3sel=any(meta.fdr_3(1:4,:)<0.05).';


overlap_rate=nnz(asel & (d3sel | d6sel))/nnz((d3sel|d6sel))
grey_overlap_rate=nnz(asel & (d3sel | d6sel) & greysel)./nnz((d3sel|d6sel) & greysel)




