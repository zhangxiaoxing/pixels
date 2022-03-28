asel=any(anovameta.anovap(:,[1,4,5,7])<0.05,2);
meta6=ephys.util.load_meta('delay',6);
meta3=ephys.util.load_meta('delay',3);
greysel=ismember(meta6.reg_tree(1,:),{'CH','BS'}).';

d6sel=any(meta6.fdr(1:4,:)<0.05).';
d3sel=any(meta3.fdr(1:4,:)<0.05).';


overlap_rate=nnz(asel & (d3sel | d6sel))/nnz((d3sel|d6sel))
grey_overlap_rate=nnz(asel & (d3sel | d6sel) & greysel)./nnz((d3sel|d6sel) & greysel)




