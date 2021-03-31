function out=get_ratio(sig_type,pair_type)
% 0=NM,1=S1 sust, 2=S1 trans, 3=S2 sust, 4=S2 trans,-1=switched

arguments
    sig_type (:,2) int32
    pair_type (:,2) int32
end
out=struct();

nm_p=nnz(all(pair_type==0,2));
if nm_p>0
    out.nm_nm=nnz(all(sig_type==0,2))/nm_p;
end

congr_p=nnz(all(ismember(pair_type,1:2),2) | all(ismember(pair_type,3:4),2));
if congr_p>0
    out.congr=nnz(all(ismember(sig_type,1:2),2) | all(ismember(sig_type,3:4),2))/congr_p;
end

incon_p=nnz(...
    (ismember(pair_type(:,1),1:2) & ismember(pair_type(:,2),3:4))...
    |(ismember(pair_type(:,1),3:4) & ismember(pair_type(:,2),1:2)));

if incon_p>0
    out.incon=nnz(...
        (ismember(sig_type(:,1),1:2) & ismember(sig_type(:,2),3:4))...
        |(ismember(sig_type(:,1),3:4) & ismember(sig_type(:,2),1:2)))/incon_p;
end

mem_nm_p=nnz(pair_type(:,1)>0 & pair_type(:,2)==0);
if mem_nm_p>0
    out.mem_nm=nnz(sig_type(:,1)>0 & sig_type(:,2)==0)/mem_nm_p;
end

nm_mem_p=nnz(pair_type(:,1)==0 & pair_type(:,2)>0);
if nm_mem_p>0
    out.nm_mem=nnz(sig_type(:,1)==0 & sig_type(:,2)>0)/nm_mem_p;
end
end