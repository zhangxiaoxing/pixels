function out=get_ratio(sig_type,pair_type,opt)
% 0=NM,1=S1 sust, 2=S1 trans, 3=S2 sust, 4=S2 trans,-1=switched

arguments
    sig_type (:,2) int32
    pair_type (:,2) int32
    opt.nm_mem (1,1) logical = false
end
out=struct();

nm_p=nnz(all(pair_type==0,2));
if nm_p>0
    sig=nnz(all(sig_type==0,2));
    [phat,pci]=binofit(sig,nm_p);
    out.nm_nm=[phat,pci,sig,nm_p];
else
    out.nm_nm=[0,0,0,0];
end

congr_p=nnz(all(ismember(pair_type,1:2),2) | all(ismember(pair_type,3:4),2));
if congr_p>0
    sig=nnz(all(ismember(sig_type,1:2),2) | all(ismember(sig_type,3:4),2));
    [phat,pci]=binofit(sig,congr_p);
    out.congr=[phat,pci,sig,congr_p];
else
    out.congr=[0,0,0,0];
end

incon_p=nnz(...
    (ismember(pair_type(:,1),1:2) & ismember(pair_type(:,2),3:4))...
    |(ismember(pair_type(:,1),3:4) & ismember(pair_type(:,2),1:2)));

if incon_p>0
    sig=nnz(...
        (ismember(sig_type(:,1),1:2) & ismember(sig_type(:,2),3:4))...
        |(ismember(sig_type(:,1),3:4) & ismember(sig_type(:,2),1:2)));
    [phat,pci]=binofit(sig,incon_p);
    out.incon=[phat,pci,sig,incon_p];
else
    out.incon=[0,0,0,0];
end

if opt.nm_mem
    mem_nm_p=nnz(pair_type(:,1)>0 & pair_type(:,2)==0);
    if mem_nm_p>0
        out.mem_nm=nnz(sig_type(:,1)>0 & sig_type(:,2)==0)/mem_nm_p;
    end
    
    nm_mem_p=nnz(pair_type(:,1)==0 & pair_type(:,2)>0);
    if nm_mem_p>0
        out.nm_mem=nnz(sig_type(:,1)==0 & sig_type(:,2)>0)/nm_mem_p;
    end
end
end