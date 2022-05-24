function out=get_cross_bin_ratio(sig_bin,pair_bin)


arguments
    sig_bin (:,2) int8
    pair_bin (:,2) int8
end

avai_bins=unique(abs(unique(sig_bin)));
out=zeros(max(avai_bins),max(avai_bins),5);
for lead_bin_id=1:max(avai_bins)
    for follow_bin_id=1:max(avai_bins)
        bin_cnt_s1=nnz(pair_bin(:,1)==lead_bin_id & pair_bin(:,2)==follow_bin_id);
        bin_cnt_s2=nnz(pair_bin(:,1)==-lead_bin_id & pair_bin(:,2)==-follow_bin_id);
        if (bin_cnt_s1+bin_cnt_s2)>0
            sig_cnt_s1=nnz(sig_bin(:,1)==lead_bin_id & sig_bin(:,2)==follow_bin_id);
            sig_cnt_s2=nnz(sig_bin(:,1)==-lead_bin_id & sig_bin(:,2)==-follow_bin_id);
            [phat,pci]=binofit(sig_cnt_s1+sig_cnt_s2,bin_cnt_s1+bin_cnt_s2);
            out(lead_bin_id,follow_bin_id,:)=[phat,pci,sig_cnt_s1+sig_cnt_s2,bin_cnt_s1+bin_cnt_s2];
        else
            out(lead_bin_id,follow_bin_id,:)=[0,0,0,0,0];
        end
    end
end
end

