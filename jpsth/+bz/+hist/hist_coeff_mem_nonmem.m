function hist_coeff_mem_nonmem(cong_range,nm_range,opt)
arguments
    cong_range (1,:) double
    nm_range (1,:) double
    opt.prefix (1,:) char = '0331'
end
sig=bz.load_sig_pair();
congrusel=find(all(ismember(sig.mem_type,1:2),2) | all(ismember(sig.mem_type,3:4),2));
if ~isempty(cong_range)
    congru_samp=congrusel(cong_range)';
    cong_spk=nan(numel(cong_range),11);
    cong_fc=nan(numel(cong_range),11);
    sidx=1;
    for i=congru_samp
        fprintf('%04d\n',sidx);
        [cong_spk(sidx,:),cong_fc(sidx,:)]=bz.hist.history_coeff(sig.sess(i),sig.suid(i,:));
        sidx=sidx+1;
    end
else
    cong_spk=[];
    cong_fc=[];
end
if ~isempty(nm_range)
    nm_sel=find(all(sig.mem_type==0,2));
    nm_samp=nm_sel(nm_range)';
    nm_spk=nan(numel(nm_range),11);
    nm_fc=nan(numel(nm_range),11);
    sidx=1;
    for i=nm_samp
        fprintf('%04d\n',sidx);
        [nm_spk(sidx,:),nm_fc(sidx,:)]=bz.hist.history_coeff(sig.sess(i),sig.suid(i,:));
        sidx=sidx+1;
    end
else
    nm_spk=[];
    nm_fc=[];
end

save(sprintf('%s_stp_%d_%d.mat',opt.prefix,min(cong_range),min(nm_range)),'cong_fc','cong_spk','nm_fc','nm_spk','cong_range','nm_range');

end