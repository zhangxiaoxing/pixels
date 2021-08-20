function rings=ring_list_bz(opt)
arguments
    opt.shufid double = []
end
%bzthres=250;  %TODO filter by spike number % not really necessary when using full-length data
if isempty(opt.shufid)
    [sig,~]=bz.load_sig_pair();
    fname='rings_bz.mat';
else
    load('bz_ring_shufs.mat','shufs')
    sig=shufs{opt.shufid};
    fname=sprintf('rings_bz_shuf_%d.mat',opt.shufid);
end

rings=cell(max(sig.sess),3);
for sess=1:max(sig.sess)
    disp(sess);
    for ring_size=3:5
        sess_sel = sig.sess==sess;
        if nnz(sess_sel)<3, continue;end
        sess_rings=bz.rings.find_rings_bz(sig.suid(sess_sel,:),ring_size);
        rings{sess,ring_size-2}=unique(bz.rings.flexsort(sess_rings),'rows');
    end
end
save(fullfile('bzdata',fname),'rings');

end