function rings=ring_list_bz_alt(opt)
arguments
    opt.shuf (1,1) logical = false
    opt.poolsize (1,1) double = 2
end
%bzthres=250;  %TODO filter by spike number % not really necessary when using full-length data

if ~opt.shuf
    [sig,~]=bz.load_sig_pair();
    rings=onerpt(sig);
    fname='rings_bz.mat';
    save(fullfile('bzdata',fname),'rings');
else
    load('bz_ring_shufs.mat','shufs')
    rings_shuf=cell(size(shufs,1),1);
    parpool(opt.poolsize);
    parfor rpt=1:size(shufs,1)
        fprintf('Shuf %d\n',rpt);
        rings_shuf{rpt}=onerpt(shufs{rpt});
    end
    fname='rings_bz_shuf.mat';
    save(fullfile('bzdata',fname),'rings_shuf');
end
end
function rings=onerpt(sig)
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
end

