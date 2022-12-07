% Preferred method as of 22.12.06

function out=ring_list_bz_alt(opt)
arguments
    opt.poolsize (1,1) double = 2
end
    [sig,~]=bz.load_sig_sums_conn_file('pair',false);
    sig_reg_sel=all(ismember(sig.reg(:,1,:),[343,567]),3);

    sig.sess=sig.sess(sig_reg_sel);
    sig.reg=sig.reg(sig_reg_sel,:,:);
    sig.suid=sig.suid(sig_reg_sel,:);
    rings=onerpt(sig);

    load('bz_ring_shufs.mat','shufs')
    rings_shuf=cell(size(shufs,1),1);
    poolo=parpool(opt.poolsize);
    parfor rpt=1:size(shufs,1)
        fprintf('Shuf %d\n',rpt);
        rings_shuf{rpt}=onerpt(shufs{rpt});
    end
    blame=vcs.blame();
    blame.tag='ring count actual vs shuffle';
    fname='rings_bz_vs_shuf.mat';
    save(fullfile('bzdata',fname),'rings_shuf','rings','blame');
    delete(poolo);
    out=struct();
    out.rings=rings;
    out.rings_shuf=rings_shuf;
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

