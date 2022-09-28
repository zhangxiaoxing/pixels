
ephys.util.dependency('buz',false);
sig=bz.load_sig_sums_conn_file();
%sig=bz.load_sig_pair();
sess=unique(sig.sess);
if isunix
    ph=parpool(16);
end
for si=1:numel(sess)
    onesess=sess(si);
    sesssel=sig.sess==onesess;
    
    if isunix
        ff(si)=parfeval(@bz.fccoding.fc_coding_one_sess,1,onesess,sig.suid(sesssel,:));
    else
        bz.fccoding.fc_coding_one_sess(onesess,sig.suid(sesssel,:));
    end
end

