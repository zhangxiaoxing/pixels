ephys.util.dependency('buz',false);
sig=bz.load_sig_sums_conn_file();
%sig=bz.load_sig_pair();
sess=unique(sig.sess);
for si=1:numel(sess)
    onesess=sess(si);
    sesssel=sig.sess==onesess;
    bz.fccoding.fc_coding_one_sess(onesess,sig.suid(sesssel,:));
end
