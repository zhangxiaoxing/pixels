% {FWD/RWD,H2L/L2H/Local, S1/S2*correct/error*3s/6s}
% congru, incongru, nonmem

ephys.util.dependency('buz',false);
sig=bz.load_sig_pair();
idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));
sess=unique(sig.sess);
parfor si=1:2%numel(sess)
    onesess=sess(si);
    sesssel=sig.sess==onesess;
    bz.fccoding.fc_coding_one_sess(onesess,sig.suid(sesssel,:));
end