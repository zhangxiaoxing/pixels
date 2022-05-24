function fc=join_fc_waveid(fc,waveid)
    meta=ephys.util.load_meta('skip_stats',true);
    usess=unique(fc.sess);
    fc.waveid=zeros(size(fc.suid));
    for ss=reshape(usess,1,[])
        metaid=int32(meta.allcid(meta.sess==ss));
        sesswave=waveid(meta.sess==ss);
        suidL=fc.suid(fc.sess==ss,1);
        [isin,idL]=ismember(suidL,metaid);
        fc.waveid(fc.sess==ss,1)=sesswave(idL);

        suidF=fc.suid(fc.sess==ss,2);
        [isin,idF]=ismember(suidF,metaid);
        fc.waveid(fc.sess==ss,2)=sesswave(idF);
    end
end