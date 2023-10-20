function fc=join_fc_waveid(fc,waveid,opt)
arguments
    fc
    waveid
    opt.pct_mat (1,1) logical = false
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
end
    meta=ephys.util.load_meta('skip_stats',true,'load_file',false,'criteria',opt.criteria);
    usess=unique(fc.sess);
    if opt.pct_mat
        fc.pct_coord=zeros(size(fc.suid,1),4)-1;
    else
        fc.waveid=zeros(size(fc.suid))-1;
    end
    for ss=reshape(usess,1,[])
        if opt.pct_mat
            metaid=int32(meta.allcid(meta.sess==ss));
            sess_coord=waveid(meta.sess==ss,:);
            suidL=fc.suid(fc.sess==ss,1);
            [~,idL]=ismember(suidL,metaid);
            fc.pct_coord(fc.sess==ss,1:2)=sess_coord(idL,:);

            suidF=fc.suid(fc.sess==ss,2);
            [~,idF]=ismember(suidF,metaid);
            fc.pct_coord(fc.sess==ss,3:4)=sess_coord(idF,:);
        else
            if ~any(meta.sess==ss)
                continue
            end
            metaid=int32(meta.allcid(meta.sess==ss));
            sesswave=waveid(meta.sess==ss);
            suidL=fc.suid(fc.sess==ss,1);
            [~,idL]=ismember(suidL,metaid);
            fc.waveid(fc.sess==ss,1)=sesswave(idL);

            suidF=fc.suid(fc.sess==ss,2);
            [~,idF]=ismember(suidF,metaid);
            fc.waveid(fc.sess==ss,2)=sesswave(idF);
        end

    end
end