function dur_sel=get_dul_sel(opt)
arguments
    opt.exclusive (1,1) logical = true
    opt.ranksum_stats (1,1) logical = true
end

persistent dur_sel_ opt_

if isempty(dur_sel_) || isequaln(opt,opt_)
    assert(opt.exclusive);
    % idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
    anovameta=wave.get_dur_waveid();
    %1:Samp,2:Dur,3:Bin,4:Samp*Dur,5:Samp*Bin,6:Dur*Bin,7:Samp*Dur*Bin
    if strcmp(dur_stats_path,'anovaonly')
        meta=ephys.util.load_meta();
        waveid=ephys.get_wave_id(meta.sess,meta.allcid);
        % waveid(waveid==0 & anovameta.dur_waveid==3)=7;
        % waveid(waveid==0 & anovameta.dur_waveid==6)=8;
        dur_sel=anovameta.dur_waveid>0 & waveid==0;
    else
        dur_sel_mix=any(anovameta.anovap(:,[2 4 6 7])<0.05,2);
        sens_sel_mix=any(anovameta.anovap(:,[1 4 5 7])<0.05,2);
        dur_sel=dur_sel_mix & ~sens_sel_mix;
    end
    opt_=opt;
    dur_sel_=dur_sel;
end

dur_sel=dur_sel_;
