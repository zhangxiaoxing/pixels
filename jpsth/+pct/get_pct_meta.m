function pct_meta=get_pct_meta(eff_meta,sens_efsz,sens_win,dur_efsz,dur_win,opt)
arguments
    eff_meta,sens_efsz,sens_win,dur_efsz,dur_win
    opt.single_mod_thresh (1,1) double {mustBeMember(opt.single_mod_thresh,[2,4,5])} = 4
end
persistent pct_meta_ opt_
if isempty(pct_meta_) || ~isequaln(opt,opt_)
    pct_meta=struct();
    pct_meta.class_id=zeros(size(sens_efsz));
    pct_meta.class_id(sens_efsz<sens_win(2) & dur_efsz<dur_win(2))=1; % none
    pct_meta.class_id(sens_efsz>=sens_win(opt.single_mod_thresh) & dur_efsz<dur_win(2))=2; % olf
    pct_meta.class_id(sens_efsz<sens_win(2) & dur_efsz>=dur_win(opt.single_mod_thresh))=3; % duration
    pct_meta.class_id(sens_efsz>=sens_win(5) & dur_efsz>=dur_win(5))=4; % both
    %% wave id >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    pct_meta.wave_id=repmat(-1,size(sens_efsz));
    %%  mixed

    sens_sel=max(abs(eff_meta.cohen_d_olf),[],2)<sens_win(2); % none
    dur_sel=max(abs(eff_meta.cohen_d_dur),[],2)<dur_win(2);
    pct_meta.wave_id(sens_sel & dur_sel)=0;

    sens_sel=max(eff_meta.cohen_d_olf,[],2)>sens_win(5); % s1 d3
    dur_sel=max(eff_meta.cohen_d_dur,[],2)>dur_win(5);
    pct_meta.wave_id(sens_sel & dur_sel)=1;

    sens_sel=max(eff_meta.cohen_d_olf,[],2)>sens_win(5); %s1 d6
    dur_sel=max(eff_meta.cohen_d_dur.*-1,[],2)>dur_win(5);
    pct_meta.wave_id(sens_sel & dur_sel)=2;

    sens_sel=max(eff_meta.cohen_d_olf.*-1,[],2)>sens_win(5); %s2 d3
    dur_sel=max(eff_meta.cohen_d_dur,[],2)>dur_win(5);
    pct_meta.wave_id(sens_sel & dur_sel)=3;

    sens_sel=max(eff_meta.cohen_d_olf.*-1,[],2)>sens_win(5); %s2 d6
    dur_sel=max(eff_meta.cohen_d_dur.*-1,[],2)>dur_win(5);
    pct_meta.wave_id(sens_sel & dur_sel)=4;

    %% olf
    sens_sel=max(eff_meta.cohen_d_olf,[],2)>sens_win(opt.single_mod_thresh); % s1
    dur_sel=max(abs(eff_meta.cohen_d_dur),[],2)<dur_win(2);
    pct_meta.wave_id(sens_sel & dur_sel)=5;

    sens_sel=max(eff_meta.cohen_d_olf.*-1,[],2)>sens_win(opt.single_mod_thresh); % s2
    dur_sel=max(abs(eff_meta.cohen_d_dur),[],2)<dur_win(2);
    pct_meta.wave_id(sens_sel & dur_sel)=6;

    %% dur
    sens_sel=max(abs(eff_meta.cohen_d_olf),[],2)<sens_win(2); % d3
    dur_sel=max(eff_meta.cohen_d_dur,[],2)>dur_win(opt.single_mod_thresh);
    pct_meta.wave_id(sens_sel & dur_sel)=7;

    sens_sel=max(abs(eff_meta.cohen_d_olf),[],2)<sens_win(2); %d6
    dur_sel=max(eff_meta.cohen_d_dur.*-1,[],2)>dur_win(opt.single_mod_thresh);
    pct_meta.wave_id(sens_sel & dur_sel)=8;
    pct_meta_=pct_meta;
    opt_=opt;
end
pct_meta=pct_meta_;
end