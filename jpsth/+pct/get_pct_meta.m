function pct_meta=get_pct_meta(eff_meta,sens_win,dur_win,opt)
arguments
    eff_meta,sens_win,dur_win
    opt.single_mod_thresh (1,1) double {mustBeMember(opt.single_mod_thresh,[4])} = 4 % option=5 removed as of 22.09.05
    opt.extended (1,1) logical = true
    opt.exlude_corner (1,1) logical = false
end

persistent pct_meta_ opt_

if isempty(pct_meta_) || ~isequaln(opt,opt_)
    pct_meta=struct();
    %% mat coord
    pct_meta.mat_coord=repmat(-1,size(eff_meta.cohen_d_dur,1),2);
    sens_efsz=max(abs(eff_meta.cohen_d_olf),[],2);
    dur_efsz=max(abs(eff_meta.cohen_d_dur),[],2);

    for olf_coord=1:(numel(sens_win)-1)
        for dur_coord=1:(numel(dur_win)-1)
            sens_sel=sens_efsz>sens_win(olf_coord) & sens_efsz<=sens_win(olf_coord+1);
            dur_sel=dur_efsz>dur_win(dur_coord) & dur_efsz<=dur_win(dur_coord+1);
            pct_meta.mat_coord(sens_sel & dur_sel,1)=olf_coord;
            pct_meta.mat_coord(sens_sel & dur_sel,2)=dur_coord;
        end
    end

    %% wave id >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    pct_meta.wave_id=repmat(-1,size(eff_meta.cohen_d_dur,1),1);
    %%  mixed
    % boundries are >= and <
    if ~opt.extended   
        sens_sel=max(abs(eff_meta.cohen_d_olf),[],2)<=sens_win(2); % none
        dur_sel=max(abs(eff_meta.cohen_d_dur),[],2)<=dur_win(2);
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
    else
        corner=false;
        if opt.exlude_corner
            sens_sel=max(abs(eff_meta.cohen_d_olf),[],2)<=sens_win(2); % none
            dur_sel=max(abs(eff_meta.cohen_d_dur),[],2)<=dur_win(2);
            corner=sens_sel & dur_sel;
        end
        sens_lsel=max(abs(eff_meta.cohen_d_olf),[],2)<=sens_win(4); % none
        dur_lsel=max(abs(eff_meta.cohen_d_dur),[],2)<=dur_win(4);
        pct_meta.wave_id(sens_lsel & dur_lsel & ~corner)=0;

        if opt.exlude_corner
            sens_sel=max(eff_meta.cohen_d_olf,[],2)>sens_win(5); % s1 d3
            dur_sel=max(eff_meta.cohen_d_dur,[],2)>dur_win(5);
            corner=sens_sel & dur_sel;
        end
        sens_lsel=max(eff_meta.cohen_d_olf,[],2)>sens_win(4); % s1 d3
        dur_lsel=max(eff_meta.cohen_d_dur,[],2)>dur_win(4);
        pct_meta.wave_id(sens_lsel & dur_lsel & ~corner)=1;
        
        if opt.exlude_corner
            sens_sel=max(eff_meta.cohen_d_olf,[],2)>sens_win(5); %s1 d6
            dur_sel=max(eff_meta.cohen_d_dur.*-1,[],2)>dur_win(5);
            corner=sens_sel & dur_sel;
        end
        sens_lsel=max(eff_meta.cohen_d_olf,[],2)>sens_win(4); %s1 d6
        dur_lsel=max(eff_meta.cohen_d_dur.*-1,[],2)>dur_win(4);
        pct_meta.wave_id(sens_lsel & dur_lsel & ~corner)=2;

        if opt.exlude_corner
            sens_sel=max(eff_meta.cohen_d_olf.*-1,[],2)>sens_win(5); %s2 d3
            dur_sel=max(eff_meta.cohen_d_dur,[],2)>dur_win(5);
            corner=sens_sel & dur_sel;
        end
        sens_lsel=max(eff_meta.cohen_d_olf.*-1,[],2)>sens_win(4); %s2 d3
        dur_lsel=max(eff_meta.cohen_d_dur,[],2)>dur_win(4);
        pct_meta.wave_id(sens_lsel & dur_lsel & ~corner)=3;

        if opt.exlude_corner
            sens_sel=max(eff_meta.cohen_d_olf.*-1,[],2)>sens_win(5); %s2 d6
            dur_sel=max(eff_meta.cohen_d_dur.*-1,[],2)>dur_win(5);
            corner=sens_sel & dur_sel;
        end
        sens_lsel=max(eff_meta.cohen_d_olf.*-1,[],2)>sens_win(4); %s2 d6
        dur_lsel=max(eff_meta.cohen_d_dur.*-1,[],2)>dur_win(4);
        pct_meta.wave_id(sens_lsel & dur_lsel & ~corner)=4;
    end
    %% olf
    if ~opt.extended
        sens_sel=max(eff_meta.cohen_d_olf,[],2)>sens_win(opt.single_mod_thresh); % s1
        dur_sel=max(abs(eff_meta.cohen_d_dur),[],2)<=dur_win(2);
        pct_meta.wave_id(sens_sel & dur_sel)=5;

        sens_sel=max(eff_meta.cohen_d_olf.*-1,[],2)>sens_win(opt.single_mod_thresh); % s2
        dur_sel=max(abs(eff_meta.cohen_d_dur),[],2)<=dur_win(2);
        pct_meta.wave_id(sens_sel & dur_sel)=6;
    else
        corner=false;
        if opt.exlude_corner
            sens_sel=max(eff_meta.cohen_d_olf,[],2)>sens_win(opt.single_mod_thresh); % s1
            dur_sel=max(abs(eff_meta.cohen_d_dur),[],2)<=dur_win(2);
            corner=sens_sel & dur_sel;
        end
        sens_lsel=max(eff_meta.cohen_d_olf,[],2)>sens_win(opt.single_mod_thresh); % s1
        dur_lsel=max(abs(eff_meta.cohen_d_dur),[],2)<=dur_win(4);
        pct_meta.wave_id(sens_lsel & dur_lsel & ~corner)=5;

        if opt.exlude_corner
            sens_sel=max(eff_meta.cohen_d_olf.*-1,[],2)>sens_win(opt.single_mod_thresh); % s2
            dur_sel=max(abs(eff_meta.cohen_d_dur),[],2)<=dur_win(2);
            corner=sens_sel & dur_sel;
        end
        sens_lsel=max(eff_meta.cohen_d_olf.*-1,[],2)>sens_win(opt.single_mod_thresh); % s2
        dur_lsel=max(abs(eff_meta.cohen_d_dur),[],2)<=dur_win(4);
        pct_meta.wave_id(sens_lsel & dur_lsel & ~corner)=6;
    end

    %% dur
    if ~opt.extended
        sens_sel=max(abs(eff_meta.cohen_d_olf),[],2)<=sens_win(2); % d3
        dur_sel=max(eff_meta.cohen_d_dur,[],2)>dur_win(opt.single_mod_thresh);
        pct_meta.wave_id(sens_sel & dur_sel)=7;

        sens_sel=max(abs(eff_meta.cohen_d_olf),[],2)<=sens_win(2); %d6
        dur_sel=max(eff_meta.cohen_d_dur.*-1,[],2)>dur_win(opt.single_mod_thresh);
        pct_meta.wave_id(sens_sel & dur_sel)=8;
    else
        corner=false;
        if opt.exlude_corner
            sens_sel=max(abs(eff_meta.cohen_d_olf),[],2)<=sens_win(2); % d3
            dur_sel=max(eff_meta.cohen_d_dur,[],2)>dur_win(opt.single_mod_thresh);
            corner=sens_sel & dur_sel;
        end
        sens_lsel=max(abs(eff_meta.cohen_d_olf),[],2)<=sens_win(4); % d3
        dur_lsel=max(eff_meta.cohen_d_dur,[],2)>dur_win(opt.single_mod_thresh);
        pct_meta.wave_id(sens_lsel & dur_lsel & ~corner)=7;

        if opt.exlude_corner
            sens_sel=max(abs(eff_meta.cohen_d_olf),[],2)<=sens_win(2); %d6
            dur_sel=max(eff_meta.cohen_d_dur.*-1,[],2)>dur_win(opt.single_mod_thresh);
            corner=sens_sel & dur_sel;
        end
        sens_lsel=max(abs(eff_meta.cohen_d_olf),[],2)<=sens_win(4); %d6
        dur_lsel=max(eff_meta.cohen_d_dur.*-1,[],2)>dur_win(opt.single_mod_thresh);
        pct_meta.wave_id(sens_lsel & dur_lsel & ~corner)=8;
    end

    pct_meta_=pct_meta;
    opt_=opt;
end
pct_meta=pct_meta_;
end