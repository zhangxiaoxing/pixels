function out_=get_wrs_mux_meta(opt)
arguments
    opt.permutation (1,1) logical = false
    opt.load_file (1,1) logical = true
    opt.save_file (1,1) logical = false
    opt.perm_repeat (1,1) double {mustBePositive,mustBeInteger} = 1000
    opt.merge_mux (1,1) logical = true
    opt.extend6s (1,1) logical = true
    opt.tag_ext6s_mem (1,1) logical = false
    opt.plot_venn (1,1) logical = false
    opt.filename (1,:) char = 'wrs_mux_meta.mat'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
    opt.fdr (1,1) logical = false
end

persistent out opt_
if isempty(out) || ~isequaln(opt,opt_)
    if opt.load_file
        load(fullfile('binary',opt.filename),'wrs_mux_meta');
        out=wrs_mux_meta;
    else
        [~,~,sessmap]=ephys.sessid2path(0,'criteria',opt.criteria);
        homedir=ephys.util.getHomedir();
        sesskeys=cell2mat(sessmap.keys());
        [out.m_pref_id,out.o_pref_id,out.d_pref_id,out.p_mux,out.p_olf,out.p_dur,out.class_fr]=deal([]);
        if opt.extend6s
            [out.o_pref_id6,out.p_olf6,out.class_fr6]=deal([]);
        end
        for sessid=sesskeys
            if rem(sessid,10)==0, disp(sessid);end
            fpath=fullfile(homedir,replace(sessmap(sessid),("\"|"/"),filesep),"FR_All_1000.hdf5");
            fr=h5read(fpath,'/FR_All');
            trials=h5read(fpath,'/Trials');
            switch opt.criteria % TODO : pending refactoring
                case 'WT'
                    s2sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,5)==8;
                    s1sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,5)==4;

                    d3sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,8)==3;
                    d6sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,8)==6;
                case 'Learning'
                    trials=behav.procPerf(trials,'criteria','Learning');
                    s2sel=trials(:,9)~=0 & trials(:,5)==8;
                    s1sel=trials(:,9)~=0 & trials(:,5)==4;

                    d3sel=trials(:,9)~=0 & trials(:,8)==3;
                    d6sel=trials(:,9)~=0 & trials(:,8)==6;
                otherwise
                    warning("Unfinished code path")
                    keyboard()
            end

            for suidx=1:size(fr,2)
                [mpref_id,opref_id,dpref_id,ppmux,ppolf,ppdur]=deal(nan(1,3));
                class_fr=nan(1,4,3);
                for bin=5:7
                    % highest condition
                    s1d3fr=fr(d3sel & s1sel,suidx,bin);
                    s1d6fr=fr(d6sel & s1sel,suidx,bin);
                    s2d3fr=fr(d3sel & s2sel,suidx,bin);
                    s2d6fr=fr(d6sel & s2sel,suidx,bin);
                    class_cell={s1d3fr,s1d6fr,s2d3fr,s2d6fr};
                    class_mm=cellfun(@(x) mean(x),class_cell);
                    class_fr(1,:,bin-4)=class_mm;
                    [~,binmpref]=max(class_mm);  % highest condition
                    
                    if mean([s1d3fr;s1d6fr])>mean([s2d3fr;s2d6fr])
                        binopref=5;
                    else
                        binopref=6;
                    end

                    if mean([s1d3fr;s2d3fr])>mean([s1d6fr;s2d6fr])
                        bindpref=7;
                    else
                        bindpref=8;
                    end

                    % determine multimodal
                    otherid=setdiff(1:4,binmpref);
                    mpref_id(bin-4)=binmpref;
                    opref_id(bin-4)=binopref;
                    dpref_id(bin-4)=bindpref;

                    ppmux(bin-4)=max(arrayfun(@(x) ranksum(class_cell{binmpref},class_cell{x}),otherid));
                    % determine monomodal
                    ppolf(bin-4)=ranksum(cell2mat(class_cell(1:2).'),cell2mat(class_cell(3:4).'));
                    ppdur(bin-4)=ranksum(cell2mat(class_cell([1,3]).'),cell2mat(class_cell([2,4]).'));
                    %                     end
                end
                out.m_pref_id=[out.m_pref_id;mpref_id];
                out.o_pref_id=[out.o_pref_id;opref_id];
                out.d_pref_id=[out.d_pref_id;dpref_id];

                out.p_mux=[out.p_mux;ppmux];
                out.p_olf=[out.p_olf;ppolf];
                out.p_dur=[out.p_dur;ppdur];
                out.class_fr=[out.class_fr;class_fr];
                if opt.extend6s
                    [pref_id6,ppolf6]=deal(nan(1,3));
                    class_fr6=nan(1,4,3);
                    for bin=8:10
                        % highest condition
                        s1d6fr=fr(d6sel & s1sel,suidx,bin);
                        s2d6fr=fr(d6sel & s2sel,suidx,bin);
                        class_cell={-1,s1d6fr,-1,s2d6fr};
                        class_mm=cellfun(@(x) mean(x),class_cell);
                        class_fr6(1,:,bin-7)=class_mm;
                        if class_mm(2)>=class_mm(4)  % highest condition
                            binopref=5;
                        else
                            binopref=6;
                        end
                        % determine monomodal
                        pref_id6(bin-7)=binopref;
                        ppolf6(bin-7)=ranksum(cell2mat(class_cell(2).'),cell2mat(class_cell(4).'));
                        %                         end
                    end
                    out.o_pref_id6=[out.o_pref_id6;pref_id6];
                    out.p_olf6=[out.p_olf6;ppolf6];
                    out.class_fr6=[out.class_fr6;class_fr6];
                end
            end
        end

        m=any(out.p_mux<0.05,2);
        o=any(out.p_olf<0.05,2);
        d=any(out.p_dur<0.05,2);

        if opt.plot_venn
            sel_mat=[m & ~o & ~d,...
                ~m & o & ~d,...
                ~m & ~o & d,...
                m & o & ~d,...
                m & ~o & d,...
                ~m & o & d,...
                m & o & d...
                ];
            figure('Position',[100,100,300,300])
            [fh,vs]=lib.venn(sum(sel_mat));
        end

        wave_id=zeros(size(m,1),1)-1;
        [~,midx]=min(out.p_mux,[],2);
        [~,oidx]=min(out.p_olf,[],2);
        [~,didx]=min(out.p_dur,[],2);

        wave_id(~o & ~m & ~d)=0;
        wave_id(m)=arrayfun(@(x) out.m_pref_id(x,midx(x)),find(m));
        wave_id(o & ~m & ~d)=arrayfun(@(x) out.o_pref_id(x,oidx(x)),find(o & ~m & ~d));
        wave_id(~o & ~m & d)=arrayfun(@(x) out.d_pref_id(x,didx(x)),find(~o & ~m & d));

        if opt.merge_mux
            wave_id(o & ~m & d)=arrayfun(@(x) out.m_pref_id(x,midx(x)),find(o & ~m & d));
        end
        
        out.wave_id=wave_id;

        if  opt.save_file
            wrs_mux_meta=out;
            blame=vcs.blame();
            save(fullfile('binary',opt.filename),'wrs_mux_meta','blame');
        end
    end

    if opt.tag_ext6s_mem
        ext_sel=out.wave_id==0 & any(out.p_olf6<0.05,2);
        out.wave_id(ext_sel)=-1;
    end
end

%% TODO: switch

opt_=opt;
out_=out;
end
