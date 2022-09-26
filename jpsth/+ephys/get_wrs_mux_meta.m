function out_=get_wrs_mux_meta(opt)
arguments
    opt.permutation (1,1) logical = false
    %     opt.merge_bin (1,1) logical = true
    opt.load_file (1,1) logical = true
    opt.save_file (1,1) logical = false
    opt.perm_repeat (1,1) double {mustBePositive,mustBeInteger} = 1000
    opt.uneven_duration (1,1) logical = false % use all 6-sec delay
    opt.merge_mux (1,1) logical = true
end

persistent out opt_
if isempty(out) || ~isequaln(opt,opt_)
    if opt.load_file
        load('wrs_mux_meta.mat','wrs_mux_meta');
        out=wrs_mux_meta;
    else
        [~,~,sessmap]=ephys.sessid2path(0);
        homedir=ephys.util.getHomedir('type','raw');
        sesskeys=cell2mat(sessmap.keys());
        [out.pref_id,out.p_mux,out.p_olf,out.p_dur,out.class_fr]=deal([]);
        for sessid=sesskeys
            disp(sessid);
            fpath=fullfile(homedir,sessmap(sessid),"FR_All_1000.hdf5");
            fr=h5read(fpath,'/FR_All');
            trials=h5read(fpath,'/Trials');

            s2sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,5)==8;
            s1sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,5)==4;

            d3sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,8)==3;
            d6sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,8)==6;

            for suidx=1:size(fr,2)
                [pref_id,ppmux,ppolf,ppdur]=deal(nan(1,3));
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
                    [~,binpref]=max(class_mm);
                    
                    if false%opt.permutation
                        wrs_p_d3(suidx,:)=arrayfun(@(x) permutation_test_1d(fr(d3sel & s1sel,suidx,x),fr(d3sel & s2sel,suidx,x),opt.perm_repeat),bin3s);
                        wrs_p_d6(suidx,:)=arrayfun(@(x) permutation_test_1d(fr(d6sel & s1sel,suidx,x),fr(d6sel & s2sel,suidx,x),opt.perm_repeat),bin6s);
                    else
                        % determine multimodal
                        otherid=setdiff(1:4,binpref);
                        pref_id(bin-4)=binpref;
                        ppmux(bin-4)=max(arrayfun(@(x) ranksum(class_cell{binpref},class_cell{x}),otherid));
                        % determine monomodal
                        ppolf(bin-4)=ranksum(cell2mat(class_cell(1:2).'),cell2mat(class_cell(3:4).'));
                        ppdur(bin-4)=ranksum(cell2mat(class_cell([1,3]).'),cell2mat(class_cell([2,4]).'));

                    end
                end
                out.pref_id=[out.pref_id;pref_id];
                out.p_mux=[out.p_mux;ppmux];
                out.p_olf=[out.p_olf;ppolf];
                out.p_dur=[out.p_dur;ppdur];
                out.class_fr=[out.class_fr;class_fr];
            end

            m=any(out.p_mux<0.05,2);
            o=any(out.p_olf<0.05,2);
            d=any(out.p_dur<0.05,2);

            sel_mat=[m & ~o & ~d,...
                ~m & o & ~d,...
                ~m & ~o & d,...
                m & o & ~d,...
                m & ~o & d,...
                ~m & o & d,...
                m & o & d...
                ];
        end
        figure('Position',[100,100,300,300])
        [fh,vs]=lib.venn(sum(sel_mat));

        wave_id=zeros(size(m,1),1)-1;
        [~,midx]=min(out.p_mux,[],2);
        [~,oidx]=min(out.p_olf,[],2);
        [~,didx]=min(out.p_dur,[],2);

        class2olf=[5 5 6 6];
        class2dur=[7 8 7 8];

        wave_id(~o & ~m & ~d)=0;
        wave_id(m)=arrayfun(@(x) out.pref_id(x,midx(x)),find(m));
        wave_id(o & ~m & ~d)=arrayfun(@(x) class2olf(out.pref_id(x,oidx(x))),find(o & ~m & ~d));
        wave_id(~o & ~m & d)=arrayfun(@(x) class2dur(out.pref_id(x,didx(x))),find(~o & ~m & d));
        
        if opt.merge_mux
            wave_id(o & ~m & d)=arrayfun(@(x) out.pref_id(x,midx(x)),find(o & ~m & d));
        end

        out.wave_id=wave_id;

        if  opt.save_file
            wrs_mux_meta=out;
            save('wrs_mux_meta.mat','wrs_mux_meta');
        end
    end
end
opt_=opt;
out_=out;
end
