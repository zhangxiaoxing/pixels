function out_=get_wrs_meta(opt)
arguments
    opt.permutation (1,1) logical = false
    %     opt.merge_bin (1,1) logical = true
    opt.load_file (1,1) logical = true
    opt.save_file (1,1) logical = false
    opt.perm_repeat (1,1) double {mustBePositive,mustBeInteger} = 1000
    opt.extend6s (1,1) logical = true
    opt.tag_ext6s_mem (1,1) logical = false
    
end

persistent out opt_
if isempty(out) || ~isequaln(opt,opt_)
    if opt.load_file
        load('wrs_meta.mat','wrs_meta');
        out=wrs_meta;
    else
        [~,~,sessmap]=ephys.sessid2path(0);
        homedir=ephys.util.getHomedir('type','raw');
        sesskeys=cell2mat(sessmap.keys());
        [out.o_pref_id,out.p_olf,out.class_fr]=deal([]);
        if opt.extend6s
            [out.m_pref_id6,out.p_olf6,out.class_fr6]=deal([]);
        end
        for sessid=sesskeys
            disp(sessid);
            fpath=fullfile(homedir,sessmap(sessid),"FR_All_1000.hdf5");
            fr=h5read(fpath,'/FR_All');
            trials=h5read(fpath,'/Trials');

            s2sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,5)==8;
            s1sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,5)==4;

            for suidx=1:size(fr,2)
                [opref_id,ppolf]=deal(nan(1,3));
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

                    % determine multimodal
                    otherid=setdiff(1:4,binmpref);
                    opref_id(bin-4)=binopref;

                    % determine monomodal
                    ppolf(bin-4)=ranksum(cell2mat(class_cell(1:2).'),cell2mat(class_cell(3:4).'));
                    %                     end
                end

                out.o_pref_id=[out.o_pref_id;opref_id];
                out.p_olf=[out.p_olf;ppolf];

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
                        [~,binmpref]=max(class_mm);  % highest condition

                        % determine monomodal
                        pref_id6(bin-7)=binmpref;
                        ppolf6(bin-7)=ranksum(cell2mat(class_cell(2).'),cell2mat(class_cell(4).'));
                        %                         end
                    end
                    out.m_pref_id6=[out.m_pref_id6;pref_id6];
                    out.p_olf6=[out.p_olf6;ppolf6];
                    out.class_fr6=[out.class_fr6;class_fr6];
                end
            end
        end
        o=any(out.p_olf<0.05,2);

        wave_id=zeros(size(m,1),1)-1;
        [~,oidx]=min(out.p_olf,[],2);

        wave_id(~o & ~m & ~d)=0;
        wave_id(o & ~m & ~d)=arrayfun(@(x) out.o_pref_id(x,oidx(x)),find(o & ~m & ~d));
        
        out.wave_id=wave_id;

        if  opt.save_file
            wrs_meta=out;
            save('wrs_meta.mat','wrs_meta');
        end
    end

    if opt.tag_ext6s_mem
        ext_sel=out.wave_id==0 & any(out.p_olf6<0.05,2);
        out.wave_id(ext_sel)=-1;
    end
end


opt_=opt;
out_=out;
end
