function out_=get_a2_meta(opt)
arguments
    opt.load_file (1,1) logical = false
    opt.save_file (1,1) logical = false
    opt.filename (1,:) char = 'a2_meta.mat'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','Naive','any'})} = 'WT'
end

persistent out opt_
if isempty(out) || ~isequaln(opt,opt_)
    if opt.load_file
        load(fullfile('binary',opt.filename),'sel_meta');
        out=sel_meta;
    else
        [~,~,sessmap]=ephys.sessid2path(0,'criteria',opt.criteria);
        sesskeys=cell2mat(sessmap.keys());
        [out.o_pref_id,out.p_olf,out.p_int,out.class_fr]=deal([]);

        for sessid=sesskeys
            if rem(sessid,10)==0, disp(sessid);end
            fpath=fullfile(replace(sessmap(sessid),("\"|"/"),filesep),"FR_All_1000.hdf5");
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
                case 'Naive'
                    trials=behav.procPerf(trials,'criteria','Learning');
                    s2sel=trials(:,5)==8;
                    s1sel=trials(:,5)==4;
                    d3sel=trials(:,8)==3;
                    d6sel=trials(:,8)==6;
                otherwise
                    warning("Unfinished code path")
                    keyboard()
            end

            for suidx=1:size(fr,2)

                s1d3fr=squeeze(fr(d3sel & s1sel,suidx,5:7));
                s1d6fr=squeeze(fr(d6sel & s1sel,suidx,5:10));
                s2d3fr=squeeze(fr(d3sel & s2sel,suidx,5:7));
                s2d6fr=squeeze(fr(d6sel & s2sel,suidx,5:10));
                %                 % anovan([s1d3fr;s1d6fr;s2d3fr;s2d6fr],...
                %     {[zeros(numel(s1d3fr)+numel(s1d6fr),1);ones(numel(s2d3fr)+numel(s2d6fr),1)],...
                %     [zeros(numel(s1d3fr),1);ones(numel(s1d6fr),1);zeros(numel(s2d3fr),1);ones(numel(s2d6fr),1)],...
                %     [reshape(repmat(1:3,numel(s1d3fr)/3,1),[],1);...
                %     reshape(repmat(1:6,numel(s1d6fr)/6,1),[],1);...
                %     reshape(repmat(1:3,numel(s2d3fr)/3,1),[],1);...
                %     reshape(repmat(1:6,numel(s2d6fr)/6,1),[],1)]},'model','linear')

                p=anovan([s1d3fr(:);s1d6fr(:);s2d3fr(:);s2d6fr(:)],...
                    {[zeros(numel(s1d3fr)+numel(s1d6fr),1);ones(numel(s2d3fr)+numel(s2d6fr),1)],...
                    [reshape(repmat(1:3,numel(s1d3fr)/3,1),[],1);...
                    reshape(repmat(4:9,numel(s1d6fr)/6,1),[],1);...
                    reshape(repmat(1:3,numel(s2d3fr)/3,1),[],1);...
                    reshape(repmat(4:9,numel(s2d6fr)/6,1),[],1)]},'model','interaction','display','off');

                out.p_olf=[out.p_olf;p(1)];
                out.p_int=[out.p_int;p(3)];

                s1m=[mean(s1d3fr),mean(s1d6fr)];
                s2m=[mean(s2d3fr),mean(s2d6fr)];

                if p(1)<0.05
                    if mean([s1d3fr(:);s1d6fr(:)])>mean([s2d3fr(:);s2d6fr(:)])
                        out.o_pref_id=[out.o_pref_id;5];    
                    else
                        out.o_pref_id=[out.o_pref_id;6];
                    end
                elseif p(3)<0.05
                    pp=[arrayfun(@(x)ranksum(s1d3fr(:,x),s2d3fr(:,x)), 1:3),...
                        arrayfun(@(x)ranksum(s1d6fr(:,x),s2d6fr(:,x)), 1:6)];
                    if all(pp>0.05)
                        % disp("unknown bin")
                        out.o_pref_id=[out.o_pref_id;-2];
                    else

                        if all(s1m(pp<0.05)>s2m(pp<0.05))
                            out.o_pref_id=[out.o_pref_id;5];
                        elseif all(s1m(pp<0.05)<s2m(pp<0.05))
                            out.o_pref_id=[out.o_pref_id;6];
                        else
                            out.o_pref_id=[out.o_pref_id;-1];
                        end
                    end
                else
                    out.o_pref_id=[out.o_pref_id;0];
                end

                out.class_fr=[out.class_fr;s1m,s2m];
            end
        end

        out.wave_id=out.o_pref_id;

        if  opt.save_file
            sel_meta=out;
            blame=vcs.blame();
            save(fullfile('binary',opt.filename),'sel_meta','blame');
        end
    end
end

opt_=opt;
out_=out;
end
