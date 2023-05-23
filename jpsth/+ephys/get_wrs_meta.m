function out_=get_wrs_meta(opt)
arguments
    opt.load_file (1,1) logical = false
    opt.save_file (1,1) logical = false
    opt.fdr (1,1) logical = false
end

% ranksum delay bin 1-3 in all trials, then ranksum delay bin 4-6 in 6s
% trials, concatenate p-values vector
% TODO: optional FDR adjust

persistent out opt_
if isempty(out) || ~isequaln(opt,opt_)
    if opt.load_file
        load('wrs_meta.mat','wrs_meta');
        out=wrs_meta;
    else
        [~,~,sessmap]=ephys.sessid2path(0);
        homedir=ephys.util.getHomedir('type','raw');
        sesskeys=cell2mat(sessmap.keys());
        [out.p_olf,out.wave_id]=deal([]);

        for sessid=sesskeys
            disp(sessid);
            fpath=fullfile(homedir,sessmap(sessid),"FR_All_1000.hdf5");
            fr=h5read(fpath,'/FR_All');
            trials=h5read(fpath,'/Trials');

            s1sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,5)==4 & ismember(trials(:,8), [3 6]);
            s2sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,5)==8 & ismember(trials(:,8), [3 6]);
            
            s1d6sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,5)==4 & trials(:,8)==6;
            s2d6sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,5)==8 & trials(:,8)==6;            
            
            for suidx=1:size(fr,2)
                % highest condition % mark switch neurons
                p_early=arrayfun(@(bin) ranksum(fr(s1sel,suidx,bin),fr(s2sel,suidx,bin)),5:7);
                p_late=arrayfun(@(bin) ranksum(fr(s1d6sel,suidx,bin),fr(s2d6sel,suidx,bin)),8:10);
                
                zero_early=arrayfun(@(bin) ~(any(fr(s1sel,suidx,bin)) || any(fr(s2sel,suidx,bin))),5:7);
                zero_late=arrayfun(@(bin) ~(any(fr(s1d6sel,suidx,bin)) || any(fr(s2d6sel,suidx,bin))),8:10);

                p_early(zero_early)=1;
                p_late(zero_late)=1;
                currp=[p_early,p_late];
                
                %optional FDR ajdust
                if opt.fdr
                    [h, ~, ~, currp]=lib.fdr_bh(currp);
                else
                    h=currp<0.05;
                end
                meanfr=[mean(squeeze(fr(s1sel,suidx,5:7))),mean(squeeze(fr(s1d6sel,suidx,8:10)));...
                    mean(squeeze(fr(s2sel,suidx,5:7))),mean(squeeze(fr(s2d6sel,suidx,8:10)))];

                if any(h)
                    currpref=diff(meanfr(:,h));
                    if any(currpref>0) && any(currpref<0) %switch
                        out.wave_id=[out.wave_id;-1];
                    elseif all(currpref>0) %s2 > s1
                        out.wave_id=[out.wave_id;6];
                    else %s1>s2
                        out.wave_id=[out.wave_id;5];
                    end
                else % not odor-selective
                    out.wave_id=[out.wave_id;0];
                end
                out.p_olf=[out.p_olf;currp];
            end
        end
        if  opt.save_file
            wrs_meta=out;
            save('wrs_meta.mat','wrs_meta');
        end
    end
end
opt_=opt;
out_=out;
end
