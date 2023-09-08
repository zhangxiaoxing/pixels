% [sschain_trl,unfound]=wave.chain_tag(chains_uf,'skip_save',true,'len_thresh',len_thresh,'odor_only',true,'extend_trial',true,'skip_ts_id',true,'DEBUG',true); % per-spk association
%
% TODO : CCG

classdef chain_tag < handle
    methods(Static)
        function [out,notfound,trials_dict]=tag(chains,len_thresh,opt)
            arguments
                chains
                len_thresh (1,1) double
                opt.ccg (1,1) logical = false
                opt.rev (1,1) logical = false
                opt.shuf_trl (1,1) logical = false
                opt.shuf_idx (1,1) double = 0
                opt.per_reg_wave (1,1) logical = true
                opt.skip_save (1,1) logical = false
                opt.odor_only (1,1) logical = false
                opt.extend_trial (1,1) logical = false % include all recording time or only preferred delay
                opt.skip_ts_id (1,1) logical = false
                % opt.DEBUG (1,1) logical = false
                opt.anti_dir (1,1) logical = false % TODO: consistent chain,reverse direction
                opt.sesses double = []
                opt.idces double = []
                opt.load_file (1,1) logical = false % reserved but not used
                opt.filename (1,:) char = 'chain_tag.mat'
                opt.poolsize (1,1) double = 2
            end

            assert(~opt.load_file,"Unfinished yet")

            if opt.shuf_trl
                assert(opt.shuf_idx>0,"index is necessary for shuffled data")
            end

            if opt.ccg
                load(fullfile("binary","sums_conn_10.mat"),'sums_conn_str');
            end

            ch_len=cellfun(@(x) numel(x),chains.cids);
            if isempty(opt.sesses)
                waveids=reshape(unique(chains.wave(ch_len>=len_thresh)),1,[]);
                sesses=reshape(unique(chains.sess(ch_len>=len_thresh)),1,[]);
            else
                sesses=opt.sesses;
                waveids=reshape(unique(chains.wave(ch_len>=len_thresh & ismember(chains.sess,opt.sesses))),1,[]);
            end

            if any(contains(waveids,'olf'))
                keyboard()
                % processed=cell2struct({[];[]},{'d3','d6'});
                % TODO: if both sxdx & olf_dx:
                % sess_indices=setdiff(sess_indices,processed.("d"+duration));
                % if isempty(sess_indices), continue;end
                % processed.("d"+duration)=[processed.("d"+duration),reshape(sess_indices,1,[])];
            end
            if opt.poolsize>1
                poolh=parpool(opt.poolsize);
            end
            F=parallel.FevalFuture.empty(0,1);
            for sessid=sesses %TODO: threads
                % if opt.DEBUG && sessid>33
                %     break
                % end
                for duration=[3 6]
                    F(end+1)=parfeval(poolh,@wave.chain_tag.per_session_func,2,chains,len_thresh,ch_len,sessid,duration,waveids,opt);
                    % [sess_out,sess_notfound]=wave.chain_tag.per_session_func(chains,len_thresh,ch_len,sessid,duration,waveids,opt);
                end
            end

            fini=[];
            notfound=cell(0,3);
            trials_dict=dictionary([],cell(0));
            while ~all([F.Read],'all')
                try
                    [ftidx,sess_out,sess_notfound]=fetchNext(F);
                    notfound=[notfound;sess_notfound];
            		if isempty(sess_out)
            			continue
            		end
                    if ~exist('out','var')
                        out=sess_out;
                    else
            		    out=[out;sess_out];
                    end

                    fini=[fini,ftidx];
                    disp(ftidx);

                    if ~opt.skip_save
                        blame=vcs.blame();
                        if opt.rev
                            save(fullfile('bzdata','chain_rev_tag.mat'),'out','notfound','blame');
                        elseif opt.shuf_trl
                            save(fullfile("bzdata","chain_shuf_tag_"+opt.shuf_idx+".mat"),'out','notfound','blame');
                        else
                            save(fullfile('binary',opt.filename),'out','notfound','blame','-v7.3')
                        end
                    end
                catch ME
                    fstate={F.State};
                    save(fullfile('binary/tagstate.mat'),'fstate','fini','ME');
                end


            end
            % if opt.ccg
            %     for dur=reshape(fieldnames(out),1,[])
            %         for wv=reshape(fieldnames(out.(dur{1})),1,[])
            %             for lp=reshape(fieldnames(out.(dur{1}).(wv{1})),1,[])
            %                 ccgs=out.(dur{1}).(wv{1}).(lp{1}).ccgs;
            %                 ccg_qual=[];
            %                 for ii=1:size(ccgs,1)
            %                     ccg_qual=[ccg_qual;bz.good_ccg(ccgs(ii,:))];
            %                     %1:Polarity 2:Time of peak 3:Noise peaks 4:FWHM 5:rising
            %                     %edge 6:falling edge
            %                 end
            %                 out.(dur{1}).(wv{1}).(lp{1}).ccg_qual=ccg_qual;
            %             end
            %         end
            %     end
            % end
            delete(poolh);
        end
        %

        function [sess_out,sess_notfound]=per_session_func(chains,len_thresh,ch_len,sessid, duration,waveids,opt)
    	    [sess_out,sess_notfound]=deal([]);
            load(fullfile('binary','trials_dict.mat'),'trials_dict');
            trials=cell2mat(trials_dict(sessid));
            for wid=waveids
                if (contains(wid,'d3') && duration==6) ...
                        || (contains(wid,'d6') && duration ==3)
                    continue
                end
                if contains(wid,'s1')
                    trial_sel=find(trials(:,5)==4 & trials(:,8)==duration & all(trials(:,9:10)>0,2));
                    % outid="olf_s1";
                elseif contains(wid,'s2')
                    trial_sel=find(trials(:,5)==8 & trials(:,8)==duration & all(trials(:,9:10)>0,2));
                    % outid="olf_s2";
                else
                    if opt.odor_only && contains(wid,'dur')
                        continue
                    else % nonmem, dur
                        trial_sel=find(trials(:,8)==duration & all(trials(:,9:10)>0,2));
                    end
                end

                outid=wid;

                if isempty(opt.idces)
                    sess_indices=reshape(find(chains.sess==sessid & strcmp(chains.wave,wid) & ch_len>=len_thresh),1,[]);
                else
                    sess_indices=opt.idces;
                end
                if ~exist('FT_SPIKE','var')
                    [spkID,spkTS,~,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true,'jagged',true);
                end
                for cc=sess_indices
                    ts_id=[];
                    cids=chains.cids{cc};
                    %                 disp({sessid,ri});
                    for in_chain_pos=1:numel(cids) % TODO, 1:chain_len
                        one_chain_sel=spkID==cids(in_chain_pos);
                        rawts=spkTS(one_chain_sel);

                        ft_sel=strcmp(FT_SPIKE.label,num2str(cids(in_chain_pos)));
                        ft_ts=FT_SPIKE.timestamp{ft_sel};
                        ft_trl_time=FT_SPIKE.time{ft_sel};
                        ft_trl=FT_SPIKE.trial{ft_sel};

                        [~,tspos]=ismember(ft_ts,rawts);
                        ext_time=repmat(-realmax,numel(rawts),1);
                        ext_time(tspos)=ft_trl_time;

                        ext_trl=repmat(-realmax,numel(rawts),1);
                        ext_trl(tspos)=ft_trl;

                        ts_id=cat(1,ts_id,[rawts,... % 1
                            repmat(cids(in_chain_pos),numel(rawts),1),...  % 2
                            ones(numel(rawts),1)*in_chain_pos,...  % 3
                            ext_time,...  % 4
                            ext_trl]); % 5
                    end
                    % optional remove non-wave spikes
                    if ~opt.extend_trial
                        ts_id=ts_id(ts_id(:,4)>=1 & ts_id(:,4)<(duration+1) & ismember(ts_id(:,5),trial_sel),:);
                    end
                    if opt.shuf_trl
                        ts_id_shuf=[];
                        for in_chain_pos=1:numel(cids)
                            nodesel=ts_id(:,3)==in_chain_pos;
                            trl_onset=arrayfun(@(x) trials(x,1),ts_id(nodesel,5));
                            delta_ts=ts_id(nodesel,1)-trl_onset;

                            shufmap=containers.Map(num2cell(trial_sel),num2cell(randsample(trial_sel,numel(trial_sel))));
                            shuf_trl=cell2mat(shufmap.values(num2cell(ts_id(nodesel,5))));
                            shuf_onset=arrayfun(@(x) trials(x,1),shuf_trl);
                            ts_id_shuf=[ts_id_shuf;sortrows([shuf_onset+delta_ts,ts_id(nodesel,2:end)],1)];
                        end

                        ts=wave.chain_tag.chain_alt(ts_id_shuf(:,[1 3]));
                    else
                        if opt.anti_dir
                            anti_ts_id=ts_id(:,[1 3]);
                            anti_ts_id(:,2)=(numel(cids)+1)-anti_ts_id(:,2);
                            ts=wave.chain_tag.chain_alt(anti_ts_id);
                        else
                            ts=wave.chain_tag.chain_alt(ts_id(:,[1 3]));
                        end
                    end
                    if size(ts,1)<2
                        ts={ts};
                    end
                    if ~isempty(ts) && ~(iscell(ts) && isempty(ts{1}))
                        outkey="s"+sessid+"c"+cc;
                        if opt.skip_ts_id
                            meta={sessid,cids,chains.cross_reg(cc)};
                            outcell={sessid,duration,outid,outkey,ts,meta,[]};
                        else
                            meta={cids,chains.tcoms(cc),sessid,chains.cross_reg(cc)};
                            ts_id=table(uint32(ts_id(:,1)),uint16(ts_id(:,2)),uint8(ts_id(:,3)),ts_id(:,4),uint16(ts_id(:,5)),'VariableNames',{'TS','CID','POS','Time','Trial'});
                            outcell={sessid,duration,outid,outkey,ts,meta,ts_id};
                        end

                        if ~isempty(sess_out)
                            sess_out=[sess_out;cell2table(outcell,"VariableNames",{ ...
                                'session','delay','wave','chain_id','ts','meta','ts_id'})];
                        else
                            sess_out=cell2table(outcell,"VariableNames",{ ...
                                'session','delay','wave','chain_id','ts','meta','ts_id'});
                        end
                        % CCG
                        if opt.ccg
                            cursess=chains.sess(cc);
                            sesspath=ephys.sessid2path(cursess);
                            strippath=regexp(sesspath,'(?<=\\).*','match');
                            sesssel=find(contains({sums_conn_str.folder},strippath));
                            ccg=sums_conn_str(sesssel).ccg_sc;
                            ccgid=sums_conn_str(sesssel).sig_con;
                            chainccg=[];
                            for ii=1:numel(cids)-1
                                sigsel=ccgid(:,1)==cids(ii) & ccgid(:,2)==cids(ii+1);
                                if nnz(sigsel)~=1
                                    keyboard()
                                end
                                chainccg=[chainccg;ccg(sigsel,:)];
                            end
                            sess_out.("d"+duration).(outid).(outkey).ccgs=chainccg;
                        end
                    else
                        sess_notfound=[sess_notfound;{sessid},{cc},{cids}];
                    end
                end
            end
        end

        function out=chain_alt(in)
            out=[];
            clen=max(in(:,2));
            persu=struct.empty(clen,0);
            for ii=1:clen
                persu(ii).tsid=in(in(:,2)==ii,:);
                persu(ii).idx=1;
            end

            curd=1;
            while true
                nxtd=curd+1;
                if persu(nxtd).idx>size(persu(nxtd).tsid,1) || persu(curd).idx>size(persu(curd).tsid,1)
                    break
                end
                if persu(nxtd).tsid(persu(nxtd).idx,1)<persu(curd).tsid(persu(curd).idx,1)+24 % ahead of window
                    persu(nxtd).idx=persu(nxtd).idx+1;
                elseif persu(nxtd).tsid(persu(nxtd).idx,1)>persu(curd).tsid(persu(curd).idx,1)+300 % beyond window
                    persu(curd).idx=persu(curd).idx+1;
                    if curd~=1
                        curd=curd-1;
                    end
                else % in window
                    if nxtd==clen
                        out=[out;arrayfun(@(x) persu(x).tsid(persu(x).idx,1),1:clen)];
                        persu(nxtd).idx=persu(nxtd).idx+1;
                    else
                        curd=nxtd;
                    end
                end
            end
        end
    end
end
