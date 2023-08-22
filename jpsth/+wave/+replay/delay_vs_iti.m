function [run_length,covered_tbl]=delay_vs_iti(chain_replay,ring_replay_tbl)
% load(fullfile('binary','rings_tag.mat'))
% [ring_replay,ring_stats,~]=wave.replay.stats(rmfield(rings_tag,"none"),'var_len',true);
% ring_replay_tbl=wave.replay.quickconvert(ring_replay);
%
% fstr=load(fullfile('binary','chain_tag_all_trl.mat'),'out','trials_dict');
% [chain_replay,chain_stats,~]=wave.replay.stats_tbl(fstr.out,fstr.trials_dict,'var_len',false);


% per session
run_length=cell2struct({[];[];[]},{'delay','iti','pre_post'});
covered_tbl=[];
for sess=reshape(unique([ring_replay_tbl.session;chain_replay.session]),1,[])
    disp(sess)
    session_tick=wave.replay.sessid2length(sess);
    covered.iti=false(ceil(session_tick/3),1);
    covered.delay=false(ceil(session_tick/3),1);
    covered.pre_post=false(ceil(session_tick/3),1);
    skip=true;
    for onewave=["s1d3","s1d6","s2d3","s2d6"]

        chain_sel=chain_replay.session==sess & chain_replay.wave==onewave;
        ring_sel=ring_replay_tbl.session==sess & ring_replay_tbl.wave==onewave;
        if nnz(chain_sel)+nnz(ring_sel)<2
            continue
        end
        skip=false;
        % chain ------------------------------------------------
        for cii=reshape(find(chain_sel),1,[])
            % per preferred trial
            trl_align=chain_replay.trl_align{cii};
            pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
            delay_cov=ceil(chain_replay.ts{cii}(pref_delay,[1 end])./3);
            for rii=1:size(delay_cov,1)
                covered.delay(delay_cov(rii,1):delay_cov(rii,2))=true;
            end

            pref_succeed_iti=all(trl_align(:,5:7)==1,2)... % WT, pref
                & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test/ 1s sample /1s test, 3s response?
                & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));

            iti_cov=ceil(chain_replay.ts{cii}(pref_succeed_iti,[1 end])./3);
            for rii=1:size(iti_cov,1)
                covered.iti(iti_cov(rii,1):iti_cov(rii,2))=true;
            end
            %  corresponding network in pre task, post task
            pre_post_motif=(trl_align(:,8)==1 & trl_align(:,9)>60) ...
                | (trl_align(:,8)<0 & trl_align(:,2)>(60+5+trl_align(:,4)));
            pre_post_cov=ceil(chain_replay.ts{cii}(pre_post_motif,[1 end])./3);
            for rii=1:size(pre_post_cov,1)
                covered.pre_post(pre_post_cov(rii,1):pre_post_cov(rii,2))=true;
            end

        end

        % loops -------------------------------------------
        for cii=reshape(find(ring_sel),1,[])
            % per preferred trial
            trl_align=ring_replay_tbl.trl_align{cii};
            pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
            delay_cov=cell2mat(cellfun(@(x) ceil(x([1,end])./3).',ring_replay_tbl.ts{cii}(pref_delay),'UniformOutput',false));
            for rii=1:size(delay_cov,1)
                covered.delay(delay_cov(rii,1):delay_cov(rii,2))=true;
            end

            pref_succeed_iti=all(trl_align(:,5:7)==1,2)... % WT, pref
                & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test/ 1s sample /1s test, 3s response?
                & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));

            iti_cov=cell2mat(cellfun(@(x) ceil(x([1,end])./3).',ring_replay_tbl.ts{cii}(pref_succeed_iti),'UniformOutput',false));
            for rii=1:size(iti_cov,1)
                covered.iti(iti_cov(rii,1):iti_cov(rii,2))=true;
            end
            % corresponding network in pre task, post task
            pre_post_motif=(trl_align(:,8)==1 & trl_align(:,9)>60) ...
                | (trl_align(:,8)<0 & trl_align(:,2)>(60+5+trl_align(:,4)));
            pre_post_cov=cell2mat(cellfun(@(x) ceil(x([1,end])./3).',ring_replay_tbl.ts{cii}(pre_post_motif),'UniformOutput',false));
            for rii=1:size(pre_post_cov,1)
                covered.pre_post(pre_post_cov(rii,1):pre_post_cov(rii,2))=true;
            end
        end
    end
    if skip
        continue
    end
    [delay_run_len,doffset,donset]=wave.replay.covered2runlength(covered.delay);
    [iti_run_len,ioffset,ionset]=wave.replay.covered2runlength(covered.iti);
    [pre_post_run_len,poffset,ponset]=wave.replay.covered2runlength(covered.pre_post);

    covered_tbl=[covered_tbl;cell2table({sess,covered.delay,covered.iti,covered.pre_post},'VariableNames',{'session','delay','iti','out_task'})];

    run_length.delay=[run_length.delay;repmat(double(sess),numel(delay_run_len),1),delay_run_len];
    run_length.iti=[run_length.iti;repmat(double(sess),numel(iti_run_len),1),iti_run_len];
    run_length.pre_post=[run_length.pre_post;repmat(double(sess),numel(pre_post_run_len),1),pre_post_run_len];
end
end

