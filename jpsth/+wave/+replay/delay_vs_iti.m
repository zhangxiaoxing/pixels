function [run_length,covered_tbl]=delay_vs_iti(chain_replay,ring_replay,opt)
arguments
    chain_replay = []
    ring_replay = []
    opt.skip_save (1,1) logical = false
end

if isempty(ring_replay)
    load(fullfile('binary','motif_replay.mat'),'ring_replay');
end

if isempty(chain_replay)
    load(fullfile('binary','motif_replay.mat'),'chain_replay');
end

% per session
run_length=cell2struct({[];[];[];[];[]},{'delay','iti','out_task','npdelay','npiti'});
covered_tbl=[];
for sess=reshape(unique([ring_replay.session;chain_replay.session]),1,[])
    disp(sess)
    session_tick=wave.replay.sessid2length(sess);
    for onewave=["s1d3","s1d6","s2d3","s2d6"]
        switch onewave
            case "s1d3"
                [samp,delay]=deal(4,3);
            case "s1d6"
                [samp,delay]=deal(4,6);
            case "s2d3"
                [samp,delay]=deal(8,3);
            case "s2d6"
                [samp,delay]=deal(8,6);
        end

        for fn=["iti","npiti","delay","npdelay","out_task"]
            covered.(fn)=false(ceil(session_tick/3),1);
        end

        chain_sel=chain_replay.session==sess & chain_replay.wave==onewave;
        ring_sel=ring_replay.session==sess & ring_replay.wave==onewave;
        if nnz(chain_sel)+nnz(ring_sel)<2
            continue
        end

        % chain ------------------------------------------------
        for cii=reshape(find(chain_sel),1,[])
            
            trl_align=chain_replay.trl_align{cii};
            % per preferred trial
            pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
            delay_cov=ceil(chain_replay.ts{cii}(pref_delay,[1 end])./3);
            for rii=1:size(delay_cov,1)
                covered.delay(delay_cov(rii,1):delay_cov(rii,2))=true;
            end

            nonpref_delay=all(trl_align(:,5:6)==1,2) & trl_align(:,3)~=samp & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
            nonpref_delay_cov=ceil(chain_replay.ts{cii}(nonpref_delay,[1 end])./3);
            for rii=1:size(nonpref_delay_cov,1)
                covered.npdelay(nonpref_delay_cov(rii,1):nonpref_delay_cov(rii,2))=true;
            end

            pref_succeed_iti=all(trl_align(:,5:7)==1,2)... % WT, pref
                & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test/ 1s sample /1s test, 3s response?
                & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));
            iti_cov=ceil(chain_replay.ts{cii}(pref_succeed_iti,[1 end])./3);
            for rii=1:size(iti_cov,1)
                covered.iti(iti_cov(rii,1):iti_cov(rii,2))=true;
            end

            nonpref_iti=all(trl_align(:,5:6)==1,2) & trl_align(:,3)~=samp...
                & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test/ 1s sample /1s test, 3s response?
                & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));
            npiti_cov=ceil(chain_replay.ts{cii}(nonpref_iti,[1 end])./3);
            for rii=1:size(npiti_cov,1)
                covered.npiti(npiti_cov(rii,1):npiti_cov(rii,2))=true;
            end


            %  corresponding network in pre task, post task
            out_task_motif=(trl_align(:,8)==1 & trl_align(:,9)>60) ...
                | (trl_align(:,8)<0 & trl_align(:,2)>(60+5+trl_align(:,4)));
            out_task_cov=ceil(chain_replay.ts{cii}(out_task_motif,[1 end])./3);
            for rii=1:size(out_task_cov,1)
                covered.out_task(out_task_cov(rii,1):out_task_cov(rii,2))=true;
            end

        end

        % loops -------------------------------------------
        for cii=reshape(find(ring_sel),1,[])
            
            trl_align=ring_replay.trl_align{cii};

            % per preferred trial
            pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
            delay_cov=cell2mat(cellfun(@(x) ceil(x([1,end])./3).',ring_replay.ts{cii}(pref_delay),'UniformOutput',false));
            for rii=1:size(delay_cov,1)
                covered.delay(delay_cov(rii,1):delay_cov(rii,2))=true;
            end

            nonpref_delay=all(trl_align(:,5:6)==1,2) & trl_align(:,3)~=samp & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
            nonpref_delay_cov=cell2mat(cellfun(@(x) ceil(x([1,end])./3).',ring_replay.ts{cii}(nonpref_delay),'UniformOutput',false));
            for rii=1:size(nonpref_delay_cov,1)
                covered.npdelay(nonpref_delay_cov(rii,1):nonpref_delay_cov(rii,2))=true;
            end

            pref_succeed_iti=all(trl_align(:,5:7)==1,2)... % WT, pref
                & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test/ 1s sample /1s test, 3s response?
                & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));
            iti_cov=cell2mat(cellfun(@(x) ceil(x([1,end])./3).',ring_replay.ts{cii}(pref_succeed_iti),'UniformOutput',false));
            for rii=1:size(iti_cov,1)
                covered.iti(iti_cov(rii,1):iti_cov(rii,2))=true;
            end

            nonpref_iti=all(trl_align(:,5:6)==1,2) & trl_align(:,3)~=samp...
                & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test/ 1s sample /1s test, 3s response?
                & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));
            npiti_cov=cell2mat(cellfun(@(x) ceil(x([1,end])./3).',ring_replay.ts{cii}(nonpref_iti),'UniformOutput',false));
            for rii=1:size(npiti_cov,1)
                covered.npiti(npiti_cov(rii,1):npiti_cov(rii,2))=true;
            end

            % corresponding network in pre task, post task
            out_task_motif=(trl_align(:,8)==1 & trl_align(:,9)>60) ...
                | (trl_align(:,8)<0 & trl_align(:,2)>(60+5+trl_align(:,4)));
            out_task_cov=cell2mat(cellfun(@(x) ceil(x([1,end])./3).',ring_replay.ts{cii}(out_task_motif),'UniformOutput',false));
            for rii=1:size(out_task_cov,1)
                covered.out_task(out_task_cov(rii,1):out_task_cov(rii,2))=true;
            end
        end
        covered_tbl=[covered_tbl;cell2table({sess,onewave,sparse(covered.delay),sparse(covered.iti),sparse(covered.out_task),sparse(covered.npdelay),sparse(covered.npiti)},'VariableNames',{'session','wave','delay','iti','out_task','npdelay','npiti'})];

        [delay_run_len,~,~]=wave.replay.covered2runlength(covered.delay);
        [iti_run_len,~,~]=wave.replay.covered2runlength(covered.iti);
        [out_task_run_len,~,~]=wave.replay.covered2runlength(covered.out_task);
        [npdelay_run_len,~,~]=wave.replay.covered2runlength(covered.npdelay);
        [npiti_run_len,~,~]=wave.replay.covered2runlength(covered.npiti);

        run_length.delay=[run_length.delay;repmat(double(sess),numel(delay_run_len),1),delay_run_len];
        run_length.iti=[run_length.iti;repmat(double(sess),numel(iti_run_len),1),iti_run_len];
        run_length.npdelay=[run_length.npdelay;repmat(double(sess),numel(npdelay_run_len),1),npdelay_run_len];
        run_length.npiti=[run_length.npiti;repmat(double(sess),numel(npiti_run_len),1),npiti_run_len];
        run_length.out_task=[run_length.out_task;repmat(double(sess),numel(out_task_run_len),1),out_task_run_len];
    end
end
if ~opt.skip_save
    blame=vcs.blame();
    save(fullfile('binary','delay_iti_runlength_covered.mat'),'run_length','covered_tbl','blame');
end
end

