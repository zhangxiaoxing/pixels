function pivot_dict=delay_vs_iti_motif_per_spike(chain_replay,ring_replay)
arguments
    chain_replay = []
    ring_replay = []
end

pivot_dict.delay=struct();
pivot_dict.iti=struct();
pivot_dict.surround=struct();
% per session
for sess=reshape(unique([ring_replay.session;chain_replay.session]),1,[])
    disp(sess)
    skip=true;
    for onewave=["s1d3","s1d6","s2d3","s2d6"]
        chain_sel=chain_replay.session==sess & chain_replay.wave==onewave;
        ring_sel=ring_replay.session==sess & ring_replay.wave==onewave;
        if nnz(chain_sel)+nnz(ring_sel)<2 % TODO: should be 1 or 2?
            continue
        end
        skip=false;
        % chain ------------------------------------------------
        for chainii=reshape(find(chain_sel),1,[])
            % per preferred trial
            trl_align=chain_replay.trl_align{chainii};
            pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);

            % per su per spk
            for suidx=1:numel(chain_replay.meta{chainii,2})
                suii=chain_replay.meta{chainii,2}(suidx);
                % disp(suii)
                if isfield(pivot_dict.delay,"S"+sess+"U"+suii)
                    for kk=reshape(chain_replay.ts{chainii}(pref_delay,suidx),1,[])
                        if pivot_dict.delay.("S"+sess+"U"+suii).isKey(kk)
                            pivot_dict.delay.("S"+sess+"U"+suii)(kk)=pivot_dict.delay.("S"+sess+"U"+suii)(kk)+1;
                        else
                            pivot_dict.delay.("S"+sess+"U"+suii)(kk)=1;
                        end
                    end
                else
                    pivot_dict.delay.("S"+sess+"U"+suii)=dictionary(chain_replay.ts{chainii}(pref_delay,suidx),1);
                end

                pref_succeed_iti=all(trl_align(:,5:7)==1,2)... % WT, pref
                    & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test/ 1s sample /1s test, 3s response?
                    & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));

                if isfield(pivot_dict.iti,"S"+sess+"U"+suii)
                    for kk=reshape(chain_replay.ts{chainii}(pref_succeed_iti,suidx),1,[])
                        if pivot_dict.iti.("S"+sess+"U"+suii).isKey(kk)
                            pivot_dict.iti.("S"+sess+"U"+suii)(kk)=pivot_dict.iti.("S"+sess+"U"+suii)(kk)+1;
                        else
                            pivot_dict.iti.("S"+sess+"U"+suii)(kk)=1;
                        end
                    end
                else
                    pivot_dict.iti.("S"+sess+"U"+suii)=dictionary(chain_replay.ts{chainii}(pref_succeed_iti,suidx),1);
                end
                %  corresponding network in pre task, post task
                pre_post_motif=(trl_align(:,8)==1 & trl_align(:,9)>60) ...
                    | (trl_align(:,8)<0 & trl_align(:,2)>(60+5+trl_align(:,4)));

                if isfield(pivot_dict.surround,"S"+sess+"U"+suii)
                    for kk=reshape(chain_replay.ts{chainii}(pre_post_motif,suidx),1,[])
                        if pivot_dict.surround.("S"+sess+"U"+suii).isKey(kk)
                            pivot_dict.surround.("S"+sess+"U"+suii)(kk)=pivot_dict.surround.("S"+sess+"U"+suii)(kk)+1;
                        else
                            pivot_dict.surround.("S"+sess+"U"+suii)(kk)=1;
                        end
                    end
                else
                    pivot_dict.surround.("S"+sess+"U"+suii)=dictionary(chain_replay.ts{chainii}(pre_post_motif,suidx),1);
                end
            end
        end

        % loops -------------------------------------------
        for cii=reshape(find(ring_sel),1,[])
            % per preferred trial
            trl_align=ring_replay.trl_align{cii};
            pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
            delay_cov=cell2mat(cellfun(@(x) ceil(x([1,end])./3).',ring_replay.ts{cii}(pref_delay),'UniformOutput',false));
            for rii=1:size(delay_cov,1)
                covered.delay(delay_cov(rii,1):delay_cov(rii,2))=true;
            end

            pref_succeed_iti=all(trl_align(:,5:7)==1,2)... % WT, pref
                & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test/ 1s sample /1s test, 3s response?
                & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));

            iti_cov=cell2mat(cellfun(@(x) ceil(x([1,end])./3).',ring_replay.ts{cii}(pref_succeed_iti),'UniformOutput',false));
            for rii=1:size(iti_cov,1)
                covered.iti(iti_cov(rii,1):iti_cov(rii,2))=true;
            end
            % corresponding network in pre task, post task
            pre_post_motif=(trl_align(:,8)==1 & trl_align(:,9)>60) ...
                | (trl_align(:,8)<0 & trl_align(:,2)>(60+5+trl_align(:,4)));
            pre_post_cov=cell2mat(cellfun(@(x) ceil(x([1,end])./3).',ring_replay.ts{cii}(pre_post_motif),'UniformOutput',false));
            for rii=1:size(pre_post_cov,1)
                covered.pre_post(pre_post_cov(rii,1):pre_post_cov(rii,2))=true;
            end
        end
    end
    if skip
        continue
    end
end
end

