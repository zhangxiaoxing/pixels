function pivot_dict=delay_vs_iti_motif_per_spike(chain_replay,ring_replay,opt)
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

pivot_dict.delay=struct();
pivot_dict.iti=struct();
pivot_dict.out_task=struct();
% per session
for sess=reshape(unique([ring_replay.session;chain_replay.session]),1,[])
    disp(sess)

    for onewave=["s1d3","s1d6","s2d3","s2d6"]
        chain_sel=chain_replay.session==sess & chain_replay.wave==onewave;
        ring_sel=ring_replay.session==sess & ring_replay.wave==onewave;
        if nnz(chain_sel)+nnz(ring_sel)<2 % TODO: should be 1 or 2?
            continue
        end
        % chain ------------------------------------------------
        for chainii=reshape(find(chain_sel),1,[])
            % per preferred trial
            trl_align=chain_replay.trl_align{chainii};
            pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
            pref_succeed_iti=all(trl_align(:,5:7)==1,2)... % WT, pref
                & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test/ 1s sample /1s test, 3s response?
                & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));
            pre_post_motif=(trl_align(:,8)==1 & trl_align(:,9)>60) ...
                | (trl_align(:,8)<0 & trl_align(:,2)>(60+5+trl_align(:,4)));

            % per su per spk
            for suidx=1:numel(chain_replay.meta{chainii,2})
                suii=chain_replay.meta{chainii,2}(suidx);
                % disp(suii)
                if isfield(pivot_dict.delay,"ses"+sess+onewave+"U"+suii)
                    for kk=reshape(chain_replay.ts{chainii}(pref_delay,suidx),1,[])
                        if pivot_dict.delay.("ses"+sess+onewave+"U"+suii).isKey(kk)
                            pivot_dict.delay.("ses"+sess+onewave+"U"+suii)(kk)=pivot_dict.delay.("ses"+sess+onewave+"U"+suii)(kk)+1;
                        else
                            pivot_dict.delay.("ses"+sess+onewave+"U"+suii)(kk)=1;
                        end
                    end
                else
                    pivot_dict.delay.("ses"+sess+onewave+"U"+suii)=dictionary(chain_replay.ts{chainii}(pref_delay,suidx),1);
                end

                if isfield(pivot_dict.iti,"ses"+sess+onewave+"U"+suii)
                    for kk=reshape(chain_replay.ts{chainii}(pref_succeed_iti,suidx),1,[])
                        if pivot_dict.iti.("ses"+sess+onewave+"U"+suii).isKey(kk)
                            pivot_dict.iti.("ses"+sess+onewave+"U"+suii)(kk)=pivot_dict.iti.("ses"+sess+onewave+"U"+suii)(kk)+1;
                        else
                            pivot_dict.iti.("ses"+sess+onewave+"U"+suii)(kk)=1;
                        end
                    end
                else
                    pivot_dict.iti.("ses"+sess+onewave+"U"+suii)=dictionary(chain_replay.ts{chainii}(pref_succeed_iti,suidx),1);
                end
                %  corresponding network in pre task, post task

                if isfield(pivot_dict.out_task,"ses"+sess+onewave+"U"+suii)
                    for kk=reshape(chain_replay.ts{chainii}(pre_post_motif,suidx),1,[])
                        if pivot_dict.out_task.("ses"+sess+onewave+"U"+suii).isKey(kk)
                            pivot_dict.out_task.("ses"+sess+onewave+"U"+suii)(kk)=pivot_dict.out_task.("ses"+sess+onewave+"U"+suii)(kk)+1;
                        else
                            pivot_dict.out_task.("ses"+sess+onewave+"U"+suii)(kk)=1;
                        end
                    end
                else
                    pivot_dict.out_task.("ses"+sess+onewave+"U"+suii)=dictionary(chain_replay.ts{chainii}(pre_post_motif,suidx),1);
                end
            end
        end

        % loops -------------------------------------------
        for rii=reshape(find(ring_sel),1,[])
            % per preferred trial
            trl_align=ring_replay.trl_align{rii};
            pref_delay=all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
            pref_succeed_iti=all(trl_align(:,5:7)==1,2)... % WT, pref
                & trl_align(:,2)>=(trl_align(:,4)+5)...  % not in delay or test/ 1s sample /1s test, 3s response?
                & (trl_align(:,8)>0|(trl_align(:,8)==-1 & trl_align(:,2)<trl_align(:,4)+1+14));
            pre_post_motif=(trl_align(:,8)==1 & trl_align(:,9)>60) ...
                | (trl_align(:,8)<0 & trl_align(:,2)>(60+5+trl_align(:,4)));

            % per su per spk
            for suidx=1:numel(ring_replay.meta{rii,2})
                suii=ring_replay.meta{rii,2}(suidx);

                spkts=[];
                for occurii=reshape(find(pref_delay),1,[])
                    spkts=[spkts;ring_replay.ts{rii}{occurii}(ring_replay.ts_seq{rii}{occurii}==suii)];
                end

                if isfield(pivot_dict.delay,"ses"+sess+onewave+"U"+suii)
                    for kk=spkts.'
                        if pivot_dict.delay.("ses"+sess+onewave+"U"+suii).isKey(kk)
                            pivot_dict.delay.("ses"+sess+onewave+"U"+suii)(kk)=pivot_dict.delay.("ses"+sess+onewave+"U"+suii)(kk)+1;
                        else
                            pivot_dict.delay.("ses"+sess+onewave+"U"+suii)(kk)=1;
                        end
                    end
                else
                    pivot_dict.delay.("ses"+sess+onewave+"U"+suii)=dictionary(ring_replay.ts{rii}(pref_delay,suidx),1);
                end

                spkts=[];
                for occurii=reshape(find(pref_succeed_iti),1,[])
                    spkts=[spkts;ring_replay.ts{rii}{occurii}(ring_replay.ts_seq{rii}{occurii}==suii)];
                end

                if isfield(pivot_dict.iti,"ses"+sess+onewave+"U"+suii)
                    for kk=spkts.'
                        if pivot_dict.iti.("ses"+sess+onewave+"U"+suii).isKey(kk)
                            pivot_dict.iti.("ses"+sess+onewave+"U"+suii)(kk)=pivot_dict.iti.("ses"+sess+onewave+"U"+suii)(kk)+1;
                        else
                            pivot_dict.iti.("ses"+sess+onewave+"U"+suii)(kk)=1;
                        end
                    end
                else
                    pivot_dict.iti.("ses"+sess+onewave+"U"+suii)=dictionary(ring_replay.ts{rii}(pref_succeed_iti,suidx),1);
                end
                %  corresponding network in pre task, post task

                spkts=[];
                for occurii=reshape(find(pre_post_motif),1,[])
                    spkts=[spkts;ring_replay.ts{rii}{occurii}(ring_replay.ts_seq{rii}{occurii}==suii)];
                end

                if isfield(pivot_dict.out_task,"ses"+sess+onewave+"U"+suii)
                    for kk=spkts.'
                        if pivot_dict.out_task.("ses"+sess+onewave+"U"+suii).isKey(kk)
                            pivot_dict.out_task.("ses"+sess+onewave+"U"+suii)(kk)=pivot_dict.out_task.("ses"+sess+onewave+"U"+suii)(kk)+1;
                        else
                            pivot_dict.out_task.("ses"+sess+onewave+"U"+suii)(kk)=1;
                        end
                    end
                else
                    pivot_dict.out_task.("ses"+sess+onewave+"U"+suii)=dictionary(ring_replay.ts{rii}(pre_post_motif,suidx),1);
                end
            end
        end
    end
end

%% TODO: SAVE FILE
if ~opt.skip_save
    blame=vcs.blame();
    save(fullfile('binary','delay_vs_iti_pivot_spike.mat'),'pivot_dict','blame')
end

end

