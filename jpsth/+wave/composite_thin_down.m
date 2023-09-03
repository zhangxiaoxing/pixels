% tag spikes in neuron with loop or chain activity
% olf, enc-both disconnect, composite proportion, bargraph
% TODO: CONTROL
classdef composite_thin_down < handle

    methods (Static)

        function per_sess_condition=merge_motif(sschain,ssloop)
            arguments
                sschain
                ssloop
            end
            usess=union(sschain.session,ssloop.session);
            %% per-session entrance
            for sessid=reshape(usess,1,[])
                for waveid=["s1d3","s1d6","s2d3","s2d6"]
                per_sess_condition.("s"+sessid+waveid)=struct();
                per_sess_condition.("s"+sessid+waveid).chains=sschain.meta(sschain.session==sessid & sschain.wave==waveid,2);
                per_sess_condition.("s"+sessid+waveid).loops=ssloop.meta(ssloop.session==sessid & ssloop.wave==waveid,2);
                end
            end
        end

        function [stats,G]=stats_remove(per_sess_condition,opt)
            arguments
                per_sess_condition
                opt.remove (1,:) char {mustBeMember(opt.remove, {'','chains','loops','UID','D5','ChainMatchFC','LoopMatchFC'})} = ''
                opt.one_net (1,:) char = ''
                opt.keep_all_subs (1,1) logical = false % for illustration only
            end

            stats=struct();
            % per_trl_nodes=cell(0);
            % complexity_sums=[];
            % degree_sums=[];

            %% systematic ablation

            switch opt.remove
                case 'chains'
                    fns=fieldnames(per_sess_condition);
                    strcell=struct2cell(per_sess_condition);
                    cellstr=cellfun(@(x) cell2struct({cell(0);x.loops},{'chains','loops'}),strcell,'UniformOutput',false);
                    per_sess_condition=cell2struct(cellstr,fns);

                case 'loops'
                    fns=fieldnames(per_sess_condition);
                    strcell=struct2cell(per_sess_condition);
                    cellstr=cellfun(@(x) cell2struct({x.chains;cell(0)},{'chains','loops'}),strcell,'UniformOutput',false);
                    per_sess_condition=cell2struct(cellstr,fns);

            end

            if isempty(opt.one_net)
                fns=fieldnames(per_sess_condition);
                G=[];
            else
                fns={opt.one_net};
            end

            for fn=reshape(fns,1,[])

                if (numel(per_sess_condition.(fn{1}).chains)+numel(per_sess_condition.(fn{1}).loops)<=1 && ~opt.keep_all_subs) ...
                        || (numel(per_sess_condition.(fn{1}).chains)+numel(per_sess_condition.(fn{1}).loops)<1 && opt.keep_all_subs)
                    continue
                end

                % build graph network
                per_sess_condition.(fn{1}).loops=cellfun(@(x) x([1:end,1]),per_sess_condition.(fn{1}).loops,'UniformOutput',false);% cyclic loops fix
                edges=categorical(unique(cell2mat(cellfun(@(x) [x(1:end-1);x(2:end)].',...
                    [per_sess_condition.(fn{1}).chains;per_sess_condition.(fn{1}).loops],...
                    'UniformOutput',false)),'rows'));
                if strcmp(opt.remove,'UID')
                    assert(~isempty(opt.UID), "missing UID input");
                    % else
                    %     error("Unfinished options")
                end


                if strcmp(opt.remove,'D10')
                    gh=graph(edges(:,1),edges(:,2));
                    gh=gh.rmnode(gh.Nodes.Name(gh.degree>=10));
                elseif strcmp(opt.remove,'D5')
                    gh=graph(edges(:,1),edges(:,2));
                    gh=gh.rmnode(gh.Nodes.Name(gh.degree>=5));
                elseif strcmp(opt.remove,'ChainMatchFC')
                    rmvedges=unique(cell2mat(cellfun(@(x) [x(1:end-1);x(2:end)].',...
                        per_sess_condition.(fn{1}).chains,'UniformOutput',false)),'rows');
                    rmv_fc_n=size(rmvedges,1);
                    keepsel=randsample(size(edges,1),size(edges,1)-rmv_fc_n);
                    gh=graph(edges(keepsel,1),edges(keepsel,2));
                elseif strcmp(opt.remove,'LoopMatchFC')
                    rmvedges=unique(cell2mat(cellfun(@(x) [x(1:end-1);x(2:end)].',...
                        per_sess_condition.(fn{1}).loops,'UniformOutput',false)),'rows');
                    rmv_fc_n=size(rmvedges,1);
                    keepsel=randsample(size(edges,1),size(edges,1)-rmv_fc_n);
                    gh=graph(edges(keepsel,1),edges(keepsel,2));
                else
                    gh=graph(edges(:,1),edges(:,2));
                end

                conncomp=gh.conncomp();
                % check module in network
                if any(conncomp~=1)
                    comps=unique(conncomp);
                    counter=zeros(numel(comps),1);
                    gnodes=cellfun(@(x) str2double(x),gh.Nodes.Name);
                    for mm=[per_sess_condition.(fn{1}).chains.',per_sess_condition.(fn{1}).loops.']
                        [isinn,nidx]=ismember(mm{1},gnodes);
                        compidx=unique(conncomp(nidx(isinn)));
                        counter(compidx)=counter(compidx)+1;
                    end
                    subs=reshape(find(counter>1),1,[]);
                else
                    subs=1;
                end
                %{subg#, subg.Node#, subg.Edge#}
                for subsidx=subs
                    if nnz(conncomp==subsidx)<2
                        continue
                    end
                    subgh=gh.subgraph(conncomp==subsidx);
                    if ~isfield(stats,fn{1}) || ~isfield(stats.(fn{1}),'subg')
                        stats.(fn{1})=cell2struct({cell(0)},{'subg'});
                    end
                    stats.(fn{1}).subg=[stats.(fn{1}).subg,{subsidx,subgh.numnodes,subgh.numedges}];
                    if false
                        complexity_sums=[complexity_sums;per_trial_motif_cid{tt,3},subgh.numnodes,subgh.numedges,subgh.numedges./nchoosek(subgh.numnodes,2),max(max(subgh.distances))];
                        per_trl_nodes=[per_trl_nodes;{per_trial_motif_cid{tt,3}(1),str2double(table2array(subgh.Nodes))}];
                        degree_sums=[degree_sums;repmat(per_trial_motif_cid{tt,3},subgh.numnodes,1),str2double(table2array(subgh.Nodes)),subgh.degree];
                    end
                end
            end
            if ~isempty(opt.one_net)
                if opt.keep_all_subs
                    if exist('gh','var')
                        G=gh;
                    else
                        G=[];
                    end
                else
                    G=gh.subgraph(conncomp==1);
                end
            end
        end

        function out=match_one(noremove,rmv)
            out=cell2struct({0;0;0;0;0;cell(0);cell(0);cell(0);cell(0);cell(0)},{'abolish','shrank','split','weaken','noeffect','weaken_net','split_net','shrank_net','noeffect_net','abolish_net'});
            nrfns=reshape(fieldnames(noremove),1,[]);
            for fn=nrfns
                if isfield(rmv,fn{1})
                    if numel(rmv.(fn{1}).subg)>numel(noremove.(fn{1}).subg)
                        out.split=out.split+1;
                        out.split_net{end+1}=fn{1};
                    elseif rmv.(fn{1}).subg{1,2}<noremove.(fn{1}).subg{1,2}
                        out.shrank=out.shrank+1;
                        out.shrank_net{end+1}=fn{1};
                    elseif rmv.(fn{1}).subg{1,3}<noremove.(fn{1}).subg{1,3}
                        out.weaken=out.weaken+1;
                        out.weaken_net{end+1}=fn{1};
                    else
                        out.noeffect=out.noeffect+1;
                        out.noeffect_net{end+1}=fn{1};
                    end
                else
                    out.abolish=out.abolish+1;
                    out.abolish_net{end+1}=fn{1};
                end
            end
        end

        function demo(sschain_trl,ssloop_trl) % This plots only the end-result bar graph
            arguments
                sschain_trl
                ssloop_trl
            end
            if isempty(sschain_trl) && isempty(ssloop_trl)
                fstr=load(fullfile('binary','motif_replay.mat'));
                sschain_trl=fstr.chain_replay;
                ssloop_trl=fstr.ring_replay;
                clear fstr
            elseif isempty(sschain_trl) || isempty(ssloop_trl)
                error("Unfinished")
            end
            
            per_sess_condition=wave.composite_thin_down.merge_motif(sschain_trl,ssloop_trl);
            noremove=wave.composite_thin_down.stats_remove(per_sess_condition);
            removechain=wave.composite_thin_down.stats_remove(per_sess_condition,'remove','chains');
            removeloops=wave.composite_thin_down.stats_remove(per_sess_condition,'remove','loops');
            removeD5=wave.composite_thin_down.stats_remove(per_sess_condition,'remove','D5');

            nochain=wave.composite_thin_down.match_one(noremove,removechain);
            noloop=wave.composite_thin_down.match_one(noremove,removeloops);
            noD5=wave.composite_thin_down.match_one(noremove,removeD5);
            nrfns=reshape(fieldnames(noremove),1,[]);

            for rpt=1:10
                removeChainMatch=wave.composite_thin_down.stats_remove(per_sess_condition,'remove','ChainMatchFC');
                noChainMatch(rpt)=wave.composite_thin_down.match_one(noremove,removeChainMatch);
                removeLoopMatch=wave.composite_thin_down.stats_remove(per_sess_condition,'remove','LoopMatchFC');
                noLoopMatch(rpt)=wave.composite_thin_down.match_one(noremove,removeLoopMatch);
            end


            bmat=[0,0,0,0,numel(nrfns);...
                nochain.abolish,nochain.split,nochain.shrank,nochain.weaken,nochain.noeffect;...
                mean([noChainMatch.abolish]),mean([noChainMatch.split]),mean([noChainMatch.shrank]),mean([noChainMatch.weaken]),mean([noChainMatch.noeffect]);...
                noloop.abolish,noloop.split,noloop.shrank,noloop.weaken,noloop.noeffect;...
                mean([noLoopMatch.abolish]),mean([noLoopMatch.split]),mean([noLoopMatch.shrank]),mean([noLoopMatch.weaken]),mean([noLoopMatch.noeffect]);...
                noD5.abolish,noD5.split,noD5.shrank,noD5.weaken,noD5.noeffect];
            figure()
            bh=bar(fliplr(bmat)./sum(bmat,2).*100,'stacked');
            legend(bh,{'Unaffected','Weaken','shranked','Splitted','Abolished'},'Location','northoutside','Orientation','horizontal')
            set(gca(),'XTick',1:6,'XTickLabel',{'Observed','Remove chains','Chains ctrl','Remove loops','Loops ctrl','Remove deg>5'})
            ylabel('Percent (%)')
            ylim([0,100])
            %% ================================================================
        end


            % per_sess_condition=wave.composite_thin_down.merge_motif(sschain_trl,ssloop_trl);
            % onebyone.ctrl=wave.composite_thin_down.one_by_one_remove(per_sess_condition,'remove','ctrl')
        % 
        % onebyone.HIP=wave.composite_thin_down.one_by_one_remove(per_sess_condition,'remove','HIP')
        function stats=one_by_one_remove(per_sess_condition,opt)
            arguments
                per_sess_condition
                opt.remove (1,:) char {mustBeMember(opt.remove, {'chains','loops','ctrl','HIP','HIPNode'})} = 'chains'
                opt.debug (1,1) logical = false
            end
            stats=struct();

            %% systematic ablation
            % TODO: remove {HIP}
            fns=fieldnames(per_sess_condition);
            dbgcnt=1;

            if ismember(opt.remove,{'HIP','HIPNode'})
                su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
                regsess=-1;
            end

            for fn=reshape(fns,1,[])
                dbgcnt=dbgcnt+1;
                if opt.debug && dbgcnt>20
                    continue
                end

                per_sess_condition.(fn{1}).loops=cellfun(@(x) x([1:end,1]),per_sess_condition.(fn{1}).loops,'UniformOutput',false);% cyclic loops fix
                chaincount=numel(per_sess_condition.(fn{1}).chains);
                loopcount=numel(per_sess_condition.(fn{1}).loops);
                if chaincount+loopcount<=1
                    continue
                end

                % build graph network
                % one by one remove chains
                switch opt.remove
                    case 'chains'
                        motifcount=chaincount;
                    case 'loops'
                        motifcount=loopcount;
                    case 'ctrl'
                        motifcount=sum(cellfun(@(x) numel(x),[per_sess_condition.(fn{1}).chains;per_sess_condition.(fn{1}).loops])-1);
                    case'HIP'
                        % enumerate FC
                        motifcell=[per_sess_condition.(fn{1}).chains;per_sess_condition.(fn{1}).loops];
                        edgesall=cell2mat(cellfun(@(x) [x(1:end-1);x(2:end)].',...
                            motifcell,'UniformOutput',false));
                        % tag FC with region
                        sessid=str2double(regexp(fn{1},'(?<=s)\d{1,3}(?=s)','match','once'));
                        if sessid~=regsess
                            id2reg=containers.Map(num2cell(su_meta.allcid(su_meta.sess==sessid)),su_meta.reg_tree(5,su_meta.sess==sessid));
                            regsess=sessid;
                        end
                        reg_all=id2reg.values(num2cell(edgesall));
                        % count region==HIP
                        HIP_FC_sel=any(strcmp(reg_all,'HIP'),2);
                        motifcount=nnz(HIP_FC_sel);
                    case 'HIPNode'
                        % enumerate FC
                        motifcell=[per_sess_condition.(fn{1}).chains;per_sess_condition.(fn{1}).loops];
                        edgesall=cell2mat(cellfun(@(x) [x(1:end-1);x(2:end)].',...
                            motifcell,'UniformOutput',false));
                        % tag FC with region
                        sessid=str2double(regexp(fn{1},'(?<=s)\d{1,3}(?=s)','match','once'));
                        if sessid~=regsess
                            id2reg=containers.Map(num2cell(su_meta.allcid(su_meta.sess==sessid)),su_meta.reg_tree(5,su_meta.sess==sessid));
                            regsess=sessid;
                        end
                        ucid=unique(edgesall);
                        ureg=id2reg.values(num2cell(ucid));
                        % count region==HIP
                        HIP_SU_sel=strcmp(ureg,'HIP');
                        motifcount=nnz(HIP_SU_sel);
                end
                if motifcount==0
                    stats.(fn{1})='NA';
                else
                    for rcn=1:motifcount
                        if nchoosek(motifcount,rcn)>100 %randsample
                            compo_list=nan(100,motifcount-rcn);
                            for rr=1:100
                                compo_list(rr,:)=randsample(1:motifcount,motifcount-rcn);
                            end
                        else % traverse
                            % 1 of 1 condition branch moved downward
                            compo_list=nchoosek(1:motifcount,motifcount-rcn);
                        end

                        % repeats of numbered removal
                        for rr=1:size(compo_list,1)
                            if ~isfield(stats,fn{1}) || ~isfield(stats.(fn{1}),string(opt.remove)+rcn) || ~isfield(stats.(fn{1}).(string(opt.remove)+rcn),"rpt"+rr) || ~isfield(stats.(fn{1}).(string(opt.remove)+rcn).("rpt"+rr),'subg')
                                stats.(fn{1}).(string(opt.remove)+rcn).("rpt"+rr).subg=cell(0);
                            end
                            % stats after removal >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                            if strcmp(opt.remove,'HIP') % expecting edges after block
                                if motifcount==1
                                    keyboard()
                                end
                                hppool=find(HIP_FC_sel);
                                HIPkeep=hppool(compo_list(rr,:));
                                edges=categorical(unique([edgesall(~HIP_FC_sel,:);edgesall(HIPkeep,:)],'rows'));
                                if isempty(edges)
                                    stats.(fn{1}).(string(opt.remove)+rcn).("rpt"+rr).subg=[stats.(fn{1}).(string(opt.remove)+rcn).("rpt"+rr).subg,{0,0,0,chaincount,loopcount,opt.remove,[rcn,motifcount]}];
                                    continue
                                end
                            elseif strcmp(opt.remove,'HIPNode') % expecting edges after block
                                nonHIP=ucid(~HIP_SU_sel);
                                HIPall=ucid(HIP_SU_sel);
                                HIPkeep=HIPall(compo_list(rr,:));
                                edges=categorical(edgesall(all(ismember(edgesall,[nonHIP;HIPkeep]),2),:));
                                if isempty(edges)
                                    stats.(fn{1}).(string(opt.remove)+rcn).("rpt"+rr).subg=[stats.(fn{1}).(string(opt.remove)+rcn).("rpt"+rr).subg,{0,0,0,chaincount,loopcount,opt.remove,[rcn,motifcount]}];
                                    continue
                                end
                            else % expecting edges after block
                                if  rcn==motifcount
                                    switch opt.remove
                                        case'chains'
                                            motifcell=per_sess_condition.(fn{1}).loops;
                                        case 'loops'
                                            motifcell=per_sess_condition.(fn{1}).chains;
                                        case 'ctrl'
                                            motifcell=[];
                                        otherwise
                                            keyboard()
                                    end
                                else
                                    if strcmp(opt.remove,'chains')
                                        motifcell=[per_sess_condition.(fn{1}).chains(compo_list(rr,:));per_sess_condition.(fn{1}).loops];
                                    elseif strcmp(opt.remove,'loops')
                                        motifcell=[per_sess_condition.(fn{1}).chains;per_sess_condition.(fn{1}).loops(compo_list(rr,:))];
                                    elseif strcmp(opt.remove,'ctrl')
                                        motifcell=[per_sess_condition.(fn{1}).chains;per_sess_condition.(fn{1}).loops];
                                    end
                                end
                                if isempty(motifcell) || numel(motifcell)==1
                                    stats.(fn{1}).(string(opt.remove)+rcn).("rpt"+rr).subg=[stats.(fn{1}).(string(opt.remove)+rcn).("rpt"+rr).subg,{0,0,0,chaincount,loopcount,opt.remove,[rcn,motifcount]}];
                                    continue
                                end

                                if strcmp(opt.remove,'ctrl')
                                    edgesall=cell2mat(cellfun(@(x) [x(1:end-1);x(2:end)].',...
                                        motifcell,'UniformOutput',false));
                                    edges=categorical(unique(edgesall(compo_list(rr,:),:),'rows'));
                                else
                                    edges=categorical(unique(cell2mat(cellfun(@(x) [x(1:end-1);x(2:end)].',...
                                        motifcell,...
                                        'UniformOutput',false)),'rows'));
                                end
                            end
                            gh=graph(edges(:,1),edges(:,2));
                            conncomp=gh.conncomp();
                            % check module in network
                            if any(conncomp~=1)
                                comps=unique(conncomp);
                                counter=zeros(numel(comps),1);
                                gnodes=cellfun(@(x) str2double(x),gh.Nodes.Name);
                                for mm=[per_sess_condition.(fn{1}).chains.',per_sess_condition.(fn{1}).loops.']
                                    [isinn,nidx]=ismember(mm{1},gnodes);
                                    compidx=unique(conncomp(nidx(isinn)));
                                    counter(compidx)=counter(compidx)+1;
                                end
                                subs=reshape(find(counter>1),1,[]); % TODO: should consider >=1
                                if isempty(subs)
                                    stats.(fn{1}).(string(opt.remove)+rcn).("rpt"+rr).subg=[stats.(fn{1}).(string(opt.remove)+rcn).("rpt"+rr).subg,{0,0,0,chaincount,loopcount,opt.remove,[rcn,motifcount]}];
                                    continue
                                end
                            else
                                subs=1;
                            end
                            %{subg#, subg.Node#, subg.Edge#}
                            for subsidx=subs
                                if nnz(conncomp==subsidx)<2 % corresponding to remove of nodes
                                    keyboard()
                                    continue
                                end
                                subgh=gh.subgraph(conncomp==subsidx);
                                stats.(fn{1}).(string(opt.remove)+rcn).("rpt"+rr).subg=[stats.(fn{1}).(string(opt.remove)+rcn).("rpt"+rr).subg,{subsidx,subgh.numnodes,subgh.numedges,chaincount,loopcount,opt.remove,[rcn,motifcount]}];
                            end
                            if isempty(stats.(fn{1}).(string(opt.remove)+rcn).("rpt"+rr).subg)
                                keyboard()
                            end
                        end
                    end
                end
            end
        end

        function [system_rmv_stats,rmvnet]=match_system(noremove,onebyone)
            % out=cell2struct({0;0;0;0;0;cell(0);cell(0);cell(0);cell(0);cell(0)},{'abolish','shrank','split','weaken','noeffect','weaken_net','split_net','shrank_net','noeffect_net','abolish_net'});
            nrfns=reshape(fieldnames(noremove),1,[]);
            system_rmv_stats=struct();
            rmvtypes=string(fieldnames(onebyone));
            for rmvtype=reshape(rmvtypes,1,[])
                rmvnet.(rmvtype)=cell2struct({cell(0);cell(0);cell(0);cell(0);cell(0);cell(0)},{'abolish','shrank','split','weaken','noeffect','WIP'});
                for currCondition=nrfns

                    if isfield(onebyone.(rmvtype),currCondition{1})
                        if isequal(onebyone.(rmvtype).(currCondition{1}),'NA')
                            disp("skipped NA condition")
                        else
                            levels=string(fieldnames(onebyone.(rmvtype).(currCondition{1})));
                            for lvl=reshape(levels,1,[])
                                out=cell2struct({0;0;0;0;0;0},{'abolish','shrank','split','weaken','noeffect','WIP'});
                                rpts=onebyone.(rmvtype).(currCondition{1}).(lvl);
                                rptcell=struct2cell(rpts);
                                
                                % {subsidx,subgh.numnodes,subgh.numedges,chaincount,loopcount,opt.remove,rcn}
                                for ii=1:numel(rptcell)
                                    if numel(noremove.(currCondition{1}).subg)>3 || noremove.(currCondition{1}).subg{1,1}>1
                                        out.WIP=out.WIP+1;
                                        if lvl==levels(end)
                                            rmvnet.(rmvtype).WIP=[rmvnet.(rmvtype).WIP,currCondition];
                                        end
                                    elseif rptcell{ii}.subg{1,2}==0
                                        out.abolish=out.abolish+1;
                                        if lvl==levels(end)
                                        rmvnet.(rmvtype).abolish=[rmvnet.(rmvtype).abolish,currCondition];
                                        end
                                    elseif (numel(rptcell{ii}.subg)/7)>(numel(noremove.(currCondition{1}).subg)/3)
                                        out.split=out.split+1;
                                        if lvl==levels(end)
                                        rmvnet.(rmvtype).split=[rmvnet.(rmvtype).split,currCondition];
                                        end
                                    elseif rptcell{ii}.subg{1,2}<noremove.(currCondition{1}).subg{1,2}
                                        out.shrank=out.shrank+1;
                                        if lvl==levels(end)
                                        rmvnet.(rmvtype).shrank=[rmvnet.(rmvtype).shrank,currCondition];
                                        end
                                    elseif rptcell{ii}.subg{1,3}<noremove.(currCondition{1}).subg{1,3}
                                        out.weaken=out.weaken+1;
                                        if lvl==levels(end)
                                        rmvnet.(rmvtype).weaken=[rmvnet.(rmvtype).weaken,currCondition];
                                        end
                                    else
                                        % if lvl==levels(end)
                                        %     keyboard()
                                        % end
                                        out.noeffect=out.noeffect+1;
                                        if lvl==levels(end)
                                        rmvnet.(rmvtype).noeffect=[rmvnet.(rmvtype).noeffect,currCondition];
                                        end
                                    end
                                end
                                % rmvlevel=str2double(regexp(lvl,'(?<=(chains|loops))\d','match','once'))
                                if contains(rmvtype,'chain')
                                    rmvlevel=floor(rptcell{1}.subg{7}(1)./rptcell{1}.subg{4}.*10);
                                elseif contains(rmvtype,'loop')
                                    rmvlevel=floor(rptcell{1}.subg{7}(1)./rptcell{1}.subg{5}.*10);
                                elseif contains(rmvtype,'ctrl') || contains(rmvtype,'HIP')
                                    rmvlevel=floor(rptcell{1}.subg{7}(1)./rptcell{1}.subg{7}(2).*10);
                                else
                                    keyboard()
                                end
                                if ~isfield(system_rmv_stats,rmvtype) || ~isfield(system_rmv_stats.(rmvtype),"L"+rmvlevel)
                                    system_rmv_stats.(rmvtype).("L"+rmvlevel).out=[];
                                end
                                if ~(numel(rptcell)==out.WIP)
                                    system_rmv_stats.(rmvtype).("L"+rmvlevel).out=...
                                        [system_rmv_stats.(rmvtype).("L"+rmvlevel).out;[out.abolish,out.split,out.shrank,out.weaken,out.noeffect]./(numel(rptcell)-out.WIP)];
                                end
                            end
                        end
                    else
                        % catch exception
                        keyboard()
                    end
                end
            end
        end

        function plot_one_by_one()
            load(fullfile("binary","onebyonermv.mat"),"system_rmv_stats");
            mm=struct();
            sem=struct();
            for motifType=["chains","loops","ctrl","HIP"]
                for lvl=0:10
                    mm.(motifType)(lvl+1,:)=mean(system_rmv_stats.(motifType).("L"+lvl).out);
                    sem.(motifType)(lvl+1,:)=std(system_rmv_stats.(motifType).("L"+lvl).out)./sqrt(size(system_rmv_stats.(motifType).("L"+lvl).out,1));
                end
            end
            hexcmap=["#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F"];
            fh=figure();
            tiledlayout(1,4)

            nexttile()
            hold on
            for cc=1:5
                fill([5:10:95,100,100,95:-10:5],[mm.ctrl(:,cc)-sem.ctrl(:,cc);flip(mm.ctrl(:,cc)+sem.ctrl(:,cc))],'k','FaceColor',hexcmap(cc),'EdgeColor','none','FaceAlpha','0.2');
                ph(cc)=plot([5:10:95,100],mm.ctrl(:,cc),'-','Color',hexcmap(cc));
            end
            
            
            set(gca,'YTick',0.2:0.2:1,'YTickLabel',20:20:100)
            xlim([0,100])
            ylim([0,1])
            xlabel('Percentage of removed random SC (%)')
            ylabel('Effect to composite loops (%)')

            nexttile()
            hold on
            for cc=1:5
                fill([5:10:95,100,100,95:-10:5],[mm.chains(:,cc)-sem.chains(:,cc);flip(mm.chains(:,cc)+sem.chains(:,cc))],'k','FaceColor',hexcmap(cc),'EdgeColor','none','FaceAlpha','0.2');
                ph(cc)=plot([5:10:95,100],mm.chains(:,cc),'-','Color',hexcmap(cc));
            end
            
            legend(ph,{'abolish','split','shrank','weakened','noeffect'},'Location','northoutside','Orientation','horizontal')
            set(gca,'YTick',0.2:0.2:1,'YTickLabel',20:20:100)
            xlim([0,100])
            ylim([0,1])
            xlabel('Percentage of removed chains (%)')
            ylabel('Effect to composite loops (%)')

            nexttile()
            hold on
            for cc=1:5
                fill([5:10:95,100,100,95:-10:5],[mm.loops(:,cc)-sem.loops(:,cc);flip(mm.loops(:,cc)+sem.loops(:,cc))],'k','FaceColor',hexcmap(cc),'EdgeColor','none','FaceAlpha','0.2');
                ph(cc)=plot([5:10:95,100],mm.loops(:,cc),'-','Color',hexcmap(cc));
            end
            set(gca,'YTick',0.2:0.2:1,'YTickLabel',20:20:100)
            xlim([0,100])
            ylim([0,1])
            xlabel('Percentage of removed loops (%)')
            ylabel('Effect to composite loops (%)')     

            nexttile()
            hold on
            for cc=1:5
                fill([5:10:95,100,100,95:-10:5],[mm.HIP(:,cc)-sem.HIP(:,cc);flip(mm.HIP(:,cc)+sem.HIP(:,cc))],'k','FaceColor',hexcmap(cc),'EdgeColor','none','FaceAlpha','0.2');
                ph(cc)=plot([5:10:95,100],mm.HIP(:,cc),'-','Color',hexcmap(cc));
            end
            
            set(gca,'YTick',0.2:0.2:1,'YTickLabel',20:20:100)
            xlim([0,100])
            ylim([0,1])
            xlabel('Percentage of removed HIP Neuron (%)')
            ylabel('Effect to composite loops (%)')
            savefig(fh,fullfile("binary","one_by_one_rmv.fig"));
        end
        function gen_one_by_one_rmv_data()
            load(fullfile("binary","motif_replay.mat"),'chain_replay','ring_replay');
            per_sess_condition=wave.composite_thin_down.merge_motif(chain_replay,ring_replay);
            noremove=wave.composite_thin_down.stats_remove(per_sess_condition);
            onebyone.ctrl=wave.composite_thin_down.one_by_one_remove(per_sess_condition,'remove','ctrl');
            onebyone.HIP=wave.composite_thin_down.one_by_one_remove(per_sess_condition,'remove','HIPNode');
            onebyone.chains=wave.composite_thin_down.one_by_one_remove(per_sess_condition,'remove','chains');
            onebyone.loops=wave.composite_thin_down.one_by_one_remove(per_sess_condition,'remove','loops');
            [system_rmv_stats,rmv_nets]=wave.composite_thin_down.match_system(noremove,onebyone);
            blame=vcs.blame();
            save(fullfile("binary","one_by_one_rmv.mat"),"onebyone","system_rmv_stats","rmv_nets","blame");
        end
    end
end