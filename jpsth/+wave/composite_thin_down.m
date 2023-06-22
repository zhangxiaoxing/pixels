% tag spikes in neuron with loop or chain activity
% olf, enc-both disconnect, composite proportion, bargraph
% TODO: overall refactoring
classdef composite_thin_down < handle

    methods (Static)
        function chains=getChains(sschain,dur,fns,sessid)
            chains=cell(0);
            for fn=fns
                sesssel=startsWith(fieldnames(sschain.out.(dur).(fn)),['s',num2str(sessid),'c']);
                if ~any(sesssel)
                    continue
                else
                    chains=[chains;...
                        cellfun(@(x) x.meta{1},...
                        subsref(struct2cell(sschain.out.(dur).(fn)),substruct('()',{sesssel})),...
                        'UniformOutput',false)...
                        ];
                end
            end
        end

        function loops=getLoops(pstats,wids,sessid)
            sesssel=startsWith(fieldnames(pstats.congru),['s',num2str(sessid),'r']);
            sess_loops=subsref(struct2cell(pstats.congru),substruct('()',{sesssel}));
            wave_sel=cellfun(@(x) all(ismember(x.rstats{4},wids)),sess_loops);
            loops=cellfun(@(x) x.rstats{3},sess_loops(wave_sel),'UniformOutput',false);
        end

        function per_sess_condition=merge_motif(sschain,pstats)
            arguments
                sschain
                pstats
            end

            %% single spike chain
            keys=[struct2cell(structfun(@(x) fieldnames(x), sschain.out.d6, 'UniformOutput', false));...
                struct2cell(structfun(@(x) fieldnames(x), sschain.out.d3, 'UniformOutput', false))];
            keys=vertcat(keys{:});
            ssc_sess=unique(str2double(regexp(keys,'(?<=s)\d{1,3}(?=c)','match','once')));

            %% single spike loop
            if isfield(pstats,'nonmem'), pstats=rmfield(pstats,"nonmem");end
            ssl_sess=unique(str2double(regexp(fieldnames(pstats.congru),'(?<=s)\d{1,3}(?=r)','match','once')));
            usess=union(ssc_sess,ssl_sess);

            %% per-session entrance
            

            for sessid=reshape(usess,1,[])
                
                %% single spike chain
                % four conditions
                
                per_sess_condition.("s"+sessid+"s1d3")=struct();
                per_sess_condition.("s"+sessid+"s1d3").chains=wave.composite_thin_down.getChains(sschain,"d3",["olf_s1","s1d3"],sessid);
                per_sess_condition.("s"+sessid+"s1d3").loops=wave.composite_thin_down.getLoops(pstats,[1 5],sessid);
                
                % d3s2
                per_sess_condition.("s"+sessid+"s1d6")=struct();
                per_sess_condition.("s"+sessid+"s1d6").chains=wave.composite_thin_down.getChains(sschain,"d6",["olf_s1","s1d6"],sessid);
                per_sess_condition.("s"+sessid+"s1d6").loops=wave.composite_thin_down.getLoops(pstats,[2 5],sessid);

                % d6s1
                per_sess_condition.("s"+sessid+"s2d3")=struct();
                per_sess_condition.("s"+sessid+"s2d3").chains=wave.composite_thin_down.getChains(sschain,"d3",["olf_s2","s2d3"],sessid);
                per_sess_condition.("s"+sessid+"s2d3").loops=wave.composite_thin_down.getLoops(pstats,[3 6],sessid);

                % d6s2
                per_sess_condition.("s"+sessid+"s2d6")=struct();
                per_sess_condition.("s"+sessid+"s2d6").chains=wave.composite_thin_down.getChains(sschain,"d6",["olf_s2","s2d6"],sessid);
                per_sess_condition.("s"+sessid+"s2d6").loops=wave.composite_thin_down.getLoops(pstats,[4 6],sessid);
            end
        end
        
        function [stats,G]=stats_remove(per_sess_condition,opt)
            arguments
                per_sess_condition
                opt.remove (1,:) char {mustBeMember(opt.remove, {'','chains','loops','UID','D10','D5','ChainEqNode','LoopEqNode'})} = ''
                opt.one_net (1,:) char = ''
                opt.keep_subs (1,1) logical = false % for illustration only
            end
            
            stats=struct();
            % per_trl_nodes=cell(0);
            % complexity_sums=[];
            % degree_sums=[];

            %% TODO: gradual removal
            % remove {all chains},{all loops},{HIP},{ORB},{OLF},{DEGREE>10}

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

                if (numel(per_sess_condition.(fn{1}).chains)+numel(per_sess_condition.(fn{1}).loops)<=1 && ~opt.keep_subs) ...
                    || (numel(per_sess_condition.(fn{1}).chains)+numel(per_sess_condition.(fn{1}).loops)<1 && opt.keep_subs)
                    continue
                end


                % build graph network
                per_sess_condition.(fn{1}).loops=cellfun(@(x) x([1:end,1]),per_sess_condition.(fn{1}).loops,'UniformOutput',false);% cyclic loops fix
                edges=categorical(unique(cell2mat(cellfun(@(x) [x(1:end-1);x(2:end)].',...
                    [per_sess_condition.(fn{1}).chains;per_sess_condition.(fn{1}).loops],...
                    'UniformOutput',false)),'rows'));
                if strcmp(opt.remove,'UID')
                    assert(~isempty(opt.UID), "missing UID input");
                end

                gh=graph(edges(:,1),edges(:,2));
                if strcmp(opt.remove,'D10')
                    gh=gh.rmnode(gh.Nodes.Name(gh.degree>=10));
                elseif strcmp(opt.remove,'D5')
                    gh=gh.rmnode(gh.Nodes.Name(gh.degree>=5));
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
                if opt.keep_subs
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
            out=cell2struct({0;0;0;0;0;cell(0);cell(0);cell(0);cell(0);cell(0)},{'abolish','shrink','split','weaken','noeffect','weaken_net','split_net','shrink_net','noeffect_net','abolish_net'});
            nrfns=reshape(fieldnames(noremove),1,[]);
            for fn=nrfns
                if isfield(rmv,fn{1})
                    if numel(rmv.(fn{1}).subg)>numel(noremove.(fn{1}).subg)
                        out.split=out.split+1;
                        out.split_net{end+1}=fn{1};
                    elseif rmv.(fn{1}).subg{1,2}<noremove.(fn{1}).subg{1,2}
                        out.shrink=out.shrink+1;
                        out.shrink_net{end+1}=fn{1};
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

        function demo(sschain,pstats)
            %% ================================================================

            per_sess_condition=wave.composite_thin_down.merge_motif(sschain,pstats);
            noremove=wave.composite_thin_down.stats_remove(per_sess_condition);
            removechain=wave.composite_thin_down.stats_remove(per_sess_condition,'remove','chains');
            removeloops=wave.composite_thin_down.stats_remove(per_sess_condition,'remove','loops');
            removeD10=wave.composite_thin_down.stats_remove(per_sess_condition,'remove','D10');
            removeD5=wave.composite_thin_down.stats_remove(per_sess_condition,'remove','D5');

            nochain=wave.composite_thin_down.match_one(noremove,removechain);
            noloop=wave.composite_thin_down.match_one(noremove,removeloops);
            noD10=wave.composite_thin_down.match_one(noremove,removeD10);
            noD5=wave.composite_thin_down.match_one(noremove,removeD5);
            nrfns=reshape(fieldnames(noremove),1,[]);
    
            bmat=[0,0,0,0,numel(nrfns);...
                nochain.abolish,nochain.split,nochain.shrink,nochain.weaken,nochain.noeffect;...
                noloop.abolish,noloop.split,noloop.shrink,noloop.weaken,noloop.noeffect;...
                noD10.abolish,noD10.split,noD10.shrink,noD10.weaken,noD10.noeffect;...
                noD5.abolish,noD5.split,noD5.shrink,noD5.weaken,noD5.noeffect];
            figure()
            bh=bar(fliplr(bmat)./sum(bmat,2).*100,'stacked');
            legend(bh,{'Unaffected','Weaken','Shrinked','Splitted','Abolished'},'Location','northoutside','Orientation','horizontal')
            set(gca(),'XTick',1:5,'XTickLabel',{'Observed','Remove chains','Remove loops','Remove deg>10','Remove deg>5'})
            ylabel('Percent (%)')
            ylim([0,100])
            %% ================================================================
        end
        

        function more_stats()

            %% stats
            num_su=cellfun(@(x) numel(x),struct2cell(per_reg));
            thresh_sel=num_su>=5;

            per_reg_mm=cellfun(@(x) mean(x),subsref(struct2cell(per_reg),substruct('()',{thresh_sel})));
            per_reg_std=cellfun(@(x) std(x),subsref(struct2cell(per_reg),substruct('()',{thresh_sel})));
            per_reg_sem=per_reg_std./sqrt(num_su(thresh_sel));

            others_mm=mean(cell2mat(subsref(struct2cell(per_reg),substruct('()',{~thresh_sel}))));
            others_sem=std(cell2mat(subsref(struct2cell(per_reg),substruct('()',{~thresh_sel}))))./sqrt(nnz(~thresh_sel));

            [plotmm,plotidx]=sort([per_reg_mm;others_mm],'descend');
            plotsem=[per_reg_sem;others_sem];
            %figure
            figure();
            hold on;
            bh=bar(plotmm,'FaceColor','w');
            xlbl=[subsref(fieldnames(per_reg),substruct('()',{thresh_sel}));'Others'];
            set(gca,'XTick',1:(nnz(thresh_sel)+1),'XTickLabel',xlbl(plotidx));
            errorbar(bh.XData,bh.YData,plotsem(plotidx),'k.')
            ylabel('Average node degree')
            title('Per region neuron degrees')

            %% overall graph stats
            ukeys=string.empty(0,1);
            uniq_net=[];
            for nidx=1:size(complexity_sums,1)
                onekey=string(sprintf('%d-',per_trl_nodes{nidx,1},sort(per_trl_nodes{nidx,2})));
                if ~ismember(onekey,ukeys)
                    ukeys=[ukeys;onekey];
                    uniq_net=[uniq_net;complexity_sums(nidx,[1,2,5:end])];
                end
            end
        end

        function before_after_loops(per_sess_condition,fn)
            %%
            % pstats=bz.rings.rings_time_constant.stats([],[],'load_file',true,'skip_save',true)
            % load(fullfile("bzdata","chain_tag.mat"),"out")
            % per_sess_condition=wave.composite_thin_down.merge_motif(sschain,pstats);
            % wave.composite_thin_down.before_after_loops(per_sess_condition,'s18s1d3')
            % ==============================
            pchains=per_sess_condition.(fn).chains;
            ploops=cellfun(@(x) x([1:end,1]),per_sess_condition.(fn).loops,'UniformOutput',false);% cyclic loops fix
            
            chainedges=unique(cell2mat(cellfun(@(x) [x(1:end-1);x(2:end)].',...
                pchains,'UniformOutput',false)),'crows');
            loopedges=unique(cell2mat(cellfun(@(x) [x(1:end-1);x(2:end)].',...
                ploops,'UniformOutput',false)),'rows');
            % ==============================
            [~,G]=wave.composite_thin_down.stats_remove(per_sess_condition,'one_net',fn,'keep_subs',false);

            gedges=str2double(G.Edges.EndNodes);
            common_edges=intersect([chainedges;loopedges],([gedges;fliplr(gedges)]),'rows');
            if ~isempty(chainedges)
                inchain=ismember(common_edges,chainedges,"rows");
            else
                inchain=false(size(common_edges,1));
            end
            if ~isempty(loopedges)
                inloop=ismember(common_edges,loopedges,"rows");
            else
                inloop=false(size(common_edges,1));
            end
            cmat=zeros(size(common_edges,1),3);
            cmat(inchain,3)=1;
            cmat(inloop,1)=1;
            DG=digraph(table(string(common_edges),cmat,'VariableNames',{'EndNodes','EdgeColor'}));
            
            figure()
            tiledlayout(1,2)
            nexttile()
            flbl=DG.plot.NodeLabel;
            fx=DG.plot.XData;
            fy=DG.plot.YData;

            plot(DG,'EdgeColor',DG.Edges.EdgeColor,'LineWidth',2,'ArrowSize',10);
            xspan=xlim();
            yspan=ylim();

            nexttile()
            [~,Gr]=wave.composite_thin_down.stats_remove(per_sess_condition,'one_net',fn,'remove','loops','keep_subs',true);
            if ~isempty(Gr)
                gredges=str2double(Gr.Edges.EndNodes);
                common_r_edges=intersect([chainedges;loopedges],([gredges;fliplr(gredges)]),'rows');
                common_r_edges=intersect(common_r_edges,[gedges;fliplr(gedges)],'rows');
                if ~isempty(chainedges)
                    inchain=ismember(common_r_edges,chainedges,"rows");
                else
                    inchain=false(size(common_r_edges,1));
                end
                if ~isempty(loopedges)
                    inloop=ismember(common_r_edges,loopedges,"rows");
                else
                    inloop=false(size(common_r_edges,1));
                end
                cmat=zeros(size(common_r_edges,1),3);
                cmat(inchain,3)=1;
                cmat(inloop,1)=1;
                DGr=digraph(table(string(common_r_edges),cmat,'VariableNames',{'EndNodes','EdgeColor'}));

                [~,rel_pos]=ismember(DGr.plot.NodeLabel,flbl);
                rx=fx(rel_pos);
                ry=fy(rel_pos);
                plot(DGr,'EdgeColor',DGr.Edges.EdgeColor,'LineWidth',2,'ArrowSize',10,...
                    'XData',rx,'YData',ry);
                xlim(xspan)
                ylim(yspan)
                set(gca(),'XTick',[],'YTick',[])
            end
            sgtitle([fn,' remove loops'])
            %%
        end

    end
end