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
                opt.remove (1,:) char {mustBeMember(opt.remove, {'','chains','loops','UID','D10','D5'})} = ''
                opt.UID = []
                opt.one_net (1,:) char = ''
            end
            
            stats=struct();
            per_trl_nodes=cell(0);
            complexity_sums=[];
            degree_sums=[];

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
                if numel(per_sess_condition.(fn{1}).chains)+numel(per_sess_condition.(fn{1}).loops)<=1
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
                %{subg#, subg.Node#}
                for subsidx=subs
                    if nnz(conncomp==subsidx)<2
                        continue
                    end
                    subgh=gh.subgraph(conncomp==subsidx);
                    if ~isfield(stats,fn{1}) || ~isfield(stats.(fn{1}),'subg')
                        stats.(fn{1})=cell2struct({cell(0)},{'subg'});
                    end
                    stats.(fn{1}).subg=[stats.(fn{1}).subg,{subsidx,subgh.numnodes}];
                    if false
                        complexity_sums=[complexity_sums;per_trial_motif_cid{tt,3},subgh.numnodes,subgh.numedges,subgh.numedges./nchoosek(subgh.numnodes,2),max(max(subgh.distances))];
                        per_trl_nodes=[per_trl_nodes;{per_trial_motif_cid{tt,3}(1),str2double(table2array(subgh.Nodes))}];
                        degree_sums=[degree_sums;repmat(per_trial_motif_cid{tt,3},subgh.numnodes,1),str2double(table2array(subgh.Nodes)),subgh.degree];
                    end
                end
            end
            if ~isempty(opt.one_net)
                G=subgh;
            end
        end

        function out=match_one(noremove,rmv)
            out=cell2struct({0;0;0;0;cell(0);cell(0);cell(0);cell(0)},{'complete','partial','split','noeffect','split_net','partial_net','noeffect_net','complete_net'});
            nrfns=reshape(fieldnames(noremove),1,[]);
            for fn=nrfns
                if isfield(rmv,fn{1})
                    if numel(rmv.(fn{1}).subg)>numel(noremove.(fn{1}).subg)
                        out.split=out.split+1;
                        out.split_net{end+1}=fn{1};
                    elseif rmv.(fn{1}).subg{1,2}<noremove.(fn{1}).subg{1,2}
                        out.partial=out.partial+1;
                        out.partial_net{end+1}=fn{1};
                    else
                        out.noeffect=out.noeffect+1;
                        out.noeffect_net{end+1}=fn{1};
                    end
                else
                    out.complete=out.complete+1;
                    out.complete_net{end+1}=fn{1};
                end
            end
        end

        function demo(sschain,pstats)
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
    
            bmat=[0,0,0,numel(nrfns);...
                nochain.complete,nochain.split,nochain.partial,nochain.noeffect;...
                noloop.complete,noloop.split,noloop.partial,noloop.noeffect;...
                noD10.complete,noD10.split,noD10.partial,noD10.noeffect;...
                noD5.complete,noD5.split,noD5.partial,noD5.noeffect];
            figure()
            bh=bar(fliplr(bmat),'stacked');
            legend(bh,{'Unaffected','Shrinked','Splitted','Abolished'},'Location','northoutside','Orientation','horizontal')
            set(gca(),'XTick',1:5,'XTickLabel',{'Observed','Remove chains','Remove loops','Remove deg>10','Remove deg>5'},'YTick',0:5:25,'YTickLabel',0:20:100)
            ylabel('Percent (%)')

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
            [~,G]=wave.composite_thin_down.stats_remove(per_sess_condition,'one_net',fn);
            figure()
            tiledlayout(1,2)
            nexttile()
            plot(G)
            flbl=G.plot.NodeLabel;
            fx=G.plot.XData;
            fy=G.plot.YData;
            xspan=xlim();
            yspan=ylim();

            nexttile()
            [~,Gr]=wave.composite_thin_down.stats_remove(per_sess_condition,'one_net',fn,'remove','loops');
            rx=fx(ismember(flbl,Gr.plot.NodeLabel));
            ry=fy(ismember(flbl,Gr.plot.NodeLabel));
            plot(Gr,'XData',rx,'YData',ry)
            xlim(xspan)
            ylim(yspan)
            sgtitle([fn,' remove loops'])
            
        end




    end
end