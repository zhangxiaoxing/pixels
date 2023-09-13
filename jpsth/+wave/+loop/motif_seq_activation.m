function motif_seq_activation(chain_replay,ring_replay,trials_dict)
arguments
    chain_replay = []
    ring_replay = []
    trials_dict = []
end
if isempty(chain_replay) || isempty (ring_replay)
    load(fullfile('binary','motif_replay.mat'),'ring_replay','chain_replay');
end
if isempty(trials_dict)
    load(fullfile('binary','trials_dict.mat'),'trials_dict');
end

load(fullfile('binary','su_meta.mat'));
load(fullfile('binary','wrs_mux_meta.mat'),'wrs_mux_meta');
reg_com_maps=wave.get_reg_com_maps(wrs_mux_meta);

usess=union(chain_replay.session,ring_replay.session);
for sessid=114%reshape(usess,1,[])
    sesscid=su_meta.allcid(su_meta.sess==sessid);
    sessreg=su_meta.reg_tree(5,su_meta.sess==sessid).';
    sess_reg_map=dictionary(sesscid(~ismissing(sessreg)),sessreg(~ismissing(sessreg)));
    trials=cell2mat(trials_dict(sessid));
    for wid="s1d6"%["s1d3","s2d3","s1d6","s2d6"]
        disp(wid+"of sess"+sessid)
        samp=str2double(regexp(wid,"(?<=s)\d(?=d)",'match','once')).*4;
        delay=str2double(regexp(wid,"(?<=d)\d$",'match','once'));
        trial_sel=find(trials(:,5)==samp & trials(:,8)==delay & all(trials(:,9:10)>0,2));
        ring_reg_sel=cellfun(@(x) all(sess_reg_map.isKey(x)),ring_replay.meta(:,2));
        chain_reg_sel=cellfun(@(x) all(sess_reg_map.isKey(x)),chain_replay.meta(:,2));
        ringsel=find(ring_replay.session==sessid & ring_replay.wave==wid & ring_reg_sel);
        chainsel=find(chain_replay.session==sessid & chain_replay.wave==wid & chain_reg_sel);
        for tt=144%reshape(trial_sel,1,[])
            tonset=trials(tt,1);
            motif_spk=cell(0,1);
            motif_id=cell(0,2);
            motif_seq=cell(0,1);
            edges=[];
            for ridx=reshape(ringsel,1,[])
                trl_align=ring_replay.trl_align{ridx};
                pref_delay=trl_align(:,1)==tt & all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
                currts=ring_replay.ts{ridx}(pref_delay);
                if ~isempty(currts)
                    motif_spk=[motif_spk;{currts}];
                    motif_id=[motif_id;{"r"+ridx,ring_replay.meta{ridx,2}}];
                    motif_seq=[motif_seq;{ring_replay.ts_seq{ridx}(pref_delay)}];
                    edges=[edges;ring_replay.meta{ridx,2}([1:end;2:end,1]).'];
                end
            end
            for cidx=reshape(chainsel,1,[])
                trl_align=chain_replay.trl_align{cidx};
                pref_delay=trl_align(:,1)==tt & all(trl_align(:,5:7)==1,2) & trl_align(:,2)>=1 & trl_align(:,2)<(trl_align(:,4)+1);
                currts=chain_replay.ts{cidx}(pref_delay,:);
                if ~isempty(currts)
                    motif_spk=[motif_spk;{mat2cell(currts.',size(currts,2),ones(size(currts,1),1)).'}];
                    motif_id=[motif_id;{"c"+cidx,chain_replay.meta{cidx,2}}];
                    motif_seq=[motif_seq;{repmat({chain_replay.meta{cidx,2}.'},nnz(pref_delay),1)}];
                    edges=[edges;chain_replay.meta{cidx,2}([1:end-1;2:end]).'];
                end
            end
            if isempty(motif_id) || ~any(startsWith([motif_id{:,1}],"r")) || ~any(startsWith([motif_id{:,1}],"c"))
                continue
            end

            if size(motif_id,1)<30 || size(motif_id,1)>100
                continue
            end

            uedges=categorical(unique(edges,'rows'));
            gh=graph(uedges(:,1),uedges(:,2));
            if any(gh.conncomp>1)
                [gc,gr]=groupcounts(gh.conncomp.');
                [~,maxid]=max(gc);
                gh=gh.subgraph(gh.conncomp==gr(maxid));
            end

            % sort by reg_tcom
            nodereg=sess_reg_map(str2double(gh.Nodes.Name));
            if ~all(reg_com_maps.("tcom"+delay+"_maps").odor_only.isKey(nodereg))
                continue
            end

            nodetcom=cell2mat(reg_com_maps.("tcom"+delay+"_maps").odor_only.values(nodereg));
            [~,sidx]=sort(nodetcom);

            if false
                csvcell=[{'Source','Target'};gh.Edges.EndNodes];
                writecell(csvcell,fullfile('binary','gephi',sprintf('SingleSpkConn4gephi_s%dt_%d.csv',sessid,tt)));

                nodecell=[{'Id','Label','Toposort','REG'};...
                    gh.Nodes.Name,gh.Nodes.Name,num2cell(nodetcom),nodereg];
                writecell(nodecell,fullfile('binary','gephi',sprintf('SingleSpkNode4gephi_s%dt_%d.csv',sessid,tt)));
            end
            % plot(gh)
            keptnodes=str2double(gh.Nodes.Name);
            gsel=cellfun(@(x) all(ismember(x,keptnodes)),motif_id(:,2));
            motif_spk=motif_spk(gsel,:);
            motif_id=motif_id(gsel,:);
            motif_seq=motif_seq(gsel,:);

            covered=false(1,3000);
            for rridx=1:size(motif_spk,1)
                for mrpt=1:numel(motif_spk{rridx})
                    raster=round((motif_spk{rridx}{mrpt}-tonset-30000)./30); % 30 @ 1ms resolution
                    for ii=2:numel(raster)
                        onset=ceil(raster(ii-1));
                        offset=ceil(raster(ii));
                        covered((onset+1):(offset+1))=true;
                    end
                end
            end
            
            edges = find(diff([0,covered,0]==1));
            eonset = edges(1:2:end-1);  % Start indices
            run_length = edges(2:2:end)-eonset;  % Consecutive ones counts

            if ~any(run_length>90)
                continue
            end
            
            [maxl,maxlid]=max(run_length);
            maxonset=eonset(maxlid);
            maxoffset=edges(2*maxlid);

            %% plot per-spike figure
            if ~exist('FT_SPIKE','var') || FT_SPIKE.sessid~=sessid
                [~,~,~,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true);
                FT_SPIKE.sessid=sessid;
            end


            figure();
            cmap=colormap('lines');
            hold on;
            % plot spikes
            for cidx=1:numel(nodetcom)
                ft_sel=strcmp(FT_SPIKE.label,gh.Nodes.Name(cidx));
                nodeid=keptnodes(cidx);
                ts=[];
                for jj=1:size(motif_spk,1)
                    mts=cell2mat(motif_spk{jj});
                    mid=cell2mat(motif_seq{jj});
                    ts=[ts;mts(mid==nodeid)];
                end
                motif_raster=unique(round((ts-tonset-30000)./30));
                plot(motif_raster,sidx(cidx),'|','MarkerEdgeColor','#FF0000','LineWidth',2);

                allts=FT_SPIKE.time{ft_sel}(FT_SPIKE.trial{ft_sel}==tt);
                other_raster=unique(round(allts*1000)-1000);
                other_raster=other_raster(other_raster>=0 & other_raster<=delay.*1000);
                other_ts=setdiff(other_raster,motif_raster);
                plot(other_ts,sidx(cidx),'|','MarkerEdgeColor',"#808080")
            end

            % plot connections
            for midx=1:size(motif_spk,1)
                if startsWith(motif_id{midx,1},"r")
                    cc=cmap(randsample([1 5 6],1),:);
                    ls='--';
                    lw=2;
                else
                    cc=cmap(randsample([2 3 4 7],1),:);
                    ls='-';
                    lw=1;
                end

                for mrpt=1:numel(motif_spk{midx})
                    mxx=round((motif_spk{midx}{mrpt}-tonset-30000)./30);
                    myy=motif_seq{midx}{mrpt};
                    for tsidx=2:numel(mxx)
                        pyy=[sidx(keptnodes==myy(tsidx-1)),...
                            sidx(keptnodes==myy(tsidx))];
                        plot(mxx([tsidx-1,tsidx]),pyy,ls,'Color',cc,'LineWidth',lw);
                    end
                end
            end
            
            xlabel('Time (s)')
            ylabel('Neuron #')
            title("Session "+num2str(sessid)+" Trial "+tt+"T"+maxonset+"~"+maxoffset);
            set(gca,'YTick',1:numel(nodereg),'YTickLabel',nodereg(sidx),'XTick',0:1000:delay*1000,'XTickLabel',0:delay)
            xlim([0,delay*1000]);
            %%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
        
            %% plot
            allspk=cell2mat(cellfun(@(x) cell2mat(x),motif_spk,'UniformOutput',false));
            allid=cell2mat(cellfun(@(x) cell2mat(x),motif_seq,'UniformOutput',false));
            uspkts=unique([allspk,allid],'rows');
            uspk_sec=(uspkts(:,1)-tonset-30000)./30000;
            sumhist=histcounts(uspk_sec,0:0.25:delay);
            %% plot motif freq
            figure()
            cmap=colormap('lines');
            tiledlayout(2,1);
            nexttile()
            hold on
            for rridx=1:size(motif_spk,1)
                raster=unique(round((cell2mat(motif_spk{rridx})-tonset-30000)./60)).*2; % 60 @ 2ms resolution
                scatter(raster./1000,rridx,'|','MarkerEdgeColor',cmap(rem(rridx,7)+1,:))
            end
            ymax=numel(motif_spk)+1;
            ylim([0,ymax])
          
            for ii=1:numel(eonset)
                % disp([eonset(ii),eonset(ii)+run_length(ii)]./1000)
                plot([eonset(ii),eonset(ii)+run_length(ii)]./1000,[ymax,ymax]-0.5,'k-','LineWidth',4);
            end
            xlim([0,delay])
            xlabel("Time (s)")
            ylabel("Motif #")
            title("S"+sessid+"T"+tt)

            nexttile()
            plot(0.125:0.25:delay,sumhist,'r-','LineWidth',1);
            xlabel("Time (s)")
            ylabel("Motif spike frequency (Hz)")
            appendfig("close",true,"fn","per_trl_motif.pdf","multi",true);
        end
    end
end
end





