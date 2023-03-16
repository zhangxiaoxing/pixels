% Plot per-spike raster showcase of chained loops
bburst=false;

if false
    global_init;
    wrs_mux_meta=ephys.get_wrs_mux_meta();
    su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
    [map_cells,pct_bar_fh]=ephys.pct_reg_bars(wrs_mux_meta,'xyscale',{'linear','linear'},'skip_plot',true); % only need map_cells for tcom-frac corr

    % single spike chain
    sschain=load('chain_tag.mat','out');
    % multi spike
    bschain=load('chain_sust_tag_600.mat','out');
    % single spike loop
    load(fullfile('bzdata','rings_spike_trial_tagged.mat'),'pstats');
    pstats=rmfield(pstats,"nonmem");
    % burst spike loop, keys only
    dbfile=fullfile("bzdata","rings_wave_burst_iter_600.db");
    % single spk chn:1, burst spk chn:2, single spk loop:4, burst spk loop:8
end

%[ 14    18    22    33    34    68   100   102   114]
sessid=100;
% extract connected components with graph tools; option: per-trial
% dynamic graph=============================
[sig,~]=bz.load_sig_sums_conn_file('pair',false);
ssel=sig.sess==sessid;
gh=graph(cellstr(int2str(sig.suid(ssel,1))),cellstr(int2str(sig.suid(ssel,2))));
[sgbin,binsize]=conncomp(gh);
bidx=binsize(sgbin)==max(binsize);
largec=subgraph(gh,bidx);
conn_cid=str2double(table2cell(largec.Nodes));
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


if ~exist('bburst','var') || bburst
    load("ChainedLoop"+num2str(sessid)+".mat",'covered','FT_SPIKE')
else
    load("SingleSpikeChainedLoop"+num2str(sessid)+".mat",'covered','FT_SPIKE')
end
edges = find(diff([0;covered;0]==1));
onset = edges(1:2:end-1);  % Start indices
run_length = edges(2:2:end)-onset;  % Consecutive ones counts

[~,maxids]=sort(run_length,'descend');
for maxid=reshape(maxids(1:10),1,[])
    maxonset=onset(maxid).*30;
    maxoffset=edges(maxid*2).*30;
    trialid=find(maxonset-FT_SPIKE.trialinfo(:,1)>0,1,'last');
    delaylen=FT_SPIKE.trialinfo(trialid,8);
    dur="d"+delaylen;

    onsetTS=FT_SPIKE.trialinfo(trialid,1)+30000;
    offsetTS=FT_SPIKE.trialinfo(trialid,1)+FT_SPIKE.trialinfo(trialid,8).*30000+30000;

    in_trial_tags=cell(0);
    in_maximum_tags=cell(0);

    %% single spike chain
    
    for wid=reshape(fieldnames(sschain.out.(dur)),1,[]) % wave-id
        for cc=reshape(fieldnames(sschain.out.(dur).(wid{1})),1,[]) % chain-id
            if ~startsWith(cc{1},['s',num2str(sessid),'c'])
                continue
            end
            onechain=sschain.out.(dur).(wid{1}).(cc{1});
            %             onechain=onechainstr.ts_id(onechainstr.ts_id(:,5)==trialid,:)
            occursel=find(all(onechain.ts>=onsetTS & onechain.ts<offsetTS,2));
            if ~isempty(occursel)
                maxsel=find(all(onechain.ts>=maxonset & onechain.ts<maxoffset,2));
                intrialsel=setdiff(occursel,maxsel);
                %                 keyboard()
                if ~isempty(intrialsel)
                    for ii=reshape(intrialsel,1,[])
                        in_trial_tags(end+1,:)={'SSC',ii,[onechain.meta{1};onechain.ts(ii,:)].'};
                    end
                end
                if ~isempty(maxsel)
                    for ii=reshape(maxsel,1,[])
                        in_maximum_tags(end+1,:)={'SSC',ii,[onechain.meta{1};onechain.ts(ii,:)].'};
                    end
                end
            end
        end
    end
    
    %% multi spike chain
    if ~exist('bburst','var') || bburst
        for wid=reshape(fieldnames(bschain.out.(dur)),1,[])
            for cc=reshape(fieldnames(bschain.out.(dur).(wid{1})),1,[])
                if ~startsWith(cc{1},['s',num2str(sessid),'c'])
                    continue
                end
                onechain=bschain.out.(dur).(wid{1}).(cc{1});
                occursel=find(cellfun(@(x) all(x(:,3)>=onsetTS,'all') && all(x(:,3)<offsetTS,'all'),onechain.ts));
                if ~isempty(occursel)
                    for ii=reshape(occursel,1,[])
                        if all(onechain.ts{ii}(:,3)>=maxonset,'all') && all(onechain.ts{ii}(:,3)<maxoffset,'all')
                            in_maximum_tags(end+1,:)={'BSC',ii,...
                                [arrayfun(@(x) onechain.meta{1}(x), onechain.ts{ii}(:,1)),onechain.ts{ii}(:,3)]};
                        else
                            in_trial_tags(end+1,:)={'BSC',ii,...
                                [arrayfun(@(x) onechain.meta{1}(x), onechain.ts{ii}(:,1)),onechain.ts{ii}(:,3)]};
                        end
                    end
                end
            end
        end
    end

    %% single spike loop
    bracket=@(x) all(x>=onsetTS,'all') && all(x<offsetTS,'all');
    for cc=reshape(fieldnames(pstats.congru),1,[])
        if ~startsWith(cc{1},['s',num2str(sessid),'r'])
            continue
        end
        onechain=pstats.congru.(cc{1});
        occursel=setdiff(arrayfun(@(x) x.*bracket(onechain.ts_id(onechain.ts_id(:,6)==x,1)),...
            setdiff(unique(onechain.ts_id(:,6)),0)),0);
        if ~isempty(occursel)
            for ii=reshape(occursel,1,[])
                curr=onechain.ts_id(onechain.ts_id(:,6)==ii,[2,1]);
                if all(curr(:,2)>=maxonset,'all') && all(curr(:,2)<maxoffset,'all')
                    in_maximum_tags(end+1,:)={'SSL',ii,curr};
                else
                    in_trial_tags(end+1,:)={'SSL',ii,curr};
                end
            end
        end
    end

    %% burst spike loop
    if ~exist('bburst','var') || bburst
        dbfile=fullfile("bzdata","rings_wave_burst_iter_600.db");
        conn=sqlite(dbfile,"readonly");
        bslkeys=table2array(conn.fetch("SELECT name FROM sqlite_master WHERE type='table'"));
        for cc=reshape(bslkeys,1,[])
            if ~contains(cc,['s',num2str(sessid),'r']) || ~endsWith(cc,'_ts')
                continue
            end
            chainmeta=table2array(conn.sqlread(replace(cc,'_ts','_meta')));
            maxrid=table2array(conn.fetch("SELECT MAX(Var1) from "+cc));
            for rid=1:1000:maxrid
                onechain=table2array(conn.fetch("SELECT * FROM "+cc+" WHERE Var1>="+num2str(rid)+ " AND Var1<"+num2str(rid+1000)));
                occursel=setdiff(arrayfun(@(x) x.*bracket(onechain(onechain(:,1)==x,4)),unique(onechain(:,1))),0);
                if ~isempty(occursel)
                    for ii=reshape(occursel,1,[])
                        curr=[arrayfun(@(x) chainmeta(x),onechain(onechain(:,1)==ii,2)),onechain(onechain(:,1)==ii,4)];
                        if all(curr(:,2)>=maxonset,'all') && all(curr(:,2)<maxoffset,'all')
                            in_maximum_tags(end+1,:)={'BSL',ii,curr};
                        else
                            in_trial_tags(end+1,:)={'BSL',ii,curr};
                        end
                    end
                end
            end
        end
    end
    if ~bburst && (all(~strcmp(in_maximum_tags(:,1),'SSC')) || all(~strcmp(in_maximum_tags(:,1),'SSL')))
        warning("No chains or loops")
        continue
    end

    %% ===========SPIKE PLOT PREPARATION===============

    sess_cid=su_meta.allcid(su_meta.sess==sessid);
    sess_reg=su_meta.reg_tree(5,su_meta.sess==sessid).';
    ucid=unique(subsref(cell2mat(in_maximum_tags(:,3)),struct(type={'()'},subs={{':',1}})));
    [~,cidmap]=ismember(ucid,sess_cid);
    ureg=sess_reg(cidmap);

    if nnz(~ismissing(ureg)) < 10
        error("Not enough region-tagged neurons")
    end

    ucid=ucid(~ismissing(ureg));
    ureg=ureg(~ismissing(ureg));

    % regional TCOM related

    if false%bburst

        reg_prop=subsref(cell2mat(map_cells.olf.values(ureg)),struct(type={'()'},subs={{':',1}}));
        [reg_prop,tcomidx]=sort(reg_prop,'descend');
        ureg=ureg(tcomidx);
        ucid=ucid(tcomidx);

    else
        chnConn=in_maximum_tags(strcmp(in_maximum_tags(:,1),'SSC'),3);
        chnConnmat=cell2mat(cellfun(@(x) [x(1:end-1,1),x(2:end,1)],chnConn,'UniformOutput',false));
        digh=digraph(cellstr(int2str(chnConnmat(:,1))),cellstr(int2str(chnConnmat(:,2))));
        if ~digh.isdag
            warning("Is not DAG")
            continue
        end
        topocid=str2double(digh.Nodes.Name(digh.toposort));
        [~,topopos]=ismember(ucid,topocid);
        topopos(topopos==0)=numel(topocid)+1;
        reg_prop=subsref(cell2mat(map_cells.olf.values(ureg)),struct(type={'()'},subs={{':',1}}));
        [reg_prop,tcomidx]=sortrows([reg_prop,-topopos],[1,2]);
        ureg=ureg(tcomidx);
        ucid=ucid(tcomidx);


        % for gephi illustration
        lopConn=in_maximum_tags(strcmp(in_maximum_tags(:,1),'SSL'),3);
        lopConnmat=cell2mat(cellfun(@(x) [x(1:end-1,1),x(2:end,1)],lopConn,'UniformOutput',false));
        
        allconns=unique([chnConnmat;lopConnmat],'rows');
        allconns(:,3)=0;
        allconns(ismember(allconns(:,1:2),chnConnmat,'rows'),3)=1;
        lopsel=ismember(allconns(:,1:2),lopConnmat,'rows');
        allconns(lopsel,3)=allconns(lopsel,3)+2;
        csvcell=[{'Source','Target','Pattern'};num2cell(allconns)];
        writecell(csvcell,fullfile('bzdata',sprintf('SingleSpkConn4gephi_s%dm_%d.csv',sessid,maxid)));
        
        % predefined layout using gephi and Ordered-Graph-Layout
        gephilayout=jsondecode(fileread(fullfile("+gephi","SS_chain_loop_230315.json")));

        nodecell=[{'Id','Label','Toposort','REG','ChnOrd'};...
            num2cell([ucid,ucid,(numel(ucid):-1:1).',reg_prop(:,1),topopos(tcomidx)])];
        writecell(nodecell,fullfile('bzdata',sprintf('SingleSpkNode4gephi_s%dm_%d.csv',sessid,maxid)));

    end


    %% unique component count======================================
    if false
    chainspks=cellfun(@(x) x(:,1),in_maximum_tags(strcmp(in_maximum_tags(:,1),'SSC'),3),'UniformOutput',false);
    for ii=5:8
        sel=cellfun(@(x) numel(x),chainspks)==ii;
        cntsu.("L"+ii)=cell2mat(chainspks(sel).').';
        ucnt.("L"+ii)=size(unique(chainsu.("L"+ii),'rows'),1);
    end

    loopspks=cellfun(@(x) unique(x(:,1)),in_maximum_tags(strcmp(in_maximum_tags(:,1),'SSL'),3),'UniformOutput',false);
    for ii=3:5
        sel=cellfun(@(x) numel(x),loopspks)==ii;
        cntsu.("C"+ii)=cell2mat(loopspks(sel).').';
        ucnt.("C"+ii)=size(unique(cntsu.("C"+ii),'rows'),1);
    end
    end

    %% ================== Actual plot ============================

    intrialmat=unique(cell2mat(in_trial_tags(:,3)),'rows');
    inmaxmat=unique(cell2mat(in_maximum_tags(:,3)),'rows');

    figure();
    hold on;
    yy=1;
    cidyymap=containers.Map('KeyType','double','ValueType','double');
    for cid=reshape(ucid,1,[])
        cidyymap(cid)=yy;
        ft_sel=strcmp(FT_SPIKE.label,num2str(cid));

        % TODO set diff tagged spike from all spikes
        ts=intrialmat(intrialmat(:,1)==cid,2);
        [~,tsidxx]=ismember(ts,FT_SPIKE.timestamp{ft_sel});
        trialxx=FT_SPIKE.time{ft_sel}(tsidxx);
        plot(trialxx,repmat(yy,1,numel(trialxx)),'|','MarkerEdgeColor','#8080FF','LineWidth',1);

        ts=inmaxmat(inmaxmat(:,1)==cid,2);
        [~,tsidxx]=ismember(ts,FT_SPIKE.timestamp{ft_sel});
        maxxx=FT_SPIKE.time{ft_sel}(tsidxx);
        plot(maxxx,repmat(yy,1,numel(maxxx)),'|','MarkerEdgeColor','#FF0000','LineWidth',2);

        xx=FT_SPIKE.time{ft_sel}(FT_SPIKE.trial{ft_sel}==trialid);
        xx=setdiff(xx,union(trialxx,maxxx));
        plot(xx,repmat(yy,1,numel(xx)),'|','MarkerEdgeColor',"#808080")
        yy=yy+1;
    end
    % xlim([0,FT_SPIKE.trialinfo(trialid,8)]);
    % xlim([0,7]);
    xlabel('Time (s)')
    ylabel('Neuron #')
    title("Session "+num2str(sessid)+" Trial "+num2str(trialid));
    set(gca,'YTick',1:numel(ureg),'YTickLabel',ureg)

    %% plot assembly
    already=[];
    xspan=[7,0];
    cc=colormap('lines');

    for ii=1:size(in_maximum_tags,1)
        curr=in_maximum_tags{ii,3};
        currfc=cell2mat(arrayfun(@(jj) [curr(jj,2),curr(jj,1),curr(jj+1,2),curr(jj+1,1)],(1:size(curr,1)-1).','UniformOutput',false));
        % skip unnecessary
        if ~isempty(already) && all(ismember(currfc,already,'rows'))
            disp("skipped "+num2str(ii))
            continue
        end
        already=unique([already;currfc],'rows');

        plotmat=[];
        for jj=1:size(curr,1)-1
            fromy=cidyymap(curr(jj,1));
            fromsel=strcmp(FT_SPIKE.label,num2str(curr(jj,1)));
            fromx=FT_SPIKE.time{fromsel}(FT_SPIKE.timestamp{fromsel}==curr(jj,2));

            toy=cidyymap(curr(jj+1,1));
            tosel=strcmp(FT_SPIKE.label,num2str(curr(jj+1,1)));
            tox=FT_SPIKE.time{tosel}(FT_SPIKE.timestamp{tosel}==curr(jj+1,2));
            plotmat=[plotmat;fromx,tox,fromy,toy];
        end
        if ismember(in_maximum_tags{ii,1},{'SSL','BSL'})  % loop
            plot(plotmat(:,1:2).',plotmat(:,3:4).','--','Color',cc(randsample([2 3 4 7],1),:));
            xspan=[min(xspan(1),min(plotmat(:,1))),max(xspan(2),max(plotmat(:,2)))];
        else
            plot(plotmat(:,1:2).',plotmat(:,3:4).','-','Color',cc(randsample([1 5 6],1),:),'LineWidth',1.5);
        end
    end
    % set(gca(),'XTick',5:0.2:6.2,'XTickLabel',0:0.2:1.2)
    if bburst
        xlim(xspan+[-0.02,0.02]);
        set(gca(),'XTick',xspan(1)-0.02:0.5:xspan(1)+2,'XTickLabel',0:0.5:2);
    else
        xlim(xspan+[-0.01,0.01]);
        set(gca(),'XTick',xspan(1)-0.01:0.05:xspan(2)+0.21,'XTickLabel',0:0.05:0.2);
    end
    ylim([0.5,yy-0.5])
    keyboard();
end





