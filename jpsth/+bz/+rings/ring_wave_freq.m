% calculate relative appearance frequence of rings vs shuffled control
function [fh,zscore,data_value,shuf_value]=ring_wave_freq(sel_meta,opt)
arguments
    sel_meta = []
    opt.burst (1,1) logical = false
    opt.repeats (1,1) double = 100
    opt.denovo (1,1) logical = true
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
    opt.odor_only (1,1) logical
    opt.poolsize (1,1) double = 2
end
%% appearance in spike sequence
if isempty(sel_meta)
    switch opt.criteria
        case 'WT'
            fstr=load(fullfile('binary','wrs_mux_meta.mat'));
            sel_meta=fstr.wrs_mux_meta;
            clear fstr
        case 'Learning'
            sel_meta=ephys.get_wrs_mux_meta('load_file',false,'save_file',false,'criteria','Learning','extend6s',true);
        case 'any'
            keyboard()
    end
end

if opt.burst
    assert(strcmp(opt.criteria,'WT'),"Unfinished");
    % nonmem
    nonmem_keys=struct();
    for bi=[150,300,600]
        nonmem_keys.("B"+bi)=[];
        dbfile=fullfile("bzdata","rings_nonmem_burst_iter_"+bi+".db");
        conn=sqlite(dbfile,'readonly');
        bikeys=table2array(conn.fetch("SELECT name FROM sqlite_master WHERE type='table'"));
        metakeys=bikeys(endsWith(bikeys,"_meta"));
        for mki=1:numel(metakeys)
            sessionstr=regexp(metakeys(mki),'s\d+(?=r)','match','once');
            cuid=sort(table2array(conn.sqlread(metakeys(mki))));
            onekey=sessionstr+"-"+string(sprintf('%d-',cuid));
            nonmem_keys.("B"+bi)=[nonmem_keys.("B"+bi);onekey];
        end
        conn.close();
    end
    % congru
    congru_keys=struct();
    for bi=[150,300,600]
        congru_keys.("B"+bi)=[];
        dbfile=fullfile("bzdata","rings_wave_burst_iter_"+bi+".db");
        conn=sqlite(dbfile,'readonly');
        bikeys=table2array(conn.fetch("SELECT name FROM sqlite_master WHERE type='table'"));
        metakeys=bikeys(endsWith(bikeys,"_meta"));
        for mki=1:numel(metakeys)
            sessionstr=regexp(metakeys(mki),'s\d+(?=r)','match','once');
            cuid=sort(table2array(conn.sqlread(metakeys(mki))));
            onekey=sessionstr+"-"+string(sprintf('%d-',cuid));
            congru_keys.("B"+bi)=[congru_keys.("B"+bi);onekey];
        end
        conn.close();
    end
    if opt.denovo
        sums_shuf=cell(opt.repeats,1);
        if opt.poolsize>1
            poolh=parpool(opt.poolsize);
            F=parallel.FevalFuture.empty(0,opt.repeats);
            for rpt=1:opt.repeats
                F(rpt)=parfeval(poolh,@bz.rings.rings_wave,1,sel_meta,'shufid',rpt,'odor_only',opt.odor_only,'criteria',opt.criteria);
            end

            while ~all([F.Read],'all')
                try
                    [pidx,oneshuf]=fetchNext(F);
                    sums_shuf{pidx}=oneshuf;
                    remain=find(~strcmp({F.State},'finished'));
                    disp(remain);
                catch ME
                    blame=vcs.blame();
                    fstate={F.State};
                    save(fullfile('binary/ring_wave_freq_err.mat'),'fstate','ME','blame','opt');
                    quit(34);
                end
            end
        else
            for rpt=1:opt.repeats
                sums_shuf{rpt}=bz.rings.rings_wave(sel_meta,'shufid',rpt,'odor_only',opt.odor_only,'criteria',opt.criteria);
            end
        end
    else
        assert(strcmp(opt.criteria,'WT'),"Unfinished");
        load(fullfile('bzdata','SS_loop_count_shuf.mat'),'sums_shuf')
    end


    data_value=struct();
    for bi=[150 300 600]
        data_value.nonmem.("B"+bi)=numel(unique(nonmem_keys.("B"+bi)));
        data_value.congru.("B"+bi)=numel(unique(congru_keys.("B"+bi)));
    end

    data_value.congru.shuf=cellfun(@(x) nnz(strcmp(x(:,4),'congru')), sums_shuf(1:opt.repeats));
    data_value.nonmem.shuf=cellfun(@(x) nnz(strcmp(x(:,4),'nonmem')), sums_shuf(1:opt.repeats));
    
    zscore=struct();
    for type={'congru','nonmem'}
%         shuffmm=mean(wavetype.(type{1}).shuf);
%         shufstd=std(wavetype.(type{1}).shuf);
            [~,bootsam]=bootstrp(1000,[],data_value.(type{1}).shuf);
            shuffmm=mean(data_value.(type{1}).shuf(bootsam));
            shufstd=std(data_value.(type{1}).shuf(bootsam));
        for bi=[150 300 600]
            zscore.(type{1}).("B"+bi)=(data_value.(type{1}).("B"+bi)-shuffmm)./shufstd;
        end
    end
    zscoremat=[zscore.congru.B150;zscore.congru.B300;zscore.congru.B600;...
        zscore.nonmem.B150;zscore.nonmem.B300;zscore.nonmem.B600];

    zscoremm=mean(zscoremat,2);
    zscorestd=std(zscoremat,0,2);
    zscoresem=zscorestd./sqrt(1000);% TODO: magic numbers everywhere

    ranksum(zscoremat(1,:),zscoremat(4,:))
    ranksum(zscoremat(2,:),zscoremat(5,:))
    ranksum(zscoremat(3,:),zscoremat(6,:))

    fh=figure('Color','w','Position',[32,32,215,215]);
    hold on
%     bh=bar(1:3,diag([zscore.nonmem,zscore.incongru,zscore.congru]),'stacked','FaceColor','k','EdgeColor','k');
    bh=bar([zscoremm(1:3),...
        zscoremm(4:6)]);
    bh(1).FaceColor='w';
    bh(2).FaceColor='k';
%     errorbar
    errorbar([bh(1).XEndPoints,bh(2).XEndPoints],zscoremm,zscoresem,'k.')
    ylabel('Z-Score')
    set(gca(),'XTick',1:3,'XTickLabel',{'200Hz','100Hz','50Hz'})
    legend(bh,{'Congruent','Non-memory'},'Location','northoutside','Orientation','horizontal')
    title('Burst rings occurrence')

else % W/o burst
    if opt.denovo
        sums_shuf=cell(opt.repeats,1);
        if opt.poolsize>1
            poolh=parpool(opt.poolsize);
            F0=parfeval(poolh,@bz.rings.rings_wave,1,sel_meta,'shufid',0,'odor_only',opt.odor_only,'criteria',opt.criteria);
            F=parallel.FevalFuture.empty(0,opt.repeats);
            for rpt=1:opt.repeats
                F(rpt)=parfeval(poolh,@bz.rings.rings_wave,1,sel_meta,'shufid',rpt,'odor_only',opt.odor_only,'criteria',opt.criteria);
            end
            sums=F0.fetchOutputs();
            while ~all([F.Read],'all')
                try
                    [pidx,oneshuf]=fetchNext(F);
                    sums_shuf{pidx}=oneshuf;
                    remain=find(~strcmp({F.State},'finished'));
                    disp(remain);
                catch ME
                    blame=vcs.blame();
                    fstate={F.State};
                    save(fullfile('binary/ring_wave_freq_err.mat'),'fstate','ME','blame','opt');
                    quit(34);
                end
            end
            delete(poolh)
        else
            sums=bz.rings.rings_wave(sel_meta,'shufid',0,'criteria',opt.criteria,'odor_only',opt.odor_only);
            for rpt=1:opt.repeats
                sums_shuf{rpt}=bz.rings.rings_wave(sel_meta,'shufid',rpt,'odor_only',opt.odor_only,'criteria',opt.criteria);
            end
        end
    else
        assert(strcmp(opt.criteria,'WT'),"Unfinished")
        % probably outdated, needs varify
        load(fullfile('bzdata','SS_loop_count_shuf.mat'),'sums_shuf','sums') 
    end
    data_value=struct();
    shuf_value=struct();
    [shuf_value.congru,shuf_value.incongru,shuf_value.nonmem]=deal([]);
    for type={'congru','incongru','nonmem'}
        data_value.(type{1})=nnz(strcmp(sums(:,4),type));
        for rpt=1:opt.repeats
            shuf_value.(type{1})=[shuf_value.(type{1}),nnz(strcmp(sums_shuf{rpt}(:,4),type))];
        end
    end

    zscore=struct();
    for type={'congru','incongru','nonmem'}
        % implementation of permutation+test
        [~,bootsam]=bootstrp(1000,[],shuf_value.(type{1}));
        shuffmm=mean(shuf_value.(type{1})(bootsam));
        shufstd=std(shuf_value.(type{1})(bootsam));
        zscore.(type{1})=(data_value.(type{1})-shuffmm)./shufstd;
    end
    zscoremat=cell2mat(struct2cell(zscore));
    zscoremm=arrayfun(@(x) mean(zscoremat(x,isfinite(zscoremat(x,:)))),1:size(zscoremat,1)); % likely under estimate due to zeros in shuffled samples
    zscorestd=arrayfun(@(x) std(zscoremat(x,isfinite(zscoremat(x,:)))),1:size(zscoremat,1));
    zscoresem=zscorestd./sqrt(1000);

    ranksum(zscoremat(1,:),zscoremat(3,:))
    ranksum(zscoremat(2,:),zscoremat(3,:))

    if false
        loopcount=cell2struct({data_value.nonmem;data_value.incongru;data_value.congru; ...
            shuf_value.nonmem;shuf_value.incongru;shuf_value.congru},...
            {'Nonmemory_observed','Incongruent_observed','Congrent_observed','Nonmemory_shuffled','Incongruent_shuffled','Congruent_shuffled'})
       fid=fopen(fullfile('binary','upload','F2L_loop_consistent_incosistent.json'),'w');
       fprintf(fid,jsonencode(loopcount));
       fclose(fid);

    end


    fh=figure('Color','w','Position',[32,32,215,215]);
    hold on
    bh=bar(1:3,diag(flip(zscoremm)),'stacked','FaceColor','k','EdgeColor','k');
    bh(2).FaceColor='b';
    bh(3).FaceColor='r';
    errorbar(3:-1:1,zscoremm,zscoresem,'k.')
    % errorbar(1:4,mm,sems,'k.')
    ylabel('Z-Score')
    set(gca(),'XTick',1:3,'XTickLabel',{'Nonmem','Incongru.','Congru.'},'XTickLabelRotation',45)
    switch opt.criteria
        case 'WT'
            savefig(fh,fullfile('binary','loop_occur_vs_shuf.fig'));
        case 'Learning'
            savefig(fh,fullfile('binary','LN_loop_occur_vs_shuf.fig'));
        otherwise
            error("Unfinished")
    end

    if false
        seltype=struct();
        seltype_shuf=struct();
        [seltype_shuf.olf,seltype_shuf.dur]=deal([]);
        for type={'olf','dur'}
            seltype.(type{1})=nnz(strcmp(sums(:,5),type));
            for rpt=1:opt.repeats
                seltype_shuf.(type{1})=[seltype_shuf.(type{1}),nnz(strcmp(sums_shuf{rpt}(:,5),type))];
            end
        end

        sel_zscore=struct();
        for type={'olf','dur'}
            shuffmm=mean(seltype_shuf.(type{1}));
            shufstd=std(seltype_shuf.(type{1}));
            sel_zscore.(type{1})=(seltype.(type{1})-shuffmm)./shufstd;
        end


        fh=figure('Color','w','Position',[32,32,215,215]);
        hold on
        bh=bar(1:2,diag([sel_zscore.olf,sel_zscore.dur]),'stacked','FaceColor','k','EdgeColor','w');
        bh(2).FaceColor='b';
        yline(1,'k:')
        % errorbar(1:4,mm,sems,'k.')
        ylabel('Z-Score')
        set(gca(),'XTick',1:2,'XTickLabel',{'Olfactory','Duration'},'XTickLabelRotation',45,'YScale','log')
    end
end
end