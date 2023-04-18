% calculate relative appearance frequence of rings vs shuffled control
function ring_wave_freq(wrs_mux_meta,opt)
arguments
    wrs_mux_meta
    opt.appearance
    opt.repeats (1,1) double = 10
end

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

nonmem_keys=struct();
for bi=[150,300,600]
    nonmem_keys.("B"+bi)=[];
    dbfile=fullfile("bzdata","rings_wave_burst_iter_"+bi+".db");
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


sums=bz.rings.rings_wave(wrs_mux_meta,'shufid',0);
sums_shuf=cell(opt.repeats,1);
for rpt=1:opt.repeats
    sums_shuf{rpt}=bz.rings.rings_wave(wrs_mux_meta,'shufid',rpt);
end
wavetype=struct();
wavetype_shuf=struct();
[wavetype_shuf.congru,wavetype_shuf.incongru,wavetype_shuf.nonmem]=deal([]);
for type={'congru','incongru','nonmem'}
    wavetype.(type{1})=nnz(strcmp(sums(:,4),type));
    for rpt=1:opt.repeats
        wavetype_shuf.(type{1})=[wavetype_shuf.(type{1}),nnz(strcmp(sums_shuf{rpt}(:,4),type))];
    end
end

zscore=struct();
for type={'congru','incongru','nonmem'}
    shuffmm=mean(wavetype_shuf.(type{1}));
    shufstd=std(wavetype_shuf.(type{1}));
    zscore.(type{1})=(wavetype.(type{1})-shuffmm)./shufstd;
end


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
bh=bar(1:3,diag([zscore.nonmem,zscore.incongru,zscore.congru]),'stacked','FaceColor','k','EdgeColor','k');
bh(2).FaceColor='b';
bh(3).FaceColor='r';
% errorbar(1:4,mm,sems,'k.')
ylabel('Z-Score')
set(gca(),'XTick',1:3,'XTickLabel',{'Nonmem','Incongru.','Congru.'},'XTickLabelRotation',45)

fh=figure('Color','w','Position',[32,32,215,215]);
hold on
bh=bar(1:2,diag([sel_zscore.olf,sel_zscore.dur]),'stacked','FaceColor','k','EdgeColor','w');
bh(2).FaceColor='b';
yline(1,'k:')
% errorbar(1:4,mm,sems,'k.')
ylabel('Z-Score')
set(gca(),'XTick',1:2,'XTickLabel',{'Olfactory','Duration'},'XTickLabelRotation',45,'YScale','log')

end