% Pending removal
% Too simple to show as figure
% Not sure if up to date
% Require test 2023.03.26

function chained_loop_extrapolation_simple()
su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
wrs_mux_meta=ephys.get_wrs_mux_meta();

figure()
hold on
for sessid=[14,18,22,33,34,68,100,102,114]
    sess_scale=nnz(su_meta.sess==sessid & wrs_mux_meta.wave_id>0);
    load("ChainedLoop"+num2str(sessid)+".mat",'covered')
    edges = find(diff([0;covered;0]==1));
    onset = edges(1:2:end-1);  % Start indices
    run_length =edges(2:2:end)-onset;  % Consecutive ones counts
    scatter(sess_scale,prctile(run_length,25),'bo')
    scatter(sess_scale,prctile(run_length,50),'ko')
    scatter(sess_scale,prctile(run_length,75),'ro')
end
end


% per session load
% calculate
% thin down neurons
% calculate again
% statistics
% plot


% interate sessions
% load single spike chain data
% load burst spike chain data
function chained_loop_extrapolation()
optmem=false;
if ispc
    [~,sysmem]=memory();
    if sysmem.PhysicalMemory.Total < 2e10
        optmem=true;
    end
end

per_sess_coverage=struct();

if false %% new
%% single spike chain
sschain=load(fullfile('bzdata','chain_tag.mat'),'out');
keys=[struct2cell(structfun(@(x) fieldnames(x), sschain.out.d6, 'UniformOutput', false));...
    struct2cell(structfun(@(x) fieldnames(x), sschain.out.d3, 'UniformOutput', false))];
keys=vertcat(keys{:});
ssc_sess=unique(str2double(regexp(keys,'(?<=s)\d{1,3}(?=c)','match','once')));
%% multi spike
bschain=load('chain_sust_tag_600.mat','out');
keys=[struct2cell(structfun(@(x) fieldnames(x), bschain.out.d6, 'UniformOutput', false));...
    struct2cell(structfun(@(x) fieldnames(x), bschain.out.d3, 'UniformOutput', false))];
keys=vertcat(keys{:});
bsc_sess=unique(str2double(regexp(keys,'(?<=s)\d{1,3}(?=c)','match','once')));
%% single spike loop
load(fullfile('bzdata','rings_spike_trial_tagged.mat'),'pstats'); % bz.rings.rings_time_constant
pstats=rmfield(pstats,"nonmem");
ssl_sess=unique(str2double(regexp(fieldnames(pstats.congru),'(?<=s)\d{1,3}(?=r)','match','once')));
%% burst spike loop, keys only
dbfile=fullfile("bzdata","rings_wave_burst_iter_600.db");
conn=sqlite(dbfile,"readonly");
bslkeys=table2array(conn.fetch("SELECT name FROM sqlite_master WHERE type='table'"));
close(conn);
bsl_sess=unique(str2double(regexp(bslkeys,'(?<=s)\d{1,3}(?=r)','match','once')));

usess=intersect(intersect(intersect(ssc_sess,bsc_sess),ssl_sess),bsl_sess);

% single spk chn:1, burst spk chn:2, single spk loop:4, burst spk loop:8
for sessid=reshape(usess,1,[])
    covered=false(1000,1); % 2hrs of milli-sec bins
    [~,~,trials,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true);
    FT_SPIKE.lc_tag=cell(size(FT_SPIKE.timestamp));
    for cidx=1:numel(FT_SPIKE.timestamp)
        FT_SPIKE.lc_tag{cidx}=zeros(size(FT_SPIKE.timestamp{cidx}));
    end

    %% single spike chain
    for dur=["d6","d3"]
        for wid=reshape(fieldnames(sschain.out.(dur)),1,[])
            for cc=reshape(fieldnames(sschain.out.(dur).(wid{1})),1,[])
                if ~startsWith(cc{1},['s',num2str(sessid),'c'])
                    continue
                end
                onechain=sschain.out.(dur).(wid{1}).(cc{1});
                for cidx=1:size(onechain.ts,2)
                    cid=onechain.meta{1}(cidx);
%                     ucid=[ucid,cid];
                    cidsel=find(strcmp(FT_SPIKE.label,num2str(cid)));
                    totag=onechain.ts(:,cidx);
                    [ism,totagidx]=ismember(totag,FT_SPIKE.timestamp{cidsel});
                    FT_SPIKE.lc_tag{cidsel}(totagidx)=FT_SPIKE.lc_tag{cidsel}(totagidx)+1;
                end
                % run length tag
                for cidx=1:size(onechain.ts,1)
                    onset=floor(onechain.ts(cidx,1)./30);
                    offset=ceil(onechain.ts(cidx,end)./30);
                    covered(onset:offset)=true;
                end                
            end
        end
    end
    if optmem
        clear sschain
    end
%% multi spike chain
    for dur=["d6","d3"]
        for wid=reshape(fieldnames(bschain.out.(dur)),1,[])
            for cc=reshape(fieldnames(bschain.out.(dur).(wid{1})),1,[])
                if ~startsWith(cc{1},['s',num2str(sessid),'c'])
                    continue
                end
                onechain=bschain.out.(dur).(wid{1}).(cc{1});
                for cidx=1:size(onechain.meta{1},2)
                    cid=onechain.meta{1}(cidx);
%                     ucid=[ucid,cid];
                    cidsel=find(strcmp(FT_SPIKE.label,num2str(cid)));
                    for bscid=1:numel(onechain.ts)
                        totag=onechain.ts{bscid}(onechain.ts{bscid}(:,1)==cidx,3);
                        [ism,totagidx]=ismember(totag,FT_SPIKE.timestamp{cidsel});
                        FT_SPIKE.lc_tag{cidsel}(totagidx)=FT_SPIKE.lc_tag{cidsel}(totagidx)+2;
                    end
                    % run length tag
                    onset=floor(onechain.ts{bscid}(1,3)./30);
                    offset=ceil(onechain.ts{bscid}(end,3)./30);
                    covered(onset:offset)=true;
                end
            end
        end
    end

    if optmem
        clear bschain
    end

    %% single spike loop
    
    for cc=reshape(fieldnames(pstats.congru),1,[])
        if ~startsWith(cc{1},['s',num2str(sessid),'r'])
            continue
        end
        onechain=pstats.congru.(cc{1});
        for cidx=1:size(onechain.rstats{3},2)
            cid=onechain.rstats{3}(cidx);
%             ucid=[ucid,cid];
            cidsel=find(strcmp(FT_SPIKE.label,num2str(cid)));
            totag=onechain.ts_id(onechain.ts_id(:,6)>0 & onechain.ts_id(:,2)==cid,1);
            [ism,totagidx]=ismember(totag,FT_SPIKE.timestamp{cidsel});
            FT_SPIKE.lc_tag{cidsel}(totagidx)=FT_SPIKE.lc_tag{cidsel}(totagidx)+4;

            % run length tag
            for tagi=reshape(setdiff(unique(onechain.ts_id(:,6)),0),1,[])
                tseq=onechain.ts_id(onechain.ts_id(:,6)==tagi,1);
                onset=floor(tseq(1)./30);
                offset=ceil(tseq(end)./30);
%                 if offset-onset<2,keyboard();end
                covered(onset:offset)=true;
            end
        end
    end
    if optmem
        clear pstats
    end
    %% burst spike loop
    dbfile=fullfile("bzdata","rings_wave_burst_iter_600.db");
    conn=sqlite(dbfile,"readonly");
    for cc=reshape(bslkeys,1,[])
        if ~contains(cc,['s',num2str(sessid),'r']) || ~endsWith(cc,'_ts')
            continue
        end
        chainmeta=table2array(conn.sqlread(replace(cc,'_ts','_meta')));
        %   Not enough memory for larger complete tables
        maxrid=table2array(conn.fetch("SELECT MAX(Var1) from "+cc));
        for rid=1:1000:maxrid
            onechain=table2array(conn.fetch("SELECT * FROM "+cc+" WHERE Var1>="+num2str(rid)+ " AND Var1<"+num2str(rid+1000)));
            %             keyboard()
            disp(rid);
            for cidx=1:numel(chainmeta)
                cid=chainmeta(cidx);
                cidsel=find(strcmp(FT_SPIKE.label,num2str(cid)));
                totag=onechain(onechain(:,2)==cidx,4);
                [ism,totagidx]=ismember(totag,FT_SPIKE.timestamp{cidsel});
                FT_SPIKE.lc_tag{cidsel}(totagidx)=FT_SPIKE.lc_tag{cidsel}(totagidx)+8;
            end
            for rrid=rid:rid+999
                if rrid>maxrid
                    break
                end
                tseq=onechain(onechain(:,1)==rrid,4);
                onset=floor(tseq(1)./30);
                offset=ceil(tseq(end)./30);
                covered(onset:offset)=true;
            end
        end
    end
    conn.close();
    blame=vcs.blame();
    save("ChainedLoop"+num2str(sessid)+".mat",'FT_SPIKE','covered','blame')
    edges = find(diff([0;covered;0]==1));
    onset = edges(1:2:end-1);  % Start indices
    run_length =[run_length; edges(2:2:end)-onset];  % Consecutive ones counts
    per_sess_coverage.("S"+num2str(sessid))=run_length;
end

end