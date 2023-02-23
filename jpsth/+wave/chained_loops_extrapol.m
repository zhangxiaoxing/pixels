su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
per_sess_coverage=struct();
if exist('dbg','var') && dbg
else
    %% single spike chain
    sschain=load('chain_tag.mat','out');
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
end

%% per session loop
rpt=10;
thinned=struct();

for sessid=reshape(usess,1,[])
    thinned.("S"+sessid)=cell(0);
    covered=false(1000,1); % 2hrs of milli-sec bins
%     [~,~,trials,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true);
    %     FT_SPIKE.lc_tag=cell(size(FT_SPIKE.timestamp));
    %     for cidx=1:numel(FT_SPIKE.timestamp)
    %         FT_SPIKE.lc_tag{cidx}=zeros(size(FT_SPIKE.timestamp{cidx}));
    %     end

    %% single spike chain
    sesscid=su_meta.allcid(su_meta.sess==sessid);
%     nRemove=round(numel(sesscid).*[0.95,0.9,0.8,0.5]);
    nRemove=round(numel(sesscid).*[0.75,0.5,0.25]);
    for thinned_num=[nRemove,0]
        rpts=cell(0);
        for rptidx=1:rpt
            disp([sessid,thinned_num,rptidx])
            if thinned_num==0 && rptidx>1
                continue
            end
            thin_down_cid=randsample(sesscid,thinned_num);
            for dur=["d6","d3"]
                for wid=reshape(fieldnames(sschain.out.(dur)),1,[])
                    for cc=reshape(fieldnames(sschain.out.(dur).(wid{1})),1,[])
                        if ~startsWith(cc{1},['s',num2str(sessid),'c'])
                            continue
                        end
                        onechain=sschain.out.(dur).(wid{1}).(cc{1});
                        if any(ismember(onechain.meta{1},thin_down_cid),'all')
                            continue
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

            %% multi spike chain
            for dur=["d6","d3"]
                for wid=reshape(fieldnames(bschain.out.(dur)),1,[])
                    for cc=reshape(fieldnames(bschain.out.(dur).(wid{1})),1,[])
                        if ~startsWith(cc{1},['s',num2str(sessid),'c'])
                            continue
                        end
                        onechain=bschain.out.(dur).(wid{1}).(cc{1});
                        if any(ismember(onechain.meta{1},thin_down_cid),'all')
                            continue
                        end
                        % run length tag
                        for bscid=1:numel(onechain.ts)
                            onset=floor(onechain.ts{bscid}(1,3)./30);
                            offset=ceil(onechain.ts{bscid}(end,3)./30);
                            covered(onset:offset)=true;
                        end
                    end
                end
            end

            %% single spike loop

            for cc=reshape(fieldnames(pstats.congru),1,[])
                if ~startsWith(cc{1},['s',num2str(sessid),'r'])
                    continue
                end
                onechain=pstats.congru.(cc{1});
                if any(ismember(onechain.rstats{3},thin_down_cid),'all')
                    continue
                end

                % run length tag
                for tagi=reshape(setdiff(unique(onechain.ts_id(:,6)),0),1,[])
                    tseq=onechain.ts_id(onechain.ts_id(:,6)==tagi,1);
                    onset=floor(tseq(1)./30);
                    offset=ceil(tseq(end)./30);
                    %                 if offset-onset<2,keyboard();end
                    covered(onset:offset)=true;
                end
            end

            %% burst spike loop
            dbfile=fullfile("bzdata","rings_wave_burst_iter_600.db");
            conn=sqlite(dbfile,"readonly");
            for cc=reshape(bslkeys,1,[])
                if ~contains(cc,['s',num2str(sessid),'r']) || ~endsWith(cc,'_ts')
                    continue
                end
                chainmeta=table2array(conn.sqlread(replace(cc,'_ts','_meta')));
%                 keyboard()
                if any(ismember(chainmeta,thin_down_cid),'all')
                    continue
                end
                %   Not enough memory for larger complete tables
                maxrid=table2array(conn.fetch("SELECT MAX(Var1) from "+cc));
                % mini batch for performance optimization
                for rid=1:1000:maxrid
                    onechain=table2array(conn.fetch("SELECT * FROM "+cc+" WHERE Var1>="+num2str(rid)+ " AND Var1<"+num2str(rid+1000)));
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
            edges = find(diff([0;covered;0]==1));
            onset = edges(1:2:end-1);  % Start indices
            run_length = edges(2:2:end)-onset;  % Consecutive ones counts
            % per thin-down level per repeat save
            rpts=[rpts;{run_length}];
        end
        thinned.("S"+sessid)=[thinned.("S"+sessid);{[numel(sesscid),thinned_num],rpts}];
    end
end
blame=vcs.blame();
save("chainned_loops_thinned.mat","thinned","blame")



function onthefly(thinned)
    sessfn=fieldnames(thinned);
%     struct2cell(thinned)
    figure()
    hold on
    for fn=reshape(sessfn,1,[])
        onesess=thinned.(fn{1});
        nsu=[];
        for ii=1:size(onesess,1)
            pct=[];
            for jj=1:size(onesess{ii,2})
                if ~isempty(onesess{ii,2}{jj})
                    pct=[pct;median(onesess{ii,2}{jj}),max(onesess{ii,2}{jj})];
                end
            end
            if isempty(pct)
                nsu=[nsu;onesess{ii,1}(1)-onesess{ii,1}(2),0,0];
            else
                nsu=[nsu;onesess{ii,1}(1)-onesess{ii,1}(2),mean(pct,1)];
            end
            plot(nsu(:,1),nsu(:,2),'--')
            plot(nsu(:,1),nsu(:,3),'-')
        end
        keyboard()
    end

    

end


function extrapolation()
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
% keyboard();
covered=struct2cell(per_sess_coverage);
% run_length=[];
for jj=1:numel(covered)
    edges = find(diff([0;covered{jj};0]==1));
    onset = edges(1:2:end-1);  % Start indices
    disp([jj,min(edges(2:2:end)-onset)]);
    %     run_length =[run_length; edges(2:2:end)-onset];  % Consecutive ones counts
    % per_sess_coverage.("S"+num2str(sessid))=run_length;
end
end