function chained_loops_extrapol()

su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
per_sess_coverage=struct();

load(fullfile('binary','motif_replay.mat'),'ring_replay');
load(fullfile('binary','motif_replay.mat'),'chain_replay');

usess=union(chain_replay.session,ring_replay.session);

if false
    %% multi spike
    bschain=load('chain_sust_tag_600.mat','out');
    keys=[struct2cell(structfun(@(x) fieldnames(x), bschain.out.d6, 'UniformOutput', false));...
        struct2cell(structfun(@(x) fieldnames(x), bschain.out.d3, 'UniformOutput', false))];
    keys=vertcat(keys{:});
    bsc_sess=unique(str2double(regexp(keys,'(?<=s)\d{1,3}(?=c)','match','once')));

    %% burst spike loop, keys only
    dbfile=fullfile("bzdata","rings_wave_burst_iter_600.db");
    conn=sqlite(dbfile,"readonly");
    bslkeys=table2array(conn.fetch("SELECT name FROM sqlite_master WHERE type='table'"));
    close(conn);
    bsl_sess=unique(str2double(regexp(bslkeys,'(?<=s)\d{1,3}(?=r)','match','once')));
    usess=intersect(intersect(intersect(ssc_sess,bsc_sess),ssl_sess),bsl_sess);
end

%% per session loop
rpt=100;
thinned=struct();


for sessid=reshape(usess,1,[])
    thinned.("S"+sessid)=cell(0);

    sesscid=su_meta.allcid(su_meta.sess==sessid);
    nRemove=round(numel(sesscid).*[0.9,0.75,0.5,0.25,0.1]);
    for thinned_num=[nRemove,0]
        rpts=cell(0);
        for rptidx=1:rpt
            covered=false(1000,1); % 2hrs of milli-sec bins
            disp([sessid,thinned_num,rptidx])
            if thinned_num==0 && rptidx>1
                continue
            end
            thin_down_cid=randsample(sesscid,thinned_num);

            for dur=["d6","d3"]
                %% single spike chain
                chainsel=find(chain_replay.session==sessid & endsWith(chain_replay.wave,dur));
                for chainidx=reshape(chainsel,1,[])
                    if any(ismember(chain_replay.meta{chainidx,2},thin_down_cid),'all')
                        continue
                    end
                    % run length tag
                    for cidx=1:size(chain_replay.ts{chainidx},1)
                        onset=ceil(chain_replay.ts{chainidx}(cidx,1)./30);
                        offset=ceil(chain_replay.ts{chainidx}(cidx,end)./30);
                        covered(onset:offset)=true;
                    end
                end
                %% single spike loop
                loopsel=find(ring_replay.session==sessid & endsWith(ring_replay.wave,dur));
                for loopidx=reshape(loopsel,1,[])
                    if any(ismember(ring_replay.meta{loopidx,2},thin_down_cid),'all')
                        continue
                    end
                    % run length tag
                    for cidx=1:size(ring_replay.ts{loopidx},1)
                        onset=ceil(ring_replay.ts{loopidx}{cidx}(1)./30);
                        offset=ceil(ring_replay.ts{loopidx}{cidx}(end)./30);
                        covered(onset:offset)=true;
                    end
                end

            end


            if false
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
            end

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
save(fullfile("binary","chainned_loops_thinned.mat"),"thinned","blame")


%%
sessfn=fieldnames(thinned);
fh=figure();
hold on
ylim([4,10000])
xlim([10,10000])
fitxx=[1:9,10:10:90,100:100:900,1000:1000:10000];
tofit=struct();
%     B=[];
for fn=reshape(sessfn,1,[])
    %         sessid=str2double(regexp(fn,'(?<=S)\d*','match','once'));
    %         sessreg=su_meta.reg_tree(5,su_meta.sess==sessid);

    onesess=thinned.(fn{1});
    nsu=[];
    for ii=1:size(onesess,1)
        pct=[];
        for jj=1:size(onesess{ii,2})
            if ~isempty(onesess{ii,2}{jj})
                pct=[pct;median(onesess{ii,2}{jj}),max(onesess{ii,2}{jj})];
            else
                pct=[pct;0,0];
            end
        end
        nsu=[nsu;onesess{ii,1}(1)-onesess{ii,1}(2),mean(pct,1)];
    end
    tofit.(fn{1}).vars=nsu;
    fp1=fit(nsu(:,1),nsu(:,3),'poly1');
    fplog=fit(log10(nsu(nsu(:,3)>0,1)),log10(nsu(nsu(:,3)>0,3)),'poly1');
    fpw2=fit(nsu(:,1),nsu(:,3),'power2');
    tofit.(fn{1}).poly1=fp1;
    tofit.(fn{1}).polylog=fplog;
    tofit.(fn{1}).pow2=fpw2;

    plot(nsu(:,1),nsu(:,3),'-','Color','#c0c0c0')
    tofit.(fn{1}).poly1fit=fp1(fitxx);
    tofit.(fn{1}).polylogfit=10.^(fplog(log10(fitxx)));
    tofit.(fn{1}).powr2fit=fpw2(fitxx);
end

poly1mat=cell2mat(cellfun(@(x) tofit.(x).poly1fit.',sessfn,'UniformOutput',false));
polylogmat=cell2mat(cellfun(@(x) tofit.(x).polylogfit.',sessfn,'UniformOutput',false));
pow2mat=cell2mat(cellfun(@(x) tofit.(x).powr2fit.',sessfn,'UniformOutput',false));

p1mm=mean(poly1mat);
p1sel=p1mm>0;
p1sem=std(poly1mat)./sqrt(size(poly1mat,1));
p1span=[p1mm(p1sel)+p1sem(p1sel),fliplr(p1mm(p1sel)-p1sem(p1sel))];
p1span(p1span<=0)=realmin();
fill([fitxx(p1sel),fliplr(fitxx(p1sel))],p1span,'g','EdgeColor','none','FaceAlpha',0.2);
poly1h=plot(fitxx(p1sel),p1mm(p1sel),'g-');


plogmm=mean(polylogmat);
plogsel=plogmm>0;
plogsem=std(polylogmat)./sqrt(size(polylogmat,1));
plogspan=[plogmm(plogsel)+plogsem(plogsel),fliplr(plogmm(plogsel)-plogsem(plogsel))];
plogspan(plogspan<=0)=realmin();
fill([fitxx(plogsel),fliplr(fitxx(plogsel))],plogspan,'r','EdgeColor','none','FaceAlpha',0.2);
loglogh=plot(fitxx(plogsel),plogmm(plogsel),'r-');


p2mm=mean(pow2mat);
p2sel=p2mm>0;
p2sem=std(pow2mat)./sqrt(size(pow2mat,1));
p2span=[p2mm(p2sel)+p2sem(p2sel),fliplr(p2mm(p2sel)-p2sem(p2sel))];
p2span(p2span<=0)=realmin();
fill([fitxx(p2sel),fliplr(fitxx(p2sel))],p2span,'b','EdgeColor','none','FaceAlpha',0.2);
pwrh=plot(fitxx(p2sel),p2mm(p2sel),'b-');

yline(3000,'k--')
yline(6000,'k--')
set(gca(),'XScale','log','YScale','log')
xlabel('Simultaneously observed neurons')
ylabel('Longest chained-loops pattern duration (msec)')
legend([loglogh,poly1h,pwrh],{'log-linear','linear','power-law'},'Location','northoutside','Orientation','horizontal')
%%
savefig(fh,fullfile('binary','chained_loops_extrapolation.fig'));
end