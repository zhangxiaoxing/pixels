% tag spikes in neuron with loop or chain activity
% olf, enc-both disconnect, composite proportion, bargraph
% TODO: overall refactoring

function [disconnected,degree_sums,complexity_sums]=module_motif_asso_composite(sschain,pstats,opt)
arguments
    sschain
    pstats
    opt.readfile=false;
    opt.skipfile=true;
end

stats=struct();
per_trl_nodes=cell(0);
complexity_sums=[];
degree_sums=[];


if true%~exist('inited','var') || ~inited  % denovo data generation

    %% single spike chain
    if opt.readfile
        sschain=load(fullfile('bzdata','chain_tag.mat'),'out');
        load(fullfile('bzdata','rings_spike_trial_tagged.mat'),'pstats'); % bz.rings.rings_time_constant
    end
    keys=[struct2cell(structfun(@(x) fieldnames(x), sschain.out.d6, 'UniformOutput', false));...
        struct2cell(structfun(@(x) fieldnames(x), sschain.out.d3, 'UniformOutput', false))];
    keys=vertcat(keys{:});
    ssc_sess=unique(str2double(regexp(keys,'(?<=s)\d{1,3}(?=c)','match','once')));

    %% single spike loop
    if isfield(pstats,'nonmem'), pstats=rmfield(pstats,"nonmem");end
    ssl_sess=unique(str2double(regexp(fieldnames(pstats.congru),'(?<=s)\d{1,3}(?=r)','match','once')));
    usess=intersect(ssc_sess,ssl_sess);

    % single spk chn:1, burst spk chn:2, single spk loop:4, burst spk loop:8

    %% per-session entrance
    per_trial_motif_cid=cell(0);
    per_trial_motif_freq=[];

    stats.chain=cell2struct({cell(0);cell(0);cell(0)},{'olf','dur','both'},1);
    stats.loop=cell2struct({cell(0);cell(0);cell(0)},{'olf','dur','both'},1);
    for sessid=reshape(usess,1,[])
        disp(sessid)
        [~,~,trials,~,~,~]=ephys.getSPKID_TS(sessid,'keep_trial',false);

        wtsel=trials(:,9)>0 & trials(:,10)>0 & ismember(trials(:,5),[4 8]) &ismember(trials(:,8),[3 6]);
        per_trial_motif_cid=[per_trial_motif_cid;cell(nnz(wtsel),3)];

        per_trial_motif_freq=[per_trial_motif_freq;...
            repmat(sessid,nnz(wtsel),1),find(wtsel),trials(wtsel,5),trials(wtsel,8),...
            zeros(nnz(wtsel),4)];

        %% single spike chain

        for dur=["d6","d3"]
            dd=str2double(replace(dur,"d",""));
            for wid=reshape(fieldnames(sschain.out.(dur)),1,[])
                for cc=reshape(fieldnames(sschain.out.(dur).(wid{1})),1,[])
                    if ~startsWith(cc{1},['s',num2str(sessid),'c'])
                        continue
                    end
                    onechain=sschain.out.(dur).(wid{1}).(cc{1});
                    % ts_id: ts, cid, pos, trial_time, trial

                    % per-trial frequency.
                    tsel=[];
                    if dd==3
                        switch wid{1}
                            case 's1d3'
                                tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==4 & per_trial_motif_freq(:,4)==3;
                                ttag='both';
                            case 's2d3'
                                tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==8 & per_trial_motif_freq(:,4)==3;
                                ttag='both';
                            case 'olf_s1'
                                tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==4 & per_trial_motif_freq(:,4)==3;
                                ttag='olf';
                            case 'olf_s2'
                                tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==8 & per_trial_motif_freq(:,4)==3;
                                ttag='olf';
                            case 'dur_d3'
                                tsel=per_trial_motif_freq(:,1)==sessid & ismember(per_trial_motif_freq(:,3),[4 8]) & per_trial_motif_freq(:,4)==3;
                                ttag='dur';
                        end
                    else
                        switch wid{1}
                            case 's1d6'
                                tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==4 & per_trial_motif_freq(:,4)==6;
                                ttag='both';
                            case 's2d6'
                                tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==8 & per_trial_motif_freq(:,4)==6;
                                ttag='both';
                            case 'olf_s1'
                                tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==4 & per_trial_motif_freq(:,4)==6;
                                ttag='olf';
                            case 'olf_s2'
                                tsel=per_trial_motif_freq(:,1)==sessid & per_trial_motif_freq(:,3)==8 & per_trial_motif_freq(:,4)==6;
                                ttag='olf';
                            case 'dur_d6'
                                tsel=per_trial_motif_freq(:,1)==sessid & ismember(per_trial_motif_freq(:,3),[4 8]) & per_trial_motif_freq(:,4)==6;
                                ttag='dur';
                        end
                    end
                    stats.chain.(ttag)=[stats.chain.(ttag);{sessid},onechain.meta(1)];
                    if ~isempty(tsel)
                        for ttt=reshape(find(tsel),1,[])
                            per_trial_motif_cid{ttt,3}=[sessid,per_trial_motif_freq(ttt,2),str2double(replace(dur,"d","")),ttt];
                            per_trial_motif_cid{ttt,1}=[per_trial_motif_cid{ttt,1},onechain.meta(1)];
                        end
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
            % .ts_id(:,6) => loop tag
            % per-trial frequency.
            [pref3,pref6]=bz.rings.preferred_trials_rings(onechain.rstats{4},trials);
            tsel=per_trial_motif_freq(:,1)==sessid & ismember(per_trial_motif_freq(:,2),[pref3;pref6]);
            % per wave
            %             if strcmp(onechain.rstats{5},'olf')
            %                 stats.loop.olf=
            %             elseif strcmp(onechain.rstats{5},'dur')
            %
            %             end
            if ismember(onechain.rstats{5},{'olf','dur','both'})
                stats.loop.(onechain.rstats{5})=[stats.loop.(onechain.rstats{5});{sessid},onechain.rstats(3)];
            else
                keyboard()
            end

            for ttt=reshape(find(tsel),1,[])
                per_trial_motif_cid{ttt,3}=[sessid,per_trial_motif_freq(ttt,2),str2double(replace(dur,"d","")),ttt];
                per_trial_motif_cid{ttt,2}=[per_trial_motif_cid{ttt,2},onechain.rstats(3)];
            end
            %check cid in largest network
        end
    end
    if ~opt.skipfile
        blame=vcs.blame();
        save(fullfile("bzdata","SingleSpikeChainedLoopCid.mat"),'per_trial_motif_cid','per_trial_motif_freq','blame','stats')
    end
else
    load(fullfile("bzdata","SingleSpikeChainedLoopCid.mat"),'per_trial_motif_cid','per_trial_motif_freq','blame')
end

processed=cell(0);
disconnected=cell(0);

for tt=1:size(per_trial_motif_cid,1)
    if isempty(per_trial_motif_cid{tt,1}) && isempty(per_trial_motif_cid{tt,2})
        continue
    end

    if numel(per_trial_motif_cid{tt,1})+numel(per_trial_motif_cid{tt,2})==1
        tkey=sprintf('%d-',per_trial_motif_freq(tt,1),cell2mat([per_trial_motif_cid{tt,1:2}]));
        if ~ismember(tkey,processed)
            disconnected=[disconnected;{per_trial_motif_freq(tt,1)},per_trial_motif_cid{tt,1:2}];
            processed=[processed;tkey];
        end
        continue
    end

    % build graph network
    edges=categorical(unique(cell2mat(cellfun(@(x) [x(1:end-1);x(2:end)].',[per_trial_motif_cid{tt,1:2}],'UniformOutput',false).'),'rows'));
    gh=graph(edges(:,1),edges(:,2));
    %   Moved to enforce composite requirement
    %   TODO: optional switch?
    %   degree_sums=[degree_sums;repmat(per_trial_motif_cid{tt,3},gh.numnodes,1),str2double(table2array(gh.Nodes)),gh.degree];
    conncomp=gh.conncomp();
    % TODO: digraph?
    % check module in network
    if any(conncomp~=1)
        comps=unique(conncomp);
        counter=zeros(numel(comps),1);
        gnodes=cellfun(@(x) str2double(x),gh.Nodes.Name);
        for mm=[per_trial_motif_cid{tt,1:2}]
            [~,nidx]=ismember(mm{1},gnodes);
            compidx=unique(conncomp(nidx));
            if numel(compidx)==1
                counter(compidx)=counter(compidx)+1;
            else % should not happen
                keyboard()
            end
        end
        if ~all(counter>1)
            for ccid=reshape(find(counter==1),1,[])
                tkey=sprintf('%d-',per_trial_motif_freq(tt,1),gnodes(conncomp==ccid));
                if ~ismember(tkey,processed)
                    disconnected=[disconnected;{per_trial_motif_freq(tt,1)},{gnodes(conncomp==ccid)}];
                    processed=[processed;tkey];
                end
            end
        end
        subs=reshape(find(counter>1),1,[]);
    else
        subs=1;
    end

    for subsidx=subs
        subgh=gh.subgraph(conncomp==subsidx);
        complexity_sums=[complexity_sums;per_trial_motif_cid{tt,3},subgh.numnodes,subgh.numedges,subgh.numedges./nchoosek(subgh.numnodes,2),max(max(subgh.distances))];
        per_trl_nodes=[per_trl_nodes;{per_trial_motif_cid{tt,3}(1),str2double(table2array(subgh.Nodes))}];
        degree_sums=[degree_sums;repmat(per_trial_motif_cid{tt,3},subgh.numnodes,1),str2double(table2array(subgh.Nodes)),subgh.degree];
    end
end

if false
    save(fullfile('bzdata','motifs_graph_stats.mat'),"disconnected","degree_sums","complexity_sums","blame")
end
disconn_count=[];
for di=1:size(disconnected,1)
    if any(cell2mat(stats.chain.olf(:,1))==disconnected{di,1} & cellfun(@(x) all(ismember(x,disconnected{di,2})),stats.chain.olf(:,2)))
        disconn_count=[disconn_count;"olfchain"];
    elseif any(cell2mat(stats.chain.both(:,1))==disconnected{di,1} & cellfun(@(x) all(ismember(x,disconnected{di,2})),stats.chain.both(:,2)))
        disconn_count=[disconn_count;"bothchain"];
    elseif any(cell2mat(stats.loop.olf(:,1))==disconnected{di,1} & cellfun(@(x) all(ismember(x,disconnected{di,2})),stats.loop.olf(:,2)))
        disconn_count=[disconn_count;"olfloop"];
    elseif any(cell2mat(stats.loop.both(:,1))==disconnected{di,1} & cellfun(@(x) all(ismember(x,disconnected{di,2})),stats.loop.both(:,2)))
        disconn_count=[disconn_count;"bothloop"];
    end
end


figure()

tiledlayout(1,2)
mcount=struct();
for wv=["olf","both"] % dur
    nexttile()
    hold on
    for motif=["chain","loop"]
        ulist=unique(arrayfun(@(kii) string(sprintf('%d-',stats.(motif).(wv){kii,1},sort(stats.(motif).(wv){kii,2}))),1:size(stats.(motif).(wv),1)));
        if isempty(disconn_count)
            disconn_n=0;
        else
            disconn_n=nnz(contains(disconn_count,wv) & contains(disconn_count,motif));
        end
        mcount.(wv).(motif)=[numel(ulist),disconn_n];
    end
    mm=[(mcount.(wv).chain(1)-mcount.(wv).chain(2))./mcount.(wv).chain(1),...
        (mcount.(wv).loop(1)-mcount.(wv).loop(2))./mcount.(wv).loop(1)];
    sem=[sqrt(mm(1).*(1-mm(1))./mcount.(wv).chain(1)),sqrt(mm(2).*(1-mm(2))./mcount.(wv).loop(1))];

    bh=bar(1:2,mm,'FaceColor','none');
    errorbar(1:2,bh.YData,sem,'k.');
    set(gca(),'XTick',1:2,'XTicklabel',{'Chain','Loop'},'YTick',0:0.5:1,'YTickLabel',0:50:100)
    xlim([0.25,2.75]);
    title(wv)
    ylabel('Associated with composite loops (%)')
    disp(wv+string(mm))
end

%% plot per region degree
su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
uuid=unique(degree_sums(:,[1,5]),'rows');
per_reg=struct();
missing_reg=[];
for uididx=1:size(uuid,1)
    % average cross trial
    su_trl_sel=(degree_sums(:,1)==uuid(uididx,1) & degree_sums(:,5)==uuid(uididx,2));
    su_mm=mean(degree_sums(su_trl_sel,6));
    % region
    metaidx=find(su_meta.sess==uuid(uididx,1) & su_meta.allcid==uuid(uididx,2));
    if ~isempty(metaidx)
        sureg=su_meta.reg_tree(5,metaidx);
    else

        error("ERROR");
    end

    % per region data add
    if ~ismissing(sureg)
        if ~isfield(per_reg,sureg{1})
            per_reg.(sureg{1})=[];
        end
        per_reg.(sureg{1})=[per_reg.(sureg{1});su_mm];
    else
        missing_reg=[missing_reg;su_mm];
    end
end
%stats
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
figure()
scatter(uniq_net(:,3),uniq_net(:,4),'ko','filled')
xlabel('Neuron (node) number')
ylabel('FC (edge) number')
set(gca,'XScale','log','YScale','log','XTick',[5:5:40])
xlim([4,50])
ylim([4,250])
title('SS composite loops neuron vs FC')


figure('Position',[100,100,250,250])
tiledlayout(1,2)
nexttile
hold on
boxplot(complexity_sums(:,5),'Colors','k','Whisker',inf,'Widths',0.5)
set(gca,'YScale','log')
ylim([1,1000])
xlim([0.5,1.5])
ylabel('Number of neurons')
nexttile
hold on
boxplot(complexity_sums(:,6),'Colors','k','Whisker',inf,'Widths',0.5)
set(gca,'YScale','log')
ylim([1,1000])
xlim([0.5,1.5])
ylabel('Number of FCs')


figure('Position',[100,100,150,400])
boxplot(uniq_net(:,5),'Colors','k','Whisker',inf,'Widths',1)
ylim([0,1])
set(gca,'XTick',[],'YTick',0:0.25:1,'YTickLabel',0:25:100)
ylabel('Network density')
end