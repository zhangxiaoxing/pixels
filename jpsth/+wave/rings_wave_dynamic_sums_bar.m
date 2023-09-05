% function  rings_wave_dynamic()
global_init();
su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
wrs_mux_meta=ephys.get_wrs_mux_meta();

skip_lowlvl_burst=true;

% Statistics for 3- 4- and 5-SU rings time constant removed on 2023.2.22
sumdata=struct();
%% individual loops ================================================
% with preselection
sumdata.ssl=[];
load(fullfile('bzdata','rings_spike_trial_tagged.mat'),'pstats','blame')
for fn=reshape(fieldnames(pstats.congru),1,[])
    ts=pstats.congru.(fn{1}).ts_id(:,[1 6]);
    for onel=setdiff(reshape(unique(ts(:,2)),1,[]),0)
        ll=ts(ts(:,2)==onel,1);
        sumdata.ssl=[sumdata.ssl;(ll(end)-ll(1))./30];
    end
end


%% burst loops ===============================
if false % skip lower level burst stats at the moment
    % TODO: update to new SQLite format
    load('rings_wave_burst_600.mat','out');
    perchaindur=struct();
    [perchaindur.d6.size,perchaindur.d6.dur,perchaindur.d3.size,perchaindur.d3.dur,perchaindur.d6.int,perchaindur.d3.int]=deal([]);
    for dur=reshape(fieldnames(out),1,[])
        %     [perchaindur.size,perchaindur.dur]=deal([]);
        for wv=reshape(fieldnames(out.(dur{1})),1,[])
            for lp=reshape(fieldnames(out.(dur{1}).(wv{1})),1,[])
                %             keyboard();
                perchaindur.(dur{1}).size=[perchaindur.(dur{1}).size,cellfun(@(x) size(x,1),out.(dur{1}).(wv{1}).(lp{1}).ts)];
                perchaindur.(dur{1}).dur=[perchaindur.(dur{1}).dur,cellfun(@(x) diff(x([1,end],3),1,1),out.(dur{1}).(wv{1}).(lp{1}).ts)./30];
                perchaindur.(dur{1}).int=[perchaindur.(dur{1}).int,cell2mat(cellfun(@(x) diff(x(:,3),1,1).',out.(dur{1}).(wv{1}).(lp{1}).ts,'UniformOutput',false))./30];
            end
        end
        statss.("d"+dur)=perchaindur;
    end
    bsl_=[perchaindur.d6.dur,perchaindur.d3.dur].';
end

%% individual FC ====================
load(fullfile("binary","sums_conn_10.mat"),'sums_conn_str');
qcmat=cell2mat({sums_conn_str.qc}.'); % from bz.sums_conn -> bz.goodccg
sumdata.FC=(qcmat(:,2)-250)./30; % #2 is peak TS, offset left half of symmatric ccg


%% wave bump width 
% skipped for now

bump3=[];
bump6=[];
com_map=wave.get_pct_com_map(wrs_mux_meta,'curve',true,'early_smooth',false);
for fn=reshape(fieldnames(com_map),1,[])
    for subfn=reshape(fieldnames(com_map.(fn{1})),1,[])
        if isfield(com_map.(fn{1}).(subfn{1}),'fwhm3')
            bump3=[bump3,cell2mat(com_map.(fn{1}).(subfn{1}).fwhm3.values())];
        end
        if isfield(com_map.(fn{1}).(subfn{1}),'fwhm6')
            bump6=[bump6,cell2mat(com_map.(fn{1}).(subfn{1}).fwhm6.values())];
        end
    end
end

bxx=0.125:0.25:6.25;
%     bump3hist=histcounts(bump3,bxx,'Normalization','probability');
%     bump6hist=histcounts(bump6,bxx,'Normalization','probability');
bumps.hist=histcounts([bump3,bump6],bxx,'Normalization','pdf');
bumps.mean3=mean(bump3).*1000;
bumps.sem3=std(bump3)./sqrt(numel(bump3)).*1000;
bumps.mean6=mean(bump6).*1000;
bumps.sem6=std(bump6)./sqrt(numel(bump6)).*1000;



%% burst chained loops
per_sess_coverage=struct();
for sessid=[14,18,22,33,34,68,100,102,114]
    load(fullfile("bzdata","ChainedLoop600S"+num2str(sessid)+".mat"),'covered')
    per_sess_coverage.("S"+num2str(sessid))=covered;
end
covered=cell2mat(struct2cell(per_sess_coverage));

edges = find(diff([0;covered;0]==1));
onset = edges(1:2:end-1);  % Start indices
sumdata.BSChLoop = edges(2:2:end)-onset;  % Consecutive ones counts

%% single spike chained loops
per_sess_coverage=struct();
for sessid=[14,18,22,33,34,68,100,102,114]
    load(fullfile("bzdata","SingleSpikeChainedLoop"+num2str(sessid)+".mat"),'covered')
    per_sess_coverage.("S"+num2str(sessid))=covered;
end
covered=cell2mat(struct2cell(per_sess_coverage));

edges = find(diff([0;covered;0]==1));
onset = edges(1:2:end-1);  % Start indices
sumdata.SSChLoop = (edges(2:2:end)-onset)./10;  % Consecutive ones counts


%% single spike composite loops
if false
    hrstats=load(fullfile('bzdata','hebbian_ring.mat'),'stats');
    C=struct2cell(hrstats.stats).';
    expd=[C{:}];
    expdd=[expd{:}];
    sumdata.SSCoL=expdd;
end

%% burst spike composite loops
% skipped for now
if false
    load(fullfile('bzdata','composite_burst_ring.mat'),'stats')
    C=struct2cell(stats).';
    expd=[C{:}];
    expdd=[expd{:}];
    burst_compo_pdf=histcounts(expdd,[0:20:200,300:300:1200],'Normalization','pdf');
end

%% single spike chain
load(fullfile('bzdata','chain_tag.mat'),'out');
perchaindur=struct();
[perchaindur.d6.size,perchaindur.d6.dur,perchaindur.d3.size,perchaindur.d3.dur,perchaindur.d6.int,perchaindur.d3.int]=deal([]);
for dur=reshape(fieldnames(out),1,[])
    perchaindur=struct();
    [perchaindur.size,perchaindur.dur]=deal([]);
    for wv=reshape(fieldnames(out.(dur{1})),1,[])
        for lp=reshape(fieldnames(out.(dur{1}).(wv{1})),1,[])
            perchaindur.size=[perchaindur.size;size(out.(dur{1}).(wv{1}).(lp{1}).ts,2)];
            perchaindur.dur=[perchaindur.dur;{diff(out.(dur{1}).(wv{1}).(lp{1}).ts(:,[1,end]),1,2)./30}];
        end
    end
    statss.("d"+dur)=perchaindur;
end
sumdata.SSchain=[cell2mat(statss.dd3.dur);cell2mat(statss.dd6.dur)];


%% burst spike chain
% skipped for now
if false
    load('chain_sust_tag_600.mat','out')
    perchaindur=struct();
    [perchaindur.d6.size,perchaindur.d6.dur,perchaindur.d3.size,perchaindur.d3.dur,perchaindur.d6.int,perchaindur.d3.int]=deal([]);
    for dur=reshape(fieldnames(out),1,[])
        %     [perchaindur.size,perchaindur.dur]=deal([]);
        for wv=reshape(fieldnames(out.(dur{1})),1,[])
            for lp=reshape(fieldnames(out.(dur{1}).(wv{1})),1,[])
                %             keyboard();
                perchaindur.(dur{1}).size=[perchaindur.(dur{1}).size,cellfun(@(x) size(x,1),out.(dur{1}).(wv{1}).(lp{1}).ts)];
                perchaindur.(dur{1}).dur=[perchaindur.(dur{1}).dur,cellfun(@(x) diff(x([1,end],3),1,1),out.(dur{1}).(wv{1}).(lp{1}).ts)./30];
                perchaindur.(dur{1}).int=[perchaindur.(dur{1}).int,cell2mat(cellfun(@(x) diff(x(:,3),1,1).',out.(dur{1}).(wv{1}).(lp{1}).ts,'UniformOutput',false))./30];
            end
        end
        statsm.("d"+dur)=perchaindur;
    end
    multichainhist=histcounts([statsm.dd3.d3.dur,statsm.dd6.d6.dur],[0:10:100,200:100:600],'Normalization','pdf');
end
%% PLOT

sumdata.FC(:,2)=1;
sumdata.SSchain(:,2)=2;
sumdata.ssl(:,2)=3;
sumdata.SSCoL(2,:)=4;
if size(sumdata.SSCoL,2)>2
    sumdata.SSCoL=sumdata.SSCoL.';
end
sumdata.SSChLoop(:,2)=5;
sumdata.BSChLoop(:,2)=6;

sumdata.sums=[sumdata.FC;sumdata.SSchain;sumdata.ssl;sumdata.SSChLoop;sumdata.BSChLoop];

figure()
hold on
boxplot(sumdata.sums(:,1),sumdata.sums(:,2),'Colors','k','Whisker',realmax)

set(gca(),'YScale','log','YTick',10.^[-1:4],'XTick',1:6,'XTickLabel',...
    {'FC','SSC','SSL','SSChL','BSChL'}); % boxplot of everything
yline([3000,6000],'k-.') % 3s and 6s constant

yline([bumps.mean3,bumps.mean6],'k:');
fill([xlim(),fliplr(xlim())],bumps.mean3+[-bumps.sem3,-bumps.sem3,bumps.sem3,bumps.sem3],'k','EdgeColor','none','FaceAlpha','0.2')
fill([xlim(),fliplr(xlim())],bumps.mean6+[-bumps.sem6,-bumps.sem6,bumps.sem6,bumps.sem6],'k','EdgeColor','none','FaceAlpha','0.2')

ylim([8e-2,1e4])
ylabel('Time (ms)')


% legend([fch,looph,bsh,sch,bch,snh,mnh,crh],{'FC','Single spike Loops','Burst spike loops','Single spike composite loops','Burst spike composite loops','Single spike chains','Burst spike chains','Chained loops'},'Location','eastoutside','Orientation','vertical')
% end
