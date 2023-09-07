% function  rings_wave_dynamic()
function rings_wave_dynamic_sums_bar()
sumdata=struct();
%% individual loops ================================================
% with preselection
sumdata.ssl=wave.loop.ssloop_time_const(skip_save=true,skip_plot=true);

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
if false
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
end

%% single spike composite loops
nestedf=load(fullfile('binary','delay_iti_runlength_covered.mat'),'run_length');
sumdata.SSChLoop=nestedf.run_length.delay(:,2);

%% single spike chain
sumdata.SSchain=wave.chain.sschain_time_const(skip_save=true,skip_plot=true);


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
sumdata.SSChLoop(:,2)=4;

sumdata.sums=[sumdata.FC;sumdata.SSchain;sumdata.ssl;sumdata.SSChLoop];

fh=figure();
hold on
boxplot(sumdata.sums(:,1)./1000,sumdata.sums(:,2),'Colors','k','Whisker',inf)

set(gca(),'YScale','log','YTick',10.^[-4:1],'XTick',1:4,'XTickLabel',...
    {'SC','Chain','Sloops','Nloops'}); % boxplot of everything
yline([3,6],'k-.') % 3s and 6s constant
if false
    yline([bumps.mean3,bumps.mean6],'k:');
    fill([xlim(),fliplr(xlim())],bumps.mean3+[-bumps.sem3,-bumps.sem3,bumps.sem3,bumps.sem3],'k','EdgeColor','none','FaceAlpha','0.2')
    fill([xlim(),fliplr(xlim())],bumps.mean6+[-bumps.sem6,-bumps.sem6,bumps.sem6,bumps.sem6],'k','EdgeColor','none','FaceAlpha','0.2')
end
ylim([1e-4,10])
ylabel('Time (s)')
savefig(fullfile('binary','hierarchy_time_const_sums.fig'))

% legend([fch,looph,bsh,sch,bch,snh,mnh,crh],{'FC','Single spike Loops','Burst spike loops','Single spike composite loops','Burst spike composite loops','Single spike chains','Burst spike chains','Chained loops'},'Location','eastoutside','Orientation','vertical')
% end
