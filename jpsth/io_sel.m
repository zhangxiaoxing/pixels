close all
load('reg_keep.mat');
plot_per_bin=false;
plot_entire=true;
plot_early_late=false;
iorepeats=cell(60,1);
% for rpt=1:60
    ioselstats=cell(1,6);
    for bin=1:6
    load(sprintf('correct_error\\0820_error_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
%     load(sprintf('correct_error\\0820_correct_resample_%03d__conn_chain_duo_6s_%d_%d.mat',rpt,bin,bin+1));
    in_out_sel=nan(length(reg_set),12);
    for reg_idx=1:length(reg_set)
        pair_fw= nnz((pair_reg(:,1)~=reg_idx) & (pair_reg(:,2)==reg_idx));
        pair_rev= nnz((pair_reg(:,1)==reg_idx) & (pair_reg(:,2)~=reg_idx));
        pair_count=pair_fw+pair_rev;

        in_conn_S1= nnz((reg_chain_S1(:,1)~=reg_idx) & (reg_chain_S1(:,2)==reg_idx));
        in_sel_S1=nnz((reg_chain_S1(:,1)~=reg_idx) & (reg_chain_S1(:,2)==reg_idx) & pref_chain_S1(:,bin)>0);

        out_conn_S1= nnz((reg_chain_S1(:,2)~=reg_idx) & (reg_chain_S1(:,1)==reg_idx));
        out_sel_S1=nnz((reg_chain_S1(:,2)~=reg_idx) & (reg_chain_S1(:,1)==reg_idx) & pref_chain_S1(:,bin+6)>0);

        auto_pair=nnz((pair_reg(:,1)==reg_idx) & (pair_reg(:,2)==reg_idx))*2;
        auto_conn_S1= nnz((reg_chain_S1(:,1)==reg_idx) & (reg_chain_S1(:,2)==reg_idx));

        in_out_sel(reg_idx,:)=[pair_count,in_conn_S1,in_conn_S1/pair_count, ...%1 2 3
            in_sel_S1,in_sel_S1/pair_count,...% 4 5
            out_conn_S1,out_conn_S1/pair_count,...% 6 7
            out_sel_S1,out_sel_S1/pair_count,...% 8 9
            auto_pair,auto_conn_S1,auto_conn_S1/auto_pair]; % 10 11 12
    end

    ioselstats{bin}=in_out_sel;
    if plot_per_bin
        plotOne(in_out_sel,reg_set,sprintf('bin %d',bin),sprintf('0714_io_selec_bin%d.png',bin));
    end
    end
    iorepeats{rpt}=ioselstats;
% end
keyboard
%% entire delay

io_entire_delay=ioselstats{1};
for bin=2:6
    io_entire_delay=io_entire_delay+ioselstats{bin};
end
io_entire_delay(:,[3 5 7 9])=io_entire_delay(:,[2 4 6 8])./io_entire_delay(:,1);
io_entire_delay(:,12)=io_entire_delay(:,11)./io_entire_delay(:,10);
if plot_entire
%     plotOne(io_entire_delay,reg_set,'sum of all bins','0714_io_congruent_bin_sum.png',[4,5,8,9]);
% plotOne(io_entire_delay,reg_set,'sum of all bins','0714_io_congruent_bin_sum.png',[2,3,11,12]);    
plotOne(io_entire_delay,reg_set,'sum of all bins','810_io_selective_in_out_bin_sum_.png',[2,3,6,7]);
end

%% early delay late delay
io_early_delay=ioselstats{1};
for bin=2:3
    io_early_delay=io_early_delay+ioselstats{bin};
end
io_early_delay(:,[3 5 7 9])=io_early_delay(:,[2 4 6 8])./io_early_delay(:,1);
io_early_delay(:,12)=io_early_delay(:,11)./io_early_delay(:,10);

io_late_delay=ioselstats{4};
for bin=5:6
    io_late_delay=io_late_delay+ioselstats{bin};
end
io_late_delay(:,[3 5 7 9])=io_late_delay(:,[2 4 6 8])./io_late_delay(:,1);
io_late_delay(:,12)=io_late_delay(:,11)./io_late_delay(:,10);

if plot_early_late
%     plotOne(io_early_delay,reg_set,'early delay','0714_io_selec_early_delay.png',[4,5,8,9]);
%     plotOne(io_late_delay,reg_set,'late delay','0810_io_selec_late_delay.png',[2 3 6 7]);
end

save('io_sel.mat','ioselstats','io_entire_delay','io_early_delay','io_late_delay','reg_set');


%% correct error
if false
    load io_sel_correct_error.mat
    error_entire=errorio{1};
    for bin=2:6
        error_entire=error_entire+errorio{bin};
    end
    error_entire(:,[3 5 7 9])=error_entire(:,[2 4 6 8])./error_entire(:,1);
    error_entire(:,12)=error_entire(:,11)./error_entire(:,10);
    
    correct_entire=zeros(140,12);
    for rpt=1:length(correct_resample)
        for bin=1:6
            correct_entire=correct_entire+correct_resample{rpt}{bin};
        end
    end
    correct_entire(:,[3 5 7 9])=correct_entire(:,[2 4 6 8])./correct_entire(:,1);
    correct_entire(:,12)=correct_entire(:,11)./correct_entire(:,10);
    
    load('reg_keep.mat');
    greymatter=cellfun(@(x) ~isempty(regexp(x,'[A-Z]','match','once')), reg_set);
    nansel=correct_entire(:,1)<100 | error_entire(:,1)<100
    sel=~nansel & greymatter;
        
    correct_WING=diff(correct_entire(sel,[3,7]),1,2)
    error_WING=diff(error_entire(sel,[3,7]),1,2)
    
    
    fh=figure('Color','w','Position',[100,100,230,230])
    hold on
    scatter(correct_WING,error_WING,10,'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.5)
    plot([-0.06,0.06],[-0.06,0.06],'k:');
    xlabel('correct trials')
    ylabel('error trials')
    xlim([-0.06,0.06])
    ylim([-0.06,0.06])
    set(gca,'XTick',-0.05:0.05:0.05,'YTick',-0.05:0.05:0.05)
    [r,p]=corr(correct_WING,error_WING);
    text(-0.05,0.05,sprintf('r=%.3f,p=%.3f',r,p))
    exportgraphics(fh,'WING_correct_error.pdf')
end

%% correct-correct resample
if false
    load('reg_keep.mat');
    greymatter=cellfun(@(x) ~isempty(regexp(x,'[A-Z]','match','once')), reg_set);
    
    load io_sel_correct_error.mat
    error_entire=errorio{1};
    for bin=2:6
        error_entire=error_entire+errorio{bin};
    end
    error_entire(:,[3 5 7 9])=error_entire(:,[2 4 6 8])./error_entire(:,1);
    error_entire(:,12)=error_entire(:,11)./error_entire(:,10);
    
    rlist=[];
    for repeat=1:100
        correct_resample=correct_resample(randperm(length(correct_resample)));
        correct_entire=zeros(140,12);
        for rpt=1
            fprintf('a %d',rpt);
            for bin=1:6
                correct_entire=correct_entire+correct_resample{rpt}{bin};
            end
        end
        correct_entire(:,[3 5 7 9])=correct_entire(:,[2 4 6 8])./correct_entire(:,1);
        correct_entire(:,12)=correct_entire(:,11)./correct_entire(:,10);
        
        corr_A=correct_entire;
        
        correct_entire=zeros(140,12);
        for rpt=2%(1:length(correct_resample)/2)+length(correct_resample)/2
            fprintf('b %d',rpt);
            for bin=1:6
                correct_entire=correct_entire+correct_resample{rpt}{bin};
            end
        end
        correct_entire(:,[3 5 7 9])=correct_entire(:,[2 4 6 8])./correct_entire(:,1);
        correct_entire(:,12)=correct_entire(:,11)./correct_entire(:,10);
        
        corr_B=correct_entire;
        
        nansel=corr_A(:,1)<100 | corr_B(:,1)<100 | error_entire(:,1)<100;
        sel=~nansel & greymatter;
        
        wingA=diff(corr_A(sel,[3,7]),1,2);
        wingB=diff(corr_B(sel,[3,7]),1,2);
        wingE=diff(error_entire(sel,[3,7]),1,2);
        [r,~]=corr(wingA,wingB);
        [re,~]=corr(wingE,wingB);
%         keyboard
        rlist(repeat,:)=[r,re];
    end
    
end

function plotOne(in_out_sel,reg_set,str_title,fname,v_idx)
    if ~exist('v_idx','var') 
        v_idx=[4,5,8,9];
    end
    t=num2cell(v_idx);
    [inCount,inFrac,outCount,outFrac]=t{:};

    greymatter=cellfun(@(x) ~isempty(regexp(x,'[A-Z]','match','once')), reg_set);
    countsel=in_out_sel(:,1)>=100 & greymatter; %126 115 104
    iosel=in_out_sel(countsel,:);
    sink_sel=false(size(iosel,1),1);
    source_sel=false(size(iosel,1),1);
    for i=1:size(iosel,1)
        currcount=iosel(i,1);
        [~,~,p]=crosstab(1:2*currcount>currcount,[1:currcount>iosel(i,inCount),1:currcount>iosel(i,outCount)]);
    %     chisqsel(i)=p*numel(chisqsel)<0.05; %bonferroni adjust
        if p<0.05 && iosel(i,outCount)>iosel(i,inCount)
            source_sel(i)=true;
        elseif p<0.05 && iosel(i,outCount)<iosel(i,inCount)
            sink_sel(i)=true;
        end
    end

    regsel=reg_set(countsel);
%     regid=zeros(size(regsel));
%     for i=1:length(regsel)
%         regid(i)=regclass.id(strcmp(regclass.reg,regsel{i})); %CTX=0,TH=1,STR=2,ELSE=3
%     end
%     
    fh=figure('Color','w','Position',[100,100,280,280]);
    hold on
    otherh=scatter(iosel(~(sink_sel|source_sel),inFrac),iosel(~(sink_sel|source_sel),outFrac),'o','MarkerEdgeColor','none','MarkerFaceColor',[0.5,0.5,0.5],'MarkerFaceAlpha',0.5);
    gainh=scatter(iosel(source_sel,inFrac),iosel(source_sel,outFrac),'o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha',0.5);
    damph=scatter(iosel(sink_sel,inFrac),iosel(sink_sel,outFrac),'o','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha',0.5);

% disp('source')
    % disp(regsel(source_sel));
    % disp('sink')
    % disp(regsel(sink_sel));
%     for i=reshape(find(source_sel | sink_sel),1,[])
%         text(iosel(i,inFrac),iosel(i,outFrac),regsel(i),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8)
%     end
    
    plot([0,1],[0,1],'k:');
    xlim([0,0.2]);
    ylim([0,0.2]);
    xlabel('in-density');
    ylabel('out-density');
    title(str_title);
    legend([gainh,damph,otherh],{'gain','damping','others'})
    
    [r,p]=corr(iosel(:,inFrac),iosel(:,outFrac));
    
    keyboard
%     print(fh,fname,'-dpng','-r300');
    exportgraphics(fh,replace(fname,'.png','.pdf'),'ContentType','vector');
    
end



function localProcessing(io_entire_delay,reg_set)

greymatter=cellfun(@(x) ~isempty(regexp(x,'[A-Z]','match','once')), reg_set);
nansel=io_entire_delay(:,1)<100 | io_entire_delay(:,10)<100;
sel=~nansel & greymatter;
selexp=io_entire_delay(:,1)>=100 & greymatter;


figure('Color','w','Position',[50,50,1680,800])
subplot(2,6,1);
scatter(io_entire_delay(selexp,3),io_entire_delay(selexp,7))
[r,p]=corr(io_entire_delay(selexp,3),io_entire_delay(selexp,7));
xlabel(sprintf('r=%.3f,p=%.3f',r,p));
title('input v.s. output')
ylabel('selective-connection')

subplot(2,6,2);
fh=figure('Color','w','Position',[100,100,280,280]);
scatter(io_entire_delay(sel,3),io_entire_delay(sel,12))
[r,p]=corr(io_entire_delay(sel,3),io_entire_delay(sel,12));
xlabel(sprintf('r=%.3f,p=%.3f',r,p));
title('input v.s. local')
%% %%%%%%%%%%%%%%%%LOCAL V.S IN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if false
    in_out_sel=io_entire_delay;
    t=num2cell([2,3,6,12]);
    [inCount,inFrac,outCount,outFrac]=t{:};
    
    greymatter=cellfun(@(x) ~isempty(regexp(x,'[A-Z]','match','once')), reg_set);
    countsel=in_out_sel(:,1)>=100 & in_out_sel(:,10)>=100 & greymatter; %126 92 115 77
    iosel=in_out_sel(countsel,:);
    sink_sel=false(size(iosel,1),1);
    source_sel=false(size(iosel,1),1);
    for i=1:size(iosel,1)
        currcount=iosel(i,1);
        [~,~,p]=crosstab(1:2*currcount>currcount,[1:currcount>iosel(i,inCount),1:currcount>iosel(i,outCount)]);
        if p<0.05 && iosel(i,outCount)>iosel(i,inCount)
            source_sel(i)=true;
        elseif p<0.05 && iosel(i,outCount)<iosel(i,inCount)
            sink_sel(i)=true;
        end
    end

    regsel=reg_set(countsel);
    fh=figure('Color','w','Position',[100,100,280,280]);
    hold on
    otherh=scatter(iosel(~(sink_sel|source_sel),inFrac),iosel(~(sink_sel|source_sel),outFrac),'o','MarkerEdgeColor','none','MarkerFaceColor',[0.5,0.5,0.5],'MarkerFaceAlpha',0.5);
    gainh=scatter(iosel(source_sel,inFrac),iosel(source_sel,outFrac),'o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha',0.5);
    damph=scatter(iosel(sink_sel,inFrac),iosel(sink_sel,outFrac),'o','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha',0.5);

    [r,p]=corr(iosel(:,inFrac),iosel(:,outFrac))
    
    plot([0,1],[0,1],'k:');
    xlim([0,0.22]);
    ylim([0,0.22]);
    set(gca,'XTick',0:0.1:0.2,'YTick',0:0.1:0.2)
    xlabel('in-density');
    ylabel('local-density');
    title('input vs local');
    legend([gainh,damph,otherh],{'gain','damping','others'})
    keyboard
    exportgraphics(gcf,'inputdden_localden_corr.pdf')
end

%% %%%%%%%%%%%%%%%%LOCAL V.S OUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,6,3);
scatter(io_entire_delay(sel,7),io_entire_delay(sel,12))
[r,p]=corr(io_entire_delay(sel,7),io_entire_delay(sel,12));
title('output v.s. local')
xlabel(sprintf('r=%.3f,p=%.3f',r,p));

if false
    in_out_sel=io_entire_delay;
    t=num2cell([2,7,6,12]);
    [inCount,inFrac,outCount,outFrac]=t{:};
    
    greymatter=cellfun(@(x) ~isempty(regexp(x,'[A-Z]','match','once')), reg_set);
    countsel=in_out_sel(:,1)>=100 & in_out_sel(:,10)>=100 & greymatter; %121 115 100
    iosel=in_out_sel(countsel,:);
    sink_sel=false(size(iosel,1),1);
    source_sel=false(size(iosel,1),1);
    for i=1:size(iosel,1)
        currcount=iosel(i,1);
        [~,~,p]=crosstab(1:2*currcount>currcount,[1:currcount>iosel(i,inCount),1:currcount>iosel(i,outCount)]);
        if p<0.05 && iosel(i,outCount)>iosel(i,inCount)
            source_sel(i)=true;
        elseif p<0.05 && iosel(i,outCount)<iosel(i,inCount)
            sink_sel(i)=true;
        end
    end

    regsel=reg_set(countsel);
    fh=figure('Color','w','Position',[100,100,280,280]);
    hold on
    otherh=scatter(iosel(~(sink_sel|source_sel),inFrac),iosel(~(sink_sel|source_sel),outFrac),'o','MarkerEdgeColor','none','MarkerFaceColor',[0.5,0.5,0.5],'MarkerFaceAlpha',0.5);
    gainh=scatter(iosel(source_sel,inFrac),iosel(source_sel,outFrac),'o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha',0.5);
    damph=scatter(iosel(sink_sel,inFrac),iosel(sink_sel,outFrac),'o','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha',0.5);

    [r,p]=corr(iosel(:,inFrac),iosel(:,outFrac))
    
    plot([0,1],[0,1],'k:');
    xlim([0,0.22]);
    ylim([0,0.22]);
    set(gca,'XTick',0:0.1:0.2,'YTick',0:0.1:0.2)
    xlabel('out-density');
    ylabel('local-density');
    title('output vs local');
    legend([gainh,damph,otherh],{'gain','damping','others'})
    keyboard
    exportgraphics(gcf,'outputden_localden_corr.pdf')
end








subplot(2,6,4);
scatter(diff(io_entire_delay(sel,[3,7]),1,2)./sum(io_entire_delay(sel,[3,7]),2),io_entire_delay(sel,12))
[r,p]=corr(diff(io_entire_delay(sel,[3,7]),1,2)./sum(io_entire_delay(sel,[3,7]),2),io_entire_delay(sel,12));
xlabel(sprintf('r=%.3f,p=%.3f',r,p));
title('gain v.s. local')

%% GLM in,local->out %%
subplot(2,6,5)
mdl=fitglm([io_entire_delay(sel,3),io_entire_delay(sel,12)],io_entire_delay(sel,7),'linear');
scatter(io_entire_delay(sel,7),mdl.Fitted.Response)
xlabel('selective-output density')
ylabel('GLM response')
legend({sprintf('input coef=%.3f, local coef=%.3f',mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(3))});
title('GLM-out');

if false
    in_out_sel=io_entire_delay;
    t=num2cell([2,7,6,12]);
    [inCount,inFrac,outCount,outFrac]=t{:};
    
    greymatter=cellfun(@(x) ~isempty(regexp(x,'[A-Z]','match','once')), reg_set);
    countsel=in_out_sel(:,1)>=100 & in_out_sel(:,10)>=100 & greymatter; %121 115 100
    iosel=in_out_sel(countsel,:);
    sink_sel=false(size(iosel,1),1);
    source_sel=false(size(iosel,1),1);
    for i=1:size(iosel,1)
        currcount=iosel(i,1);
        [~,~,p]=crosstab(1:2*currcount>currcount,[1:currcount>iosel(i,inCount),1:currcount>iosel(i,outCount)]);
        if p<0.05 && iosel(i,outCount)>iosel(i,inCount)
            source_sel(i)=true;
        elseif p<0.05 && iosel(i,outCount)<iosel(i,inCount)
            sink_sel(i)=true;
        end
    end

    regsel=reg_set(countsel);
    
    mdl=fitglm(io_entire_delay(countsel,[3 12]),io_entire_delay(countsel,7),'linear');
    scatter(io_entire_delay(countsel,7),mdl.Fitted.Response)
    
    fh=figure('Color','w','Position',[100,100,280,280]);
    hold on
    otherh=scatter(iosel(~(sink_sel|source_sel),7),mdl.Fitted.Response(~(sink_sel|source_sel)),'o','MarkerEdgeColor','none','MarkerFaceColor',[0.5,0.5,0.5],'MarkerFaceAlpha',0.5);
    gainh=scatter(iosel(source_sel,7),mdl.Fitted.Response(source_sel),'o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha',0.5);
    damph=scatter(iosel(sink_sel,7),mdl.Fitted.Response(sink_sel),'o','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha',0.5);
   
    plot([0,1],[0,1],'k:');
    xlim([0,0.2]);
    ylim([0,0.2]);
    set(gca,'XTick',0:0.1:0.2,'YTick',0:0.1:0.2)
    xlabel('out-density');
    ylabel('fitted response');
    title('output vs local');
    legend([gainh,damph,otherh],{'gain','damping','others'})
    keyboard
    exportgraphics(gcf,'outputden_fitted.pdf')
end
    
subplot(2,6,6);
mdl=fitglm([io_entire_delay(sel,3),io_entire_delay(sel,12)],diff(io_entire_delay(sel,[3,7]),1,2)./sum(io_entire_delay(sel,[3,7]),2),'linear');
scatter(diff(io_entire_delay(sel,[3,7]),1,2)./sum(io_entire_delay(sel,[3,7]),2),mdl.Fitted.Response)
xlabel('selective-output density')
ylabel('GLM response')
legend({sprintf('input coef=%.3f, local coef=%.3f',mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(3))});
title('GLM-gain');



subplot(2,6,7);
scatter(io_entire_delay(sel,5),io_entire_delay(sel,9))
[r,p]=corr(io_entire_delay(sel,5),io_entire_delay(sel,9));
xlabel(sprintf('r=%.3f,p=%.3f',r,p));
title('input v.s. output')
ylabel('congruent coding connection')

subplot(2,6,8);
scatter(io_entire_delay(:,5),io_entire_delay(:,12))
[r,p]=corr(io_entire_delay(sel,5),io_entire_delay(sel,12));
xlabel(sprintf('r=%.3f,p=%.3f',r,p));
title('input v.s. local')

subplot(2,6,9);
scatter(io_entire_delay(:,9),io_entire_delay(:,12))
[r,p]=corr(io_entire_delay(sel,9),io_entire_delay(sel,12));
title('output v.s. local')
xlabel(sprintf('r=%.3f,p=%.3f',r,p));

subplot(2,6,10);
scatter(diff(io_entire_delay(sel,[5,9]),1,2)./sum(io_entire_delay(sel,[5,9]),2),io_entire_delay(sel,12))
[r,p]=corr(diff(io_entire_delay(sel,[5,9]),1,2)./sum(io_entire_delay(sel,[5,9]),2),io_entire_delay(sel,12));
title('gain v.s. local')
xlabel(sprintf('r=%.3f,p=%.3f',r,p));


subplot(2,6,11);
mdl=fitglm([io_entire_delay(sel,5),io_entire_delay(sel,12)],io_entire_delay(sel,9),'linear');
scatter(io_entire_delay(sel,9),mdl.Fitted.Response)
xlabel('selective-output density')
ylabel('GLM response')
legend({sprintf('input coef=%.3f, local coef=%.3f',mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(3))});
title('GLM-out density');


subplot(2,6,12);
mdl=fitglm([io_entire_delay(sel,5),io_entire_delay(sel,12)],diff(io_entire_delay(sel,[5,9]),1,2)./sum(io_entire_delay(sel,[5,9]),2),'linear');
scatter(diff(io_entire_delay(sel,[5,9]),1,2)./sum(io_entire_delay(sel,[5,9]),2),mdl.Fitted.Response)
xlabel('selective-output density')
ylabel('GLM response')
legend({sprintf('input coef=%.3f, local coef=%.3f',mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(3))});
title('GLM-gain');


exportgraphics(gcf,'i_o_local_corr.pdf','ContentType','vector');
end

function localProcessOne(io_entire_delay,id1,id2)
fh=figure();
scatter(io_entire_delay(:,id1),io_entire_delay(:,id2))
nansel=isnan(io_entire_delay(:,1)) | isnan(io_entire_delay(:,10));
% corr(io_entire_delay(:,id1)
pause()
close(fh)
end


