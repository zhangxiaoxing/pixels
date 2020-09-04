% [pWing,to_test]=baseline_delay_WING()
% return

% type='correct_resample';
type='selec';
close all
load('reg_keep.mat');
plot_per_bin=false;
plot_entire=false;
plot_early_late=false;

if strcmp(type,'correct_resample')
    rptN=100;
    iorepeats=cell(rptN,1);
else
    rptN=1;
end
for rpt=1:rptN
    ioselstats=cell(1,6);
    for bin=1:6
        if strcmp(type,'error')
            load(sprintf('correct_error\\0820_error_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
        elseif strcmp(type,'correct_resample') && exist('rpt','var')
            load(sprintf('correct_error\\0820_correct_resample_%03d__conn_chain_duo_6s_%d_%d.mat',rpt,bin,bin+1));
        elseif strcmp(type,'baseline')
            if bin>1
                break;
            end
            load('0826_selec_conn_chain_duo_6s_-2_-1.mat');
        else
            load(sprintf('0831_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
        end
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
    if strcmp(type,'correct_resample')
        iorepeats{rpt}=ioselstats;
    end
end
% keyboard
if strcmp(type,'baseline')
    io_baseline=in_out_sel;
    save('io_sel_baseline.mat','io_baseline');
    return
end





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

return
%% per bin wing stability
load io_sel_baseline.mat
load io_sel.mat

greymatter=cellfun(@(x) ~isempty(regexp(x,'[A-Z]','match','once')), reg_set);
any(cell2mat(cellfun(@(x) x(:,1)<100,ioselstats,'UniformOutput',false)))
nansel=any(cell2mat(cellfun(@(x) x(:,1)<100,ioselstats,'UniformOutput',false)),2) | io_baseline(:,1)<100;
sel=~nansel & greymatter;
to_test=find(sel)';

base_WING=diff(io_baseline(sel,[3,7]),1,2);
delay_WING=cell2mat(cellfun(@(x) diff(x(sel,[3,7]),1,2),ioselstats,'UniformOutput',false));

corrmat=[base_WING,delay_WING];
rmat=nan(6,6);
for i=1:6
    for j=i+1:7
        [r,~]=corr(corrmat(:,i),corrmat(:,j));
        rmat(i,j-1)=r;
    end
end
figure('Color','w')
imagesc(rmat,[-0.5,0.5])
colormap('jet')
xlabel('delay bin (sec)')
set(gca,'YTickLabel',{'baseline','1 sec','2 sec','3 sec','4 sec','5 sec'})
colorbar()

%% per region wing per bin

base_WING(:,1)=diff(io_baseline(sel,[3,7]),1,2);
early_WING(:,1)=diff(io_early_delay(sel,[3,7]),1,2);
late_WING(:,1)=diff(io_late_delay(sel,[3,7]),1,2);

base_WING(:,2)=io_baseline(sel,3);
early_WING(:,2)=io_early_delay(sel,3);
late_WING(:,2)=io_late_delay(sel,3);

base_WING(:,3)=io_baseline(sel,7);
early_WING(:,3)=io_early_delay(sel,7);
late_WING(:,3)=io_late_delay(sel,7);


figure('Color','w','Position',[100,100,750,250])
for k=1:3
    corrmat=[base_WING(:,k),early_WING(:,k),late_WING(:,k)];
    rmat=nan(2,2);
    pmat=nan(2,2);
    for i=1:2
        for j=i+1:3
            [r,p]=corr(corrmat(:,i),corrmat(:,j));
            rmat(i,j-1)=r;
            pmat(i,j-1)=p;
        end
    end
    subplot(1,3,1)
    scatter(early_WING,late_WING,'MarkerFaceColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none')
    title('early vs late');
    xlim([-0.075,0.075]);
    ylim([-0.075,0.075]);
    hold on
    plot([-0.075,0.075],[-0.075,0.075],'r:')
    subplot(1,3,2)
    scatter(base_WING,early_WING,'MarkerFaceColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none')
    title('baseline vs early')
    xlim([-0.075,0.075]);
    ylim([-0.075,0.075]);
    hold on
    plot([-0.075,0.075],[-0.075,0.075],'r:')
    subplot(1,3,3)
    scatter(base_WING,late_WING,'MarkerFaceColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none')
    title('baseline vs late')
    xlim([-0.075,0.075]);
    ylim([-0.075,0.075]);
    hold on
    plot([-0.075,0.075],[-0.075,0.075],'r:')
end


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
    
    
    
   
    error_WING=[];
    correct_WING=[];
    z_WING=[];
    to_test=find(sel)';
    for i=to_test
        %%reshape repeats
        one_correct_WING=diff(correct_entire(i,[3,7]),1,2);
        one_error_WING=diff(error_entire(i,[3,7]),1,2);
        rpts=cell2mat(cellfun(@(x) sum(cell2mat(cellfun(@(y) y(i,:),x,'UniformOutput',false)')),correct_resample,'UniformOutput',false));
        wing_rpt=rpts(:,6)./rpts(:,1)-rpts(:,2)./rpts(:,1);
        %%mean, std
        mm=nanmean(wing_rpt);
        stdv=nanstd(wing_rpt);
        %%z-score
        one_z=(one_error_WING-mm)./stdv;
        
        correct_WING=[correct_WING;one_correct_WING];
        error_WING=[error_WING;one_error_WING];
        z_WING=[z_WING;one_z];
    end
    
    incsel=(z_WING>1.96);
    decsel=(z_WING<-1.96);
    
    disp(reg_set(to_test(find(incsel))));
    disp(reg_set(to_test(find(decsel))));
    
    insig=~(incsel|decsel);
    
    close all
    fh=figure('Color','w','Position',[100,100,130,240])
    hold on
    randx=rand(size(error_WING))*0.1;
    plot(1+randx,correct_WING,'ko')
    plot(2+randx,error_WING,'ko','MarkerFaceColor','k')
    plot(repmat([1;2],1,nnz(insig))+randx(insig)',[correct_WING(insig)';error_WING(insig)'],'k-');
    hi=plot(repmat([1;2],1,nnz(incsel))+randx(incsel)',[correct_WING(incsel)';error_WING(incsel)'],'r-','LineWidth',2);
    hd=plot(repmat([1;2],1,nnz(decsel))+randx(decsel)',[correct_WING(decsel)';error_WING(decsel)'],'b-','LineWidth',2);
    xlim([0.5,2.5])
    ylim([-0.06,0.06])
    set(gca,'YTick',-0.05:0.05:0.05)
    ylabel('WING')
    legend([hi(1),hd(1)],{'increased','decreased'},'Location','northoutside');
    set(gca,'XTick',1:2,'XTickLabel',{'correct','error'},'XTickLabelRotation',90)
    exportgraphics(fh,'WING_correct_error.pdf')
%     
%     fh=figure('Color','w','Position',[100,100,230,230])
%     hold on
%     scatter(correct_WING,error_WING,10,'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.5)
%     plot([-0.06,0.06],[-0.06,0.06],'k:');
%     xlabel('correct trials')
%     ylabel('error trials')
%     xlim([-0.06,0.06])
%     ylim([-0.06,0.06])
%     set(gca,'XTick',-0.05:0.05:0.05,'YTick',-0.05:0.05:0.05)
%     [r,p]=corr(correct_WING,error_WING);
%     text(-0.05,0.05,sprintf('r=%.3f,p=%.3f',r,p))
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
    for repeat=1:1000
        correct_resample=correct_resample(randperm(length(correct_resample)));
        correct_entire=zeros(140,12);
        for rpt=1
%             fprintf('a %d',rpt);
            for bin=1:6
                correct_entire=correct_entire+correct_resample{rpt}{bin};
            end
        end
        correct_entire(:,[3 5 7 9])=correct_entire(:,[2 4 6 8])./correct_entire(:,1);
        correct_entire(:,12)=correct_entire(:,11)./correct_entire(:,10);
        
        corr_A=correct_entire;
        
        correct_entire=zeros(140,12);
        for rpt=2%(1:length(correct_resample)/2)+length(correct_resample)/2
%             fprintf('b %d',rpt);
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
    
    fh=figure('Color','w','Position',[100,100,80,240])
    hold on
    mm=mean(rlist);
    ci=bootci(1000,@(x) mean(x), rlist);
    bh=bar(1,mm(1),'FaceColor','k','EdgeColor','k','LineWidth',1);
    bk=bar(2,mm(2),'FaceColor','w','EdgeColor','k','LineWidth',1);
    errorbar(1:2,mm,ci(1,:)-mm,ci(2,:)-mm,'k.','Color',[0.25,0.25,0.25])
    xlim([0.5,2.5])
    ylim([-0.05,0.25])
    set(gca,'YTick',0:0.1:0.2,'XTick',1:2,'XTickLabel',{'corr.-corr.','error-corr.'},'XTickLabelRotation',90)
    ylabel('correlation coefficient')
    exportgraphics(fh,'WING_correct_error_corr.pdf')
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

function p=WING_permutation_test(statsA,statsB,rpt)
     poolLen=statsA(1)+statsB(1);
     inCount=statsA(2)+statsB(2);
     outCount=statsA(6)+statsB(6);
     permDiff=nan(rpt,1);
     for i=1:rpt
        fullPool=randperm(poolLen);
        inA=nnz(fullPool(1:statsA(1))<=inCount);
        inB=inCount-inA;
        
        outA=nnz(fullPool(1:statsA(1))>inCount & fullPool(1:statsA(1))<=inCount+outCount);
        outB=outCount-outA;
        
        permWingA=(outA-inA)/statsA(1);
        permWingB=(outB-inB)/statsB(1);
        perfDiff(i)=permWingB-permWingA;
     end
     actualDiff=(statsB(6)-statsB(2))/statsB(1)-(statsA(6)-statsA(2))/statsA(1);
     z=-abs(actualDiff-mean(perfDiff))/std(perfDiff);
     p=2*normcdf(z);

end


%% baseline delay
function [p_WING,to_test]=baseline_delay_WING()
    load io_sel_baseline.mat
    load io_sel.mat
    
    greymatter=cellfun(@(x) ~isempty(regexp(x,'[A-Z]','match','once')), reg_set);
    nansel=io_entire_delay(:,1)<100 | io_baseline(:,1)<100;
    sel=~nansel & greymatter;
    to_test=find(sel)';
    
    
    base_WING=diff(io_baseline(sel,[3,7]),1,2);
    delay_WING=diff(io_entire_delay(sel,[3,7]),1,2);
    p_WING=[];

    for i=to_test
%         catBase=zeros(1,io_baseline(i,1));
%         catBase(1:io_baseline(i,2))=1;
%         catBase(io_baseline(i,2)+1:sum(io_baseline(i,[2,6])))=2;
%         
%         catDelay=zeros(1,io_entire_delay(i,1));
%         catDelay(1:io_entire_delay(i,2))=1;
%         catDelay(io_entire_delay(i,2)+1:sum(io_entire_delay(i,[2,6])))=2;
%         
%         [~,~,p]=crosstab(1:io_baseline(i,1)+io_entire_delay(i,1)>io_baseline(i,1),...
%             [catBase,catDelay]);
        p_WING(end+1)=WING_permutation_test(io_baseline(i,:),io_entire_delay(i,:),1000);
    end
    
    incsel=(p_WING'<0.05 & delay_WING>base_WING);
    decsel=(p_WING'<0.05 & delay_WING<base_WING);
    
    disp(reg_set(to_test(find(incsel))));
    disp(reg_set(to_test(find(decsel))));
    
    insig=~(incsel|decsel);
    
    
    fh=figure('Color','w','Position',[100,100,130,240])
    hold on
    randx=rand(size(delay_WING))*0.1;
    plot(1+randx,base_WING,'ko')
    plot(2+randx,delay_WING,'ko','MarkerFaceColor','k')
    plot(repmat([1;2],1,nnz(insig))+randx(insig)',[base_WING(insig)';delay_WING(insig)'],'k-');
    hi=plot(repmat([1;2],1,nnz(incsel))+randx(incsel)',[base_WING(incsel)';delay_WING(incsel)'],'r-','LineWidth',2);
    hd=plot(repmat([1;2],1,nnz(decsel))+randx(decsel)',[base_WING(decsel)';delay_WING(decsel)'],'b-','LineWidth',2);
    xlim([0.5,2.5])
    ylim([-0.08,0.08])
    set(gca,'YTick',-0.05:0.05:0.05)
    ylabel('WING')
    legend([hi(1),hd(1)],{'increased','decreased'},'Location','northoutside');
    set(gca,'XTick',1:2,'XTickLabel',{'baseline','delay'},'XTickLabelRotation',90)
    exportgraphics(fh,'WING_base_delay.pdf')
    
    
    
    if false % scatter plot
        [r,p]=corr(base_WING,delay_WING);
        figure('Color','w')
        hold on
        sh=scatter(base_WING,delay_WING)
        xlim([-0.08,0.08])
        ylim([-0.08,0.08])
        plot([-0.08,0.08],[-0.08,0.08],'k:')
        xlabel('baseline WING')
        ylabel('delay WING')
        legend(sh,sprintf('r = %.3f, p = %.3f',r,p));
    end
   
end