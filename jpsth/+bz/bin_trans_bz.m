
%% one side
close all
% corr_coactive_wing()
% return
all_reg=true;
if all_reg
    [local_hat,local_ci]=coactive_all(true);
    [inter_hat,inter_ci]=coactive_all(false);
    plot_coactivate_all(local_hat,local_ci,inter_hat,inter_ci);
    return
end

if true
    per_regl=coactive_per_reg('local');
    fhl=plot_coactive_per_reg(per_regl,'local');
    
    per_regi=coactive_per_reg('inter_input');
    fhi=plot_coactive_per_reg(per_regi,'inter region input');
    
    per_rego=coactive_per_reg('inter_output');
    fho=plot_coactive_per_reg(per_rego,'inter region output');
    
    keyboard
    save('bin_trans_per_region.mat','per_regi','per_rego','per_regl')
    
    exportgraphics(fhl,'bin_trans_per_region_local.pdf')
    exportgraphics(fhi,'bin_trans_per_region_input.pdf')
    exportgraphics(fho,'bin_trans_per_region_output.pdf')
end



function corr_coactive_wing()
load('io_sel.mat','ioselstats','io_entire_delay');
load('bin_trans_per_region.mat','per_regi','per_rego','per_regl');
load('reg_keep.mat');
greymatter=find(cellfun(@(x) ~isempty(regexp(x,'[A-Z]','match','once')), reg_set));

%             in_out_sel(reg_idx,:)=[pair_count,in_conn_S1,in_conn_S1/pair_count, ...%1 2 3
%                 in_sel_S1,in_sel_S1/pair_count,...% 4 5
%                 out_conn_S1,out_conn_S1/pair_count,...% 6 7
%                 out_sel_S1,out_sel_S1/pair_count,...% 8 9
%                 auto_pair,auto_conn_S1,auto_conn_S1/auto_pair]; % 10 11 12


%     per_reg(onereg,:)=[nnz(curr_code),numel(curr_code),... %1 2
%         nnz(curr_none),numel(curr_none),...% 3 4
%         nnz(prev_code),numel(prev_code),...% 5 6
%         nnz(prev_none),numel(prev_none),...% 7 8
%         nnz(prev2_code),numel(prev2_code),...% 9 10
%         nnz(prev2_none),numel(prev2_none)];% 11 12
per_reg_list={per_regl,per_regi,per_rego};
titles={'local','input','output'};
for type=1:3
    fh=figure('Color','w','Position',[100,100,750,250]);
    per_reg=per_reg_list{type};
    for bin=1:3
        corrpair=[];
        for onereg=greymatter'
            if io_entire_delay(onereg,1)>=100 && min(per_reg(onereg,[6 8]))>=100
                corrpair(end+1,:)=[diff(io_entire_delay(onereg,[3 7]),1,2),...
                  (per_reg(onereg,bin*4-3)/per_reg(onereg,bin*4-2))/(per_reg(onereg,bin*4-1)/per_reg(onereg,bin*4))];
            end
        end
        
        [r,p]=corr(corrpair(:,1),corrpair(:,2));
        subplot(1,3,bin);
        scatter(corrpair(:,1),corrpair(:,2),12,'Marker','o','MarkerFaceAlpha',0.5,'MarkerFaceColor','k','MarkerEdgeColor','none');
        legend(sprintf('r=%.3f, p=%.3f',r,p));
        xlabel('WING')
        ylabel(sprintf('select. enhance t+%d',bin-1));
    end
    sgtitle(titles{type});
    exportgraphics(fh,sprintf('WING_corr_%s_enhance.pdf',titles{type}))
end
end



function plot_coactivate_all(local_hat,local_ci,inter_hat,inter_ci)
fh=figure('Color','w','Position',[100,100,550,235]);
subplot(1,2,1)
hold on
hln=bar((1:3)-0.2,local_hat([1 3 5]),0.4,'FaceAlpha',1,'FaceColor','k');
hlc=bar((1:3)+0.2,local_hat([2 4 6]),0.4,'FaceAlpha',0.5,'FaceColor','w');
errorbar((1:3)-0.2,local_hat([1 3 5]),local_ci([1 3 5],1)'-local_hat([1 3 5]),local_ci([1 3 5],2)'-local_hat([1 3 5]),'.','Color',[0.2,0.2,0.2])
errorbar((1:3)+0.2,local_hat([2 4 6]),local_ci([2 4 6],1)'-local_hat([2 4 6]),local_ci([2 4 6],2)'-local_hat([2 4 6]),'k.')
xlim([0.4,3.6])
ylim([0,0.5])
ylabel('fraction of selective post-unit')
set(gca,'XTick',1:3,'XTickLabel',{'t','t+1sec','t+2sec'},'YTick',0:0.2:0.4)
title('Local')

subplot(1,2,2)
hold on;
hin=bar((1:3)-0.2,inter_hat([1 3 5]),0.4,'FaceAlpha',1,'FaceColor','k');
hic=bar((1:3)+0.2,inter_hat([2 4 6]),0.4,'FaceAlpha',0.5,'FaceColor','w');
errorbar((1:3)-0.2,inter_hat([1 3 5]),inter_ci([1 3 5],1)'-inter_hat([1 3 5]),inter_ci([1 3 5],2)'-inter_hat([1 3 5]),'.','Color',[0.2,0.2,0.2])
errorbar((1:3)+0.2,inter_hat([2 4 6]),inter_ci([2 4 6],1)'-inter_hat([2 4 6]),inter_ci([2 4 6],2)'-inter_hat([2 4 6]),'k.')
xlim([0.4,3.6])
ylim([0,0.5])
ylabel('fraction of selective post-unit')
set(gca,'XTick',1:3,'XTickLabel',{'t','t+1sec','t+2sec'},'YTick',0:0.2:0.4)
title('Inter-region')
disp('export pdf?')
keyboard
% print('local_inter_coactivate_all.pdf','-dpdf')
exportgraphics(fh,'bin_trans_all_reg_bz.pdf');
end

function [phat,pci]=coactive_all(local)
load('reg_keep.mat');
greymatter=find(cellfun(@(x) ~isempty(regexp(x,'[A-Z]','match','once')), reg_set));
per_reg=zeros(0,12);

curr_code=[];
curr_none=[];
prev_code=[];
prev_none=[];

prev2_code=[];
prev2_none=[];

for bin=1:6
    bzthres=250;
    load(sprintf('0831_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
    s1sel=all(bz_spk_count_S1>bzthres,2);
    reg_chain_S1=bz_conn_reg_S1(s1sel,:);
    pref_chain_S1=bz_pref_S1(s1sel,:);
    
    localselS1=(reg_chain_S1(:,1)==reg_chain_S1(:,2)) & ismember(reg_chain_S1(:,1),greymatter);
    if local
        pref_chain_S1=pref_chain_S1(localselS1,:);
    else
        pref_chain_S1=pref_chain_S1(~localselS1,:);
    end
        
    for i=1:size(pref_chain_S1,1)
        %T+0
        if pref_chain_S1(i,bin)>0
            curr_code(end+1)=pref_chain_S1(i,bin+6);
        else
            curr_none(end+1)=pref_chain_S1(i,bin+6);
        end
        %T+1
        if bin<6 && pref_chain_S1(i,bin)>0
            prev_code(end+1)=pref_chain_S1(i,bin+7);
        elseif bin<6 && pref_chain_S1(i,bin)==0
            prev_none(end+1)=pref_chain_S1(i,bin+7);
        end
        
        %T+2
        if bin<5 && pref_chain_S1(i,bin)>0
            prev2_code(end+1)=pref_chain_S1(i,bin+8);
        elseif bin<5 && pref_chain_S1(i,bin)==0
            prev2_none(end+1)=pref_chain_S1(i,bin+8);
        end
    end
    
end %per bin end

[phat(1),pci(1,:)]=binofit(nnz(curr_none>0),numel(curr_none));
[phat(2),pci(2,:)]=binofit(nnz(curr_code>0),numel(curr_code));
[phat(3),pci(3,:)]=binofit(nnz(prev_none>0),numel(prev_none));
[phat(4),pci(4,:)]=binofit(nnz(prev_code>0),numel(prev_code));
[phat(5),pci(5,:)]=binofit(nnz(prev2_none>0),numel(prev2_none));
[phat(6),pci(6,:)]=binofit(nnz(prev2_code>0),numel(prev2_code));

[tbl,chi,p]=crosstab((1:numel(curr_code)+numel(curr_none))>numel(curr_code),[curr_code>0,curr_none>0])
[tbl,chi,p]=crosstab((1:numel(prev_code)+numel(prev_none))>numel(prev_code),[prev_code>0,prev_none>0])
[tbl,chi,p]=crosstab((1:numel(prev2_code)+numel(prev2_none))>numel(prev2_code),[prev2_code>0,prev2_none>0])
end

function fh=plot_coactive_per_reg(per_reg,figTitle)
load('reg_keep.mat');
greymatter=find(cellfun(@(x) ~isempty(regexp(x,'[A-Z]','match','once')), reg_set));

fh=figure('Color','w','Position',[100,100,900,220]);

subplot(1,3,1)
hold on
plot([0,1],[0,1],'r:')
for i=greymatter'%size(per_reg,1) %
    if per_reg(i,12)>=50 && per_reg(i,10)>=50
        [~,~,p]=crosstab(1:per_reg(i,2)+per_reg(i,4)>per_reg(i,2),[1:per_reg(i,2)>per_reg(i,1),1:per_reg(i,4)>per_reg(i,3)]);
        if p<0.05
            mkColor='r';
        else
            mkColor='k';
        end
        scatter(per_reg(i,3)./per_reg(i,4),per_reg(i,1)./per_reg(i,2),12,'MarkerEdgeColor','none','MarkerFaceColor',mkColor,'MarkerFaceAlpha',0.5)
    end
end
xlabel('t, pre-unit non-selective');
ylabel('t, pre-unit selective');
xlim([0,0.7])
ylim([0,0.7])
set(gca,'XTick',[0,0.5],'YTick',[0,0.5]);

subplot(1,3,2)
hold on
plot([0,1],[0,1],'r:')
for i=1:size(per_reg,1)
    if per_reg(i,12)>=50 && per_reg(i,10)>=50
        [~,~,p]=crosstab(1:per_reg(i,6)+per_reg(i,8)>per_reg(i,6),[1:per_reg(i,6)>per_reg(i,5),1:per_reg(i,8)>per_reg(i,7)]);
        if p<0.05
            mkColor='r';
        else
            mkColor='k';
        end
        scatter(per_reg(i,7)./per_reg(i,8),per_reg(i,5)./per_reg(i,6),12,'MarkerEdgeColor','none','MarkerFaceColor',mkColor,'MarkerFaceAlpha',0.5)
    end
end
xlabel('t+1, pre-unit non-selective');
ylabel('t+1, pre-unit selective');
xlim([0,0.7])
ylim([0,0.7])
set(gca,'XTick',[0,0.5],'YTick',[0,0.5]);

subplot(1,3,3)
hold on
plot([0,1],[0,1],'r:')
for i=1:size(per_reg,1)
    if per_reg(i,12)>=50 && per_reg(i,10)>=50
        [~,~,p]=crosstab(1:per_reg(i,10)+per_reg(i,12)>per_reg(i,10),[1:per_reg(i,10)>per_reg(i,9),1:per_reg(i,12)>per_reg(i,11)]);
        if p<0.05
            mkColor='r';
        else
            mkColor='k';
        end
        scatter(per_reg(i,11)./per_reg(i,12),per_reg(i,9)./per_reg(i,10),12,'MarkerEdgeColor','none','MarkerFaceColor',mkColor,'MarkerFaceAlpha',0.5)
    end
end
xlabel('t+2, pre-unit non-selective');
ylabel('t+2, pre-unit selective');
xlim([0,0.7])
ylim([0,0.7])
set(gca,'XTick',[0,0.5],'YTick',[0,0.5]);

sgtitle(figTitle)
end

function per_reg=coactive_per_reg(type)
load('reg_keep.mat');
% greymatter=find(cellfun(@(x) ~isempty(regexp(x,'[A-Z]','match','once')), reg_set));
per_reg=zeros(0,12);
for onereg=1:length(reg_set)
    curr_code=[];
    curr_none=[];
    prev_code=[];
    prev_none=[];
    
    prev2_code=[];
    prev2_none=[];
    
    for bin=1:6
        load(sprintf('0831_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
        if strcmp(type,'local')
            localselS1=(reg_chain_S1(:,1)==reg_chain_S1(:,2)) & reg_chain_S1(:,1)==onereg;
        elseif strcmp(type,'inter_input')
            localselS1=(reg_chain_S1(:,1)~=reg_chain_S1(:,2)) & reg_chain_S1(:,2)==onereg;
        elseif strcmp(type,'inter_output')
            localselS1=(reg_chain_S1(:,1)~=reg_chain_S1(:,2)) & reg_chain_S1(:,1)==onereg;
        end
        
        pref_chain_S1=pref_chain_S1(localselS1,:);
        
        for i=1:size(pref_chain_S1,1)
            %T+0
            if pref_chain_S1(i,bin)>0
                curr_code(end+1)=pref_chain_S1(i,bin+6);
            else
                curr_none(end+1)=pref_chain_S1(i,bin+6);
            end
            %T+1
            if bin<6 && pref_chain_S1(i,bin)>0
                prev_code(end+1)=pref_chain_S1(i,bin+7);
            elseif bin<6 && pref_chain_S1(i,bin)==0
                prev_none(end+1)=pref_chain_S1(i,bin+7);
            end
            
            %T+2
            if bin<5 && pref_chain_S1(i,bin)>0
                prev2_code(end+1)=pref_chain_S1(i,bin+8);
            elseif bin<5 && pref_chain_S1(i,bin)==0
                prev2_none(end+1)=pref_chain_S1(i,bin+8);
            end
        end
        
    end %per bin end
    per_reg(onereg,:)=[nnz(curr_code),numel(curr_code),...
        nnz(curr_none),numel(curr_none),...
        nnz(prev_code),numel(prev_code),...
        nnz(prev_none),numel(prev_none),...
        nnz(prev2_code),numel(prev2_code),...
        nnz(prev2_none),numel(prev2_none)];
end %per region end
end

%%


function otherstats()


if false
    %% test last current next bin selectivity
    para_pre_bin_S1=[];
    para_ctrl_bin_S1=[];
    
    para_pre_bin_S2=[];
    para_ctrl_bin_S2=[];
    
    para_pre_bin_both=[];
    para_ctrl_bin_both=[];
    for bin=2:5
        disp(bin);
        load(sprintf('0831_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
        for i=1:length(pref_chain_S1)
            if pref_chain_S1(i,bin)>0 && pref_chain_S1(i,bin-1)==0 && pref_chain_S1(i,bin+1)==0
                para_pre_bin_S1(end+1,:)=pref_chain_S1(i,(bin+5):(bin+7));
            elseif pref_chain_S1(i,bin)==0 && pref_chain_S1(i,bin-1)==0 && pref_chain_S1(i,bin+1)==0
                para_ctrl_bin_S1(end+1,:)=pref_chain_S1(i,(bin+5):(bin+7));
            end
        end
        for i=1:length(pref_chain_S2)
            if pref_chain_S2(i,bin)>0 && pref_chain_S2(i,bin-1)==0 && pref_chain_S2(i,bin+1)==0
                para_pre_bin_S2(end+1,:)=pref_chain_S2(i,(bin+5):(bin+7));
            elseif pref_chain_S2(i,bin)==0 && pref_chain_S2(i,bin-1)==0 && pref_chain_S2(i,bin+1)==0
                para_ctrl_bin_S2(end+1,:)=pref_chain_S2(i,(bin+5):(bin+7));
            end
        end
        
        for i=1:length(pref_chain_both)
            if pref_chain_both(i,bin)>0 && pref_chain_both(i,bin-1)==0 && pref_chain_both(i,bin+1)==0
                para_pre_bin_both(end+1,:)=pref_chain_both(i,(bin+5):(bin+7));
            elseif pref_chain_both(i,bin)==0 && pref_chain_both(i,bin-1)==0 && pref_chain_both(i,bin+1)==0
                para_ctrl_bin_both(end+1,:)=pref_chain_both(i,(bin+5):(bin+7));
            end
        end
    end
    figure()
    subplot(1,3,1)
    bar([mean(para_pre_bin_S1>0)',mean(para_ctrl_bin_S1>0)'])
    
    subplot(1,3,2)
    bar([mean(para_pre_bin_S2>0)',mean(para_ctrl_bin_S2>0)'])
    
    subplot(1,3,3)
    bar([mean(para_pre_bin_both>0)',mean(para_ctrl_bin_both >0)'])
end
%% two consec bins


para_pre_bin_S1=[];
para_ctrl_bin_S1=[];

para_pre_bin_S2=[];
para_ctrl_bin_S2=[];

para_pre_bin_both=[];
para_ctrl_bin_both=[];
for bin=2:4
    disp(bin);
    load(sprintf('0831_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
    for i=1:length(pref_chain_S1)
        if all(pref_chain_S1(i,bin:bin+1)>0) && all(pref_chain_S1(i,[bin-1,bin+2])==0)
            para_pre_bin_S1(end+1,:)=pref_chain_S1(i,(bin+5):(bin+8));
        elseif all(pref_chain_S1(i,bin-1:bin+2)==0)
            para_ctrl_bin_S1(end+1,:)=pref_chain_S1(i,(bin+5):(bin+8));
        end
    end
    for i=1:length(pref_chain_S2)
        if all(pref_chain_S2(i,bin:bin+1)>0) && all(pref_chain_S2(i,[bin-1,bin+2])==0)
            para_pre_bin_S2(end+1,:)=pref_chain_S2(i,(bin+5):(bin+8));
        elseif all(pref_chain_S2(i,bin-1:bin+2)==0)
            para_ctrl_bin_S2(end+1,:)=pref_chain_S2(i,(bin+5):(bin+8));
        end
    end
    
    for i=1:length(pref_chain_both)
        if all(pref_chain_both(i,bin:bin+1)>0) && all(pref_chain_both(i,[bin-1,bin+2])==0)
            para_pre_bin_both(end+1,:)=pref_chain_both(i,(bin+5):(bin+8));
        elseif all(pref_chain_both(i,bin-1:bin+2)==0)
            para_ctrl_bin_both(end+1,:)=pref_chain_both(i,(bin+5):(bin+8));
        end
    end
end

pbin=nan(1,4);
for bin=1:4
    [~,~,pbin(bin)]=crosstab([zeros(length(para_pre_bin_S1),1);ones(length(para_ctrl_bin_S1),1)],...
        [para_pre_bin_S1(:,bin)>0;para_ctrl_bin_S1(:,bin)>0]);
end
disp(pbin);

figure('Color','w','Position',[100,100,230,260])
bh=bar([mean(para_pre_bin_S1>0)',mean(para_ctrl_bin_S1>0)'],0.75);
bh(1).FaceColor='k';
bh(1).EdgeColor='k';
bh(2).FaceColor='w';
bh(2).EdgeColor='k';
ylim([0,0.6])
set(gca(),'YTick',0:0.2:0.6,'XTick',1:4,...
    'XTickLabel',{'prev','select.','select.','next'})
ylabel('')


% set(gca(),'YTick',0:0.2:0.6,'XTick',[])
legend([bh(1),bh(2)],{'with selective input','with non-select. input'})
ylabel('proportion of selective neuron');
exportgraphics(gcf(),'with_sel_input.pdf','Resolution',300)

return
figure()
bar([mean(para_pre_bin_S1>0)',mean(para_ctrl_bin_S1>0)'])

subplot(1,3,2)
bar([mean(para_pre_bin_S2>0)',mean(para_ctrl_bin_S2>0)'])

subplot(1,3,3)
bar([mean(para_pre_bin_both>0)',mean(para_ctrl_bin_both >0)'])


%% 3 consec bins
if false
    
    para_pre_bin_S1=[];
    para_ctrl_bin_S1=[];
    
    para_pre_bin_S2=[];
    para_ctrl_bin_S2=[];
    
    para_pre_bin_both=[];
    para_ctrl_bin_both=[];
    for bin=2:3
        disp(bin);
        load(sprintf('0629_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
        for i=1:length(pref_chain_S1)
            if all(pref_chain_S1(i,bin:bin+2)>0) && all(pref_chain_S1(i,[bin-1,bin+3])==0)
                para_pre_bin_S1(end+1,:)=pref_chain_S1(i,(bin+5):(bin+9));
            elseif all(pref_chain_S1(i,bin-1:bin+3)==0)
                para_ctrl_bin_S1(end+1,:)=pref_chain_S1(i,(bin+5):(bin+9));
            end
        end
        for i=1:length(pref_chain_S2)
            if all(pref_chain_S2(i,bin:bin+2)>0) && all(pref_chain_S2(i,[bin-1,bin+3])==0)
                para_pre_bin_S2(end+1,:)=pref_chain_S2(i,(bin+5):(bin+9));
            elseif all(pref_chain_S2(i,bin-1:bin+3)==0)
                para_ctrl_bin_S2(end+1,:)=pref_chain_S2(i,(bin+5):(bin+9));
            end
        end
        
        for i=1:length(pref_chain_both)
            if all(pref_chain_both(i,bin:bin+2)>0) && all(pref_chain_both(i,[bin-1,bin+3])==0)
                para_pre_bin_both(end+1,:)=pref_chain_both(i,(bin+5):(bin+9));
            elseif all(pref_chain_both(i,bin-1:bin+3)==0)
                para_ctrl_bin_both(end+1,:)=pref_chain_both(i,(bin+5):(bin+9));
            end
        end
    end
    figure('Color','w','Position',[100,100,230,280])
    bar([mean(para_pre_bin_S1>0)',mean(para_ctrl_bin_S1>0)'])
    
    figure()
    subplot(1,3,2)
    bar([mean(para_pre_bin_S2>0)',mean(para_ctrl_bin_S2>0)'])
    
    subplot(1,3,3)
    bar([mean(para_pre_bin_both>0)',mean(para_ctrl_bin_both >0)'])
end
%%








pre_post_bin=[];
for i=1:length(pref_chain_S1)
    pre=(pref_chain_S1(i,1:6)>0)*(1:6)';  %TODO should do both pref >0 and ==1 test
    post=(pref_chain_S1(i,7:12)>0)*(1:6)';
    pre_post_bin(end+1,:)=[pre,post];
end
% doesn't really work
% scatter(pre_post_bin(:,1),pre_post_bin(:,2),20,'o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.01);

dbin=diff(pre_post_bin,1,2);
histogram(dbin)




pre_post_bin=[];
for i=1:length(pref_chain_S2)
    %     pre=(pref_chain_S1(i,1:6)>0).*(1:6);  %TODO should do both pref >0 and ==1 test
    %     post=(pref_chain_S1(i,7:12)>0).*(1:6);
    pre1st=find(pref_chain_S2(i,1:6)>0,1);
    post1st=find(pref_chain_S2(i,7:12)>0,1);
    pre_post_bin(end+1,:)=[pre1st,post1st];
end
% doesn't really work
% scatter(pre_post_bin(:,1),pre_post_bin(:,2),20,'o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.01);

dbin=diff(pre_post_bin,1,2);
histogram(dbin)

end

%% old bar graph

% figure('Color','w','Position',[100,100,235,260])
% hold on
% bar(1,phat(1),0.8,'EdgeColor','k','FaceColor','k')
% bar(4,phat(3),0.8,'EdgeColor','k','FaceColor','k')
% bar(7,phat(5),0.8,'EdgeColor','k','FaceColor','k')
%
% bar(2,phat(2),0.8,'EdgeColor','k','FaceColor','w')
% bar(5,phat(4),0.8,'EdgeColor','k','FaceColor','w')
% bar(8,phat(6),0.8,'EdgeColor','k','FaceColor','w')
%
%
%
% errorbar([1 2 4 5 7 8],phat,pci(:,1)'-phat,pci(:,2)'-phat,'.','Color',[0.4,0.4,0.4],'LineWidth',1,'CapSize',10)
% xlim([0.25,8.75])
% set(gca,'YTick',0:0.2:0.4,'XTick',[1 2 4 5 7 8],'XTickLabel',{'t->t','t->t','t->t+1','t->t+1','t->t+2','t->t+2'},'XTickLabelRotation',90);
% ylabel('post-synaptic coding fraction');

%     exportgraphics(gcf,'bin_transfer_coding.pdf','ContentType','vector');