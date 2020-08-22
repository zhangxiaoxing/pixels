
%% one side
local=true;

if true
    load('reg_keep.mat');
    greymatter=find(cellfun(@(x) ~isempty(regexp(x,'[A-Z]','match','once')), reg_set));
    
    
    curr_code=[];
    curr_none=[];
    prev_code=[];
    prev_none=[];
    
    prev2_code=[];
    prev2_none=[];
    
    for bin=1:6
        disp(bin);
        load(sprintf('0814_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
        localselS1=(reg_chain_S1(:,1)==reg_chain_S1(:,2)) & ismember(reg_chain_S1(:,1),greymatter);
%         localselS2=(reg_chain_S2(:,1)==reg_chain_S2(:,2)) & ismember(reg_chain_S2(:,1),greymatter);
        if local
            pref_chain_S1=pref_chain_S1(localselS1,:);
%             pref_chain_S2=pref_chain_S2(localselS2,:);

        else
            pref_chain_S1=pref_chain_S1(~localselS1,:);
%             pref_chain_S2=pref_chain_S2(~localselS2,:);
        end
%             pref_chain_S1=[pref_chain_S1;pref_chain_S2];
        
        for i=1:length(pref_chain_S1)
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

    end

    [phat(1),pci(1,:)]=binofit(nnz(curr_none>0),numel(curr_none));
    [phat(2),pci(2,:)]=binofit(nnz(curr_code>0),numel(curr_code));
    [phat(3),pci(3,:)]=binofit(nnz(prev_none>0),numel(prev_none));
    [phat(4),pci(4,:)]=binofit(nnz(prev_code>0),numel(prev_code));
    [phat(5),pci(5,:)]=binofit(nnz(prev2_none>0),numel(prev2_none));
    [phat(6),pci(6,:)]=binofit(nnz(prev2_code>0),numel(prev2_code));
    
    
    figure('Color','w','Position',[100,100,235,260])
    hold on
    bar(1,phat(1),0.8,'EdgeColor','k','FaceColor','k')
    bar(4,phat(3),0.8,'EdgeColor','k','FaceColor','k')
    bar(7,phat(5),0.8,'EdgeColor','k','FaceColor','k')
    
    bar(2,phat(2),0.8,'EdgeColor','k','FaceColor','w')
    bar(5,phat(4),0.8,'EdgeColor','k','FaceColor','w')
    bar(8,phat(6),0.8,'EdgeColor','k','FaceColor','w')
    
    
    
    errorbar([1 2 4 5 7 8],phat,pci(:,1)'-phat,pci(:,2)'-phat,'.','Color',[0.4,0.4,0.4],'LineWidth',1,'CapSize',10)
    xlim([0.25,8.75])
    set(gca,'YTick',0:0.2:0.4,'XTick',[1 2 4 5 7 8],'XTickLabel',{'PlaceHolder','PlaceHolder','PlaceHolder','PlaceHolder','PlaceHolder','PlaceHolder'},'XTickLabelRotation',90);
    ylabel('post-synaptic coding fraction');
    exportgraphics(gcf,'bin_transfer_coding.pdf','ContentType','vector');


    [tbl,chi,p]=crosstab((1:numel(curr_code)+numel(curr_none))>numel(curr_code),[curr_code>0,curr_none>0])
    [tbl,chi,p]=crosstab((1:numel(prev_code)+numel(prev_none))>numel(prev_code),[prev_code>0,prev_none>0])
    [tbl,chi,p]=crosstab((1:numel(prev2_code)+numel(prev2_none))>numel(prev2_code),[prev2_code>0,prev2_none>0])
    return
end

%%





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
        load(sprintf('0629_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
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
    load(sprintf('0629_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
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