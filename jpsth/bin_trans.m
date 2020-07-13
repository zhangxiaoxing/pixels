
%load 0629_selec_conn_chain_duo_6s_1_2.mat

% pre_post_bin=[];
% for i=1:length(pref_chain_S1)
% %     pre=(pref_chain_S1(i,1:6)>0).*(1:6);  %TODO should do both pref >0 and ==1 test
% %     post=(pref_chain_S1(i,7:12)>0).*(1:6);
%       pre1st=find(pref_chain_S1(i,1:6)>0,1);
%       post1st=find(pref_chain_S1(i,7:12)>0,1);
%       pre_post_bin(end+1,:)=[pre1st,post1st];
% end
% % doesn't really work
% % scatter(pre_post_bin(:,1),pre_post_bin(:,2),20,'o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.01);
% 
% dbin=diff(pre_post_bin,1,2);
% histogram(dbin)

%% test of first selective bin. result negative
% pre_post_bin=[];
% for i=1:length(pref_chain_S1)
% %     pre=(pref_chain_S1(i,1:6)>0).*(1:6);  %TODO should do both pref >0 and ==1 test
% %     post=(pref_chain_S1(i,7:12)>0).*(1:6);
%       pre1st=find(pref_chain_S1(i,1:6)>0,1);
%       post1st=find(pref_chain_S1(i,7:12)>0,1);
% %       disp(pref_chain_S1(i,1:6));
% %       disp(pref_chain_S1(i,7:12));
%       
% %       keyboard
%       if ~isempty(pre1st) && ~isempty(post1st)
%         pre_post_bin(end+1,:)=[pre1st,post1st];
% %         disp(post1st-pre1st);
%       end
% end
% 
% dbin=diff(pre_post_bin,1,2);
% histogram(dbin)
% 
% return

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