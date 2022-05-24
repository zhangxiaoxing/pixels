% assumes all stats file in memory for data generation

if false%~exist('su_var_t.mat','file')
    su_var_t_S1=[];
    su_var_t_S2=[];
    for bin=1:6
        disp(bin);
        load(sprintf('0116_memory_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
        %% TODO: brain region
        for i=1:length(conn_chain_S1)
            key=conn_chain_S1(i,1)*100000+rem(conn_chain_S1(i,2),100000);
            if ~isempty(su_var_t_S1)
                idx=find(su_var_t_S1(:,1)==key);
            else
                idx=[];
            end
            if pref_chain_S1(i,bin) == pref_chain_S1(i,bin+6) && pref_chain_S1(i,bin)>0
                selType=3;
            elseif max(pref_chain_S1(i,1:6)) == max(pref_chain_S1(i,7:12)) && max(pref_chain_S1(i,1:6))>0
                selType=2;
            else
                selType=1;
            end
            if isempty(idx)
                su_var_t_S1(end+1,:)=zeros(1,7);
                su_var_t_S1(end,1)=key;
                su_var_t_S1(end,bin+1)=selType;
            else
                su_var_t_S1(idx,bin+1)=selType;
            end
        end
        
        for i=1:length(conn_chain_S2)
            key=conn_chain_S2(i,1)*100000+rem(conn_chain_S2(i,2),100000);
            if ~isempty(su_var_t_S2)
                idx=find(su_var_t_S2(:,1)==key);
            else
                idx=[];
            end
            if pref_chain_S2(i,bin) == pref_chain_S2(i,bin+6) && pref_chain_S2(i,bin)>0
                selType=3;
            elseif max(pref_chain_S2(i,1:6)) == max(pref_chain_S2(i,7:12)) && max(pref_chain_S2(i,1:6))>0
                selType=2;
            else
                selType=1;
            end
            if isempty(idx)
                su_var_t_S2(end+1,:)=zeros(1,7);
                su_var_t_S2(end,1)=key;
                su_var_t_S2(end,bin+1)=selType;
            else
                su_var_t_S2(idx,bin+1)=selType;
            end
        end
        
        
    end
    save('su_var_t_0116.mat','su_var_t_S1','su_var_t_S2');
else
    load('su_var_t.mat','su_var_t_S1','su_var_t_S2');
end

su_var_t_all=[su_var_t_S1;su_var_t_S2];

for bin=2:6
    d_conn_pair(bin)=nnz((su_var_t_all(:,bin)>0) ~= (su_var_t_all(:,bin+1)>0));
    d_conn_num(bin)=(nnz(su_var_t_all(:,bin+1)>0)-nnz(su_var_t_all(:,bin)>0));
    d_conn_mat(:,bin)=(su_var_t_all(:,bin+1)>0)-(su_var_t_all(:,bin)>0);
    
    d_sel_pair(bin)=nnz((su_var_t_all(:,bin)>1) ~= (su_var_t_all(:,bin+1)>1));
    d_sel_num(bin)=(nnz(su_var_t_all(:,bin+1)>1)-nnz(su_var_t_all(:,bin)>1));
end

regain=false(size(su_var_t_all,1),1);
for bin=2:5
    regain(d_conn_mat(bin,:)<0 & max(d_conn_mat((bin+1):end,:),[],2)>0)=true;
end
disp(mean(regain))

%% Delta Conn
fh=figure('Color','w','Position',[100,100,200,150]);
hold on;
[phat,pci]=binofit(sum(d_conn_mat(:,2:end)>0),size(su_var_t_all,1));
phat=phat*100;
pci=pci*100;
fill([2:6,6:-1:2]',[pci(:,1);flip(pci(:,2))],'m','EdgeColor','none','FaceAlpha',0.2)
ph1=plot(2:6,phat,'m-','LineWidth',1);

[phat,pci]=binofit(sum(d_conn_mat(:,2:end)<0),size(su_var_t_all,1));
phat=phat*100;
pci=pci*100;
fill([2:6,6:-1:2]',[pci(:,1);flip(pci(:,2))],'c','EdgeColor','none','FaceAlpha',0.2)
ph2=plot(2:6,phat,'c-','LineWidth',1);
legend([ph1,ph2],{'Gained coupling','Lost coupling'});
xlim([0,6])
ylim([0,25])
ylabel('Change in coupled pairs(%)')
xlabel('Time (1-sec bin)')
set(gca(),'XTick',0:6)
exportgraphics(fh,'Func.coup.gain.loss.pdf')
%%


to_plot=true;
if to_plot

    %assumes su_var_t etc in memory
    close all
    if true
        for s=1:2
            if s==1
                su_var_t=su_var_t_S1;
            else
                su_var_t=su_var_t_S2;
            end
            fh=figure('Color','w','Position',[20,20,215,1000]);
            %% TODO: sort 
            sort_mat=su_var_t(:,2:end)>0;
            sel_mat=su_var_t(:,2:end)>1;

            sortsum=sort_mat*(2.^(11:-1:6)')+sel_mat*(2.^(5:-1:0)');
            [~,sumI]=sort(sortsum,'descend');

            imagesc(su_var_t(sumI,2:end));
            cmap=[1 1 1;0 0 1;0 1 1;1 0 0];
            colormap(cmap)
            xlabel('Time (1-sec bin)')
            ylabel('Coupled memory Neuron #')
            set(gca(),'YTick',0:25000:75000,'YTickLabel',{'0','25K','50K','75K'},'XTick',1:6);
            ylim([-0.5,max(ylim())]);
            % title('DP -> TTd')
            exportgraphics(fh,sprintf('conn_per_bin_change_%d.pdf',s),'ContentType','vector')
            clear su_var_t
        end
    end
    %% Histogram
    hc=histcounts(sum(su_var_t_all(:,2:end)>0,2),1:7)./length(su_var_t_all)*100;
    fh=figure('Color','w','Position',[100,100,215,150]);
    bar(1:6,hc,0.8,'FaceColor','W','EdgeColor','k');
    set(gca,'YScale','log','YTick',[0.1,1,10,100]);
    ylim([0.05,100])
    ylabel('Fraction of coupled pairs(%)')
    xlabel('Coupled bins')
    exportgraphics(fh,'func.coup.bins.pdf')

    %%
    
    
    if true
        fh=figure('Color','w','Position',[100,100,200,150]);
    %     subplot(3,1,2);
        bin=1;
        load(sprintf('0114_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
        pair_candi=length(pair_chain)*2;
        hold on
        
        [phat,pci]=binofit(sum(su_var_t_all(:,2:end)>0),pair_candi);
        phat=phat*100;
        pci=pci*100;
        fill([1:6,6:-1:1]',[pci(:,1);flip(pci(:,2))],'b','EdgeColor','none','FaceAlpha',0.2)
        ph1=plot(phat,'b-','LineWidth',1);
        ylim([0,20])
        ylabel('Coupling pair / all pair (%)');
%         set(gca(),'YTick',0:5:20)
        
        [phat2,pci2]=binofit(sum(su_var_t_all(:,2:end)>1),pair_candi);
        phat2=phat2*100;
        pci2=pci2*100;
 
        
        
        [phat,pci]=binofit(sum(su_var_t_all(:,2:end)>2),pair_candi);
        phat=phat*100;
        pci=pci*100;
        
        pci2=pci2-pci;
        fill([1:6,6:-1:1]',[pci2(:,1);flip(pci2(:,2))],'c','EdgeColor','none','FaceAlpha',0.2)
        fill([1:6,6:-1:1]',[pci(:,1);flip(pci(:,2))],'r','EdgeColor','none','FaceAlpha',0.2)
        
        ph2=plot(phat2-phat,'c-','LineWidth',1);
        ph3=plot(phat,'r-','LineWidth',1);
        
        ylabel('Functional coupling density (%)');
        xlim([0.5,6.5])
%         ylim([0,4])
%         set(gca(),'YTick',0:0.02:0.04)
        xlabel('Delay time(1-sec bin)')
%         ylabel('selective pair ratio')
%         legend([ph1,ph2],{'connective','selective'});
        ax=gca();
%         ax.YAxis(1).Color=[0,0,1];
%         ax.YAxis(2).Color=[1,0,0];
        ax.XTick=1:6;
        
%         legend([ph1,ph2,ph3],{'All func. coupling','Congruent func. coupling','Congruent active func. coupling'})
        keyboard
        exportgraphics(fh,'conn_per_bin_change_curve.pdf','ContentType','vector')
    end
    
    figure('Color','w','Position',[100,100,220,250])
%     subplot(3,2,4)
    hold on
    [phat,pci]=binofit(abs(d_conn_num(2:end)),length(su_var_t));
    negsel=d_conn_num(2:end)<0;
    phat(negsel)=phat(negsel)*-1*100;
    pci(negsel,:)=pci(negsel,:)*-1*100;
    fill([2:6,6:-1:2],[pci(:,1);flip(pci(:,2))],'b','FaceAlpha',0.2,'EdgeColor','none');
    phn=plot(2:6,phat,'-b','LineWidth',1.5);
    
    
    [phat,pci]=binofit(d_conn_pair(2:end),length(su_var_t));
    phat=phat*100;
    pci=pci*100;
    fill([2:6,6:-1:2],[pci(:,1);flip(pci(:,2))],'b','FaceAlpha',0.2,'EdgeColor','none');
    php=plot(2:6,phat,':b','LineWidth',1.5);
    yline(0,'k:');
    legend([phn,php],{'In coupling density','In individual pairs'});
    xlabel('Time (1-sec bin)')
    ylabel('Change in fraction(%)')
    xlim([0,6])
%     ylim([-0.1,0.7])
    set(gca,'XTick',0:6)

    keyboard
    exportgraphics(gcf(),'per_bin_change_number_n_pair.pdf','ContentType','vector')

    figure()
%     subplot(3,2,6)
    hold on
    phn=plot(2:6,d_sel_num(2:end),'-r','LineWidth',1.5);
    php=plot(2:6,d_sel_pair(2:end),'--r','LineWidth',1.5);
    legend([phn,php],{'delta global selectivity','delta pair selectivity'})
    xlabel('delay bin (s)')
    ylabel('fraction of all pairs')

    sgtitle('brain-wide')
    print('conn_dyn_brain-wide.png','-dpng')

    % colorbar()


    key=su_var_t(:,1);
    source=idivide(int32(key),int32(100));
    target=rem(key,100);
    stable=all(su_var_t(:,2:end)');

    stable_sel=all(su_var_t(:,2:end)'>1);


    unisource=unique(source);

    stable_count=zeros(0,2);
    transient_count=zeros(0,2);
    for s=unisource'
        stable_count(end+1,:)=[s,nnz(source==s & stable')];
        transient_count(end+1,:)=[s,nnz(source==s & ~stable')];
    end


    [shist,edges]=histcounts(stable_count(:,2),1:2:60);
    thist=histcounts(transient_count(:,2),1:2:60);
    figure('Color','w','Position',[100,100,450,300])
    bh=bar(2:2:59, [thist',shist'],'stacked');
    legend({'transient','sustained'})
    xlabel('number of target neuron from same source neuron')
    ylabel('number of source neuron')
    print('Neuron-neuron-hist.png','-dpng')


%     figure()
%     histogram(stable_count(:,2),1:2:80)
%     figure()
%     histogram(transient_count(:,2),1:2:80)
end



