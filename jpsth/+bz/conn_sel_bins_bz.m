% assumes all stats file in memory for data generation

if ~exist('su_var_t_bz.mat','file')
    su_var_t=[];
    for bin=1:6
        disp(bin);
        bzthres=250;
        load(sprintf('0831_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
        s1sel=all(bz_spk_count_S1>bzthres,2);
%         s2sel=all(bz_spk_count_S2>bzthres,2);
        bz_conn_chain_S1=bz_conn_chain_S1(s1sel,:);
        bz_pref_S1=bz_pref_S1(s1sel,:);
        %% TODO: brain region
        for i=1:length(bz_conn_chain_S1)
            key=bz_conn_chain_S1(i,1)*100000+rem(bz_conn_chain_S1(i,2),100000);
            if ~isempty(su_var_t)
                idx=find(su_var_t(:,1)==key);
            else
                idx=[];
            end
            if bz_pref_S1(i,bin)>0 && bz_pref_S1(i,bin+6)>0
                selType=2;
            else
                selType=1;
            end
            if isempty(idx)
                su_var_t(end+1,:)=zeros(1,7);
                su_var_t(end,1)=key;
                su_var_t(end,bin+1)=selType;
            else
                su_var_t(idx,bin+1)=selType;
            end
        end
    end
    save('su_var_t_bz.mat','su_var_t');
else
    load('su_var_t_bz.mat','su_var_t');
end


for bin=2:6
    d_conn_pair(bin)=nnz((su_var_t(:,bin)>0) ~= (su_var_t(:,bin+1)>0));
    d_conn_num(bin)=(nnz(su_var_t(:,bin+1)>0)-nnz(su_var_t(:,bin)>0));

    d_sel_pair(bin)=nnz((su_var_t(:,bin)>1) ~= (su_var_t(:,bin+1)>1));
    d_sel_num(bin)=(nnz(su_var_t(:,bin+1)>1)-nnz(su_var_t(:,bin)>1));
end


to_plot=true;
if to_plot

    %assumes su_var_t etc in memory
    close all
    if true
        figure('Color','w','Position',[20,20,280,700])
        %% TODO: sort 
        sort_mat=su_var_t(:,2:end)>0;
        sel_mat=su_var_t(:,2:end)>1;

        sortsum=sort_mat*(2.^(11:-1:6)')+sel_mat*(2.^(5:-1:0)');
        [~,sumI]=sort(sortsum,'descend');

        imagesc(su_var_t(sumI,2:end));
        cmap=[1 1 1;0 0 1;1 0 0];
        colormap(cmap)
        xlabel('delay time (1-sec bin)')
        ylabel('pair #')
        set(gca(),'YTick',[0,100000,200000],'YTickLabel',0:100000:200000);
        ylim([-0.5,max(ylim())]);
        % title('DP -> TTd')

        keyboard
        exportgraphics(gcf(),'conn_per_bin_change_bz.pdf','Resolution',300)
    end
    
    if true
        figure('Color','w','Position',[100,100,220,250])
    %     subplot(3,1,2);
        bin=1;
        load(sprintf('0831_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
        magicN=length(pair_chain);
        hold on
        yyaxis left
        [phat,pci]=binofit(sum(su_var_t(:,2:end)>0),magicN);
        fill([1:6,6:-1:1]',[pci(:,1);flip(pci(:,2))],'b','EdgeColor','none','FaceAlpha',0.2)
        ph1=plot(phat,'b-','LineWidth',1.5);
        ylim([0,0.002])
        ylabel('connective pair ratio');
        set(gca(),'YTick',0:0.05:0.2)
        yyaxis right
        [phat,pci]=binofit(sum(su_var_t(:,2:end)>1),magicN);
        fill([1:6,6:-1:1]',[pci(:,1);flip(pci(:,2))],'r','EdgeColor','none','FaceAlpha',0.2)
        ph2=plot(phat,'r-','LineWidth',1.5);
        ylabel('selective ratio');
        ylim([0,max(ylim())])
        xlim([0.5,6.5])
        ylim([0,0.0006])
        set(gca(),'YTick',0:0.02:0.04)
        xlabel('delay time(1-sec bin)')
        ylabel('selective pair ratio')
        legend([ph1,ph2],{'connective','selective'});

        keyboard
        exportgraphics(gcf(),'conn_per_bin_change_curve_bz.pdf','Resolution',300)
    end
    
    figure('Color','w','Position',[100,100,220,250])
%     subplot(3,2,4)
    hold on
    [phat,pci]=binofit(abs(d_conn_num(2:end)),length(su_var_t));
    negsel=d_conn_num(2:end)<0;
    phat(negsel)=phat(negsel)*-1;
    pci(negsel,:)=pci(negsel,:)*-1;
    fill([1.5:5.5,5.5:-1:1.5],[pci(:,1);flip(pci(:,2))],'b','FaceAlpha',0.2,'EdgeColor','none');
    phn=plot(1.5:5.5,phat,'-b','LineWidth',1.5);
    
    
    [phat,pci]=binofit(d_conn_pair(2:end),length(su_var_t));
    fill([1.5:5.5,5.5:-1:1.5],[pci(:,1);flip(pci(:,2))],'b','FaceAlpha',0.2,'EdgeColor','none');
    php=plot(1.5:5.5,phat,':b','LineWidth',1.5);
    
    legend([phn,php],{'change in total number','change in underlying pairs'})
    xlabel('delay time (1-sec bin)')
    ylabel('fraction of all pairs')
    xlim([0.5,6.5])
    ylim([-0.1,0.7])
    set(gca,'YTick',0:0.2:0.6)

    keyboard
    exportgraphics(gcf(),'per_bin_change_number_n_pair_bz.pdf','Resolution',300)
    return
    
    %% other plots
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



