function ridge_fit_plot(tasks,featRange, effs, lambda,use_mean, featFile, plot_cv_curve,subset)


%% plot
close all
fl=ls('Ridge_frac*nphr*.mat');
plot_cv_curve=false;
plot_table=true;

full_features_tag={'','selective early','selective late','conn dens early','conn dens late'};

for i=1:size(fl,1)
    disp(fl(i,:))
    load(strtrim(fl(i,:)));
    int_result=int_result(strcmp(int_result(:,5),'linear'),:);
    features_tag=full_features_tag(2:end);
%     fh=figure('Color','w','Position',[50,50,960,720],'DefaultAxesFontSize',10);
%     subplot(2,3,1);
%     hold on;
%     [~,I]=sort(cell2mat(int_result(:,3)));
% 
%     yyaxis right;
%     if plot_cv_curve
%         h_cv_rsq=plot(cell2mat(int_result(I,9)),'b-');
%     end
%     hrsq=plot(cell2mat(int_result(I,3)),'r-','LineWidth',2);
%     ylabel('r-squared');
%     set(gca(),'YColor','r','XColor','k','FontSize',10)
%     
%     yyaxis left;
%     hAIC=plot(cell2mat(int_result(I,2)),'k-');
%     ylabel('Akaike Information Criterion');
%     set(gca(),'YColor','k')
%     if plot_cv_curve
%         legend([hAIC,hrsq,h_cv_rsq],{'AIC','r-squared','CV rsq.'})
%     else
%         legend([hAIC,hrsq],{'AIC','r-squared'})
%     end
%     
%     xlabel('Model #');


    
%     [~,Iaic]=min(cell2mat(int_result(:,2)));
    sigIdces=find(cell2mat(int_result(:,4))<0.05);
    if isempty(sigIdces)
        continue
    end
    sigIdces=reshape(sigIdces,1,[]);
    for Iaic=15
        figure('Color','w','Position',[100,100,1200,400])
        subplot(1,3,1);
%         [~,Iaic]=max(cell2mat(int_result(:,3)));
        int_result{Iaic,7}=int_result{Iaic,7}-mean(int_result{Iaic,7});
        plot(int_result{Iaic,7}(:,2),int_result{Iaic,7}(:,1),'r.','MarkerSize',10);
        arrayfun(@(x) text(int_result{Iaic,7}(x,2),...
            int_result{Iaic,7}(x,1)-0.02*diff(ylim()),...
            int2str(x), 'Color','k','HorizontalAlignment','center','VerticalAlignment','middle'),...
        1:size(int_result{Iaic,7},1))

        text(min(xlim())+0.15*diff(xlim()),min(ylim())+0.85*diff(ylim()),sprintf('rsq = %.3f, p = %.3f',int_result{Iaic,3},int_result{Iaic,4}),'HorizontalAlignment','left')
        xlim(diff(xlim())*0.05*[-1,1]+xlim());
        ylim(diff(ylim())*0.1*[-1,1]+ylim());
        xlabel('Ephys GLM prediction');
        ylabel('optogenetic effect size')
        set(gca(),'FontSize',10)

    %     cv_results=int_result{Iaic,8};
    %     r=int_result{Iaic,9};
    %     p=int_result{Iaic,10};

    %% we don't have cv at the moment
    %     subplot(2,3,5)
    %     cv_results=cv_results-mean(cv_results);
    %     plot(cv_results(:,2),cv_results(:,1),'b.','MarkerSize',10);
    %     arrayfun(@(x) text(cv_results(x,2),...
    %         cv_results(x,1)-0.02*diff(ylim()),...
    %         int2str(x), 'Color','k','HorizontalAlignment','center','VerticalAlignment','middle'),...
    %     1:size(int_result{Iaic,7},1))
    %     text(min(xlim())+0.15*diff(xlim()),min(ylim())+0.85*diff(ylim()),sprintf('rsq = %.3f, p = %.3f',r*r,p),'HorizontalAlignment','left')
    %     xlabel('Leave-one-region-out prediction');
    %     ylabel('optogenetic effect size')
    %     set(gca(),'FontSize',10)

    %%    

        def=zeros(1,numel(int_result{Iaic,1}));
        for j=1:numel(int_result{Iaic,1})
            if strcmp(int_result{Iaic,5},'interaction') && j>=4
                def(j)=int_result{Iaic,1}(j)*0.01;
            else
                def(j)=int_result{Iaic,1}(j)*0.1;
            end
        end
        [~,I]=sort(abs(def),'descend');



        subplot(1,3,2)
        lh=[];
        hold on
        colors={'r','m','b','c','k'};
        for coeIdx=1:min(4,numel(I))
            if strcmp(int_result{Iaic,5},'interaction') && j>=4
                lh(coeIdx)=plot(0:0.01:0.1,((0:0.01:0.1).^2)*int_result{Iaic,1}(I(coeIdx)),'-','Color',colors{coeIdx});
            else
                lh(coeIdx)=plot(0:0.01:0.1,(0:0.01:0.1)*int_result{Iaic,1}(I(coeIdx)),'-','Color',colors{coeIdx});
            end
        end
        legend(lh,features_tag{int_result{Iaic,6}(I(1:min(4,numel(I))))-1},'Interpreter','none','FontSize',10,'Location','southwest')
        ylabel('Delta effect size')
        xlabel('Delta fraction of neurons')
        set(gca(),'FontSize',10) 
    %     print(sprintf('GLM_explain%d.eps',i),'-depsc','-painters','-r300')

    %      if contains(fl(i,:),'EFF')
    %         tiStr=strjoin({'Effect size',...
    %             regexp(fl(i,:),'(Mean|Median)','match','once'),...
    %             replace(regexp(fl(i,:),'lambda_[0-9\.]*_','match','once'),'_',' '),...
    %             replace(regexp(fl(i,:),'(?<=EFF_)(ED|LD|DM).*(?=\.mat)','match','once'),'-',' ')});
    %     else
    %         tiStr=strjoin({'Optogenetic specificity',...
    %             regexp(fl{i},'(Mean|Median)','match','once'),...
    %             replace(regexp(fl{i},'lambda_[0-9\.]*_(ED|LD|DM).*(?=\.mat)','match','once'),'_',' ')});
    %     end
    %     


            tiStr=strjoin({'Effect size',...
                regexp(fl(i,:),'(Mean|Median)','match','once'),...
                replace(regexp(fl(i,:),'lambda_[0-9\.]*_','match','once'),'_',' '),...
                replace(regexp(fl(i,:),'(?<=EFF_).*(?=\.mat)','match','once'),'-',' '),num2str(Iaic)});

        if plot_table

            uit=uitable();
            t=cell(size(regions,1),1);

            t(:,1)=regions(:,2);
            uit.Data=t;
            uit.RowName=num2cell((1:size(regions,1))');
            uit.Position=[840,20,320,360];
            uit.FontSize=10;
    %         uit.ColumnWidth={50,50};
        end
    sgtitle(replace(tiStr,'_',' '));
        print(strcat(replace(tiStr,' ','_'),'.png'),'-dpng','-painters','-r200');
        close all
    end
end


end

