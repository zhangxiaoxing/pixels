function ridge_fit_plot(tasks,featRange, effs, lambda,use_mean)
early_late='early';

%% plot
close all
fl=ls('Ridge_frac*vgat*early.mat');
plot_cv_curve=false;
plot_table=false;
if strcmp(early_late,'early')
    full_features_tag={'','selective fraction early','conn density early','decoding early','out-in density'};
elseif strcmp(early_late,'late')
    full_features_tag={'','selective fraction late','conn density late','decoding late','out-in density'};
else
end
for i=1:size(fl,1)
    disp(fl(i,:))
    if contains(fl(i,:),'vgat')
        opsin='vgat';
    else
        opsin='nphr';
    end
    load(strtrim(fl(i,:)));
    int_result=int_result(strcmp(int_result(:,5),'linear'),:);
    features_tag=full_features_tag(2:end);
    sigIdces=find(cell2mat(int_result(:,4))<0.05);
    if isempty(sigIdces)
        continue
    end
    sigIdces=reshape(sigIdces,1,[]);
    for Iaic=13
        fh=figure('Color','w','Position',[100,100,240,500]);
        subplot(2,1,1);
        %         [~,Iaic]=max(cell2mat(int_result(:,3)));
        int_result{Iaic,7}=int_result{Iaic,7}-mean(int_result{Iaic,7});
        plot(int_result{Iaic,7}(:,2),int_result{Iaic,7}(:,1),'r.','MarkerSize',10);
        arrayfun(@(x) text(int_result{Iaic,7}(x,2),...
            int_result{Iaic,7}(x,1)-0.02*diff(ylim()),...
            regions(x), 'Color','k','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8),...
            1:size(int_result{Iaic,7},1))
        
        text(min(xlim())+0.15*diff(xlim()),min(ylim())+0.85*diff(ylim()),sprintf('r=%.3f, p=%.3f',int_result{Iaic,3},int_result{Iaic,4}),'HorizontalAlignment','left')
        xlim(diff(xlim())*0.05*[-1,1]+xlim());
        ylim(diff(ylim())*0.1*[-1,1]+ylim());
        xlabel('Ephys GLM prediction');
        ylabel('optogenetic effect size')
        set(gca(),'FontSize',10)
        
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
        
        
        
        subplot(2,1,2)
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
        
        tiStr=strjoin({opsin, 'Effect size',...
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
%         print(strcat(replace(tiStr,' ','_'),'.png'),'-dpng','-painters','-r200');
        exportgraphics(fh,strcat(replace(tiStr,' ','_'),'.pdf'),'ContentType','vector');
        close all
    end
end


end

