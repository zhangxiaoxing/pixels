to_save=false;
close all;
gen_ratio_map=true;
if gen_ratio_map
    load('reg_keep.mat');
    connkeep=false(length(reg_set),1);
    pair_list=cell(6,1);
    conn_list=cell(6,1);
    ratio_list=cell(6,1);
    for bin=1:6
        load(sprintf('0626_selec_pair_mat_duo_6s_%d_%d.mat',bin,bin+1));
        pair_list{bin}=pair_mat;
        for i=1:length(pair_mat)
            if nnz(pair_mat(i,:)>=40)+nnz(pair_mat(:,i)>=40)>=20
                connkeep(i)=true;
            end
        end
    end
    grayMatter=cellfun(@(x) ~isempty(regexp(x,'[A-Z]','match','once')), reg_set);
    ckeep=connkeep & grayMatter;
    reg_keep=reg_set();
    nnz(ckeep) %likely around 27
    save('reg_keep.mat','reg_keep','reg_set','grayMatter','connkeep','ckeep')
    ratio_sum=zeros(nnz(ckeep),1);
    count_all=zeros(length(reg_set),1);
    for bin =1:6
        load(sprintf('0626_selec_conn_mat_duo_6s_%d_%d.mat',bin,bin+1));
        conn_mat=conn_mat_S1;
        conn_list{bin}=conn_mat;
        conn_mat_k=conn_mat(ckeep,ckeep);
        pair_mat_k=pair_list{bin}(ckeep,ckeep);
        ratio_mat=conn_mat_k./pair_mat_k;
        ratio_list{bin}=ratio_mat;
        for i=1:length(ratio_mat)
            ratio_sum(i)=ratio_sum(i)+nansum(ratio_mat(i,:))+nansum(ratio_mat(:,i))-2*ratio_mat(i,i);
        end
        for i=1:length(reg_set)
            count_all(i)=count_all(i)+sum(conn_mat(i,:))+sum(conn_mat(:,i))-2*conn_mat(i,i);
        end
    end
%     keyboard
    [~,ratioIdx]=sort(ratio_sum);
    for bin=1:6
        figure('Color','w','Position',[100,100,450,400])
        h=imagesc(ratio_list{bin}(ratioIdx,ratioIdx),[0,0.2]);
        
        set(h,'alphadata',~isnan(ratio_list{bin}(ratioIdx,ratioIdx))); 
        ax=gca();
        ax.YTick=(1:nnz(ckeep));
        ax.YTickLabel=reg_keep(ratioIdx);
      
        ax.XTick=(1:nnz(ckeep));
        ax.XTickLabel=reg_keep(ratioIdx);
        ax.XTickLabelRotation=90;
        
        ax.YDir='normal';
        ax.Color=[0.5,0.5,0.5];
        ylabel('target');
        xlabel('source');
        colorbar;
        if to_save
            print(sprintf('ratio_map_%d_%d.pdf',bin,bin+1),'-dpdf','-painters','-r300')
            print(sprintf('ratio_map_%d_%d.png',bin,bin+1),'-dpng','-painters','-r300')
        end
    end
    
    for bin=1:6
        figure('Color','w','Position',[100,100,450,400])
        conn_mat_k=conn_list{bin}(ckeep,ckeep);
        h=imagesc(conn_mat_k(ratioIdx,ratioIdx),[0,100]);
        
        set(h,'alphadata',~isnan(ratio_list{bin}(ratioIdx,ratioIdx))); 
        ax=gca();
        ax.YTick=(1:nnz(ckeep));
        ax.YTickLabel=reg_keep(ratioIdx);
      
        ax.XTick=(1:nnz(ckeep));
        ax.XTickLabel=reg_keep(ratioIdx);
        ax.XTickLabelRotation=90;
        
        ax.YDir='normal';
        ax.Color=[0.5,0.5,0.5];
        ylabel('target');
        xlabel('source');
        colormap('cool')
        colorbar;
        if to_save
            print(sprintf('count_map_%d_%d.pdf',bin,bin+1),'-dpdf','-painters','-r300')
            print(sprintf('count_map_%d_%d.png',bin,bin+1),'-dpng','-painters','-r300')
        end
    end
    keyboard
%     save('reg_keep.mat','reg_keep','keep');
    
        
return    
%     for bin=1:6
%         figure('Color','w','Position',[100,100,450,400])
%         h=imagesc(conn_list{bin}(ratioIdx,ratioIdx),[0,50]);
%         
%         set(h,'alphadata',~isnan(ratio_list{bin}(ratioIdx,ratioIdx))); 
%         ax=gca();
%         ax.YTick=(1:nnz(keep));
%         ax.YTickLabel=reg_keep(ratioIdx);
%       
%         ax.XTick=(1:nnz(keep));
%         ax.XTickLabel=reg_keep(ratioIdx);
%         ax.XTickLabelRotation=90;
%         
%         ax.YDir='normal';
%         ax.Color=[0.5,0.5,0.5];
%         ylabel('target');
%         xlabel('source');
%         colormap('jet')
%         colorbar;
%         
%         print(sprintf('count_map_%d_%d.pdf',bin,bin+1),'-dpdf','-painters','-r300')
%         print(sprintf('count_map_%d_%d.png',bin,bin+1),'-dpng','-painters','-r300')
%     end
end
