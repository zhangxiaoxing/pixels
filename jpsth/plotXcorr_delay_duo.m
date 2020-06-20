% assume 'sums' is loaded in workspace. Otherwise load corresponding
% XCORR_delay_bin.mat file first
currbin=1;
close('all')
prefix='0613_sel_';
to_plot=false;
to_save=true;
prepare_stats_file=false;
if prepare_stats_file
    %debugging
    %fs=dir('0604_nonsel_XCORR_duo_*delay_6_1_*.mat');
    for bin=currbin %:6 % debuging
        % debuging
        %load(fullfile(fs(bin).folder,fs(bin).name));
        
        stats=cell(0);
        thresh=norminv(0.995); %0.05 bonferroni corrected
        bin_range=[bin,bin+1];
        for sidx=1:size(sums,1)
            fprintf('%d of %d\n',sidx, size(sums,1));
            xc_s1=sums{sidx,5};
            xshuf_s1=sums{sidx,6};
            xc_s2=sums{sidx,7};
            xshuf_s2=sums{sidx,8};
%             sustCount=numel(sums{sidx,3});useless?
            for si=1:(size(xc_s1.xcorr,1)-1)
                su1id=str2double(xc_s1.label{si,1});
                if ismember(su1id,sums{sidx,3})
                    su1='sust';
                elseif ismember(su1id,sums{sidx,4})
                    su1='transient';
                else
                    su1='non_selective';
                end
                for sj=(si+1):size(xc_s1.xcorr,2)
                    su2id=str2double(xc_s1.label{sj,1});
                    if ismember(su2id,sums{sidx,3})
                        su2='sust';
                    elseif ismember(su2id,sums{sidx,4})
                        su2='transient';
                    else
                        su2='non_selective';
                    end
                    totalCount=nansum(squeeze(xc_s1.xcorr(si,sj,:)));
                    if numel(xc_s1.cfg.trials)<20 ||  numel(xc_s2.cfg.trials)<20 || totalCount<1000
                        onepair=struct();
                        onepair.fileidx=sidx;
                        onepair.su1_label_idx=si;
                        onepair.su2_label_idx=sj;
                        onepair.su1_sel_type=su1;
                        onepair.su2_sel_type=su2;
                        onepair.totalcount=totalCount;
                        onepair.su1_clusterid=su1id;
                        onepair.su2_clusterid=su2id;
                        onepair.wf_stats_su1=xc_s1.label{si,2}; %cluster_id, FR, err_type,vollay-peak, fwhm
                        onepair.wf_stats_su2=xc_s1.label{sj,2};
                        onepair.wf_su1=xc_s1.label{si,3};
                        onepair.wf_su2=xc_s1.label{sj,3};
                        onepair.prefered_sample_su1=xc_s1.label{si,4};
                        onepair.prefered_sample_su2=xc_s1.label{sj,4};
                        onepair.reg_su1=xc_s1.label{si,5};
                        onepair.reg_su2=xc_s1.label{sj,5};
                        onepair.s1_trials=numel(xc_s1.cfg.trials);
                        onepair.s2_trials=numel(xc_s2.cfg.trials);
                        stats{end+1}=onepair;
                        clear onepair
                        continue
                    end

                    %sample 1 stats
                    hists1=squeeze(xc_s1.xcorr(si,sj,:));
                    shufs1=squeeze(xshuf_s1.shiftpredictor(si,sj,:));
                    diffs1=hists1-smooth(shufs1);
                    diffs1(50:51)=0;
                    stds1=std(shufs1);
                    scores1=diffs1(46:55)./stds1;
                    %any score > thresh
                    if any(scores1>thresh)
                        bincounts1=diffs1(46:55);
                        bincounts1(scores1<=thresh)=0;
                        sumdiffs1=(sum(bincounts1(1:5))-sum(bincounts1(6:10)));
                        if sumdiffs1==0
                            AIs1=0;
                        else
                            AIs1=sumdiffs1/(sum(bincounts1(1:5))+sum(bincounts1(6:end)));
                        end
                    else
                        AIs1=0;
                    end


                    %sample 2 stats, duplicate of sample 1. TODO:DRY
                    hists2=squeeze(xc_s2.xcorr(si,sj,:));
                    shufs2=squeeze(xshuf_s2.shiftpredictor(si,sj,:));
                    diffs2=hists2-smooth(shufs2);
                    diffs2(50:51)=0;
                    stds2=std(shufs2);
                    scores2=diffs2(46:55)./stds2;
                    %any score > thresh
                    if any(scores2>thresh)
                        bincounts2=diffs2(46:55);
                        bincounts2(scores2<=thresh)=0;
                        sumdiffs2=(sum(bincounts2(1:5))-sum(bincounts2(6:10)));
                        if sumdiffs2==0
                            AIs2=0;
                        else
                            AIs2=sumdiffs2/(sum(bincounts2(1:5))+sum(bincounts2(6:end)));
                        end
                    else
                        AIs2=0;
                    end

                    onepair=struct();
                    onepair.fileidx=sidx;
                    onepair.su1_label_idx=si;
                    onepair.su2_label_idx=sj;
                    onepair.su1_sel_type=su1;
                    onepair.su2_sel_type=su2;
                    onepair.totalcount=totalCount;
                    onepair.su1_clusterid=su1id;
                    onepair.su2_clusterid=su2id;
                    onepair.hists1=hists1;
                    onepair.hists2=hists2;
                    onepair.shufs1=shufs1;
                    onepair.shufs2=shufs2;
                    onepair.diffs1=diffs1;
                    onepair.diffs2=diffs2;
                    onepair.s1_peak_significant=any(scores1>thresh);
                    onepair.s2_peak_significant=any(scores2>thresh);
                    onepair.AIs1=AIs1;
                    onepair.AIs2=AIs2;
                    onepair.wf_stats_su1=xc_s1.label{si,2}; %cluster_id, FR, err_type,vollay-peak, fwhm
                    onepair.wf_stats_su2=xc_s1.label{sj,2};
                    onepair.wf_su1=xc_s1.label{si,3};
                    onepair.wf_su2=xc_s1.label{sj,3};
                    onepair.prefered_sample_su1=xc_s1.label{si,4};
                    onepair.prefered_sample_su2=xc_s1.label{sj,4};
                    onepair.reg_su1=xc_s1.label{si,5};
                    onepair.reg_su2=xc_s1.label{sj,5};
                    onepair.s1_trials=numel(xc_s1.cfg.trials);
                    onepair.s2_trials=numel(xc_s2.cfg.trials);
                    stats{end+1}=onepair;
                    
                    %TODO: plot
                    if to_plot && abs(onepair.AIs1)>0.9 && onepair.totalcount>2000 && onepair.wf_stats_su2(4)>390 && onepair.wf_stats_su1(4)>390 && nnz(scores1>thresh)>2  && onepair.prefered_sample_su1(currbin+1)>0 && onepair.prefered_sample_su2(currbin+1)>0 && onepair.prefered_sample_su1(currbin+1)==onepair.prefered_sample_su2(currbin+1)
                        fh=figure('Color','w','Position',[100,100,600,800]);
                        subplot(3,2,1);
                        hold on
                        stem(xc_s1.time(50:51)*1000,onepair.hists1(50:51),'Marker','none','LineWidth',15,'Color',[0.8,0.8,0.8])
                        stem(xc_s1.time([46:49,52:55])*1000,onepair.hists1([46:49,52:55]),'Marker','none','LineWidth',15)
                        title('cross-correlogram')
                        xlabel('time (ms)')
                        ylabel('spike-pair count')
%                         print(fh,sprintf('%s_%s_%d_%d_%d.png',su1,su2,sidx,si,sj),'-dpng')
                        subplot(3,2,3);
                        hold on;
                        
                        stem(xc_s1.time(50:51)*1000,onepair.shufs1(50:51),'Marker','none','LineWidth',15,'Color',[0.8,0.8,0.8])
                        stem(xc_s1.time([46:49,52:55])*1000,onepair.shufs1([46:49,52:55]),'Marker','none','LineWidth',15)
%                         stem(xc_s1.time(46:55)*1000,onepair.shufs1(46:55),'Marker','none','LineWidth',15)
                        ssf=smooth(onepair.shufs1);
                        plot(xc_s1.time(46:55)*1000,ssf(46:55),'-r','LineWidth',1.5)
                        title('shift-predictor')
                        xlabel('time (ms)')
                        ylabel('spike-pair count')
                        subplot(3,2,5);
                        hold on
                        sigbin=find((onepair.diffs1(46:55)./stds1)>thresh);
                        insigbin=find((onepair.diffs1(46:55)./stds1)<=thresh);
                        stem(xc_s1.time(insigbin+45)*1000,onepair.diffs1(45+insigbin)./stds1,'Marker','none','LineWidth',15,'Color',[0.8,0.8,0.8]);
                        stem(xc_s1.time(sigbin+45)*1000,onepair.diffs1(45+sigbin)./stds1,'Marker','none','LineWidth',15)
                        yline(thresh,'--r');
                        title('significant peak')
                        xlabel('time (ms)')
                        ylabel('z-score')
                        subplot(3,4,3)
                        plot(onepair.wf_su1)
                        set(gca,'XTick',[31,91],'XTickLabel',[0,2]);
                        xline(31+12,':k');
                        xline(31,':k');
                        xlabel(sprintf('time (ms), v-p=%d',onepair.wf_stats_su1(4)));
                        ylabel('uV');
                        subplot(3,4,4)
                        plot(onepair.wf_su1)
                        xline(31+12,':k');
                        xline(31,':k');
                        set(gca,'XTick',[31,91],'XTickLabel',[0,2])
                        
                        xlabel(sprintf('time (ms), v-p=%d',onepair.wf_stats_su2(4)));
                        ylabel('uV');
                        
                        subplot(3,2,4)
                        hold on
                        text(0,0.6,sums{sidx,2});
                        text(0,0.3,sprintf('%d, %d',onepair.su1_clusterid,onepair.su2_clusterid));
                        print(sprintf('xcorr_showcase_%d_%d_%d.pdf',sidx,onepair.su1_clusterid,onepair.su2_clusterid),'-dpdf','-painters');
                        print(sprintf('xcorr_showcase_%d_%d_%d.png',sidx,onepair.su1_clusterid,onepair.su2_clusterid),'-dpng','-painters');
                        close(fh)
                    end 
                    
                    clear onepair

                end
            end
        end
        if to_save
            disp(sprintf('%s_XCORR_stats_delay_6_%d_%d_2msbin.mat',prefix, bin_range(1),bin_range(2)));
%             keyboard
            save(sprintf('%s_XCORR_stats_delay_6_%d_%d_2msbin.mat',prefix, bin_range(1),bin_range(2)),'stats','bin_range','-v7.3')
        end
    end
% return
end

gen_join_set=false;
if gen_join_set
    join_reg_list=cell(0);
    for sidx=1:length(stats)
        s=stats{sidx};
        join_reg_list{end+1}=s.reg_su1;
        join_reg_list{end+1}=s.reg_su2;
    end
    join_reg_set=unique(join_reg_list);
    keyboard
    save(fullfile('..',sprintf('join_reg_set_%d_%d.mat',bin_range(1),bin_range(2))),'join_reg_set');
return
end

gen_pair_mat=false;
if gen_pair_mat
    if ~exist('join_reg_set','var')
        load(fullfile('..','join_reg_set.mat'));
        reg_set=join_reg_set;
    end
    
    reg_set=reg_set(~strcmp(reg_set,'Unlabeled'));
    reg_set=reg_set(~strcmp(reg_set,'root'));
    
    pair_mat=zeros(length(reg_set),length(reg_set));
    pair_sel_mat=zeros(length(reg_set),length(reg_set));
    for sidx=1:length(stats)
        s=stats{sidx};
        
        if s.s1_trials<20 || s.s2_trials<20 
            continue
        end

        su1reg_idx=find(strcmp(s.reg_su1,reg_set));
        su2reg_idx=find(strcmp(s.reg_su2,reg_set));
        if isempty(su1reg_idx) || isempty(su2reg_idx)
            fprintf('%s, %s\n',s.reg_su1,s.reg_su2);
            %keyboard
            continue
        end
        
        pair_mat(su1reg_idx,su2reg_idx)=pair_mat(su1reg_idx,su2reg_idx)+1;
        pair_mat(su2reg_idx,su1reg_idx)=pair_mat(su2reg_idx,su1reg_idx)+1;

        if s.prefered_sample_su1(currbin+1) && s.prefered_sample_su2(currbin+1) && s.prefered_sample_su1(currbin+1)==s.prefered_sample_su2(currbin+1)  
            pair_sel_mat(su1reg_idx,su2reg_idx)=pair_sel_mat(su1reg_idx,su2reg_idx)+1;
            pair_sel_mat(su2reg_idx,su1reg_idx)=pair_sel_mat(su2reg_idx,su1reg_idx)+1;
        end
    end
% keyboard   
save(sprintf('%s_pair_mat_duo_6s_%d_%d.mat',prefix,bin_range(1),bin_range(2)),'pair_mat','pair_sel_mat');
% return
end


test_selectivity=false;
if test_selectivity
    for sidx=1:length(stats)
        s=stats{sidx};
        if ~any(s.prefered_sample_su1(2:end)) || ~any(s.prefered_sample_su1(2:end))
            keyboard
        end
    end
end



gen_conn_mat=false;
if gen_conn_mat
    conn_mat_all=cell(0);
    if ~exist('join_reg_set','var')
        load(fullfile('..','join_reg_set.mat'));
        reg_set=join_reg_set;
    end
    
    reg_set=reg_set(~strcmp(reg_set,'Unlabeled'));
    reg_set=reg_set(~strcmp(reg_set,'root'));
    
    for bin=currbin
%         load(sprintf('XCORR_stats_delay_6_%d_%d_2msbin.mat',bin,bin+1));
        conn_mat=zeros(length(reg_set),length(reg_set));
        conn_sel_mat=zeros(length(reg_set),length(reg_set));
        tbin=bin;
        for pidx=1:length(stats)
            s=stats{pidx};
            if s.totalcount<1000 || s.s1_trials<20 || s.s2_trials<20 || strcmp(s.reg_su1,'Unlabeled') || strcmp(s.reg_su2,'Unlabeled') || strcmp(s.reg_su1,'root') || strcmp(s.reg_su2,'root')
                continue
            end
            
            su1reg_idx=find(strcmp(s.reg_su1,reg_set));
            su2reg_idx=find(strcmp(s.reg_su2,reg_set));
            if isempty(su1reg_idx) || isempty(su2reg_idx)
                fprintf('%s, %s\n', s.reg_su1, s.reg_su2);
                %keyboard
                continue
            end
            
%             not helpful in selective or non-selective type of datasets
%             if ~any(s.prefered_sample_su1(2:end)) || ~any(s.prefered_sample_su2(2:end))
%                 continue
%             end
            
            if s.prefered_sample_su1(currbin+1) && s.prefered_sample_su2(currbin+1) && s.prefered_sample_su1(currbin+1)==s.prefered_sample_su2(currbin+1)  
                sel_flag=true;
            else
                sel_flag=false;
            end

            if s.s1_peak_significant && s.s2_peak_significant
                if (s.AIs1>0 && s.AIs2>=0) || (s.AIs1>=0 && s.AIs2>0) %2 to 1
                    conn_mat(su1reg_idx,su2reg_idx)=conn_mat(su1reg_idx,su2reg_idx)+1;
                    if sel_flag
                        conn_sel_mat(su1reg_idx,su2reg_idx)=conn_sel_mat(su1reg_idx,su2reg_idx)+1;
                    end
                elseif (s.AIs1<0 && s.AIs2<=0) || (s.AIs1<=0 && s.AIs2<0)
                    conn_mat(su2reg_idx,su1reg_idx)=conn_mat(su2reg_idx,su1reg_idx)+1;
                    if sel_flag
                        conn_sel_mat(su2reg_idx,su1reg_idx)=conn_sel_mat(su2reg_idx,su1reg_idx)+1;
                    end                    
                elseif (s.AIs1<0 && s.AIs2>0) || (s.AIs1>0 && s.AIs2<0)
                    conn_mat(su2reg_idx,su1reg_idx)=conn_mat(su2reg_idx,su1reg_idx)+1;
                    conn_mat(su1reg_idx,su2reg_idx)=conn_mat(su1reg_idx,su2reg_idx)+1;
                    if sel_flag
                        conn_sel_mat(su1reg_idx,su2reg_idx)=conn_sel_mat(su1reg_idx,su2reg_idx)+1;
                        conn_sel_mat(su2reg_idx,su1reg_idx)=conn_sel_mat(su2reg_idx,su1reg_idx)+1;
                    end
                end
            elseif s.s1_peak_significant
                if s.AIs1>0 % su2 to su1
                    conn_mat(su1reg_idx,su2reg_idx)=conn_mat(su1reg_idx,su2reg_idx)+1;
                    if sel_flag
                        conn_sel_mat(su1reg_idx,su2reg_idx)=conn_sel_mat(su1reg_idx,su2reg_idx)+1;
                    end                    
                elseif s.AIs1<0  %=0 is possible
                    conn_mat(su2reg_idx,su1reg_idx)=conn_mat(su2reg_idx,su1reg_idx)+1;
                    if sel_flag
                        conn_sel_mat(su2reg_idx,su1reg_idx)=conn_sel_mat(su2reg_idx,su1reg_idx)+1;
                    end
                end
                
            elseif s.s2_peak_significant
                if s.AIs2>0 % su2 to su1
                    conn_mat(su1reg_idx,su2reg_idx)=conn_mat(su1reg_idx,su2reg_idx)+1;
                    if sel_flag
                        conn_sel_mat(su1reg_idx,su2reg_idx)=conn_sel_mat(su1reg_idx,su2reg_idx)+1;
                    end
                elseif s.AIs2<0  %=0 is possible
                    conn_mat(su2reg_idx,su1reg_idx)=conn_mat(su2reg_idx,su1reg_idx)+1;
                    if sel_flag
                        conn_sel_mat(su2reg_idx,su1reg_idx)=conn_sel_mat(su2reg_idx,su1reg_idx)+1;
                    end
                end
            end
                
        end
    end
disp('check file name')
% keyboard
save(sprintf('%s_conn_mat_duo_6s_%d_%d.mat',prefix,bin_range(1),bin_range(2)),'conn_mat','conn_sel_mat')    
% return
end

gen_ratio_map=false;
if gen_ratio_map
%     load('pair_mat_duo_6s_1_2.mat','pair_mat');
%     load('conn_mat_duo_6s_1_2.mat','conn_mat');
%     pair_mat(pair_mat<=20)=0;
load(fullfile('..','join_reg_set.mat'));
reg_set=join_reg_set;
reg_set=reg_set(~strcmp(reg_set,'Unlabeled'));
reg_set=reg_set(~strcmp(reg_set,'root'));
for bin=currbin
    load(sprintf('conn_mat_duo_6s_%d_%d.mat',bin,bin+1));
    load(sprintf('pair_mat_duo_6s_%d_%d.mat',bin,bin+1));

    pair_mat(pair_mat<10)=0;
    keep=false(length(pair_mat),1);
    for i=1:length(pair_mat)
        if nnz(pair_mat(i,:)>=10)+nnz(pair_mat(:,i)>=10)>=64
            keep(i)=true;
        end
    end
%     keyboard
    pair_mat_k=pair_mat(keep,keep);
    conn_mat_k=conn_mat(keep,keep);
    reg_keep=reg_set(keep);

    ratio_mat=conn_mat_k./pair_mat_k;
    sum_ratio=zeros(length(ratio_mat),1);
    for i=1:length(ratio_mat)
        sum_ratio(i)=nansum(ratio_mat(i,:))+nansum(ratio_mat(:,i));
    end
%     
    sort_mat=ratio_mat;
    sort_mat(isnan(sort_mat))=nanmean(ratio_mat(:));
    
%     T=clusterdata(sort_mat,'criterion','distance','linkage','complete','maxclust',3);
%     T = kmeans(sort_mat,3,'Distance','sqeuclidean','Display','off','MaxIter',100000); 
%     I=T*1000+(1:length(sort_mat))';
    [~,indices]=sort(sum_ratio);
%     [~,indices]=sort(I);
    figure('Color','w','Position',[100,100,450,400])
    h=imagesc(ratio_mat(indices,indices),[0,0.3]);

    set(h,'alphadata',~isnan(ratio_mat(indices,indices))); 
    ax=gca();
    ax.YTick=(1:nnz(keep));
    ax.YTickLabel=reg_keep(indices);
    
    
    ax.XTick=(1:nnz(keep));
    ax.XTickLabel=reg_keep(indices);
    ax.XTickLabelRotation=90;
    
    ax.YDir='normal';
    ax.Color=[0.5,0.5,0.5];
%     colormap('jet')  
    colorbar;

%     print(sprintf('ratio_map_%d_%d.pdf',bin,bin+1),'-dpdf','-painters','-r300')
%     print(sprintf('ratio_map_%d_%d.png',bin,bin+1),'-dpng','-painters','-r300')
end
return
end

return



% join_reg_set=unique([reg_set_all{1};reg_set_all{2};reg_set_all{3};reg_set_all{4};reg_set_all{5};reg_set_all{6}]);

% member_mat=zeros(length(join_reg_set),6);
% for i=1:6
%     member_mat(:,i)=ismember(join_reg_set,reg_set_all{i});
% end

% min_set=join_reg_set(sum(member_mat')==6);
keep=false(length(join_reg_set),1);
for i=1:length(join_reg_set)
%     reg_sum=0;
%     for j=1:6
%         reg_sum=reg_sum+sum([conn_mat_all{j}(i,:),conn_mat_all{j}(:,i)']);
%     end
%     if reg_sum>500
    if sum([conn_mat(i,:),conn_mat(:,i)'])>100
        keep(i)=true;
    end
end

for i=0:6
    bin_range=[i,i+1];
    conn_mat=conn_mat_all{i+1}
    
    ff=figure('Color','w');
    imagesc(conn_mat)
    colormap('jet')
    set(gca(),'XTick',1:length(reg_set),'XTickLabel',reg_set,'XTickLabelRotation',90,...
        'YTick',1:length(reg_set),'YTickLabel',reg_set)
    xlabel('target')
    ylabel('source')
    xlim([-0.5,length(reg_set)+0.5])
    ylim([-0.5,length(reg_set)+0.5])
    colorbar
    print(sprintf('fullmap_%d_%d.png',bin_range(1),bin_range(2)),'-dpng')
    close(ff)
    coremap=conn_mat(keep,keep);
    figure('Color','w','Position',[100,100,450,375])
    imagesc(coremap,[0,500])
    colormap('jet')
    set(gca(),'XTick',1:length(coremap),'XTickLabel',reg_set(keep),'XTickLabelRotation',90,...
        'YTick',1:length(coremap),'YTickLabel',reg_set(keep))
    xlabel('target')
    ylabel('source')
    xlim([0.5,length(coremap)+0.5])
    ylim([0.5,length(coremap)+0.5])
    colorbar
    print(sprintf('error_core_%d_%d.png',bin_range(1),bin_range(2)),'-dpng')
    print(sprintf('error_core_%d_%d.eps',bin_range(1),bin_range(2)),'-depsc')
end

coefmat=zeros(12,12)
all_mat=[conn_mat_all(:);error_conn_mat_all(:)];
for i=1:12
    for j=1:12
        coef=corrcoef(all_mat{i}(keep,keep),all_mat{j}(keep,keep));
        coefmat(i,j)=min(coef(:));
    end
end

figure('Color','w','Position',[100,100,400,300])
imagesc(coefmat(1:end,1:end),[0.5,1])
colormap(jet)
xlabel('bin #')
ylabel('bin #')
colorbar()
print('cross_bin_corr_coef.pdf','-dpdf')


