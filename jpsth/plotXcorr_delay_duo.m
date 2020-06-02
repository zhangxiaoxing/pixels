% assume 'sums' is loaded in workspace. Otherwise load corresponding
% XCORR_delay_bin.mat file first


prepare_stats_file=false;
if prepare_stats_file
    fs=dir('selec_duo_XCORR_de*');
    for bin=1:6
        load(fullfile(fs(bin).folder,fs(bin).name));
        to_plot=false;
        stats=cell(0);
        thresh=norminv(0.995);
        bin_range=[bin,bin+1];
        for sidx=1:size(sums,1)
            xc_s1=sums{sidx,5};
            xshuf_s1=sums{sidx,6};
            xc_s2=sums{sidx,7};
            xshuf_s2=sums{sidx,8};
            sustCount=numel(sums{sidx,3});
            for si=1:(size(xc_s1.xcorr,1)-1)
                su1id=str2double(xc_s1.label{si,1});
                if ismember(su1id,sums{sidx,3})
                    su1='sust';
                else
                    su1='transient';
                end
                for sj=(si+1):size(xc_s1.xcorr,2)
                    su2id=str2double(xc_s1.label{sj,1});
                    if ismember(su2id,sums{sidx,3})
                        su2='sust';
                    else
                        su2='transient';
                    end
                    totalCount=nansum(squeeze(xc_s1.xcorr(si,sj,:)));
                    if totalCount<100
        %                 stats(end+1,:)={sidx,si,sj,su1,su2,totalCount,NaN}; % TODO: match updated data format
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
                            AIs1=sumdiffs1/(sum(bincounts1(1:10))+sum(bincounts1(11:end)));
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
                            AIs2=sumdiffs2/(sum(bincounts2(1:10))+sum(bincounts2(11:end)));
                        end
                    else
                        AIs2=0;
                    end


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

                    stats{end+1}=onepair;

                    %TODO: plot
        %             if to_plot
        %                 fh=figure('Color','w','Position',[100,100,400,300]);
        %                 stem(xc_s1.time*1000,squeeze(xc_s1.xcorr(si,sj,:)),'Marker','none','LineWidth',10)
        %                 title(sprintf('%s-%s,%d,%d,%d',su1,su2,sidx,si,sj))
        %                 xlabel('time (ms)')
        %                 ylabel('spike-pair count')
        %                 print(fh,sprintf('%s_%s_%d_%d_%d.png',su1,su2,sidx,si,sj),'-dpng')
        % 
        % %                 pause
        %                 close(fh)
        %             end
                end
            end
        end
        save(sprintf('error_XCORR_stats_delay_6_%d_%d_2msbin.mat',bin_range(1),bin_range(2)),'stats','bin_range')
    end
end

gen_conn_mat=true;
if gen_conn_mat
    conn_mat_all=cell(0);
    if ~exist('join_reg_set','var')
        load('join_reg_set.mat')
    end
    for bin=1
        others=0;

%         load(sprintf('XCORR_stats_delay_6_%d_%d_2msbin.mat',bin,bin+1));
        conn_dirs=cell(0,5);
        tbin=bin;
        for pidx=1:length(stats)
            onepair=stats{pidx};
            if onepair.totalcount<100
                others=others+1;
                continue
            end
%             if onepair.prefered_sample_su1(tbin+1)>0 && onepair.prefered_sample_su2(tbin+1)>0 &&...
%                 onepair.s1_peak_significant && ~strcmp(onepair.reg_su1,onepair.reg_su2) && ...
            if onepair.s1_peak_significant && any(onepair.prefered_sample_su1(2:end)) && any(onepair.prefered_sample_su2(2:end)) &&...
                    ~strcmp(onepair.reg_su1,'Unlabeled') && ~strcmp(onepair.reg_su2,'Unlabeled')
                if onepair.AIs1>0
                    conn_dirs(end+1,:)={onepair.reg_su1,onepair.reg_su2,onepair.AIs1,onepair.totalcount,pidx};
                elseif onepair.AIs1<0  %=0 is possible
                    conn_dirs(end+1,:)={onepair.reg_su2,onepair.reg_su1,-onepair.AIs1,onepair.totalcount,pidx};
                end
            else
                others=others+1;
            end
        end


    %     reg_set=unique([conn_dirs(:,1);conn_dirs(:,2)]);
        reg_set=join_reg_set;
        conn_mat=zeros(length(reg_set),length(reg_set));

        for i=1:length(conn_dirs)
            conn_from_idx=find(strcmp(conn_dirs{i,1},reg_set));
            conn_to_idx=find(strcmp(conn_dirs{i,2},reg_set));
            conn_mat(conn_from_idx,conn_to_idx)=conn_mat(conn_from_idx,conn_to_idx)+1;
        end
        conn_mat_all{end+1}=conn_mat;
    %     reg_set_all{end+1}=reg_set;
    end
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