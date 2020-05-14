% assume 'sums' is loaded in workspace. Otherwise load corresponding
% XCORR_delay_bin.mat file first
to_plot=false;
stats=cell(0);
thresh=norminv(0.995);
bin_range=[1,2];
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
save(sprintf('XCORR_stats_delay_6_%d_%d_2msbin.mat',bin_range(1),bin_range(2)),'stats','bin_range')
% autoCount=sum(strcmp(stats(:,5),'auto-corr'));
% stCount=sum(strcmp(stats(:,4),'sust') & strcmp(stats(:,5),'transient'));
% stCount50=sum(strcmp(stats(:,4),'sust') & strcmp(stats(:,5),'transient') & cell2mat(stats(:,6))>=50);
% % ssCount=sum(strcmp(stats(:,4),'sust') & strcmp(stats(:,5),'sust'));
% ttCount=sum(strcmp(stats(:,4),'transient') & strcmp(stats(:,5),'transient'))/2;
% ttCount50=sum(strcmp(stats(:,4),'transient') & strcmp(stats(:,5),'transient')& cell2mat(stats(:,6))>=50)/2;

% conn_dirs=cell(0);
% for pidx=1:length(stats)
%     onepair=stats{pidx};
%     if onepair.s1_peak_significant && ~strcmp(onepair.reg_su1,onepair.reg_su2) && ...
%             ~strcmp(onepair.reg_su1,'Unlabeled') && ~strcmp(onepair.reg_su2,'Unlabeled')
%         if onepair.AIs1>0
%             conn_dirs{end+1}=sprintf('%s->%s',onepair.reg_su1,onepair.reg_su2);
%         elseif onepair.AIs1<0  %=0 is possible
%             conn_dirs{end+1}=sprintf('%s->%s',onepair.reg_su2,onepair.reg_su1);
%         end
%     end
%     
% end
% 
% conn_sum=[];
% cidx=[];
% conn_set=unique(conn_dirs);
% for cidx=1:length(conn_set)
%     conn_sum(cidx)=nnz(strcmp(conn_dirs,conn_set{cidx}));
% end
% 
% [~,conn_idx]=sort(conn_sum,'descend');
% 
% figure()
% hold on
% for i=1:100
%     bar(i,conn_sum(conn_idx(i)))
% end
% set(gca(),'XTick',1:100,'XTickLabel',conn_set(conn_idx(1:100)),'XTickLabelRotation',90)



conn_dirs=cell(0,4);
tbin=1;
for pidx=1:length(stats)
    onepair=stats{pidx};
    if onepair.prefered_sample_su1(tbin+1)>0 && onepair.prefered_sample_su2(tbin+1)>0 &&...
        onepair.s1_peak_significant && ~strcmp(onepair.reg_su1,onepair.reg_su2) && ...
            ~strcmp(onepair.reg_su1,'Unlabeled') && ~strcmp(onepair.reg_su2,'Unlabeled')
        if onepair.AIs1>0
            conn_dirs(end+1,:)={onepair.reg_su1,onepair.reg_su2,onepair.AIs1,onepair.totalcount};
        elseif onepair.AIs1<0  %=0 is possible
            conn_dirs(end+1,:)={onepair.reg_su2,onepair.reg_su1,-onepair.AIs1,onepair.totalcount};
        end
    end
end


reg_set=unique([conn_dirs(:,1);conn_dirs(:,2)]);
conn_mat=zeros(length(reg_set),length(reg_set));

for i=1:length(conn_dirs)
    conn_from_idx=find(strcmp(conn_dirs{i,1},reg_set));
    conn_to_idx=find(strcmp(conn_dirs{i,2},reg_set));
    conn_mat(conn_from_idx,conn_to_idx)=conn_mat(conn_from_idx,conn_to_idx)+1;
end

figure('Color','w')
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

coremap=conn_mat;
keep=[];
for i=1:length(reg_set)
    if sum([conn_mat(i,:),conn_mat(:,i)'])>100
        keep(end+1)=i;
    end
end
coremap=conn_mat(keep,keep)

figure('Color','w')
imagesc(coremap,[0,100])
colormap('jet')
set(gca(),'XTick',1:length(coremap),'XTickLabel',reg_set(keep),'XTickLabelRotation',90,...
    'YTick',1:length(coremap),'YTickLabel',reg_set(keep))
xlabel('target')
ylabel('source')
xlim([0.5,length(coremap)+0.5])
ylim([0.5,length(coremap)+0.5])
colorbar
print(sprintf('core_%d_%d.png',bin_range(1),bin_range(2)),'-dpng')


