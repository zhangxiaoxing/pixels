%  A peak at a negative lag for stat.xcorr(chan1,chan2,:) means that chan1 is leading
%  chan2. Thus, a negative lag represents a spike in the second dimension of
%  stat.xcorr before the channel in the third dimension of stat.stat.


% assume 'sums' is loaded in workspace. Otherwise load corresponding
% XCORR_delay_bin.mat file first

fs=dir('selec*delay_6_1_2*');
sums=cell(0);
for i=1:size(fs,1)
    fprintf('%d of %d\n',i,size(fs,1));
    fstr=load(fullfile(fs(i).folder,fs(i).name));
    sums(end+1,:)=fstr.sums(1,:);
end



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
%             continue;
        elseif ismember(su1id,sums{sidx,4})
            su1='transient';
%             continue;
        else
            su1='non-selective';
            continue;
        end
        for sj=(si+1):size(xc_s1.xcorr,2)
            su2id=str2double(xc_s1.label{sj,1});
            if ismember(su2id,sums{sidx,3})
                su2='sust';
%                 continue
            elseif ismember(su1id,sums{sidx,4})
                su2='transient';
%                 continue
            else
                su2='non-selective';
                continue;
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
                    AIs1=sumdiffs1/(sum(bincounts1(1:5))+sum(bincounts1(6:10)));
                end
            else
                AIs1=0;
            end
            
            %  A peak at a negative lag for stat.xcorr(chan1,chan2,:) means that chan1 is leading
            %  chan2. Thus, a negative lag represents a spike in the second dimension of
            %  stat.xcorr before the channel in the third dimension of stat.stat.

            
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
                    AIs2=sumdiffs2/(sum(bincounts2(1:5))+sum(bincounts2(6:10)));
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
            onepair.s1_trialNum=numel(sums{sidx,5}.cfg.trials);
            onepair.s2_trialNum=numel(sums{sidx,7}.cfg.trials);
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
save(sprintf('non-selective_XCORR_stats_delay_6_%d_%d_2msbin.mat',bin_range(1),bin_range(2)),'stats','bin_range')
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



c_conn_dirs=cell(0,4);%{to, from, AI,count}
tbin=1;
for pidx=1:length(cfstr.stats)
    onepair=cfstr.stats{pidx};
    if onepair.prefered_sample_su1(tbin+1)>0 && onepair.prefered_sample_su2(tbin+1)>0 &&...
        onepair.s1_peak_significant && ~strcmp(onepair.reg_su1,onepair.reg_su2) && ...
            ~strcmp(onepair.reg_su1,'Unlabeled') && ~strcmp(onepair.reg_su2,'Unlabeled')
        if onepair.AIs1>0  %<0 is 1 leading 2 >0 is 2 leading 1, =0 is possible
            c_conn_dirs(end+1,:)={onepair.reg_su1,onepair.reg_su2,onepair.AIs1,onepair.totalcount};
        elseif onepair.AIs1<0  %
            c_conn_dirs(end+1,:)={onepair.reg_su2,onepair.reg_su1,-onepair.AIs1,onepair.totalcount};
        end
    end
end

e_conn_dirs=cell(0,4); %{to, from, AI,count}
tbin=1;
for pidx=1:length(efstr.stats)
    onepair=efstr.stats{pidx};
    if onepair.prefered_sample_su1(tbin+1)>0 && onepair.prefered_sample_su2(tbin+1)>0 &&...
        onepair.s1_peak_significant && ~strcmp(onepair.reg_su1,onepair.reg_su2) && ...
            ~strcmp(onepair.reg_su1,'Unlabeled') && ~strcmp(onepair.reg_su2,'Unlabeled')
        if onepair.AIs1>0 %<0 is 1 leading 2 >0 is 2 leading 1, =0 is possible
            e_conn_dirs(end+1,:)={onepair.reg_su1,onepair.reg_su2,onepair.AIs1,onepair.totalcount};
        elseif onepair.AIs1<0  
            e_conn_dirs(end+1,:)={onepair.reg_su2,onepair.reg_su1,-onepair.AIs1,onepair.totalcount};
        end
    end
end


reg_set=unique([c_conn_dirs(:,1);c_conn_dirs(:,2)]);
c_conn_mat=zeros(length(reg_set),length(reg_set));
e_conn_mat=zeros(length(reg_set),length(reg_set));

for i=1:length(c_conn_dirs)
    conn_to_idx=find(strcmp(c_conn_dirs{i,1},reg_set)); 
    conn_from_idx=find(strcmp(c_conn_dirs{i,2},reg_set));
    c_conn_mat(conn_to_idx,conn_from_idx)=c_conn_mat(conn_to_idx,conn_from_idx)+1;
end

c_conn_list=cell(0);
for i=1:length(c_conn_mat)
    for j=1:length(c_conn_mat)
        c_conn_list(end+1,:)={reg_set{i},reg_set{j},c_conn_mat(i,j)};
    end
end
keep=[];
for i=1:length(reg_set)
    if sum([c_conn_mat(i,:),c_conn_mat(:,i)'])>100
        keep(end+1)=i;
    end
end


for i=1:length(e_conn_dirs)
    conn_to_idx=find(strcmp(c_conn_dirs{i,1},reg_set));
    conn_from_idx=find(strcmp(c_conn_dirs{i,2},reg_set));
    e_conn_mat(conn_to_idx,conn_from_idx)=e_conn_mat(conn_to_idx,conn_from_idx)+1;
end

%%%% connection graph
G = digraph(c_conn_mat(keep,keep),reg_set(keep),'omitselfloops');
GE = digraph(e_conn_mat(keep,keep),reg_set(keep),'omitselfloops');


fh=figure('Color','w','Position',[100,100,900,600]);
Gh=plot(G, 'LineWidth', G.Edges.Weight/20,'EdgeAlpha',0.5,'ArrowSize',5,'XData',Gref.XData,'YData',Gref.YData);

set(fh,'visible','off')
set(fh,'PaperSize',[15,10])
print(fh,sprintf('connection_layer_%d_%d.pdf',bin_range(1),bin_range(2)),'-painters','-dpdf','-r300')
set(fh, 'PaperPosition', [0 0 12 12])
print(fh,sprintf('connection_layer_%d_%d.png',bin_range(1),bin_range(2)),'-dpng','-r300')



fh=figure('Color','w');
plot(G, 'LineWidth', G.Edges.Weight/20,'EdgeAlpha',0.5,'ArrowSize',6,'Layout','force')
set(fh,'visible','off')
set(fh,'PaperSize',[15,10])
print(fh,sprintf('connection_force_%d_%d.pdf',bin_range(1),bin_range(2)),'-dpdf','-r300')
set(fh, 'PaperPosition', [0 0 12 12])
print(fh,sprintf('connection_force_%d_%d.png',bin_range(1),bin_range(2)),'-dpng','-r300')


fh=figure('Color','w','Position',[100,100,900,600]);
GEh=plot(GE, 'LineWidth', GE.Edges.Weight/20,'EdgeAlpha',0.5,'ArrowSize',5,'XData',Gh.XData,'YData',Gh.YData);
title('error trials');
set(fh,'visible','off')
set(fh,'PaperSize',[15,10])
print(fh,sprintf('connection_error_layer_%d_%d.eps',bin_range(1),bin_range(2)),'-depsc','-r300')
set(fh, 'PaperPosition', [0 0 12 12])
print(fh,sprintf('connection_error_layer_%d_%d.png',bin_range(1),bin_range(2)),'-dpng','-r300')



fh=figure('Color','w');
plot(GE, 'LineWidth', GE.Edges.Weight/20,'EdgeAlpha',0.25,'ArrowSize',10,'Layout','force')
set(fh,'visible','off')
set(fh,'PaperSize',[12,12])
print(fh,sprintf('connection_force_%d_%d.pdf',bin_range(1),bin_range(2)),'-dpdf','-r300')
set(fh, 'PaperPosition', [0 0 12 12])
print(fh,sprintf('connection_force_%d_%d.png',bin_range(1),bin_range(2)),'-dpng','-r300')

%%

fh=figure('Color','w');
% subplot(1,2,1)
imagesc(c_conn_mat)
colormap('jet')
set(gca(),'XTick',1:length(reg_set),'XTickLabel',reg_set,'XTickLabelRotation',90,...
    'YTick',1:length(reg_set),'YTickLabel',reg_set)
xlabel('source')
ylabel('target')
xlim([0.5,length(reg_set)+0.5])
ylim([0.5,length(reg_set)+0.5])
title('correct trials');

subplot(1,2,2)
imagesc(e_conn_mat)
colormap('jet')
set(gca(),'XTick',1:length(reg_set),'XTickLabel',reg_set,'XTickLabelRotation',90,...
    'YTick',1:length(reg_set),'YTickLabel',reg_set)
ylabel('target')
xlabel('source')
xlim([-0.5,length(reg_set)+0.5])
ylim([-0.5,length(reg_set)+0.5])
title('error trials');
colorbar
set(fh,'visible','off')
set(fh,'PaperSize',[21,20])
print(fh,sprintf('all_reg_non_sel_%d_%d.pdf',bin_range(1),bin_range(2)),'-dpdf','-r300')
set(fh, 'PaperPosition', [0 0 21 20])
print(fh,sprintf('all_reg_non_sel_%d_%d.png',bin_range(1),bin_range(2)),'-dpng','-r300')



keep=[];
for i=1:length(reg_set)
    if sum([c_conn_mat(i,:),c_conn_mat(:,i)'])>100
        keep(end+1)=i;
    end
end
c_coremap=c_conn_mat(keep,keep);
e_coremap=e_conn_mat(keep,keep);

diff_coremap=(c_coremap-e_coremap)./(c_coremap+e_coremap);
diff_coremap(c_coremap==e_coremap)=0;

fh=figure('Color','w');
subplot(1,3,1);
imagesc(c_coremap,[0,100])
colormap('jet')
set(gca(),'XTick',1:length(c_coremap),'XTickLabel',reg_set(keep),'XTickLabelRotation',90,...
    'YTick',1:length(c_coremap),'YTickLabel',reg_set(keep))
ylabel('target')
xlabel('source')
xlim([0.5,length(c_coremap)+0.5])
ylim([0.5,length(c_coremap)+0.5])
title('correct trials')
colorbar

subplot(1,3,2);
imagesc(e_coremap,[0,100])
colormap('jet')
set(gca(),'XTick',1:length(e_coremap),'XTickLabel',reg_set(keep),'XTickLabelRotation',90,...
    'YTick',1:length(e_coremap),'YTickLabel',reg_set(keep))
ylabel('target')
xlabel('source')
xlim([0.5,length(e_coremap)+0.5])
ylim([0.5,length(e_coremap)+0.5])
title('error trials')
colorbar

subplot(1,3,3);
imagesc(diff_coremap,[-1,1])
colormap('jet')
set(gca(),'XTick',1:length(e_coremap),'XTickLabel',reg_set(keep),'XTickLabelRotation',90,...
    'YTick',1:length(e_coremap),'YTickLabel',reg_set(keep))
ylabel('target')
xlabel('source')
xlim([0.5,length(e_coremap)+0.5])
ylim([0.5,length(e_coremap)+0.5])
title('correct-error/correct+error')
colorbar

set(fh,'visible','off')
set(fh,'PaperSize',[4.5,4])
print(fh,sprintf('core_all_su_%d_%d.pdf',bin_range(1),bin_range(2)),'-dpdf','-r300')
set(fh, 'PaperPosition', [0 0 4.5 4])
print(fh,sprintf('corr_all_su_%d_%d.png',bin_range(1),bin_range(2)),'-dpng','-r300')


