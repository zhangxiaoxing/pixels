close all
load('pv_sst.mat');
cellSums=cell(0,2);
for i=4:height(s2)
    if ismember(combineSubRegion(s2{i,1}),cellSums(:,1))
%         cidx=find(strcmp(combineSubRegion(s2{i,1}),cellSums(:,1)));
%         cellSums{cidx,2}=cellSums{cidx,2}+round(table2array(s2(i,4:end)));
    else
        cellSums(end+1,:)={combineSubRegion(s2{i,1}),round(table2array(s2(i,4:end)))};
    end
end
load('reg_keep.mat');
load('su_region_recip.mat')
load('io_sel.mat')

sus_trans=h5read('../transient_6.hdf5','/sus_trans');%4+2+7
reg_all=deblank(h5read('../transient_6.hdf5','/reg'));


wing=diff(ioselstats{1}(:,[3 7]),1,2);
per_reg_recip_idx=nanmean(recip_count./input_to_count,2);
per_reg_shuf=nanmean(recip_count_shuf./input_to_count_shuf,2);


per_reg_recip_idx=[per_reg_recip_idx,(1:140)',wing];
recip_idx=per_reg_recip_idx(sum(input_to_count,2)>120,:);
recip_shuf=squeeze(per_reg_shuf(sum(input_to_count,2)>120,:,:));

mmshuf=nanmean(recip_shuf,2);
stdshuf=nanstd(recip_shuf,0,2);
z=(recip_idx(:,1)-mmshuf)./stdshuf;
%%figure 3
fh=figure('Color','w','Position',[100,100,235,235]);
hold on;
sh=scatter(mmshuf(z>1.96),recip_idx(z>1.96,1),'MarkerFaceAlpha',0.5,'MarkerFaceColor','r','MarkerEdgeColor','none');
sh=scatter(mmshuf(z<-1.96),recip_idx(z<-1.96,1),'MarkerFaceAlpha',0.5,'MarkerFaceColor','b','MarkerEdgeColor','none');
sh=scatter(mmshuf(z>=-1.96 & z<=1.96),recip_idx(z>=-1.96 & z<=1.96,1),'MarkerFaceAlpha',0.5,'MarkerFaceColor','k','MarkerEdgeColor','none');

plot([0,1],[0,1],'k:')
xlim([0,1]);
ylim([0,1]);
set(gca,'YTick',0:0.5:1)
xlabel('shuffled reciprocal fraction')
ylabel('recorded reciprocal fraction')
exportgraphics(fh,'recip_neuron_region.pdf')

corr_list=[];
su_count_list=[];
reg_list=cell(0);
for i=1:length(recip_idx)
    reg=reg_set{recip_idx(i,2)};
%     if strcmp(reg,'GPi')
%         keyboard
%     end
%     disp(reg);
    pvsstIdx=find(strcmp(reg,cellSums(:,1)));
    
    all_su_count=nnz(strcmp(reg_all,reg));
    mem_su_count=nnz(strcmp(reg_all,reg) & (sus_trans(:,1)|sus_trans(:,2)));
    
    if all_su_count==0
        all_su_count=nnz(startsWith(reg_all,reg));
        mem_su_count=nnz(startsWith(reg_all,reg) & (sus_trans(:,1)|sus_trans(:,2)));
    end
    
    if ~isempty(pvsstIdx)
        %recip_index,pv/pv+sst,wing
        corr_list=[corr_list;recip_idx(i,1),cellSums{pvsstIdx,2}(1)/sum(cellSums{pvsstIdx,2}([1 3])),recip_idx(i,3)];
        reg_list=[reg_list;reg];
        su_count_list=[su_count_list;all_su_count,mem_su_count];
    else
        disp(reg);
    end
end

fh=figure('Color','w','Position',[100,100,400,300]);
s=scatter(corr_list(:,1),corr_list(:,2),'MarkerFaceAlpha',0.5,'MarkerFaceColor','k','MarkerEdgeColor','none');
xlabel('reciprocal fraction')
ylabel('PV/(PV+SST)')
[r,p]=corr(corr_list(:,1),corr_list(:,2));
legend(s,sprintf('r=%.3f,p=%.3f',r,p));
keyboard
exportgraphics(fh,'recip_frac_PVSST.pdf')
% print(fh,'-dpng','recip_frac_PVSST.png','-r300')

fh=figure('Color','w','Position',[100,100,400,300]);
s=scatter(corr_list(:,1),corr_list(:,3),'MarkerFaceAlpha',0.5,'MarkerFaceColor','k','MarkerEdgeColor','none');
[r,p]=corr(corr_list(:,1),corr_list(:,3));
xlabel('reciprocal fraction');
ylabel('WING');
legend(s,sprintf('r=%.3f,p=%.3f',r,p));
if ~verLessThan('matlab','9.8')
    exportgraphics(fh,'recip_frac_WING.pdf')
end


sus_trans=h5read('../transient_6.hdf5','/sus_trans');%4+2+7
reg_all=h5read('../transient_6.hdf5','/reg');

fh=figure('Color','w','Position',[100,100,240,240]);
% subplot(1,2,1)
% s=scatter(corr_list(:,1),su_count_list(:,1),'MarkerFaceAlpha',0.5,'MarkerFaceColor','k','MarkerEdgeColor','none');
% [r,p]=corr(corr_list(:,1),su_count_list(:,1));
% xlabel('reciprocal fraction')
% ylabel('unit count')
% legend(s,sprintf('r=%.3f,p=%.3f',r,p));

% subplot(1,2,2)
s=scatter(su_count_list(:,1),corr_list(:,2),'MarkerFaceAlpha',0.5,'MarkerFaceColor','k','MarkerEdgeColor','none');
[r,p]=corr(corr_list(:,1),su_count_list(:,2));
ylabel('reciprocal fraction')
xlabel('memory unit count')
legend(s,sprintf('r=%.3f,p=%.3f',r,p));
set(gca(),'YScale','linear','XScale','log')
exportgraphics(fh,'su_count_recip_frac.pdf');


figure();
s=scatter(su_count_list(:,1),su_count_list(:,2));
[r,p]=corr(su_count_list(:,1),su_count_list(:,2));
legend(s,sprintf('r=%.3f,p=%.3f',r,p));
title('selective fraction')



function out=combineSubRegion(r)
    if ~isempty(regexp(r,'ACB.*','once'))
        out='ACB';
        return
    end

    if ~isempty(regexp(r,'CA[13]','once'))
        [t,~]=regexp(r,'(CA[13])','tokens','match');
        out=t{1}{1};
        return
    end
    if ~isempty(regexp(r,'([A-Za-z]+)[1-6/-]{0,3}[a-z]{0,3}','once'))
        [t,~] = regexp(r,'([A-Za-z]+)[1-6/-]{0,3}[a-z]{0,3}','tokens','match');
        out=t{1}{1};
        return
    end
end