%% TODO merge with ring_span_plot.m
%% TODO converge to new data structure

keyboard()
%% within
figure('Color','w','Position',[32,32,750,235])
% plot_one(true,within_3,within_3_shuf,1);
plot_one(true,ring_meta.congru.within_4,ring_meta.congru.within_4_shuf,2);
% plot_one(true,within_5,within_5_shuf,3);
exportgraphics(gcf,'loops_within_region_norm.pdf')

figure('Color','w','Position',[32,32,750,235])
plot_one(true,ring_meta.nonmem.within_4,ring_meta.nonmem.within_4_shuf,2);
ylim([-10,62])
exportgraphics(gcf,'loops_within_region_norm_nonmem.pdf')



%% cross
figure('Color','w','Position',[32,32,750,235])
% plot_one(false,cross_3,cross_3_shuf,1);
plot_one(false,ring_meta.congru.cross_4,ring_meta.congru.cross_4_shuf,2);
% [c5r,c5n]=plot_one(false,cross_5,cross_5_shuf,3);
exportgraphics(gcf,'loops_cross_region_norm.pdf')

figure('Color','w','Position',[32,32,750,235])
% plot_one(false,ring_meta.nonmem.cross_3,ring_meta.nonmem.cross_3_shuf,1);
plot_one(false,ring_meta.nonmem.cross_4,ring_meta.nonmem.cross_4_shuf,2,'ylim',[0,50]);
% plot_one(false,ring_meta.nonmem.cross_5,ring_meta.nonmem.cross_5_shuf,3);
exportgraphics(gcf,'loops_cross_region_norm_nonmem.pdf')

%% dependency
function [ureg,normw]=norm_within(rreal,shuf)
wshuf=cellfun(@(x) x.reg,shuf,'UniformOutput',false);
ureg=unique([vertcat(wshuf{:});rreal.reg]);
[B,BG]=groupcounts(rreal.reg);
normw=zeros(numel(ureg),3);
for ri=1:numel(ureg)
    if ismember(ureg(ri),BG)
        normw(ri,1)=B(strcmp(BG,ureg(ri)));
    end
    shufvec=cellfun(@(x) nnz(strcmp(x,ureg(ri))),wshuf);
    normw(ri,2)=mean(shufvec);
    normw(ri,3)=std(shufvec);
end
normw(:,4)=diff(normw(:,[2,1]),1,2)./normw(:,3);
end

function [ureg,normw]=norm_cross(rreal,shuf)
wshuf=cellfun(@(x) x.reg,shuf,'UniformOutput',false);
ureg=unique(cellfun(@(x) x,[vertcat(wshuf{:});rreal.reg]));
[B,BG]=groupcounts(reshape(cellfun(@(x) x,rreal.reg),[],1));
normw=zeros(numel(ureg),3);
for ri=1:numel(ureg)
    if ismember(ureg(ri),BG)
        normw(ri,1)=B(strcmp(BG,ureg(ri)));
    end
    shufvec=cellfun(@(x) nnz(strcmp(x,ureg(ri))),cellfun(@(x) cellfun(@(y) y,x),wshuf,'UniformOutput',false));
    normw(ri,2)=mean(shufvec);
    normw(ri,3)=std(shufvec);
end
normw(:,4)=diff(normw(:,[2,1]),1,2)./normw(:,3);
end



function [ureg,normw]=plot_one(within,ddata,sshuf,spidx,opt)
arguments
    within
    ddata
    sshuf
    spidx
    opt.ylim (1,:) double = []
end
subplot(1,3,spidx)
hold on
if within
    [ureg,normw]=norm_within(ddata,sshuf);
else
    [ureg,normw]=norm_cross(ddata,sshuf);
end
[~,~,ratiomap]=ref.get_pv_sst();
ratios=cellfun(@(x) ratiomap(x),ureg);
[ratios,sidx]=sort(ratios);
ureg=ureg(sidx);
normw=normw(sidx,:);

bh=bar(normw(:,4),'FaceColor','w','EdgeColor','k');
for ri=1:size(normw,1)
    if normw(ri,4)>6
        text(ri,5,sprintf('%.1f',normw(ri,4)),'HorizontalAlignment','center','VerticalAlignment','middle','Rotation',90);
    end
end
set(gca(),'XTick',1:numel(ureg),'XTickLabel',ureg,'XTickLabelRotation',90)
ylabel('Normalized number of within-region loops')
if within
    ylim([-2,20])
    ylabel('Normalized number of within-region loops')
else
    ylim([-2,60])
    ylabel('Normalized number of cross-region loops')
end
if ~isempty(opt.ylim)
    ylim(opt.ylim)
end
title(sprintf('%d-neuron loops',spidx+2));
end

