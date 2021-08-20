keyboard()
[cross_3,within_3]=bz.rings.rings_span('ring_size',3,'memtype','congru');
[cross_4,within_4]=bz.rings.rings_span('ring_size',4,'memtype','congru');
[cross_5,within_5]=bz.rings.rings_span('ring_size',5,'memtype','congru');

[cross_3_shuf,cross_4_shuf,cross_5_shuf,within_3_shuf,within_4_shuf,within_5_shuf]=deal(cell(100,1));
for ri=1:100
    disp(ri);
    [cross_3_shuf{ri},within_3_shuf{ri}]=bz.rings.rings_span('ring_size',3,'memtype','congru','shufid',ri);
    [cross_4_shuf{ri},within_4_shuf{ri}]=bz.rings.rings_span('ring_size',4,'memtype','congru','shufid',ri);
    [cross_5_shuf{ri},within_5_shuf{ri}]=bz.rings.rings_span('ring_size',5,'memtype','congru','shufid',ri);
end

%% within
figure('Color','w','Position',[32,32,750,235])
% plot_one(true,within_3,within_3_shuf,1);
plot_one(true,within_4,within_4_shuf,2);
% plot_one(true,within_5,within_5_shuf,3);
% exportgraphics(gcf,'loops_within_region_norm.pdf')

%% cross
figure('Color','w','Position',[32,32,750,235])
% plot_one(false,cross_3,cross_3_shuf,1);
plot_one(false,cross_4,cross_4_shuf,2);
% [c5r,c5n]=plot_one(false,cross_5,cross_5_shuf,3);
% exportgraphics(gcf,'loops_cross_region_norm.pdf')

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



function [ureg,normw]=plot_one(within,ddata,sshuf,spidx)
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

bh=bar(normw(:,4));
for ri=1:size(normw,1)
    if normw(ri,4)>6
        text(ri,5,sprintf('%.1f',normw(ri,4)),'HorizontalAlignment','center','VerticalAlignment','middle','Rotation',90);
    end
end
set(gca(),'XTick',1:numel(ureg),'XTickLabel',ureg,'XTickLabelRotation',90)
ylabel('Normalized number of within-region loops')
if within
    ylim([-4,4])
    ylabel('Normalized number of within-region loops')
else
    ylim([-1,8])
    ylabel('Normalized number of cross-region loops')
end
title(sprintf('%d-neuron loops',spidx+2));
end

