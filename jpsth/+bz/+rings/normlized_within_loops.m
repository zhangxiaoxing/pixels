%% TODO merge with ring_span_plot.m
%% TODO converge to new data structure
if ~exist('ring_meta','var')
    ring_meta=bz.rings.get_ring_meta('loadfile',true);
end

keyboard()
%% within
plot_two(true,ring_meta.congru.within_4,...
    ring_meta.nonmem.within_4,...
    ring_meta.congru.within_4_shuf,...
    ring_meta.nonmem.within_4_shuf,'ylim',[-10,65]);
exportgraphics(gcf,'loops_within_region_norm.pdf')



%% cross
plot_two(false,ring_meta.congru.cross_4,...
    ring_meta.nonmem.cross_4,...
    ring_meta.congru.cross_4_shuf,...
    ring_meta.nonmem.cross_4_shuf,'ylim',[-20,65]);
exportgraphics(gcf,'loops_cross_region_norm.pdf')


plot_two(false,ring_meta.congru.cross_3,...
    ring_meta.nonmem.cross_3,...
    ring_meta.congru.cross_3_shuf,...
    ring_meta.nonmem.cross_3_shuf,'ylim',[-20,65]);

% plot_one(false,ring_meta.congru.cross_4,ring_meta.congru.cross_4_shuf,2);
% exportgraphics(gcf,'loops_cross_region_norm.pdf')
% 
% plot_one(false,ring_meta.nonmem.cross_4,ring_meta.nonmem.cross_4_shuf,2,'ylim',[0,50]);
% exportgraphics(gcf,'loops_cross_region_norm_nonmem.pdf')

%% dependency
function [ureg,normw]=norm_within(rreal,shuf)
wshuf=cellfun(@(x) x.reg,shuf,'UniformOutput',false);
ureg=unique([vertcat(wshuf{:});rreal.reg]);
[B,BG]=groupcounts(rreal.reg);
normw=zeros(numel(ureg),6);
for ri=1:numel(ureg)
    if ismember(ureg(ri),BG)
        normw(ri,1)=B(strcmp(BG,ureg(ri)));
    end
    shufvec=cellfun(@(x) nnz(strcmp(x,ureg(ri))),wshuf);
    ci=bootci(500,@(x) mean(x),shufvec);
    normw(ri,2)=mean(shufvec);
    normw(ri,3:4)=reshape(ci,1,[]);
    normw(ri,5)=std(shufvec);
end
normw(:,6)=diff(normw(:,[2,1]),1,2)./normw(:,5);
normw(:,7)=diff(normw(:,[3,1]),1,2)./normw(:,5);
normw(:,8)=diff(normw(:,[4,1]),1,2)./normw(:,5);
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
    ci=bootci(500,@(x) mean(x),shufvec);
    normw(ri,2)=mean(shufvec);
    normw(ri,3:4)=reshape(ci,1,[]);
    normw(ri,5)=std(shufvec);
end
normw(:,6)=diff(normw(:,[2,1]),1,2)./normw(:,5);
normw(:,7)=diff(normw(:,[3,1]),1,2)./normw(:,5);
normw(:,8)=diff(normw(:,[4,1]),1,2)./normw(:,5);
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



function [ureg,normw]=plot_two(within,datac,datan,shufc,shufn,opt)
arguments
    within (1,1) logical
    datac
    datan
    shufc
    shufn
    opt.ylim (1,:) double = []
end

if within
    [uregc,normwc]=norm_within(datac,shufc); %reg,{real-mean, shuf-mean,ci(1),ci(2),shuf-std, z-score,zci(1),zci(2)}
    [uregn,normwn]=norm_within(datan,shufn);
else
    [uregc,normwc]=norm_cross(datac,shufc);
    [uregn,normwn]=norm_cross(datan,shufn);
end
ureg=unique([uregc;uregn]);
[~,~,ratiomap]=ref.get_pv_sst();
ratios=cellfun(@(x) ratiomap(x),unique(ureg));
[ratios,sidx]=sort(ratios);
ureg=ureg(sidx);

bardata=[];
for ri=1:numel(ureg)
    %congru and nonmem
    if any(strcmp(uregc,ureg(ri)))
        bardata=[bardata;normwc(strcmp(uregc,ureg(ri)),6:8),normwn(strcmp(uregn,ureg(ri)),6:8),0];
    %only nonmem
    else
        bardata=[bardata;0,0,0,normwn(strcmp(uregn,ureg(ri)),6:8),1];
    end
end


figure('Color','w','Position',[32,32,200,175])
hold on
bh=bar(bardata(:,[1,4]),1,'FaceColor','w','EdgeColor','k');
bh(1).FaceColor='k';
for ri=1:size(bardata,1)
    if any(bardata(ri,1:2)>62)
        text(ri,5,sprintf('%.1f,%.1f',bardata(ri,1),bardata(ri,4)),'HorizontalAlignment','center','VerticalAlignment','middle','Rotation',90);
    end
    
    if bardata(ri,7)
        plot(bh(1).XEndPoints(ri),5,'kv','MarkerSize',3,'MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor','none');
    end
end
errorbar(bh(1).XEndPoints,bardata(:,1),diff(bardata(:,1:2),1,2),diff(bardata(:,1:2:3),1,2),'k','CapSize',2,'Linestyle','none');
errorbar(bh(2).XEndPoints,bardata(:,4),diff(bardata(:,4:5),1,2),diff(bardata(:,4:2:6),1,2),'k','CapSize',2,'Linestyle','none');
set(gca(),'XTick',1:numel(ureg),'XTickLabel',ureg,'XTickLabelRotation',90)
ylabel('Number of neurons in loops (Z-score)')
if ~isempty(opt.ylim)
    ylim(opt.ylim)
end
% title(sprintf('%d-neuron loops',spidx+2));
end



