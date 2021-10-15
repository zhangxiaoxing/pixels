%% TODO merge with ring_span_plot.m
%% TODO converge to new data structure
if ~exist('ring_meta','var')
    ring_meta=bz.rings.get_ring_meta('loadfile',true);
end

keyboard()
%% within
% if false
%     w=join_within(ring_meta);
%     plot_two(true,w.congru.within,...
%         w.nonmem.within,...
%         w.congru.within_shuf,...
%         w.nonmem.within_shuf,'ylim',[-10,70]);
%     exportgraphics(gcf,'loops_within_region_norm.pdf')
% end
plot_two(ring_meta.congru.within_4,...
    ring_meta.nonmem.within_4,...
    ring_meta.congru.within_4_shuf,...
    ring_meta.nonmem.within_4_shuf,...
    ring_meta.congru.cross_4,...
    ring_meta.nonmem.cross_4,...
    ring_meta.congru.cross_4_shuf,...
    ring_meta.nonmem.cross_4_shuf,'ylim',[-10,90]);
exportgraphics(gcf,'loops_per_region_norm.pdf')


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


function [ureg,normw]=plot_two(datacw,datanw,shufcw,shufnw,datacc,datanc,shufcc,shufnc,opt)
arguments
    datacw
    datanw
    shufcw
    shufnw
    datacc
    datanc
    shufcc
    shufnc
    opt.ylim (1,:) double = []
end

[uregcw,normwcw]=norm_within(datacw,shufcw); %reg,{real-mean, shuf-mean,ci(1),ci(2),shuf-std, z-score,zci(1),zci(2)}
[uregnw,normwnw]=norm_within(datanw,shufnw);
[uregcc,normwcc]=norm_cross(datacc,shufcc);
[uregnc,normwnc]=norm_cross(datanc,shufnc);

ureg=unique([uregcw;uregnw;uregcc;uregnc]);
% [~,~,ratiomap]=ref.get_pv_sst();
load('OBM1map.mat','OBM1map');
msel=cellfun(@(x) OBM1map.isKey(x), ureg);
ureg=ureg(msel);
trans_idces=cellfun(@(x) OBM1map(x),ureg);
[trans_idces,sidx]=sort(trans_idces);
ureg=ureg(sidx);

bardata=[];
for ri=1:numel(ureg)
    linedata=[];
    %congru and nonmem
    if any(strcmp(uregcw,ureg(ri)))
        linedata=[1,normwcw(strcmp(uregcw,ureg(ri)),6:8)];
    else
        linedata=[0,0,0,0];
    end
    
    if any(strcmp(uregnw,ureg(ri)))
        linedata=[linedata,1,normwnw(strcmp(uregnw,ureg(ri)),6:8)];
    else
        linedata=[linedata,0,0,0,0];
    end    

    if any(strcmp(uregcc,ureg(ri)))
        linedata=[linedata,1,normwcc(strcmp(uregcc,ureg(ri)),6:8)];
    else
        linedata=[linedata,0,0,0,0];
    end
    
    if any(strcmp(uregnc,ureg(ri)))
        linedata=[linedata,1,normwnc(strcmp(uregnc,ureg(ri)),6:8)];
    else
        linedata=[linedata,0,0,0,0];
    end    
    bardata=[bardata;linedata];
end


figure('Color','w','Position',[32,32,700,175])
hold on
bh=bar(bardata(:,2:4:16),1,'FaceColor','w','EdgeColor','k');
[bh(1).FaceColor,bh(2).FaceColor,bh(3).FaceColor,bh(4).FaceColor]=deal('r','k','r','k');
[bh(1).FaceAlpha,bh(2).FaceAlpha]=deal(0.5);
for ri=1:size(bardata,1)
    if ~bardata(ri,1)
        plot(bh(1).XEndPoints(ri),5,'kv','MarkerSize',3,'MarkerFaceColor',[1,0.5,0.5],'MarkerEdgeColor','none');
    end
    if ~bardata(ri,5)
        plot(bh(2).XEndPoints(ri),5,'kv','MarkerSize',3,'MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor','none');
    end
    if ~bardata(ri,9)
        plot(bh(3).XEndPoints(ri),5,'kv','MarkerSize',3,'MarkerFaceColor','r','MarkerEdgeColor','none');
    end
end
errorbar(bh(1).XEndPoints,bardata(:,2),diff(bardata(:,2:3),1,2),diff(bardata(:,2:2:4),1,2),'k','CapSize',3,'Linestyle','none');
errorbar(bh(2).XEndPoints,bardata(:,6),diff(bardata(:,6:7),1,2),diff(bardata(:,6:2:8),1,2),'k','CapSize',3,'Linestyle','none');
errorbar(bh(3).XEndPoints,bardata(:,10),diff(bardata(:,10:11),1,2),diff(bardata(:,10:2:12),1,2),'k','CapSize',3,'Linestyle','none');
errorbar(bh(4).XEndPoints,bardata(:,14),diff(bardata(:,15:15),1,2),diff(bardata(:,14:2:16),1,2),'k','CapSize',3,'Linestyle','none');
set(gca(),'XTick',1:numel(ureg),'XTickLabel',ureg,'XTickLabelRotation',90)
ylabel('Number of neurons in loops (Z-score)')
if ~isempty(opt.ylim)
    ylim(opt.ylim)
end
% title(sprintf('%d-neuron loops',spidx+2));
end


function within=join_within(ring_meta)
within=struct();

 for mtype=["congru","nonmem"]
     reg=cell(0);
     shuf=cell(0);
%      meta=[];
     for rsize=3:5
         reg=[reg;ring_meta.(mtype).(sprintf('within_%d',rsize)).reg];
%          meta=[meta;ring_meta.(mtype).(sprintf('within_%d',rsize)).meta,...
%              nan(numel(ring_meta.(mtype).(sprintf('within_%d',rsize)).reg),5-rsize)];
         shuf=[shuf;ring_meta.(mtype).(sprintf('within_%d_shuf',rsize))];
     end
     within.(mtype).within.reg=reg;
     within.(mtype).within_shuf=shuf;
%      within.(mtype).meta=meta;
 end
end
