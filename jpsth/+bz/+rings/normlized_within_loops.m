%% TODO merge with ring_span_plot.m
%% TODO converge to new data structure
if ~exist('ring_meta','var')
    ring_meta=bz.rings.get_ring_meta('loadfile',true);
    jwithin=join_within(ring_meta);
    jcross=join_cross(ring_meta);
end

% keyboard()
plot_two(jwithin.congru.within,...
    jwithin.nonmem.within,...
    jwithin.congru.within_shuf,...
    jwithin.nonmem.within_shuf,...
    jcross.congru.cross,...
    jcross.nonmem.cross,...
    jcross.congru.cross_shuf,...
    jcross.nonmem.cross_shuf);
exportgraphics(gcf,'loops_per_region_norm.pdf')
ylim([0.01,100])
set(gca(),'YScale','log')
yline(1.96,'--k')
exportgraphics(gcf,'loops_per_region_norm_log.pdf')

%% dependency
function [ureg,normw]=norm_reg(rreal,shuf)
% wshuf=cellfun(@(x) x.reg,shuf,'UniformOutput',false);
ureg=unique([rreal;shuf]);
[B,BG]=groupcounts(rreal);
normw=zeros(numel(ureg),6);
[BS,BGS]=groupcounts(shuf);
for ri=1:numel(ureg)
    if ismember(ureg(ri),BG)
        normw(ri,1)=B(strcmp(BG,ureg(ri)));
    end
    [phat,pci]=binofit(BS(strcmp(BGS,ureg(ri))),numel(shuf));
    normw(ri,2)=phat;
    normw(ri,3:4)=reshape(pci,1,[]);
    normw(ri,5)=sqrt(numel(shuf)*phat*(1-phat));
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

[uregcw,normwcw]=norm_reg(datacw,shufcw); %reg,{real-mean, shuf-mean,ci(1),ci(2),shuf-std, z-score,zci(1),zci(2)}
[uregnw,normwnw]=norm_reg(datanw,shufnw);
[uregcc,normwcc]=norm_reg(datacc,shufcc);
[uregnc,normwnc]=norm_reg(datanc,shufnc);

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
    if any(strcmp(uregnw,ureg(ri)))
        linedata=[linedata,1,normwnw(strcmp(uregnw,ureg(ri)),6:8)];
    else
        linedata=[linedata,0,0,0,0];
    end
    if any(strcmp(uregnc,ureg(ri)))
        linedata=[linedata,1,normwnc(strcmp(uregnc,ureg(ri)),6:8)];
    else
        linedata=[linedata,0,0,0,0];
    end
    if any(strcmp(uregcw,ureg(ri)))
        linedata=[linedata,1,normwcw(strcmp(uregcw,ureg(ri)),6:8)];
    else
        linedata=[linedata,0,0,0,0];
    end
    
    if any(strcmp(uregcc,ureg(ri)))
        linedata=[linedata,1,normwcc(strcmp(uregcc,ureg(ri)),6:8)];
    else
        linedata=[linedata,0,0,0,0];
    end
    bardata=[bardata;linedata];
end


figure('Color','w','Position',[32,32,700,175])
hold on
bh=bar(bardata(:,2:4:16),1,'FaceColor','w','EdgeColor','k');
[bh(1).FaceColor,bh(2).FaceColor,bh(3).FaceColor,bh(4).FaceColor]=deal('k','k','r','r');
[bh(1).FaceAlpha,bh(3).FaceAlpha]=deal(0.5);
[bh(1).EdgeColor,bh(2).EdgeColor,bh(3).EdgeColor,bh(4).EdgeColor]=deal('none');
for ri=1:size(bardata,1)
    if ~bardata(ri,1)
        plot(bh(1).XEndPoints(ri),1,'kv','MarkerSize',3,'MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor','none');
    end
    if ~bardata(ri,5)
        plot(bh(2).XEndPoints(ri),1,'kv','MarkerSize',3,'MarkerFaceColor','k','MarkerEdgeColor','none');
    end
    if ~bardata(ri,9)
        plot(bh(3).XEndPoints(ri),1,'kv','MarkerSize',3,'MarkerFaceColor',[1,0.5,0.5]','MarkerEdgeColor','none');
    end
    if ~bardata(ri,13)
        plot(bh(4).XEndPoints(ri),1,'kv','MarkerSize',3,'MarkerFaceColor','r','MarkerEdgeColor','none');
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
        shufvec=cellfun(@(x) x.reg,ring_meta.(mtype).(sprintf('within_%d_shuf',rsize)),'UniformOutput',false);
        shuf=[shuf;cat(1,shufvec{:})];
    end
    within.(mtype).within=reg;
    within.(mtype).within_shuf=shuf;
    %      within.(mtype).meta=meta;
end
end


function cross=join_cross(ring_meta)
persistent flat
if isempty(flat)
    flat=@(x) cellfun(@(y) y{1},x,'UniformOutput',false);
end
cross=struct();
for mtype=["congru","nonmem"]
    reg=cell(0);
    shuf=cell(0);
    %      meta=[];
    for rsize=3:5
        reg=[reg;reshape(ring_meta.(mtype).(sprintf('cross_%d',rsize)).reg,[],1)];
        %          meta=[meta;ring_meta.(mtype).(sprintf('within_%d',rsize)).meta,...
        %              nan(numel(ring_meta.(mtype).(sprintf('within_%d',rsize)).reg),5-rsize)];
        shufvec=cellfun(@(x) reshape(flat(x.reg),[],1),ring_meta.(mtype).(sprintf('cross_%d_shuf',rsize)),'UniformOutput',false);
        shuf=[shuf;cat(1,shufvec{:})];
    end
    cross.(mtype).cross=flat(reg);
    cross.(mtype).cross_shuf=shuf;
end
end


