function plot_one_stp_curve(ref_data,cmp_data,title_,opt)
arguments
    ref_data (1,1) struct {mustBeNonempty}
    cmp_data (1,1) struct {mustBeNonempty}
    title_ (1,:) char
    opt.plot_dual (1,1) logical = false
    opt.ref_label (1,:) char = 'ref'
    opt.cmp_label (1,:) char = 'cmp'
end
hdl=[];
lgd=cell(0);

% fh=figure('Color','w');
hold on
if isfield(ref_data,'congru')
    if opt.plot_dual
        shadow(fliplr(cmp_data.congru(:,2:end)),200:200:2000,'r');
        hlnc=plot(200:200:2000,fliplr(mean(cmp_data.congru(:,2:end)))*100,'--r');
        hdl=cat(2,hdl,hlnc);
        lgd=cat(2,lgd,['Congruent ',opt.cmp_label]);
    end
    shadow(fliplr(ref_data.congru(:,2:end)),200:200:2000,'r');
    hwtc=plot(200:200:2000,fliplr(mean(ref_data.congru(:,2:end)))*100,'-r');
    hdl=cat(2,hdl,hwtc);
    lgd=cat(2,lgd,['Congruent ',opt.ref_label]);
end
if isfield(ref_data,'nonmem')
    if opt.plot_dual
        shadow(fliplr(cmp_data.nonmem(:,2:end)),200:200:2000,'k');
        hlnn=plot(200:200:2000,fliplr(mean(cmp_data.nonmem(:,2:end)))*100,'--k');
        hdl=cat(2,hdl,hlnn);
        lgd=cat(2,lgd,['NonMem. ',opt.cmp_label]);
    end
    shadow(fliplr(ref_data.nonmem(:,2:end)),200:200:2000,'k');
    hwtn=plot(200:200:2000,fliplr(mean(ref_data.nonmem(:,2:end)))*100,'-k');
    hdl=cat(2,hdl,hwtn);
    lgd=cat(2,lgd,['NonMem. ',opt.ref_label]);
end
if isfield(ref_data,'incong')
    if opt.plot_dual
    shadow(fliplr(cmp_data.incong(:,2:end)),200:200:2000,'b');
    hlni=plot(200:200:2000,fliplr(mean(cmp_data.incong(:,2:end)))*100,'--b');
    hdl=cat(2,hdl,hlni);
    lgd=cat(2,lgd,['Incongru. ',opt.cmp_label]);
    end
    shadow(fliplr(ref_data.incong(:,2:end)),200:200:2000,'b');
    
    hwti=plot(200:200:2000,fliplr(mean(ref_data.incong(:,2:end)))*100,'-b');
    hdl=cat(2,hdl,hwti);
    lgd=cat(2,lgd,['Incongru. ',opt.ref_label]);
end
% if ~isempty(hdl),legend(hdl,lgd);end
title(title_);
ylabel('Post spike gain (%)')
ylim([-1,5]);
xlim([0,2000])
set(gca,'XTick',500:500:2000)
xlabel('Time lag (ms)')
% grid on
nstr=cellfun(@(x) sprintf('%s,%d',x,size(ref_data.(x),1)),fieldnames(ref_data),'UniformOutput',false);
text(max(xlim()),max(ylim()),nstr,'HorizontalAlignment','right','VerticalAlignment','top');
end


function shadow(data,xx,color)
boots=bootstrp(1000,@(x) mean(x),data);
ci=[prctile(boots,2.5);prctile(boots,97.5)].*100;
fill([xx,fliplr(xx)],[ci(1,:),fliplr(ci(2,:))],color,'EdgeColor','none','FaceAlpha',0.1);
end