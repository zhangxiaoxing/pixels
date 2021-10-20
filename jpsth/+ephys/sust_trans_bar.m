function sust_trans_bar(opt)
arguments
    opt.good_wf (1,1) logical = false;
    opt.delay (1,1) double {mustBeMember(opt.delay,[3,6])} = 6
end
if opt.delay==6
    warning('Delay set to default 6')
end

meta=ephys.util.load_meta('delay',opt.delay);
upath=unique(meta.allpath);
for ii=1:numel(upath)
    if opt.good_wf
        sesssel=strcmp(meta.allpath,upath{ii}) & meta.good_waveform;
    else
        sesssel=strcmp(meta.allpath,upath{ii});
    end
    sustfrac(ii)=nnz(ismember(meta.mem_type(sesssel),[1,3]))./nnz(sesssel);
    transfrac(ii)=nnz(ismember(meta.mem_type(sesssel),[2,4]))./nnz(sesssel);
end
sustci=bootci(1000,@(x) mean(x), sustfrac).*100;
transci=bootci(1000,@(x) mean(x), transfrac).*100;
sustmm=mean(sustfrac.*100);
transmm=mean(transfrac.*100);
fh=figure('Color','w','Position',[100,100,150,235]);
hold on

bar(1,sustmm,'FaceColor','w','EdgeColor','k');
bar(2,transmm,'FaceColor','w','EdgeColor','k');
errorbar(1,sustmm,sustci(1)-sustmm,sustci(2)-sustmm,'k.');
errorbar(2,transmm,transci(1)-transmm,transci(2)-transmm,'k.');
ylabel('Fraction of all neurons (%)')
set(gca(),'XTick',1:2,'XTickLabel',{'Sustained','Transient'},'XTickLabelRotation',90,'FontSize',10)
if opt.good_wf
    text(1,mean(ylim()),num2str(nnz(ismember(meta.mem_type,[1 3]) & meta.good_waveform.')),'Rotation',90,'FontSize',10)
    text(2,mean(ylim()),num2str(nnz(ismember(meta.mem_type,[2 4]) & meta.good_waveform.')),'Rotation',90,'FontSize',10)
    text(max(xlim()),max(ylim()),num2str(nnz(meta.good_waveform)),'HorizontalAlignment','right','VerticalAlignment','top')
else
    text(1,mean(ylim()),num2str(nnz(ismember(meta.mem_type,[1 3]))),'Rotation',90,'FontSize',10)
    text(2,mean(ylim()),num2str(nnz(ismember(meta.mem_type,[2 4]))),'Rotation',90,'FontSize',10)
    text(max(xlim()),max(ylim()),num2str(numel(meta.mem_type)),'HorizontalAlignment','right','VerticalAlignment','top')
end
exportgraphics(fh,sprintf('sust_trans_fraction_%d.pdf',opt.delay));