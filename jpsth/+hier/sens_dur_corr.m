function fh=sens_dur_corr(dur_only_map,sens_only_map,opt)
arguments
    dur_only_map
    sens_only_map
    opt.range (1,:) char {mustBeMember(opt.range,{'grey','CH','CTX'})} = 'grey'
end
grey_regs=ephys.getGreyRegs('range',opt.range);
idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
ureg=intersect(intersect(dur_only_map.keys(),sens_only_map.keys()),grey_regs);

% dur_only_map K:\code\jpsth\+ephys\duration_reg_bars.m
% sens_only_map K:\code\jpsth\+ephys\Both_either_reg_bars.m
fh=figure('Color','w','Position',[32,32,320,320]);
hold on
corrmat=[];
for reg=reshape(ureg,1,[])
    c=ephys.getRegColor(reg{1},'large_area',true);
    sens=sens_only_map(reg{1});
    dur=dur_only_map(reg{1});
    scatter(sens(1),dur(1),4,c,'filled','o')
    text(sens(1),dur(1),reg{1},'HorizontalAlignment','center','VerticalAlignment','top','Color',c);
    corrmat=[corrmat;idmap.reg2ccfid(reg{1}),sens(1),dur(1)];
end
[r,p]=corr(corrmat(:,2),corrmat(:,3),'type','Spearman');
xlabel('Sensory neuron proportion (%)');
ylabel('Duration neuron proportion (%)');
set(gca,'XScale','log','YScale','log')
title(sprintf(' r = %.3f, p = %.3f',r,p));
% xlim([0.003,0.4])
% if max(ylim())>0.3
%     ylim([0.1,0.35])
%     set(gca(),'YTick',[0.1,0.2,0.3],'YTickLabel',[0.1,0.2,0.3].*100,'XTick',[0.02,0.1,0.2],'XTickLabel',[0.02,0.1,0.2].*100)
% else
%     ylim([0.02,0.2])
%     set(gca(),'YTick',[0.02,0.1,0.2],'YTickLabel',[0.02,0.1,0.2].*100,'XTick',[0.02,0.1,0.2],'XTickLabel',[0.02,0.1,0.2].*100)
% end

% exportgraphics(fh,'sens_dur_corr.pdf','ContentType','vector')
%             close(fh)
end