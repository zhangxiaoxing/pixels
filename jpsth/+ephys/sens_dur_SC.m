function fh=sens_dur_SC(metaidx,sens_meta,dur_meta,opt)
arguments
    metaidx
    sens_meta
    dur_meta
    opt.anova_meta
end
close all
% sens_meta=ephys.util.load_meta();
% dur_meta=ephys.get_dur_meta();
set(groot,'defaultTextFontSize',10);
% sens_meta=ephys.util.load_meta();
homedir=ephys.util.getHomedir('type','raw');
for ii=reshape(metaidx,1,[])
    
    fpath=fullfile(homedir,sens_meta.allpath{ii},'FR_All_ 250.hdf5');
    sessid=sens_meta.sess(ii);
    trials=h5read(fpath,'/Trials');
    suid=h5read(fpath,'/SU_id');
    fr=h5read(fpath,'/FR_All');
    
%     sesspath=regexp(fpath,'(?<=\\)[0-9a-zA-Z_\-]*\\[0-9a-zA-Z_\-]*(?=\\FR_All)','match','once');
%     fidx=ephys.path2sessid(sesspath);

    s1d3t=find(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==4 & trials(:,8)==3);
    s2d3t=find(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==8 & trials(:,8)==3);
    s1d6t=find(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==4 & trials(:,8)==6);
    s2d6t=find(trials(:,9)~=0 & trials(:,10)~=0 &trials(:,5)==8 & trials(:,8)==6);
    if min([numel(s1d3t),numel(s2d3t),numel(s1d6t),numel(s2d6t)])<2
        continue
    end

    fh=figure('Color','w','Position',[32,180,800,800]);
    %%
    ax3=subplot(2,2,1);
    hold on
    s1d3c=squeeze(fr(s1d3t,suid==sens_meta.allcid(ii),5:44));
    s2d3c=squeeze(fr(s2d3t,suid==sens_meta.allcid(ii),5:44));
    
%     mean(s2d3c)
    s1d3ci=[-1;1]*std(s1d3c)./sqrt(size(s1d3c,1)).'+mean(s1d3c);
    s2d3ci=[-1;1]*std(s2d3c)./sqrt(size(s2d3c,1)).'+mean(s2d3c);
    
    fill([1:40,fliplr(1:40)],[smooth(s1d3ci(1,:),5);flip(smooth(s1d3ci(2,:),5))],'r','EdgeColor','none','FaceAlpha',0.1);
    fill([1:40,fliplr(1:40)],[smooth(s2d3ci(1,:),5);flip(smooth(s2d3ci(2,:),5))],'b','EdgeColor','none','FaceAlpha',0.1);
    
    plot(smooth(mean(s1d3c),5),'--r')
    plot(smooth(mean(s2d3c),5),'--b')
    arrayfun(@(x) xline(x,'--k'),[8.5 12.5 24.5 28.5]);
    title('s1 v. s2, d3')
    set(gca,'XTick',12.5:4:40.5,'XTickLabel',0:6)
    xlim([4.5,36.5]);
    xlabel('Time (s)');
    ylabel('F.R.');

    arrayfun(@(x) text(10.5+x*4,max(ylim()),formatp(sens_meta.fdr_d3(ii,x+1)),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9),1:3)

    %%
    ax6=subplot(2,2,2);
    hold on
    s1d6c=squeeze(fr(s1d6t,suid==sens_meta.allcid(ii),5:44));
    s2d6c=squeeze(fr(s2d6t,suid==sens_meta.allcid(ii),5:44));
    
    s1d6ci=[-1;1]*std(s1d6c)./sqrt(size(s1d6c,1)).'+mean(s1d6c);
    s2d6ci=[-1;1]*std(s2d6c)./sqrt(size(s2d6c,1)).'+mean(s2d6c);
    
    fill([1:40,fliplr(1:40)],[smooth(s1d6ci(1,:),5);flip(smooth(s1d6ci(2,:),5))],'r','EdgeColor','none','FaceAlpha',0.1);
    fill([1:40,fliplr(1:40)],[smooth(s2d6ci(1,:),5);flip(smooth(s2d6ci(2,:),5))],'b','EdgeColor','none','FaceAlpha',0.1);

    plot(smooth(mean(s1d6c),5),'-r')
    plot(smooth(mean(s2d6c),5),'-b')
    
    xline(8.5,'--k')
    xline(12.5,'--k')
    title('s1 v. s2, d6')
    set(gca,'XTick',12.5:4:40.5,'XTickLabel',0:6)
    xlim([4.5,36.5]);
    xlabel('Time (s)');
    ylabel('F.R.');

    arrayfun(@(x) text(10.5+x*4,max(ylim()),formatp(sens_meta.fdr_d6(ii,x+1)),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9),1:3)


    %%
    axs1=subplot(2,2,3);
    hold on
    s1d3c=squeeze(fr(s1d3t,suid==sens_meta.allcid(ii),5:44));
    s1d6c=squeeze(fr(s1d6t,suid==sens_meta.allcid(ii),5:44));
    
    s1d3ci=[-1;1]*std(s1d3c)./sqrt(size(s1d3c,1)).'+mean(s1d3c);
    s1d6ci=[-1;1]*std(s1d6c)./sqrt(size(s1d6c,1)).'+mean(s1d6c);
    
    fill([1:40,fliplr(1:40)],[smooth(s1d3ci(1,:),5);flip(smooth(s1d3ci(2,:),5))],'r','EdgeColor','none','FaceAlpha',0.1);
    fill([1:40,fliplr(1:40)],[smooth(s1d6ci(1,:),5);flip(smooth(s1d6ci(2,:),5))],'r','EdgeColor','none','FaceAlpha',0.1);
    
    plot(smooth(mean(s1d3c),5),'--r')
    plot(smooth(mean(s1d6c),5),'-r')
    arrayfun(@(x) xline(x,'--k'),[8.5 12.5 24.5 28.5]);
    title('d3 v. d6 ,s1')
    set(gca,'XTick',12.5:4:40.5,'XTickLabel',0:6)
    xlim([4.5,36.5]);
    xlabel('Time (s)');
    ylabel('F.R.');
    arrayfun(@(x) text(10.5+x*4,max(ylim()),formatp(dur_meta.fdr_s1(ii,x+1)),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9),1:3)

    %%
    axs2=subplot(2,2,4);
    hold on
    s2d3c=squeeze(fr(s2d3t,suid==sens_meta.allcid(ii),5:44));
    s2d6c=squeeze(fr(s2d6t,suid==sens_meta.allcid(ii),5:44));
    
    s2d3ci=[-1;1]*std(s2d3c)./sqrt(size(s2d3c,1)).'+mean(s2d3c);
    s2d6ci=[-1;1]*std(s2d6c)./sqrt(size(s2d6c,1)).'+mean(s2d6c);
    
    fill([1:40,fliplr(1:40)],[smooth(s2d3ci(1,:),5);flip(smooth(s2d3ci(2,:),5))],'b','EdgeColor','none','FaceAlpha',0.1);
    fill([1:40,fliplr(1:40)],[smooth(s2d6ci(1,:),5);flip(smooth(s2d6ci(2,:),5))],'b','EdgeColor','none','FaceAlpha',0.1);

    plot(smooth(mean(s2d3c),5),'--b')
    plot(smooth(mean(s2d6c),5),'-b')
    
    arrayfun(@(x) xline(x,'--k'),[8.5 12.5 24.5 28.5]);
    title('d3 v. d6 ,s2')
    set(gca,'XTick',12.5:4:40.5,'XTickLabel',0:6)
    xlim([4.5,36.5]);
    xlabel('Time (s)');
    ylabel('F.R.');
    arrayfun(@(x) text(10.5+x*4,max(ylim()),formatp(dur_meta.fdr_s2(ii,x+1)),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9),1:3)
    arrayfun(@(x) set(x,'YLim',[min([ax3.YLim,ax6.YLim].'),max([ax3.YLim,ax6.YLim]).']),[ax3,ax6,axs1,axs2]);
    
    
    %%table
%     uitable(fh,'Data',rand(2,2),'Position',[0.55,0.1,0.4,0.8])
    
    
    if exist('opt','var') && isfield(opt,'anova_meta')
        ttstr=sprintf('#%d,ANOVA SEN%.3f,DUR%.3f,INT%.3f',ii,opt.anova_meta.anovap(ii,1),opt.anova_meta.anovap(ii,2),opt.anova_meta.anovap(ii,3));
        sgtitle(ttstr);
    else
        sgtitle(ii);
    end
    % exportgraphics(fh,'correct_error_sc.pdf','ContentType','vector')
end
end
function s=formatp(v)
if v>=0.1
    s='ns';
elseif v<0.05
    s=sprintf('%0.3f',v);
    s=s(2:end);
else
    s=sprintf('%0.3f',v);
    s=['#',s(3:end)];
end
end