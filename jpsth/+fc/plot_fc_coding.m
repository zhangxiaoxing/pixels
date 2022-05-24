% {pre,post,fc} x pair x bin

%                 onepair(sbin,:,2)=[mean(s61fc),std(s61fc),numel(s61fc),...
%                 mean(s62fc),std(s62fc),numel(s62fc),...
%                 mean(s61pre),std(s61pre),numel(s61pre),...
%                 mean(s62pre),std(s62pre),numel(s62pre),...
%                 mean(s61post),std(s61post),numel(s61post),...
%                 mean(s62post),std(s62post),numel(s62post),...
%                 p6,p6pre,p6post];
if exist('denovo','var') && denovo
    stats=struct();
    stats.d3=nan(0,14,3);
    stats.d6=nan(0,14,3);
    for fidx=1:163
        fpath=fullfile('K:','code','jpsth','fcdata',sprintf('fc_coding_f%d.mat',fidx));
        if ~isfile(fpath)
            continue
        end
        load(fpath)
        for i=1:size(sums,1)
            stats.d3(end+1,:,:)=sums{i,3}(:,19:21,1);
            stats.d6(end+1,:,:)=sums{i,3}(:,19:21,2);
        end
    end
end
expd=@(x) sprintf('d%d',x);
fcdata=struct();
fcdata.d3=struct();
fcdata.d6=struct();
pp=0.05/14;
for bin=1:14
    for d=[3 6]
        fcdata.(expd(d)).all(bin)=nnz(stats.(expd(d))(:,bin,1)<pp);
        fcdata.(expd(d)).non(bin)=nnz(stats.(expd(d))(:,bin,1)<pp ...
            & stats.(expd(d))(:,bin,2)>=pp...
            & stats.(expd(d))(:,bin,3)>=pp);
           
        fcdata.(expd(d)).pre(bin)=nnz(stats.(expd(d))(:,bin,1)<pp ...
            & stats.(expd(d))(:,bin,2)<pp...
            & stats.(expd(d))(:,bin,3)>=pp);
            
        fcdata.(expd(d)).post(bin)=nnz(stats.(expd(d))(:,bin,1)<pp ...
            & stats.(expd(d))(:,bin,2)>=pp...
            & stats.(expd(d))(:,bin,3)<pp);
        fcdata.(expd(d)).both(bin)=nnz(stats.(expd(d))(:,bin,1)<pp ...
            & stats.(expd(d))(:,bin,2)<pp...
            & stats.(expd(d))(:,bin,3)<pp);
        
        
        fcdata.(expd(d)).alln(bin)=size(stats.(expd(d)),1);
        fcdata.(expd(d)).nonn(bin)=nnz(stats.(expd(d))(:,bin,2)>=pp...
            & stats.(expd(d))(:,bin,3)>=pp);
        fcdata.(expd(d)).pren(bin)=nnz(stats.(expd(d))(:,bin,2)<pp...
            & stats.(expd(d))(:,bin,3)>=pp);
        fcdata.(expd(d)).postn(bin)=nnz(stats.(expd(d))(:,bin,2)>=pp...
            & stats.(expd(d))(:,bin,3)<pp);
        fcdata.(expd(d)).bothn(bin)=nnz(stats.(expd(d))(:,bin,2)<pp...
            & stats.(expd(d))(:,bin,3)<pp);
    end
end

posd=@(x) (x-3)/3*4;
fh=figure('Position',[100,100,1200,800],'Color','w');
for d=[3 6]
    h(d/3)=subplot(2,4,posd(d)+1);
    hold on
    pall=plot((-3:10)+0.5,fcdata.(expd(d)).all./fcdata.(expd(d)).alln*100,'k-');
    pnon=plot((-3:10)+0.5,fcdata.(expd(d)).non./fcdata.(expd(d)).nonn*100,'r-');
    ppre=plot((-3:10)+0.5,fcdata.(expd(d)).pre./fcdata.(expd(d)).pren*100,'b-');
    ppost=plot((-3:10)+0.5,fcdata.(expd(d)).post./fcdata.(expd(d)).postn*100,'m-');
    pboth=plot((-3:10)+0.5,fcdata.(expd(d)).both./fcdata.(expd(d)).bothn*100,'c-');
    xlim([-3,11]);
    ylim([0,50]);
    arrayfun(@(x) xline(x,'k:'),[0,1,d+1,d+2]);
    title([num2str(d) 's delay, n=65591 pairs'])
    xlabel('Time (s)');
    ylabel('Fraction of selective coupling (%)');
    legend([pall,pnon,ppre,ppost,pboth],{'All pairs',...
        'non-rate-selective',...
        'pre-rate-selective',...
        'post-rate-selective',...
        'both-rate-selective'});
    
    
    h(d/3)=subplot(2,4,posd(d)+2);
    hold on
    plot((-3:10)+0.5,fcdata.(expd(d)).all./fcdata.(expd(d)).alln*100,'k-');
    plot((-3:10)+0.5,fcdata.(expd(d)).non./fcdata.(expd(d)).nonn*100,'r-');
    xlim([-3,11]);
    ylim([0,1]);
    arrayfun(@(x) xline(x,'k:'),[0,1,d+1,d+2]);
    title([num2str(d) 's delay, n=65591 pairs'])
    xlabel('Time (s)');
    ylabel('Fraction of selective coupling (%)');
    
    h(d/3)=subplot(2,4,posd(d)+3);
    hold on
    plot((-3:10)+0.5,fcdata.(expd(d)).nonn./fcdata.(expd(d)).alln*100,'r-');
    plot((-3:10)+0.5,fcdata.(expd(d)).pren./fcdata.(expd(d)).alln*100,'b-');
    plot((-3:10)+0.5,fcdata.(expd(d)).postn./fcdata.(expd(d)).alln*100,'m-');
    plot((-3:10)+0.5,fcdata.(expd(d)).bothn./fcdata.(expd(d)).alln*100,'c-');
    xlim([-3,11]);
    ylim([0,100]);
    arrayfun(@(x) xline(x,'k:'),[0,1,d+1,d+2]);
    title([num2str(d) 's delay, n=65591 pairs'])
    xlabel('Time (s)');
    ylabel('Fraction of all pairs (%)');
    legend([pall,pnon,ppre,ppost,pboth],{'All pairs',...
        'non-rate-selective',...
        'pre-rate-selective',...
        'post-rate-selective',...
        'both-rate-selective'});
    
    
    h(d/3)=subplot(2,4,posd(d)+4);
    hold on
%     pall=plot((-3:10)+0.5,fcdata.(expd(d)).all./fcdata.(expd(d)).alln*100,'k-');
    plot((-3:10)+0.5,fcdata.(expd(d)).non./fcdata.(expd(d)).alln*100,'r-');
    plot((-3:10)+0.5,fcdata.(expd(d)).pre./fcdata.(expd(d)).alln*100,'b-');
    plot((-3:10)+0.5,fcdata.(expd(d)).post./fcdata.(expd(d)).alln*100,'m-');
    plot((-3:10)+0.5,fcdata.(expd(d)).both./fcdata.(expd(d)).alln*100,'c-');
    xlim([-3,11]);
    ylim([0,1]);
    arrayfun(@(x) xline(x,'k:'),[0,1,d+1,d+2]);
    title([num2str(d) 's delay, n=65591 pairs'])
    xlabel('Time (s)');
    ylabel('Fraction of selective coupling (%)');
    legend([pall,pnon,ppre,ppost,pboth],{'All pairs',...
        'non-rate-selective',...
        'pre-rate-selective',...
        'post-rate-selective',...
        'both-rate-selective'});
    
    
    exportgraphics(fh,'fc_coding.pdf')
    
    
end

return


% chisq6=ones(11,1);
% for bin=4:14
%     [~,~,chisq6(bin-3)]=crosstab([zeros(nnz(stats6(:,1,2)>pp & stats6(:,1,3)>pp),1);...
%         zeros(nnz(stats6(:,2,2)>pp & stats6(:,2,3)>pp),1);...
%         zeros(nnz(stats6(:,3,2)>pp & stats6(:,3,3)>pp),1);...
%         ones(nnz(stats6(:,bin,2)>pp & stats6(:,bin,3)>pp),1)],...
%         [1:nnz(stats6(:,1,2)>pp & stats6(:,1,3)>pp)>nnz(stats6(:,1,1)<pp & stats6(:,1,2)>pp & stats6(:,1,3)>pp),...
%         1:nnz(stats6(:,2,2)>pp & stats6(:,2,3)>pp)>nnz(stats6(:,2,1)<pp & stats6(:,2,2)>pp & stats6(:,2,3)>pp),...
%         1:nnz(stats6(:,3,2)>pp & stats6(:,3,3)>pp)>nnz(stats6(:,3,1)<pp & stats6(:,3,2)>pp & stats6(:,3,3)>pp),...
%         1:nnz(stats6(:,bin,2)>pp & stats6(:,bin,3)>pp)>nnz(stats6(:,bin,1)<pp & stats6(:,bin,2)>pp & stats6(:,bin,3)>pp)]');
%      if chisq6(bin-3)<(0.001/11)
%          text(bin-3.5,0.4,'***','HorizontalAlignment','center','Color','r')
%      elseif chisq6(bin-3)<(0.01/11)
%          text(bin-3.5,0.4,'**','HorizontalAlignment','center','Color','r')
%      elseif chisq6(bin-3)<(0.05/11)
%          text(bin-3.5,0.4,'*','HorizontalAlignment','center','Color','r')
%      else
%          text(bin-3.5,0.4,'ns','HorizontalAlignment','center','Color','r')
%      end
% end
