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
        fcdata.(expd(d)).all(bin)=nnz(stats.(expd(d))(:,bin,1)<pp)./size(stats.(expd(d)),1);
        fcdata.(expd(d)).non(bin)=nnz(stats.(expd(d))(:,bin,1)<pp ...
            & stats.(expd(d))(:,bin,2)>=pp...
            & stats.(expd(d))(:,bin,3)>=pp)...
            ./nnz(stats.(expd(d))(:,bin,2)>=pp...
            & stats.(expd(d))(:,bin,3)>=pp);
        fcdata.(expd(d)).pre(bin)=nnz(stats.(expd(d))(:,bin,1)<pp ...
            & stats.(expd(d))(:,bin,2)<pp...
            & stats.(expd(d))(:,bin,3)>=pp)...
            ./nnz(stats.(expd(d))(:,bin,2)<pp...
            & stats.(expd(d))(:,bin,3)>=pp);
        fcdata.(expd(d)).post(bin)=nnz(stats.(expd(d))(:,bin,1)<pp ...
            & stats.(expd(d))(:,bin,2)>=pp...
            & stats.(expd(d))(:,bin,3)<pp)...
            ./nnz(stats.(expd(d))(:,bin,2)>=pp...
            & stats.(expd(d))(:,bin,3)<pp);
        fcdata.(expd(d)).both(bin)=nnz(stats.(expd(d))(:,bin,1)<pp ...
            & stats.(expd(d))(:,bin,2)<pp...
            & stats.(expd(d))(:,bin,3)<pp)...
            ./nnz(stats.(expd(d))(:,bin,2)<pp...
            & stats.(expd(d))(:,bin,3)<pp);
    end
end

fh=figure('Position',[100,100,400,300],'Color','w');
for d=[3 6]
    subplot(1,2,d/3)
    hold on
    pall=plot((-3:10)+0.5,fcdata.(expd(d)).all*100,'k-');
    pnon=plot((-3:10)+0.5,fcdata.(expd(d)).non*100,'r-');
    ppre=plot((-3:10)+0.5,fcdata.(expd(d)).pre*100,'b-');
    ppost=plot((-3:10)+0.5,fcdata.(expd(d)).post*100,'m-');
    pboth=plot((-3:10)+0.5,fcdata.(expd(d)).post*100,'c-');
    xlim([-3,11]);
    ylim([0,50]);
    arrayfun(@(x) xline(x,'k:'),[0,1,7,8]);
    title('6s delay, n=65591 pairs')
    xlabel('Time (s)');
    ylabel('Fraction of selective coupling (%)');
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
    legend([pall,pnon,ppre,ppost,pboth],{'All pairs',...
        'non-rate-selective',...
        'pre-rate-selective',...
        'post-rate-selective',...
        'both-rate-selective'});
    % exportgraphics(fh,'fc_coding_6s.pdf')
end

return

for fidx=1:163
    fpath=fullfile('K:','code','jpsth','fcdata',sprintf('fc_coding_f%d.mat',fidx));
    if ~isfile(fpath)
        continue
    end
    load(fpath)
    for i=1:size(sums,1)
        stats3(end+1,:,:)=sums{i,3}(:,19:21,1);
    end
end

pp=0.05/14;
for bin=1:14

end
fh=figure('Position',[100,100,400,300],'Color','w');
hold on
pall=plot((-3:10)+0.5,fcFrac3*100,'b-');
pnonsel=plot((-3:10)+0.5,fcOnly3*100,'r-');
xlim([-3,11]);
ylim([0,1]);
arrayfun(@(x) xline(x,'k:'),[0,1,4,5]);
title('3s delay, n=65591 pairs')

xlabel('Time (s)');
ylabel('Fraction of selective coupling (%)');

chisq3=ones(11,1);
for bin=4:14
    [~,~,chisq3(bin-3)]=crosstab([zeros(nnz(stats3(:,1,2)>pp & stats3(:,1,3)>pp),1);...
        zeros(nnz(stats3(:,2,2)>pp & stats3(:,2,3)>pp),1);...
        zeros(nnz(stats3(:,3,2)>pp & stats3(:,3,3)>pp),1);...
        ones(nnz(stats3(:,bin,2)>pp & stats3(:,bin,3)>pp),1)],...
        [1:nnz(stats3(:,1,2)>pp & stats3(:,1,3)>pp)>nnz(stats3(:,1,1)<pp & stats3(:,1,2)>pp & stats3(:,1,3)>pp),...
        1:nnz(stats3(:,2,2)>pp & stats3(:,2,3)>pp)>nnz(stats3(:,2,1)<pp & stats3(:,2,2)>pp & stats3(:,2,3)>pp),...
        1:nnz(stats3(:,3,2)>pp & stats3(:,3,3)>pp)>nnz(stats3(:,3,1)<pp & stats3(:,3,2)>pp & stats3(:,3,3)>pp),...
        1:nnz(stats3(:,bin,2)>pp & stats3(:,bin,3)>pp)>nnz(stats3(:,bin,1)<pp & stats3(:,bin,2)>pp & stats3(:,bin,3)>pp)]');
     if chisq3(bin-3)<(0.001/11)
         text(bin-3.5,0.4,'***','HorizontalAlignment','center','Color','r')
     elseif chisq3(bin-3)<(0.01/11)
         text(bin-3.5,0.4,'**','HorizontalAlignment','center','Color','r')
     elseif chisq3(bin-3)<(0.05/11)
         text(bin-3.5,0.4,'*','HorizontalAlignment','center','Color','r')
     else
         text(bin-3.5,0.4,'ns','HorizontalAlignment','center','Color','r')
     end
end
legend([pall,pnonsel],{'All pairs','non-rate-selective pairs'});
exportgraphics(fh,'fc_coding_3s.pdf')



for bin=1:14
    fconlypair3(bin)=nnz(stats3(:,bin,2)>pp & stats3(:,bin,3)>pp);
    fconlypair6(bin)=nnz(stats6(:,bin,2)>pp & stats6(:,bin,3)>pp);
    
    fconlycnt3(bin)=nnz(stats3(:,bin,1)<pp & stats3(:,bin,2)>pp & stats3(:,bin,3)>pp);
    fconlycnt6(bin)=nnz(stats6(:,bin,1)<pp & stats6(:,bin,2)>pp & stats6(:,bin,3)>pp);
    
end

figure('Position',[100,100,400,300]);
hold on
plot(fconlypair6,'--b');
plot(fconlypair3,'--r');
figure('Position',[100,100,400,300]);
hold on
plot(fconlycnt6,'-b');
plot(fconlycnt3,'-r');
