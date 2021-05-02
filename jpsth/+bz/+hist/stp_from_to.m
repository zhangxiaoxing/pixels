function fh=stp_from_to(stats,opt)
arguments
    stats (1,1) struct
    opt.treedepth (1,1) double {mustBeInteger,mustBeInRange(opt.treedepth,1,7)} = 5;
    opt.title (1,:) char = []
    opt.minpair (1,1) double = 30
end

regs=containers.Map('KeyType','double','ValueType','any');
for dep=opt.treedepth
    regs(dep)=cell(0);
    for mt=fieldnames(stats)
        C=cellfun(@(x) horzcat(x{1}(dep),x{2}(dep)),stats.(mt{1}).reg(stats.(mt{1}).same_reg(:,dep) |stats.(mt{1}).diff_reg(:,dep)),'UniformOutput',false);
        regs(dep)=unique(horzcat(regs(dep),(cat(2,C{:}))));
    end
end
%maxiter->[SPK,FC_EFF,FC_PROB]

perdep=containers.Map('KeyType','double','ValueType','any');
for dep=opt.treedepth
    perreg=struct();
    for creg=regs(dep)
        perreg.(creg{1})=struct();
        for mt=fieldnames(stats).'
            perreg.(creg{1}).(mt{1})=struct();
            fromsel=cellfun(@(x) strcmp(x{1}{dep},creg{1}),stats.(mt{1}).reg);
            tosel=cellfun(@(x) strcmp(x{2}{dep},creg{1}),stats.(mt{1}).reg);
            itersel=stats.(mt{1}).maxiter(:,1)~=1;
            perreg.(creg{1}).(mt{1}).samedata=sum(stats.(mt{1}).postspk(itersel & fromsel & stats.(mt{1}).same_reg(:,dep),2:end),2);
            perreg.(creg{1}).(mt{1}).fromdata=sum(stats.(mt{1}).postspk(itersel & fromsel & stats.(mt{1}).diff_reg(:,dep),2:end),2);
            perreg.(creg{1}).(mt{1}).todata=sum(stats.(mt{1}).postspk(itersel & tosel & stats.(mt{1}).diff_reg(:,dep),2:end),2);
        end
    end
    perdep(dep)=perreg;
end

dcolor={'m','b','k'};

mtypemat={'congru','nonmem','Congruent','Non-memory';...
    'congru','incongru','Congruent','Incongruent';...
    'incongru','nonmem','Incongruent','Non-memory'};

corrmat={'fromdata','todata','From region','To region';...
    'fromdata','samedata','From region','Within region';...
    'todata','samedata','To region','Within region'};
for dep=opt.treedepth
    fh=figure('Color','w','Position',[1441,31,1280,917]);
    for j=1:size(mtypemat,1)
        for i=1:size(corrmat,1)
            %     allcorr=[];
            
            subplot(3,3,(j-1)*3+i);
            hold on;
            regs=fieldnames(perdep(dep));
            count=cell2mat(cellfun(@(x) [nnz(~isnan(perdep(dep).(x).(mtypemat{j,1}).fromdata)),nnz(~isnan(perdep(dep).(x).(mtypemat{j,2}).fromdata))],regs,'UniformOutput',false));
            regs=regs(~any(count<opt.minpair,2));
            xx=cell2mat(cellfun(@(x) nanmean(perdep(dep).(x).(mtypemat{j,1}).(corrmat{i,1}))-nanmean(perdep(dep).(x).(mtypemat{j,2}).(corrmat{i,1})),regs,'UniformOutput',false));
            yy=cell2mat(cellfun(@(x) nanmean(perdep(dep).(x).(mtypemat{j,1}).(corrmat{i,2}))-nanmean(perdep(dep).(x).(mtypemat{j,2}).(corrmat{i,2})),regs,'UniformOutput',false));
            scatter(xx,yy,((6-dep)*5).^2,dcolor{dep-2},'filled','MarkerFaceAlpha',0.4)
            text(xx,yy,regs,'FontSize',7)
            %         nsel=~any(isnan([xx,yy]),2);
            [r,p]=corrcoef(xx,yy);
            text(max(xlim()),max(ylim()),sprintf('r = %0.3f\np = %0.3f',r(1,2),p(1,2)),'VerticalAlignment','top','HorizontalAlignment','right','FontSize',8);
            xlabel(corrmat{i,3});
            ylabel(corrmat{i,4});
            title(sprintf('%s - %s Area',mtypemat{j,3},mtypemat{j,4}));
        end
    end
end