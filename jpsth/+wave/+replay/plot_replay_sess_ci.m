function fhb=plot_replay_sess_ci(per_sess_mat,xlbl,opt)
arguments
    per_sess_mat double
    xlbl char
    opt.title (1,:) char = 'chains'
    opt.ref_line (1,1) logical = false
    opt.ref_p_value (1,1) logical = true
    opt.median_value (1,1) logical = false
    opt.ratio_block (1,1) double = 0
    opt.stats char {mustBeMember(opt.stats,{'mean','median'})} = 'median'
end

if opt.ratio_block>0
    if opt.ratio_block==2
        ratios=per_sess_mat(:,2:2:size(per_sess_mat,2)) ...
            ./per_sess_mat(:,1:2:size(per_sess_mat,2));
    elseif opt.ratio_block==3
        ratios=per_sess_mat(:,[2:3:size(per_sess_mat,2),3:3:size(per_sess_mat,2)]) ...
            ./repmat(per_sess_mat(:,1:3:size(per_sess_mat,2)),1,2);
    end
    if strcmp(opt.stats,'median')
        mm=nanmedian(ratios(:,[1 4 2 5 3 6])); % TODO: magic numbers?
        ci=bootci(1000,@(x) nanmedian(x),ratios(:,[1 4 2 5 3 6]));
    else
        keyboard(); % ratio might not work

        mm=nanmean(ratios(:,[1 4 2 5 3 6]));
        ci=nanstd(ratios(:,[1 4 2 5 3 6]))./sqrt(sum(isfinite(ratios)));
        ci=[ci;-ci]+mm;
    end
else
    cmatv=reshape(per_sess_mat,[],1);
    cmatg=reshape(repmat(1:size(per_sess_mat,2),size(per_sess_mat,1),1),[],1);
    if strcmp(opt.stats,'median')
        mm=nanmedian(per_sess_mat);
        ci=bootci(1000,@(x) nanmedian(x),per_sess_mat);
    else
        mm=nanmean(per_sess_mat);
        ci=nanstd(per_sess_mat)./sqrt(sum(isfinite(per_sess_mat)));
        ci=[ci;-ci]+mm;
    end
end

fhb=figure('Position',[100,100,600,400]);
hold on
bar(mm.','grouped','FaceColor','none','EdgeColor','k')
errorbar(1:numel(mm),mm,ci(1,:)-mm,ci(2,:)-mm,'k.');
if opt.ref_line
    if opt.ratio_block>0
        yline(1,'--r');
    else
        yline(median(per_sess_mat(:,1)),'--r');
    end
end
set(gca(),'XTick',1:size(per_sess_mat,2),'XTickLabel',xlbl,'XTickLabelRotation',90,'YLim',[0,4])
if opt.ref_p_value
    for jj=2:size(per_sess_mat,2)
        finisel=isfinite(per_sess_mat(:,jj));
        pp=signrank(per_sess_mat(finisel,1),per_sess_mat(finisel,jj));
        text(jj,0.01,sprintf('%.3f',pp),'VerticalAlignment','bottom','HorizontalAlignment','center')
    end
end

if opt.median_value
    for jj=1:size(per_sess_mat,2)
        % mm=median(cmatv(cmatg==jj & isfinite(cmatv)));
        text(jj,1,sprintf('%.1f',mm(jj)),'VerticalAlignment','bottom','HorizontalAlignment','center')
    end
end

title(opt.title)
ylabel('Motif spike frequency (Hz)')

end
