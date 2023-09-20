function [fh,out]=region_replay(motif_replay,opt)
arguments
    motif_replay
    opt.reg="HIP"
end
out=cell2struct({[];[]},{'REG','Others'});
su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);

for midx=1:size(motif_replay,1)
    sessid=motif_replay.session(midx);
    cids=motif_replay.meta{midx,2};
    motif_reg=arrayfun(@(x) string(su_meta.reg_tree{5,su_meta.sess==sessid & su_meta.allcid==x}),cids);
    onefreq=[motif_replay.freqstats{midx}.pref_delay_correct,...
        motif_replay.freqstats{midx}.nonpref_delay_correct,...
        motif_replay.freqstats{midx}.pref_succeed_ITI,...
        motif_replay.freqstats{midx}.before_session,...
        motif_replay.freqstats{midx}.after_session];

    if any(motif_reg==opt.reg)
        out.REG=[out.REG;onefreq(1:2:9)./onefreq(2:2:10)];
    else
        out.Others=[out.Others;onefreq(1:2:9)./onefreq(2:2:10)];
    end
end
dd=reshape([out.REG;out.Others],[],1);
gg=reshape([repmat([1:2:10],size(out.REG,1),1);...
    repmat([2:2:10],size(out.Others,1),1)],[],1);

mm=arrayfun(@(x) median(dd(gg==x & isfinite(dd))),unique(gg));
rci=cell2mat(arrayfun(@(x) bootci(100,@(x) median(x), dd(gg==x & isfinite(dd))),unique(gg),'UniformOutput',false).');

fh=figure();
hold on
bh=bar(mm([1:2:9;2:2:10].'),'grouped','FaceColor','none','EdgeColor','k');
bh(1).FaceColor='k';
errorbar([bh.XEndPoints],[bh.YEndPoints],rci(1,[1:2:9,2:2:10])-mm([1:2:9,2:2:10]).',rci(2,[1:2:9,2:2:10])-mm([1:2:9,2:2:10]).','k.');
set(gca,'XTick',1:5,'XTickLabel',{'Delay','N.P.Delay','ITI','Before task','After task'})
% ylim([0,2])

for ii=1:5
    pp=ranksum(dd(gg==ii), dd(gg==ii+1));
    text(ii,max(ylim()),sprintf('%.3f',pp),'VerticalAlignment','top','HorizontalAlignment','center');
end

end

