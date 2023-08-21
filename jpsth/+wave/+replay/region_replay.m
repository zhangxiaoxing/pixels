
function [fh,out]=region_replay(motif_replay,opt)
arguments
    motif_replay
    opt.reg="HIP"
end
out=cell2struct({[];[]},{'REG','Others'});
su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
for dd=reshape(string(fieldnames(motif_replay)),1,[])
    for samp=reshape(string(fieldnames(motif_replay.(dd))),1,[])
        for condi=reshape(string(fieldnames(motif_replay.(dd).(samp))),1,[])
            % sessid=str2double(regexp(condi,'(?<=^s)\d{1,3}(?=(c|r))','match','once'));
            onemotif=motif_replay.(dd).(samp).(condi);
            sessid=onemotif.meta{1};
            cids=onemotif.meta{2};
            motif_reg=arrayfun(@(x) string(su_meta.reg_tree{5,su_meta.sess==sessid & su_meta.allcid==x}),cids);
            onefreq=[onemotif.freqstats.pref_delay_correct,...
                onemotif.freqstats.pref_succeed_ITI,...
                onemotif.freqstats.before_session,...
                onemotif.freqstats.after_session];

            if any(motif_reg==opt.reg)
                out.REG=[out.REG;onefreq(1:2:7)./onefreq(2:2:8)];
            else
                out.Others=[out.Others;onefreq(1:2:7)./onefreq(2:2:8)];
            end
        end
    end
end
dd=reshape([out.REG;out.Others],[],1);
gg=reshape([repmat([1:2:8],size(out.REG,1),1);...
    repmat([2:2:8],size(out.Others,1),1)],[],1);

mm=arrayfun(@(x) median(dd(gg==x & isfinite(dd))),unique(gg));
rci=cell2mat(arrayfun(@(x) bootci(100,@(x) median(x), dd(gg==x & isfinite(dd))),unique(gg),'UniformOutput',false).');

fh=figure();
hold on
bar(mm.','grouped','FaceColor','none','EdgeColor','k')
errorbar(1:numel(mm),mm,rci(1,:)-mm.',rci(2,:)-mm.','k.');
set(gca,'XTick',1.5:2:8,'XTickLabel',{'Delay','ITI','Before task','After task'})
% ylim([0,2])

for ii=1:2:8
    pp=ranksum(dd(gg==ii), dd(gg==ii+1));
    text(ii+0.5,max(ylim()),sprintf('%.3f',pp),'VerticalAlignment','top');
end

end

