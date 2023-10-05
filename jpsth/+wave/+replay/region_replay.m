function [fh,out]=region_replay(motif_replay,opt)
arguments
    motif_replay
    opt.reg="HIP"
    opt.iti (1,:) char {mustBeMember(opt.iti,{'nonpref_precede_ITI','pref_succeed_ITI'})}='pref_succeed_ITI'
    opt.normalize (1,1) logical = false
end
out=cell2struct({[];[]},{'REG','Others'});
su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);

for midx=1:size(motif_replay,1)
    sessid=motif_replay.session(midx);
    cids=motif_replay.meta{midx,2};
    motif_reg=arrayfun(@(x) string(su_meta.reg_tree{5,su_meta.sess==sessid & su_meta.allcid==x}),cids);
    onefreq=[motif_replay.freqstats{midx}.pref_delay_correct,...
        motif_replay.freqstats{midx}.nonpref_delay_correct,...
        motif_replay.freqstats{midx}.(opt.iti),...
        motif_replay.freqstats{midx}.before_session,...
        motif_replay.freqstats{midx}.after_session];

    if any(motif_reg==opt.reg)
        out.REG=[out.REG;onefreq(1:2:9)./onefreq(2:2:10)];
    else
        out.Others=[out.Others;onefreq(1:2:9)./onefreq(2:2:10)];
    end
end

if false
    dataout.HIP_preferred_delay=out.REG(:,1);
    dataout.HIP_nonpreferred_delay=out.REG(:,2);
    dataout.HIP_ITI=out.REG(:,3);
    dataout.HIP_before_task=out.REG(:,4);
    dataout.HIP_after_task=out.REG(:,5);
    dataout.Others_preferred_delay=out.Others(:,1);
    dataout.Others_nonpreferred_delay=out.Others(:,2);
    dataout.Others_ITI=out.Others(:,3);
    dataout.Others_before_task=out.Others(:,4);
    dataout.Others_after_task=out.Others(:,5);
    if false
        fid=fopen(fullfile('binary','upload','F3I_Chain_spike_frequency_with_without_HIP_neuron.json'),'w');
    else
        fid=fopen(fullfile('binary','upload','F3J_Loop_spike_frequency_with_without_HIP_neuron.json'),'w');
    end
    fprintf(fid,jsonencode(dataout));
    fclose(fid);

end

if opt.normalize
    out.Others=out.Others./out.Others(:,1);
    out.REG=out.REG./out.REG(:,1);
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

for ii=1:2:10
    pp=ranksum(dd(gg==ii), dd(gg==ii+1));
    text((ii+1)./2,max(ylim()),sprintf('%.3f',pp),'VerticalAlignment','top','HorizontalAlignment','center');
end

end

