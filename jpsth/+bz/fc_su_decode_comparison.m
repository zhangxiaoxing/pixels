
if false % denovo
    for nsu=[5 10 50 100 200 300 400 500]
        odor.("n"+nsu)=pct.pct_decoding_correct_error(wrs_mux_meta,5:6,'lblidx',5,'n_su',nsu);% odor
        both_odor.("n"+nsu)=pct.pct_decoding_correct_error(wrs_mux_meta,1:4,'lblidx',5,'n_su',nsu);% odor
        both_dur.("n"+nsu)=pct.pct_decoding_correct_error(wrs_mux_meta,1:4,'lblidx',8,'n_su',nsu);% dur
    end
    fc_dec=bz.fccoding.plot_coding(wrs_mux_meta,'nfc_grp',[5 10 50 100 200 300 400 500],'skip_plot',true);
end

odormat=cell2mat(cellfun(@(x) subsref(struct2cell(x.olf),substruct('{}',{1})),struct2cell(odor),'UniformOutput',false).');
summ=mean(odormat);
susem=std(odormat)./sqrt(size(odormat,1));


fcmat=cell2mat(cellfun(@(x) x.cvcorr{1},struct2cell(fc_dec),'UniformOutput',false).');
fcmm=mean(fcmat);
fsem=std(fcmat)./sqrt(size(fcmat,1));

figure()
hold on
fill([nsu,fliplr(nsu)],[summ+susem,fliplr(summ-susem)],'k','EdgeColor','none','FaceAlpha',0.2);
fill([nsu,fliplr(nsu)],[fcmm+fsem,fliplr(fcmm-fsem)],'k','EdgeColor','none','FaceAlpha',0.2);

hsu=plot([5 10 50 100 200 300 400 500],summ,'k-');
hfc=plot([5 10 50 100 200 300 400 500],fcmm,'k--');
ylim([0.5 1.0])
xlim([3,1000])
set(gca,'XScale','log','YTick',0.5:0.25:1,'YTickLabel',50:25:100)
ylabel('Classification accuracy (%)')
xlabel('Number of entities')
title('Decode for sample odor')
legend([hsu,hfc],{'Neuron FR','FCSP rate'})
