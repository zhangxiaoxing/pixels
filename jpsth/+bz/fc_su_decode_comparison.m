function fc_su_decode_comparison()
load(fullfile('binary','wrs_mux_meta.mat'),'wrs_mux_meta');
sel_meta=wrs_mux_meta;
load(fullfile('binary','su_meta.mat'));

for nsu=[5 10 50 100 200 300 400 500]
    o3=pct.pct_decoding_correct_error(su_meta,sel_meta,[1 3 5 6],lblidx=5,rpt=25,n_su=nsu,delay=3);
    o6=pct.pct_decoding_correct_error(su_meta,sel_meta,[2 4 5 6],lblidx=5,rpt=25,n_su=nsu,delay=6);
    odor.("n"+nsu).olf=cell2struct({[o3.olf.("c_result_"+nsu+"su");o6.olf.("c_result_"+nsu+"su")];...
        [o3.olf.("e_result_"+nsu+"su");o6.olf.("e_result_"+nsu+"su")]},...
        {sprintf('c_result_%dsu',nsu),sprintf('e_result_%dsu',nsu)});
    % both_odor.("n"+nsu)=pct.pct_decoding_correct_error(wrs_mux_meta,1:4,'lblidx',5,'n_su',nsu);% odor
    % both_dur.("n"+nsu)=pct.pct_decoding_correct_error(wrs_mux_meta,1:4,'lblidx',8,'n_su',nsu);% dur
end
fc_dec=bz.fccoding.plot_coding(wrs_mux_meta,'nfc_grp',[5 10 50 100 200 300 400 500],'skip_plot',true);
if false
    for nsu=[5 10 50 100 200 300 400 500]
        dec.CSP.("N"+nsu)=fc_dec.("FC"+nsu).cvcorr{1};
        dec.neuron.("N"+nsu)=odor.("n"+nsu).olf.("c_result_"+nsu+"su");
    end
    fid=fopen(fullfile('binary','upload','SF4C_decoding_with_CSP_or_neuron.json'),'w');
    fprintf(fid,jsonencode(dec));
    fclose(fid);
end



odormat=cell2mat(cellfun(@(x) subsref(struct2cell(x.olf),substruct('{}',{1})),struct2cell(odor),'UniformOutput',false).');
summ=mean(odormat);
susem=std(odormat)./sqrt(size(odormat,1));

fcmat=cell2mat(cellfun(@(x) x.cvcorr{1},struct2cell(fc_dec),'UniformOutput',false).');
fcmm=mean(fcmat);
fsem=std(fcmat)./sqrt(size(fcmat,1));

nsu=[5 10 50 100 200 300 400 500];

fh=figure();
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
savefig(fullfile('binary','FC_SU_decode_comparison.fig'));