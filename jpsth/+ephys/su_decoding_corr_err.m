function fh=su_decoding_corr_err(su_meta,sel_meta,opt)
arguments
    su_meta = []
    sel_meta = []
    opt.skip_save (1,1) logical = false
end
if isempty(su_meta)
    load(fullfile('binary','su_meta.mat'))
end
if isempty(sel_meta)
    fstr=load(fullfile('binary','wrs_mux_meta.mat'));
    sel_meta=fstr.wrs_mux_meta;
    clear fstr
end
o3=pct.pct_decoding_correct_error(su_meta,sel_meta,[1 3 5 6],lblidx=5,rpt=25,n_su=50,delay=3);
o6=pct.pct_decoding_correct_error(su_meta,sel_meta,[2 4 5 6],lblidx=5,rpt=25,n_su=50,delay=6);
if ~opt.skip_save
    blame=vcs.blame();
    save(fullfile("binary","su_decoding.mat"),"o3","o6","blame");
end

mm=[mean([o3.olf.c_result_50su;o6.olf.c_result_50su]),mean([o3.olf.e_result_50su;o6.olf.e_result_50su])...
    ];
sem=sqrt(mm.*(1-mm)./numel([o3.olf.c_result_50su;o6.olf.c_result_50su]));


fh=figure('Color','w','Position',[100,100,300,240]);
hold on
bh=bar(mm);
errorbar(1:2,mm,sem,'k.');
ylim([0.5,1]);
set(gca(),'XTick',1:2,...
    'XTickLabel',{'Correct','Error'},...
    'YTick',0.5:0.25:1,'YTickLabel',50:25:100)
ylabel('Classification accuracy');
title('Odor only su dec');

if ~opt.skip_save
    savefig(fh,fullfile("binary","su_decoding.fig"));
end
end