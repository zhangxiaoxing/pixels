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

rptdata=[o3.olf.c_result_50su,o3.olf.e_result_50su;o6.olf.c_result_50su,o6.olf.e_result_50su];
mm=mean(rptdata);
sem=sqrt(mm.*(1-mm)./size(rptdata,1));

[~,~,p]=crosstab([zeros(size(rptdata,1),1);ones(size(rptdata,1),1)],[rptdata(:,1);rptdata(:,2)]);


fh=figure('Color','w','Position',[100,100,300,240]);
hold on
bh=bar(mm);
errorbar(1:2,mm,sem,'k.','CapSize',12);
ylim([0.5,1]);
set(gca(),'XTick',1:2,...
    'XTickLabel',{'Correct','Error'},...
    'YTick',0.5:0.25:1,'YTickLabel',50:25:100)
ylabel('Classification accuracy');

title(sprintf('Odor chisq %.4f',p));

if ~opt.skip_save
    savefig(fh,fullfile("binary","su_decoding.fig"));
end
end