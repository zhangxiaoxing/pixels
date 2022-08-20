function [fc_com_pvsst_stats,fh]=fc_com_pct(com_map,pct_meta,opt)
arguments
    com_map
    pct_meta
    opt.level (1,1) double {mustBeInteger} = 5 %allen ccf structure detail level
    opt.decision (1,1) logical = false % return statistics of decision period, default is delay period
    opt.pair_type (1,:) char {mustBeMember(opt.pair_type,{'congru','incong'})} = 'congru' %TCOM for NONMEM undefined
    opt.reg_min_su (1,1) double = 100
    opt.tbl_title (1,:) char = []
    opt.hier_map
    opt.descend
end
[sig,~]=bz.load_sig_sums_conn_file();

sig=bz.join_fc_waveid(sig,pct_meta.wave_id);
% same order as file-list, different to sessionid
% TODO:
if isfield(opt,'hier_map')
    [~,is_same_unsort,h2l_unsort,l2h_unsort]=bz.util.diff_at_level(sig.reg,'hierarchy',true,'range','grey','hiermap',opt.hier_map,'mincount',opt.reg_min_su,'descend',opt.descend);
else
    [is_diff_unsort,is_same_unsort]=bz.util.diff_at_level(sig.reg,'hierarchy',false,'range','grey','mincount',opt.reg_min_su);
end
disp(nnz(is_same_unsort | is_diff_unsort));
fc_com_pvsst_stats=[];
usess=unique(sig.sess);
[is_same,is_diff]=deal([]);

for sii=reshape(usess,1,[]) %iterate through sessions
    sesssel=sig.sess==sii;
    is_same=[is_same;is_same_unsort(sesssel,:)];
    is_diff=[is_diff;is_diff_unsort(sesssel,:)];
    
    suid=sig.suid(sesssel,:);
    waveid=sig.waveid(sesssel,:);
    regsess=squeeze(sig.reg(sesssel,5,:));

    com_sess_pct=nan(size(suid));
    for ff=["s1d3","s1d6","s2d3","s2d6","olf_s1","olf_s2","dur_d3","dur_d6"]
        sukeys=com_map.(['s',num2str(sii)]).(ff).com.keys(); % prefered SUid
        susel=ismember(suid,int32(cell2mat(sukeys)));% same dim as suid
        com_sess_pct(susel)=cell2mat(com_map.(['s',num2str(sii)]).(ff).com.values(num2cell(suid(susel)))); % out put is nx2 in dim
    end
    fc_com_pvsst_stats=[fc_com_pvsst_stats;double(sii).*ones(size(suid(:,1))),double(suid),com_sess_pct,nan(size(suid)),double(regsess),double(waveid)];
    %==================================================sess=====================suid=======COM_context1=====COM_context2=====ccfid===========waveid======
end

switch opt.pair_type
%     case 'nonmem'
%         congrusel=pct.su_pairs.get_nonmem(fc_com_pvsst_stats(:,10:11));
    case 'congru'
        congrusel=pct.su_pairs.get_congru(fc_com_pvsst_stats(:,10:11));
    case 'incong'
        congrusel=pct.su_pairs.get_incongru(fc_com_pvsst_stats(:,10:11));
end

fini_sel=@(x) x(all(isfinite(x),2),:);


% >>>>>>>>>>>>>>>>>>>>>population FC-TCOM comparison>>>>>>>>>
    congru_tcom=fini_sel(fc_com_pvsst_stats(congrusel,4:5));
    ps_total=signrank(congru_tcom(:,1),congru_tcom(:,2));
    [~,pt_total]=ttest(congru_tcom(:,1),congru_tcom(:,2));
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


% COM vs. hier
comsame=(fc_com_pvsst_stats(is_same(:,opt.level) & congrusel,4:5));
comdiff=(fc_com_pvsst_stats(is_diff(:,opt.level) & congrusel,4:5));
% coml2h=(fc_com_pvsst_stats(l2h(:,opt.level)  & congrusel,4:7));
% comh2l=(fc_com_pvsst_stats(h2l(:,opt.level)  & congrusel,4:7));

%>>>>>>> per subgroup TCOM comparison>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ps_same=signrank(diff(fini_sel(comsame),1,2));
[~,pt_same]=ttest(diff(fini_sel(comsame),1,2));

ps_diff=signrank(diff(fini_sel(comdiff),1,2));
[~,pt_diff]=ttest(diff(fini_sel(comdiff),1,2));

% ps_l2h=signrank(diff(fini_sel(fork_join(coml2h)),1,2));
% [~,pt_l2h]=ttest(diff(fini_sel(fork_join(coml2h)),1,2));
% ps_h2l=signrank(diff(fini_sel(fork_join(comh2l)),1,2));
% [~,pt_h2l]=ttest(diff(fini_sel(fork_join(comh2l)),1,2));
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%>>>>>>> per subgroup ration comparison <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
[mmsame,sameci]=binofit( ...
    nnz(diff(fini_sel(comsame),1,2)>0), ...
    size(fini_sel(comsame),1));

[mmdiff,diffci]=binofit( ...
    nnz(diff(fini_sel(comdiff),1,2)>0), ...
    size(fini_sel(comdiff),1));


% [mml2h,l2hci]=binofit( ...
%     nnz(diff(fini_sel(fork_join(coml2h)),1,2)>0), ...
%     size(fini_sel(fork_join(coml2h)),1));
% [mmh2l,h2lci]=binofit( ...
%     nnz(diff(fini_sel(fork_join(comh2l)),1,2)>0), ...
%     size(fini_sel(fork_join(comh2l)),1));
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


fh=figure('Color','w','Position',[100,100,720,320]);
tiledlayout(1,3)
nexttile(1);
hold on;
bh=bar([mmsame,1-mmsame;mmdiff,1-mmdiff]);
% bh=bar([mmsame,1-mmsame;mml2h,1-mml2h;mmh2l,1-mmh2l]);
% bh=bar([mml2h,1-mml2h;mmh2l,1-mmh2l]);

[bh(1).FaceColor,bh(1).EdgeColor,bh(2).EdgeColor]=deal('k');
bh(2).FaceColor='w';
errorbar(bh(1).XEndPoints,...
    bh(1).YEndPoints,...
    [sameci(1),diffci(1)]-bh(1).YEndPoints,...
    [sameci(2),diffci(2)]-bh(1).YEndPoints,...
    'k.')
errorbar(bh(2).XEndPoints,...
    bh(2).YEndPoints,...
    1-[sameci(1),diffci(1)]-bh(2).YEndPoints,...
    1-[sameci(2),diffci(2)]-bh(2).YEndPoints,...
    'k.')

legend({'Lead neuron active earlier','Follow neuron activate earlier'},'Location','northoutside');
ylabel('Proportion of all F.C. (%)');
set(gca(),'XTick',1:3,'XTickLabel',{'Same','Diff'},'XTickLabelRotation',0);

th=nexttile(2,[1,2]);

s=([nnz(comsame(:,2)>comsame(:,1)),nnz(all(isfinite(comsame(:,1:2)),2))]);
d=([nnz(comdiff(:,2)>comdiff(:,1)),nnz(all(isfinite(comdiff(:,1:2)),2))]);
% l=([nnz([coml2h(:,2);coml2h(:,4)]>[coml2h(:,1);coml2h(:,3)]),nnz(all(isfinite([coml2h(:,1:2);coml2h(:,3:4)]),2))]);
% h=([nnz([comh2l(:,2);comh2l(:,4)]>[comh2l(:,1);comh2l(:,3)]),nnz(all(isfinite([comh2l(:,1:2);comh2l(:,3:4)]),2))]);
% [tbl,chi2,p]=crosstab([zeros(s(2),1);ones(l(2),1);2*ones(h(2),1)],[(1:s(2))>s(1),(1:l(2))>l(1),(1:h(2))>h(1)].');

[tbl,chi2,p]=crosstab([zeros(s(2),1);ones(d(2),1)],[(1:s(2))>s(1),(1:d(2))>d(1)].');

% binocdfp=[2*binocdf(min(tbl(1,:)),sum(tbl(1,:)),0.5),2*binocdf(min(tbl(2,:)),sum(tbl(2,:)),0.5),2*binocdf(min(tbl(3,:)),sum(tbl(3,:)),0.5),];
binocdfp=[2*binocdf(min(tbl(1,:)),sum(tbl(1,:)),0.5),2*binocdf(min(tbl(2,:)),sum(tbl(2,:)),0.5)];

% binocdf(min(135+144),sum([135 279 144 208]),0.5)

binocdfptotal=(2*binocdf(min(sum(tbl)),sum(tbl,'all'),0.5));

statstbl={['signrank p total ', num2str(ps_total),' ttest p total ',num2str(pt_total)]; ...
    ['signrank p same ',num2str(ps_same),' ttest p same ',num2str(pt_same)]; ...
    ['signrank p diff ',num2str(ps_diff),' ttest p diff ',num2str(pt_diff)]; ...
%     ['signrank p high2low ',num2str(ps_h2l),' ttest p high2low ',num2str(pt_h2l)]; ...
    ['same count ',num2str(tbl(1,:))]; ...
    ['diff count ',num2str(tbl(2,:))]; ...
%     ['h2l count ',num2str(tbl(3,:))]; ...
    ['chisq ',num2str(chi2),' chisq p ',num2str(p)];...
    ['bino cdf p total ', num2str(binocdfptotal)];...
    'bino cdf p same/diff';...
    num2str(binocdfp) ....
    };

ephys.util.figtable(fh,th,statstbl,'title',opt.tbl_title)

% exportgraphics(fh,'fc_prog_regres_bars_hier.pdf')
end

