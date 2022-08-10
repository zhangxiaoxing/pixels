function [fc_com_pvsst_stats,fh]=fc_com_pvsst(com_map_1half,com_map_2half,selmeta,opt)
arguments
    com_map_1half
    com_map_2half
    selmeta
    opt.level (1,1) double {mustBeInteger} = 5 %allen ccf structure detail level
    opt.decision (1,1) logical = false % return statistics of decision period, default is delay period
    opt.hiermap (1,:) char 
    opt.mem_type (1,:) char {mustBeMember(opt.mem_type,{'congru','incong','any','indep','mixed'})} = 'congru' %TCOM for NONMEM undefined
    opt.reg_min_su (1,1) double = 100
    opt.tbl_title (1,:) char = []
end
[sig,~]=bz.load_sig_sums_conn_file();
% sig=bz.load_sig_pair('type','neupix','prefix','BZWT','criteria','WT','load_waveid',false);
sig=bz.join_fc_waveid(sig,selmeta.wave_id);
% same sort as file-list, not same to sessionid
[~,is_same_unsort,h2l_unsort,l2h_unsort]=bz.util.diff_at_level(sig.reg,'hierarchy',true,'range','grey','hiermap',opt.hiermap,'mincount',opt.reg_min_su,'descend',true);
disp(nnz(is_same_unsort | l2h_unsort | h2l_unsort))
fc_com_pvsst_stats=[];
usess=unique(sig.sess);
[is_same,h2l,l2h]=deal([]);

for sii=reshape(usess,1,[]) %iterate through sessions
    sesssel=sig.sess==sii;
    is_same=[is_same;is_same_unsort(sesssel,:)];
    h2l=[h2l;h2l_unsort(sesssel,:)];
    l2h=[l2h;l2h_unsort(sesssel,:)];

    suid=sig.suid(sesssel,:);
    waveid=sig.waveid(sesssel,:);
    com_sess_2half=nan(size(suid));
    if isfield(com_map_2half,['s',num2str(sii)])
        s1keys=int32(cell2mat(com_map_2half.(['s',num2str(sii)]).c1.keys())); % s1 prefered SUid
        s1sel=ismember(suid,s1keys);% same dim as suid
        com_sess_2half(s1sel)=cell2mat(com_map_2half.(['s',num2str(sii)]).c1.values(num2cell(suid(s1sel)))); % out put is nx2 in dim
        s2keys=int32(cell2mat(com_map_2half.(['s',num2str(sii)]).c2.keys())); % differernt SUs, e.g. prefer S2
        s2sel=ismember(suid,s2keys);
        com_sess_2half(s2sel)=cell2mat(com_map_2half.(['s',num2str(sii)]).c2.values(num2cell(suid(s2sel))));
    end
    com_sess_1half=nan(size(suid));
    if isfield(com_map_1half,['s',num2str(sii)])
        s1keys=int32(cell2mat(com_map_1half.(['s',num2str(sii)]).c1.keys()));
        s1sel=ismember(suid,s1keys);
        com_sess_1half(s1sel)=cell2mat(com_map_1half.(['s',num2str(sii)]).c1.values(num2cell(suid(s1sel))));
        s2keys=int32(cell2mat(com_map_1half.(['s',num2str(sii)]).c2.keys()));
        s2sel=ismember(suid,s2keys);
        com_sess_1half(s2sel)=cell2mat(com_map_1half.(['s',num2str(sii)]).c2.values(num2cell(suid(s2sel))));
    end

    regsess=squeeze(sig.reg(sesssel,5,:));
    fc_com_pvsst_stats=[fc_com_pvsst_stats;double(sii).*ones(size(suid(:,1))),double(suid),com_sess_2half,com_sess_1half,double(regsess),double(waveid)];
    %==================================================sess=====================suid=======COM_context1=====COM_context2=====ccfid===========waveid======
end
% save('fc_com_pvsst_stats.mat','fc_com_pvsst_stats');

switch opt.mem_type
    case 'indep'
        congrusel=all(fc_com_pvsst_stats(:,10:11)==5,2) ...
        | all(fc_com_pvsst_stats(:,10:11)==6,2);

    case 'mixed'
        congrusel=all(fc_com_pvsst_stats(:,10:11)==1,2) ...
            | all(fc_com_pvsst_stats(:,10:11)==2,2) ...
            | all(fc_com_pvsst_stats(:,10:11)==3,2) ...
            | all(fc_com_pvsst_stats(:,10:11)==4,2);

    case 'congru'
        congrusel=all(ismember(fc_com_pvsst_stats(:,10:11),[1 5]),2) ...
        | all(ismember(fc_com_pvsst_stats(:,10:11),[3 5]),2) ...
        | all(ismember(fc_com_pvsst_stats(:,10:11),[2 6]),2) ...
        | all(ismember(fc_com_pvsst_stats(:,10:11),[4 6]),2);
    case 'incong'
        congrusel=any(ismember(fc_com_pvsst_stats(:,10:11),[1 3 5]),2) & any(ismember(fc_com_pvsst_stats(:,10:11),[2 4 6]),2);
end

fini_sel=@(x) x(all(isfinite(x),2),:);
fork_join=@(x) [x(:,1:2);x(:,3:4)];

% >>>>>>>>>>>>>>>>>>>>>population FC-TCOM comparison>>>>>>>>>
    congru_tcom=fini_sel(fork_join(fc_com_pvsst_stats(congrusel,4:7)));
    ps_total=signrank(congru_tcom(:,1),congru_tcom(:,2));
    [~,pt_total]=ttest(congru_tcom(:,1),congru_tcom(:,2));
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


% COM vs. hier
comsame=(fc_com_pvsst_stats(is_same(:,opt.level) & congrusel,4:7));
coml2h=(fc_com_pvsst_stats(l2h(:,opt.level)  & congrusel,4:7));
comh2l=(fc_com_pvsst_stats(h2l(:,opt.level)  & congrusel,4:7));

%>>>>>>> per subgroup TCOM comparison>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ps_same=signrank(diff(fini_sel(fork_join(comsame)),1,2));
[~,pt_same]=ttest(diff(fini_sel(fork_join(comsame)),1,2));
ps_l2h=signrank(diff(fini_sel(fork_join(coml2h)),1,2));
[~,pt_l2h]=ttest(diff(fini_sel(fork_join(coml2h)),1,2));
ps_h2l=signrank(diff(fini_sel(fork_join(comh2l)),1,2));
[~,pt_h2l]=ttest(diff(fini_sel(fork_join(comh2l)),1,2));
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%>>>>>>> per subgroup ration comparison <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
[mmsame,sameci]=binofit( ...
    nnz(diff(fini_sel(fork_join(comsame)),1,2)>0), ...
    size(fini_sel(fork_join(comsame)),1));
[mml2h,l2hci]=binofit( ...
    nnz(diff(fini_sel(fork_join(coml2h)),1,2)>0), ...
    size(fini_sel(fork_join(coml2h)),1));
[mmh2l,h2lci]=binofit( ...
    nnz(diff(fini_sel(fork_join(comh2l)),1,2)>0), ...
    size(fini_sel(fork_join(comh2l)),1));
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


fh=figure('Color','w','Position',[100,100,720,320]);
tiledlayout(1,3)
nexttile(1);
hold on;
bh=bar([mmsame,1-mmsame;mml2h,1-mml2h;mmh2l,1-mmh2l]);
% bh=bar([mml2h,1-mml2h;mmh2l,1-mmh2l]);

[bh(1).FaceColor,bh(1).EdgeColor,bh(2).EdgeColor]=deal('k');
bh(2).FaceColor='w';
errorbar(bh(1).XEndPoints,...
    bh(1).YEndPoints,...
    [sameci(1),l2hci(1),h2lci(1)]-bh(1).YEndPoints,...
    [sameci(2),l2hci(2),h2lci(2)]-bh(1).YEndPoints,...
    'k.')
errorbar(bh(2).XEndPoints,...
    bh(2).YEndPoints,...
    1-[sameci(1),l2hci(1),h2lci(1)]-bh(2).YEndPoints,...
    1-[sameci(2),l2hci(2),h2lci(2)]-bh(2).YEndPoints,...
    'k.')

legend({'Lead neuron active earlier','Follow neuron activate earlier'},'Location','northoutside');
ylabel('Proportion of all F.C. (%)');
set(gca(),'XTick',1:3,'XTickLabel',{'Same','Along','Against'},'XTickLabelRotation',0);

th=nexttile(2,[1,2]);

s=([nnz([comsame(:,2);comsame(:,4)]>[comsame(:,1);comsame(:,3)]),nnz(all(isfinite([comsame(:,1:2);comsame(:,3:4)]),2))]);
l=([nnz([coml2h(:,2);coml2h(:,4)]>[coml2h(:,1);coml2h(:,3)]),nnz(all(isfinite([coml2h(:,1:2);coml2h(:,3:4)]),2))]);
h=([nnz([comh2l(:,2);comh2l(:,4)]>[comh2l(:,1);comh2l(:,3)]),nnz(all(isfinite([comh2l(:,1:2);comh2l(:,3:4)]),2))]);

[tbl,chi2,p]=crosstab([zeros(s(2),1);ones(l(2),1);2*ones(h(2),1)],[(1:s(2))>s(1),(1:l(2))>l(1),(1:h(2))>h(1)].');
binocdfp=[2*binocdf(min(tbl(1,:)),sum(tbl(1,:)),0.5),2*binocdf(min(tbl(2,:)),sum(tbl(2,:)),0.5),2*binocdf(min(tbl(3,:)),sum(tbl(3,:)),0.5),];

% binocdf(min(135+144),sum([135 279 144 208]),0.5)



binocdfptotal=(2*binocdf(min(sum(tbl)),sum(tbl,'all'),0.5));

statstbl={['signrank p total ', num2str(ps_total),' ttest p total ',num2str(pt_total)]; ...
    ['signrank p same ',num2str(ps_same),' ttest p same ',num2str(pt_same)]; ...
    ['signrank p low2high ',num2str(ps_l2h),' ttest p low2high ',num2str(pt_l2h)]; ...
    ['signrank p high2low ',num2str(ps_h2l),' ttest p high2low ',num2str(pt_h2l)]; ...
    ['same count ',num2str(tbl(1,:))]; ...
    ['l2h count ',num2str(tbl(2,:))]; ...
    ['h2l count ',num2str(tbl(3,:))]; ...
    ['chisq ',num2str(chi2),' chisq p ',num2str(p)];...
    ['bino cdf p total ', num2str(binocdfptotal)];...
    'bino cdf p same/l2h/h2l';...
    num2str(binocdfp) ....
    };

ephys.util.figtable(fh,th,statstbl,'title',opt.tbl_title)

% exportgraphics(fh,'fc_prog_regres_bars_hier.pdf')
end

function com_map=list2map(com_meta)
com_map=struct();
for ii=1:size(com_meta,1)
    if ~isfield(com_map,['s',num2str(com_meta{ii,1})])
        com_map.(['s',num2str(com_meta{ii,1})])=containers.Map('KeyType','int32','ValueType','any');
    end
    com_map.(['s',num2str(com_meta{ii,1})])(com_meta{ii,2})=com_meta(ii,3:end);
end
end