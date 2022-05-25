function [fc_com_pvsst_stats,fh]=fc_com_pvsst(com_map_3,com_map_6,selmeta,opt)
arguments
    com_map_3
    com_map_6
    selmeta
    opt.level (1,1) double {mustBeInteger} = 5 %allen ccf structure detail level
    opt.decision (1,1) logical = false % return statistics of decision period, default is delay period
    opt.hiermap (1,:) char {mustBeMember(opt.hiermap,{'pvsst','OBM1','AON','CP'})} = 'AON'
    opt.mem_type (1,:) char {mustBeMember(opt.mem_type,{'congru','incong','any'})} = 'congru' %TCOM for NONMEM undefined
    opt.reg_min_su (1,1) double = 100
end


sig=bz.load_sig_pair('type','neupix','prefix','BZWT','criteria','WT','load_waveid',false);
sig=bz.join_fc_waveid(sig,selmeta.wave_id);
[~,is_same,h2l,l2h]=bz.util.diff_at_level(sig.reg,'hierarchy',true,'range','grey','hiermap',opt.hiermap,'mincount',opt.reg_min_su);
disp(nnz(is_same | l2h | h2l))
fc_com_pvsst_stats=[];
usess=unique(sig.sess);

for sii=reshape(usess,1,[]) %iterate through sessions
    sesssel=sig.sess==sii;
    suid=sig.suid(sesssel,:);

    waveid=sig.waveid(sesssel,:);
    comsess_6=nan(size(suid));
    if isfield(com_map_6,['s',num2str(sii)])
        s1keys=int32(cell2mat(com_map_6.(['s',num2str(sii)]).c1.keys())); % s1 prefered SUid
        s1sel=ismember(suid,s1keys);% same dim as suid
        comsess_6(s1sel)=cell2mat(com_map_6.(['s',num2str(sii)]).c1.values(num2cell(suid(s1sel)))); % out put is nx2 in dim
        s2keys=int32(cell2mat(com_map_6.(['s',num2str(sii)]).c2.keys())); % differernt SUs, e.g. prefer S2
        s2sel=ismember(suid,s2keys);
        comsess_6(s2sel)=cell2mat(com_map_6.(['s',num2str(sii)]).c2.values(num2cell(suid(s2sel))));
    end
    comsess_3=nan(size(suid));
    if isfield(com_map_3,['s',num2str(sii)])
        s1keys=int32(cell2mat(com_map_3.(['s',num2str(sii)]).c1.keys()));
        s1sel=ismember(suid,s1keys);
        comsess_3(s1sel)=cell2mat(com_map_3.(['s',num2str(sii)]).c1.values(num2cell(suid(s1sel))));
        s2keys=int32(cell2mat(com_map_3.(['s',num2str(sii)]).c2.keys()));
        s2sel=ismember(suid,s2keys);
        comsess_3(s2sel)=cell2mat(com_map_3.(['s',num2str(sii)]).c2.values(num2cell(suid(s2sel))));
    end

    regsess=squeeze(sig.reg(sesssel,5,:));
    fc_com_pvsst_stats=[fc_com_pvsst_stats;double(sii).*ones(size(suid(:,1))),double(suid),comsess_6,comsess_3,double(regsess),double(waveid)];
    %==================================================sess=====================suid=======COM_6s=====COM_3s=====ccfid===========waveid?======
end
save('fc_com_pvsst_stats.mat','fc_com_pvsst_stats');
switch opt.mem_type
    case 'congru'
        congrusel=all(ismember(fc_com_pvsst_stats(:,10:11),[1 3 5]),2) |all(ismember(fc_com_pvsst_stats(:,10:11),[2 4 6]),2);
    case 'incong'
        congrusel=any(ismember(fc_com_pvsst_stats(:,10:11),[1 3 5]),2) & any(ismember(fc_com_pvsst_stats(:,10:11),[2 4 6]),2);
end

% COM vs. hier
comsame=(fc_com_pvsst_stats(is_same(:,opt.level) & congrusel,4:7));
coml2h=(fc_com_pvsst_stats(l2h(:,opt.level)  & congrusel,4:7));
comh2l=(fc_com_pvsst_stats(h2l(:,opt.level)  & congrusel,4:7));

[mmsame,sameci]=binofit(nnz([comsame(:,2);comsame(:,4)]>[comsame(:,1);comsame(:,3)]),nnz(all(isfinite([comsame(:,1:2);comsame(:,3:4)]),2)));
[mml2h,l2hci]=binofit(nnz([coml2h(:,2);coml2h(:,4)]>[coml2h(:,1);coml2h(:,3)]),nnz(all(isfinite([coml2h(:,1:2);coml2h(:,3:4)]),2)));
[mmh2l,h2lci]=binofit(nnz([comh2l(:,2);comh2l(:,4)]>[comh2l(:,1);comh2l(:,3)]),nnz(all(isfinite([comh2l(:,1:2);comh2l(:,3:4)]),2)));

fh=figure('Color','w','Position',[100,100,225,235]);
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

legend({'Consistent','Inconsistent'},'Location','northoutside');
ylabel('Proportion of all F.C. (%)');
set(gca(),'XTick',1:3,'XTickLabel',{'Same','From-OLF','To-OLF'},'XTickLabelRotation',0);

s=([nnz([comsame(:,2);comsame(:,4)]>[comsame(:,1);comsame(:,3)]),nnz(all(isfinite([comsame(:,1:2);comsame(:,3:4)]),2))]);
l=([nnz([coml2h(:,2);coml2h(:,4)]>[coml2h(:,1);coml2h(:,3)]),nnz(all(isfinite([coml2h(:,1:2);coml2h(:,3:4)]),2))]);
h=([nnz([comh2l(:,2);comh2l(:,4)]>[comh2l(:,1);comh2l(:,3)]),nnz(all(isfinite([comh2l(:,1:2);comh2l(:,3:4)]),2))]);

[tbl,chi2,p]=crosstab([zeros(s(2),1);ones(l(2),1);2*ones(h(2),1)],[(1:s(2))>s(1),(1:l(2))>l(1),(1:h(2))>h(1)].')
disp([2*binocdf(min(tbl(1,:)),sum(tbl(1,:)),0.5),2*binocdf(min(tbl(2,:)),sum(tbl(2,:)),0.5),2*binocdf(min(tbl(3,:)),sum(tbl(3,:)),0.5),]);

disp(2*binocdf(min(sum(tbl)),sum(tbl,'all'),0.5))

% 
% import scipy.stats as stats
% stats.binom_test(745,1318,0.5) %binocdf
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