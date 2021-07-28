function [comdiff_stats,com_pair]=fc_com_pvsst(opt)
arguments
    opt.level (1,1) double {mustBeInteger} = 5 %allen ccf structure detail level
    opt.decision (1,1) logical = false % return statistics of decision period, default is delay period
end

[~,com_meta]=wave.per_region_COM('decision',opt.decision);
com_map=list2map(com_meta);
[~,~,ratiomap]=ref.get_pv_sst();
sig=bz.load_sig_pair('type','neupix','prefix','BZWT','criteria','WT');
[is_diff,is_same]=bz.util.diff_at_level(sig.reg);
fc_com_pvsst_stats=[];

for ii=1:size(sig.sess,1)
    if ~all(ismember(sig.mem_type,[2 4]),2), continue;end %only concerns congruent pairs
    if is_same(ii,opt.level) || is_diff(ii,opt.level)%same
        sess=sig.sess(ii);
        suid=sig.suid(ii,:);
        if ~all(arrayfun(@(x) com_map.(['s',num2str(sess)]).isKey(x),suid))
%             warning('Missing keys in COM meta data, session %d, suid %d, %d',sess,suid(1),suid(2));
            continue
        end
        com_meta_one{1}=com_map.(['s',num2str(sess)])(suid(1));
        com_meta_one{2}=com_map.(['s',num2str(sess)])(suid(2));
        [com1,com2]=deal(com_meta_one{1}{1},com_meta_one{2}{1});
        [reg1,reg2]=deal(com_meta_one{1}{6},com_meta_one{2}{6});
        if ~ratiomap.isKey(reg1)||~ratiomap.isKey(reg2), continue;end
        [ratio1,ratio2]=deal(ratiomap(reg1),ratiomap(reg2));
        fc_com_pvsst_stats=[fc_com_pvsst_stats;{is_diff(ii,opt.level),sess,suid(1),suid(2),com1,com2,reg1,reg2,ratio1,ratio2}];
    end
end
% save('fc_com_pvsst_stats.mat','fc_com_pvsst_stats');

samesel=cell2mat(fc_com_pvsst_stats(:,1))==0;
low2high=cell2mat(fc_com_pvsst_stats(:,1))==1 & diff(cell2mat(fc_com_pvsst_stats(:,9:10)),1,2)<0;
high2low=cell2mat(fc_com_pvsst_stats(:,1))==1 & diff(cell2mat(fc_com_pvsst_stats(:,9:10)),1,2)>0;

comsame=cell2mat(fc_com_pvsst_stats(samesel,5:6));
coml2h=cell2mat(fc_com_pvsst_stats(low2high,5:6));
comh2l=cell2mat(fc_com_pvsst_stats(high2low,5:6));

mmsame=mean(diff(comsame,1,2)>0).*100;
mmh2l=mean(diff(comh2l,1,2)>0).*100;
mml2h=mean(diff(coml2h,1,2)>0).*100;
sameci=bootci(1000, @(x) nnz(diff(x,1,2)>0)./size(x,1), comsame).*100;
l2hci=bootci(1000, @(x) nnz(diff(x,1,2)>0)./size(x,1), coml2h).*100;
h2lci=bootci(1000, @(x) nnz(diff(x,1,2)>0)./size(x,1), comh2l).*100;

fh=figure('Color','w','Position',[100,100,225,235]);
hold on;

bh=bar([100-mmsame,mmsame;100-mmh2l,mmh2l;100-mml2h,mml2h]);
[bh(1).FaceColor,bh(1).EdgeColor,bh(2).EdgeColor]=deal('k');
bh(2).FaceColor='w';
errorbar(bh(1).XEndPoints,...
    bh(1).YEndPoints,...
    [100-sameci(1)-bh(1).YEndPoints(1),100-h2lci(1)-bh(1).YEndPoints(2),100-l2hci(1)-bh(1).YEndPoints(3)],...
    [100-sameci(2)-bh(1).YEndPoints(1),100-h2lci(2)-bh(1).YEndPoints(2),100-l2hci(2)-bh(1).YEndPoints(3)],...
    'k.')
errorbar(bh(2).XEndPoints,...
    bh(2).YEndPoints,...
    [sameci(1)-bh(2).YEndPoints(1),h2lci(1)-bh(2).YEndPoints(2),l2hci(1)-bh(2).YEndPoints(3)],...
    [sameci(2)-bh(2).YEndPoints(1),h2lci(2)-bh(2).YEndPoints(2),l2hci(2)-bh(2).YEndPoints(3)],...
    'k.')
legend({'Backward','Forward'},'Location','northoutside');
ylabel('Proportion of all F.C. (%)');
set(gca(),'XTick',1:3,'XTickLabel',{'Within reg.','High to low','Low to high'},'XTickLabelRotation',45);
 
keyboard()
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