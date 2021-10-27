function [comdiff_stats,com_pair]=fc_com_pvsst(opt)
arguments
    opt.level (1,1) double {mustBeInteger} = 5 %allen ccf structure detail level
    opt.decision (1,1) logical = false % return statistics of decision period, default is delay period
    opt.pvsst (1,1) logical = false % use MOB-MOp axis index since 20210928
    opt.mem_type (1,:) char {mustBeMember(opt.mem_type,{'congru','incong'})} = 'incong' %TCOM for NONMEM undefined
end

[~,com_meta]=wave.per_region_COM('decision',opt.decision);
com_map=list2map(com_meta);
if opt.pvsst
    [~,~,ratiomap]=ref.get_pv_sst();
else
    load('OBM1map.mat','OBM1map');
end
sig=bz.load_sig_pair('type','neupix','prefix','BZWT','criteria','WT');
[~,is_same,h2l,l2h]=bz.util.diff_at_level(sig.reg,'hierarchy',true);

fc_com_pvsst_stats=[];
for ii=1:numel(sig.sess) %iterate through all pairs
    if is_same(ii,opt.level) || h2l(ii,opt.level) || l2h(ii,opt.level)%same
        sess=sig.sess(ii);
        suid=sig.suid(ii,:);
        mem_type=sig.mem_type(ii,:);
        if strcmp(opt.mem_type,'congru') && ~all(ismember(mem_type,1:2)) && ~all(ismember(mem_type,3:4))
            continue
        end
        if strcmp(opt.mem_type,'incong') && (~any(ismember(mem_type,1:2)) || ~any(ismember(mem_type,3:4)))
            continue
        end

        if ~all(arrayfun(@(x) com_map.(['s',num2str(sess)]).isKey(x),suid))
%             warning('Missing keys in COM meta data, session %d, suid %d, %d',sess,suid(1) ,suid(2));
            continue
        end
%         keyboard()
        com_meta_one{1}=com_map.(['s',num2str(sess)])(suid(1));
        com_meta_one{2}=com_map.(['s',num2str(sess)])(suid(2));
        [com1,com2]=deal(com_meta_one{1}{1},com_meta_one{2}{1});
        [reg1,reg2]=deal(com_meta_one{1}{6},com_meta_one{2}{6});
        if opt.pvsst
            if (~ratiomap.isKey(reg1)||~ratiomap.isKey(reg2)), continue;end
            [ratio1,ratio2]=deal(ratiomap(reg1),ratiomap(reg2));
        else
            if (~OBM1map.isKey(reg1)||~OBM1map.isKey(reg2)), continue;end
            [ratio1,ratio2]=deal(OBM1map(reg1),OBM1map(reg2));
        end
        fc_com_pvsst_stats=[fc_com_pvsst_stats;{is_same(ii,opt.level),l2h(ii,opt.level),h2l(ii,opt.level),sess,suid(1),suid(2),com1,com2,reg1,reg2,ratio1,ratio2}];
    end
end
save('fc_com_pvsst_stats.mat','fc_com_pvsst_stats');

samesel=cell2mat(fc_com_pvsst_stats(:,1))==1;
low2high=cell2mat(fc_com_pvsst_stats(:,2))==1;
high2low=cell2mat(fc_com_pvsst_stats(:,3))==1;

comsame=cell2mat(fc_com_pvsst_stats(samesel,7:8));
coml2h=cell2mat(fc_com_pvsst_stats(low2high,7:8));
comh2l=cell2mat(fc_com_pvsst_stats(high2low,7:8));

mmsame=mean(diff(comsame,1,2)>0).*100;
mmh2l=mean(diff(comh2l,1,2)>0).*100;
mml2h=mean(diff(coml2h,1,2)>0).*100;
sameci=bootci(1000, @(x) nnz(diff(x,1,2)>0)./size(x,1), comsame).*100;
l2hci=bootci(1000, @(x) nnz(diff(x,1,2)>0)./size(x,1), coml2h).*100;
h2lci=bootci(1000, @(x) nnz(diff(x,1,2)>0)./size(x,1), comh2l).*100;

fh=figure('Color','w','Position',[100,100,225,235]);
hold on;
bh=bar([mmsame,100-mmsame;mml2h,100-mml2h;mmh2l,100-mmh2l]);
[bh(1).FaceColor,bh(1).EdgeColor,bh(2).EdgeColor]=deal('k');
bh(2).FaceColor='w';
errorbar(bh(1).XEndPoints,...
    bh(1).YEndPoints,...
    [sameci(1)-bh(1).YEndPoints(1),l2hci(1)-bh(1).YEndPoints(2),h2lci(1)-bh(1).YEndPoints(3)],...
    [sameci(2)-bh(1).YEndPoints(1),l2hci(2)-bh(1).YEndPoints(2),h2lci(2)-bh(1).YEndPoints(3)],...
    'k.')
errorbar(bh(2).XEndPoints,...
    bh(2).YEndPoints,...
    [100-sameci(1)-bh(2).YEndPoints(1),100-l2hci(1)-bh(2).YEndPoints(2),100-h2lci(1)-bh(2).YEndPoints(3)],...
    [100-sameci(2)-bh(2).YEndPoints(1),100-l2hci(2)-bh(2).YEndPoints(2),100-h2lci(2)-bh(2).YEndPoints(3)],...
    'k.')
legend({'Progressive','Regressive'},'Location','northoutside');
ylabel('Proportion of all F.C. (%)');
set(gca(),'XTick',1:3,'XTickLabel',{'Within reg.','Low to high','High to Low'},'XTickLabelRotation',30);

disp([nnz(diff(comsame,1,2)>0),size(comsame)])
disp([nnz(diff(coml2h,1,2)>0),size(coml2h)])
disp([nnz(diff(comh2l,1,2)>0),size(comh2l)])
keyboard()
% 
% import scipy.stats as stats
% stats.binom_test(745,1318,0.5)
exportgraphics(fh,'fc_prog_regres_bars_hier.pdf')
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