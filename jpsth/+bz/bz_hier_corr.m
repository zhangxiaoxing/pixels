keyboard()
%V.S. Frac, V.S. COM V.S. PV/SST

non_stats=onestat('nonmem',5); % 4:5,PV-SST; 6:7,COM; 8:9, frac;
congru_stats=onestat('congru',5);

congru_exp={'Source','Target','Weight'};
non_exp=congru_exp;
for fi=1:size(congru_stats.meta,1)
    pre=congru_stats.reg{fi,1};
    post=congru_stats.reg{fi,2};
    congru_exp=[congru_exp;...
        {pre},{post},{congru_stats.meta(fi,1)/sum(congru_stats.meta(fi,1:2))};...
        {post},{pre},{congru_stats.meta(fi,2)/sum(congru_stats.meta(fi,1:2))}];
end
writecell(congru_exp,'Congru_asym_fc.csv')
ureg=unique(congru_exp(2:end,1:2));
ffrac.collection=ephys.per_region_fraction('memtype','any'); % *100
[~,~,ratiomap]=ref.get_pv_sst();
node_exp={'Id','Label','XX','YY'};
for ni=1:numel(ureg)
    node_exp=[node_exp;...
        ureg(ni),ureg(ni),{ratiomap(ureg{ni})},1-ffrac.collection{strcmp(ffrac.collection(:,2),ureg(ni)),1}];
end
writecell(node_exp,'Congru_asym_node.csv')

% corrmat=[];
% figure('Color','w','Position',[32,32,400,800])
% hold on;
% for fi=1:size(congru_stats.meta,1)
%     pre=congru_stats.reg{fi,1};
%     post=congru_stats.reg{fi,2};
%     fwdsel=strcmp(non_stats.reg(:,1),pre)&strcmp(non_stats.reg(:,2),post);
%     revsel=strcmp(non_stats.reg(:,2),pre)&strcmp(non_stats.reg(:,1),post);
%     if any(fwdsel)
%         corrmat=[corrmat;congru_stats.meta(fi,3),non_stats.meta(fwdsel,3)];
%         scatter(congru_stats.meta(fi,3),non_stats.meta(fwdsel,3),4,'ro');
%         text(congru_stats.meta(fi,3),non_stats.meta(fwdsel,3),strjoin(congru_stats.reg(fi,:),'-'),'HorizontalAlignment','center','VerticalAlignment','top');
%     elseif any(revsel)
%         corrmat=[corrmat;congru_stats.meta(fi,3),-non_stats.meta(revsel,3)];
%         scatter(congru_stats.meta(fi,3),-non_stats.meta(revsel,3),4,'ro');
%         text(congru_stats.meta(fi,3),-non_stats.meta(revsel,3),strjoin(congru_stats.reg(fi,:),'-'),'HorizontalAlignment','center','VerticalAlignment','top');
%     end
% end
% plot([0,1],[0,1],'--k')
% plot([0,1],[0,0],'--k')
% plot([0,1],[0,-1],'--k')
% xlim([0,1])
% ylim([-1,1])
% xlabel('Congruent direction index')
% ylabel('Non-memory direction index')
% [r,p]=corr(corrmat(:,1),corrmat(:,2));

fh1=plot_one(congru_stats,4:5,'PV/(PV+SST)',100);
fh2=plot_one(congru_stats,6:7,'Firing rate C.O.M.',0.25);
fh3=plot_one(congru_stats,8:9,'Memory fraction',100);
exportgraphics(fh1,'congru_fc_corr_pvsst.pdf')
exportgraphics(fh2,'congru_fc_corr_COM.pdf')
exportgraphics(fh3,'congru_fc_corr_Frac.pdf')

fh1n=plot_one(non_stats,4:5,'PV/(PV+SST)',100);
fh2n=plot_one(non_stats,6:7,'Firing rate C.O.M.',0.25);
fh3n=plot_one(non_stats,8:9,'Memory fraction',100);

function fh=plot_one(stats,idx,ylbl,k)
arguments
    stats (1,1) struct
    idx (1,2) double
    ylbl (1,:) char
    k (1,1) double
end
fh=figure('Color','w','Position',[32,32,135,235]);
hold on
plot(repmat([1;2],1,size(stats.meta,1)),...
    stats.meta(:,idx).'.*k,...
    '-k.')
mm=mean(stats.meta(:,idx)).*k;
ci=bootci(1000,@(x) mean(x).*k,stats.meta(:,idx));
plot([0.5,2.5],mm,'ko');
errorbar([0.5,2.5],mm,ci(1,:)-mm,ci(2,:)-mm,'k.')
p=signrank(stats.meta(:,idx(1)),stats.meta(:,idx(2)));
text(1.5,max(ylim()),sprintf('%.3f',p),'HorizontalAlignment','center','VerticalAlignment','bottom');
xlim([0,3])
ylabel(ylbl)
set(gca(),'XTick',1:2,'XTickLabel',{'Pre-region','Post-region'},'XTickLabelRotation',90)

end


function stats=onestat(type,minconn)
[~,~,ratiomap]=ref.get_pv_sst();
[sig,pair]=bz.load_sig_pair('pair',true);
idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));
fcom=load('per_region_com_collection.mat','collection'); % /4
ffrac.collection=ephys.per_region_fraction('memtype','any'); % *100
if strcmp(type,'congru')
    mem_sel=all(ismember(sig.mem_type,1:2),2) | all(ismember(sig.mem_type,3:4),2);
    pmem_sel=all(ismember(pair.mem_type,1:2),2) | all(ismember(pair.mem_type,3:4),2);
%     ftitle='Congruent memory';
elseif strcmp(type,'nonmem')
    mem_sel=all(sig.mem_type==0,2);
    pmem_sel=all(pair.mem_type==0,2);
%     ftitle='Non-memory';
end
diff_sel=(sig.reg(:,5,1)~=sig.reg(:,5,2));
pdiff_sel=(pair.reg(:,5,1)~=pair.reg(:,5,2));
ctx_sel=all(squeeze(sig.reg(:,2,:)==688),2) & all(squeeze(sig.reg(:,5,:)>0),2);
pctx_sel=all(squeeze(pair.reg(:,2,:)==688),2) & all(squeeze(pair.reg(:,5,:)>0),2);

hier_sel=diff_sel & ctx_sel & mem_sel;
phier_sel=pdiff_sel & pctx_sel & pmem_sel;

ureg=unique(sig.reg(hier_sel,5,:));
cross_reg=nchoosek(ureg,2);
hreg=squeeze(sig.reg(hier_sel,5,:));
preg=squeeze(pair.reg(phier_sel,5,:));

stats=struct();
stats.reg=cell(0);
stats.meta=[];
for ri=1:size(cross_reg,1)
    nfwd=nnz(hreg(:,1)==cross_reg(ri,1) & hreg(:,2)==cross_reg(ri,2));
    nbck=nnz(hreg(:,1)==cross_reg(ri,2) & hreg(:,2)==cross_reg(ri,1));
    pcnt=nnz(preg(:,1)==cross_reg(ri,1) & preg(:,2)==cross_reg(ri,2))*2;
    if pcnt>100 && (nfwd + nbck)>minconn
        reg1=char(idmap.ccfid2reg(cross_reg(ri,1)));
        reg2=char(idmap.ccfid2reg(cross_reg(ri,2)));
        if ~any(strcmp(fcom.collection(:,2),reg1)) || ~any(strcmp(fcom.collection(:,2),reg2)) ...
                || ~any(strcmp(ffrac.collection(:,2),reg1)) || ~any(strcmp(ffrac.collection(:,2),reg2))
            warning('Missing meta data :%s, %s',reg1,reg2);
            continue
        end
            
        sidx=(nfwd-nbck)./(nfwd+nbck);
        if sidx>=0
            stats.reg=[stats.reg;{reg1},{reg2}];
            %FRAC,COM,PV
            stats.meta=[stats.meta;...
                nfwd,nbck,sidx,...
                ratiomap(reg1),ratiomap(reg2),...
                fcom.collection{strcmp(fcom.collection(:,2),reg1),1},...
                fcom.collection{strcmp(fcom.collection(:,2),reg2),1},...
                ffrac.collection{strcmp(ffrac.collection(:,2),reg1),1},...
                ffrac.collection{strcmp(ffrac.collection(:,2),reg2),1},...
                pcnt];
        else
            stats.reg=[stats.reg;{reg2},{reg1}];
            %FRAC,COM,PV
            stats.meta=[stats.meta;...
                nbck,nfwd,-sidx,...
                ratiomap(reg2),ratiomap(reg1),...
                fcom.collection{strcmp(fcom.collection(:,2),reg2),1},...
                fcom.collection{strcmp(fcom.collection(:,2),reg1),1}...
                ffrac.collection{strcmp(ffrac.collection(:,2),reg2),1},...
                ffrac.collection{strcmp(ffrac.collection(:,2),reg1),1},...
                pcnt];
        end
    end
end
end
