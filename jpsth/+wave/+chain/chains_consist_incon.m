%%
function fh=chains_consist_incon(su_meta,sel_meta,opt)
arguments
    su_meta = []
    sel_meta = []
    opt.skip_save (1,1) logical = true
    opt.shuf (1,1) logical = false
    opt.non_mem (1,1) logical = false

end
cross_only=false;
if isempty(su_meta)
    load(fullfile('binary','su_meta.mat'))
end
if isempty(sel_meta)
    fstr=load(fullfile('binary','wrs_mux_meta.mat'));
    sel_meta=fstr.wrs_mux_meta;
    clear fstr
end

global_init;
reg_com_maps=wave.get_reg_com_maps(sel_meta);

chains_uf_all=wave.COM_chain_reg(su_meta,sel_meta,reg_com_maps,'cross_only',cross_only,'non_mem',opt.non_mem);
chains_uf_rev_all=wave.COM_chain_reg(su_meta,sel_meta,reg_com_maps,'reverse',true,'cross_only',cross_only,'non_mem',opt.non_mem);

usess=unique([chains_uf_all.sess;chains_uf_rev_all.sess]);
sums=[];

for ss=reshape(double(usess),1,[])
    fwdcnt=nnz(chains_uf_all.sess==ss & chains_uf_all.cross_reg);
    wtncnt=nnz(chains_uf_all.sess==ss & ~chains_uf_all.cross_reg);
    revcnt=nnz(chains_uf_rev_all.sess==ss & chains_uf_rev_all.cross_reg);

    ssr=[fwdcnt,revcnt,wtncnt]./sum([wtncnt;revcnt;fwdcnt]);
    sums=[sums;ss,fwdcnt,revcnt,wtncnt,ssr];
end

mm=mean(sums(:,[7,5,6])).*100;
sem=std(sums(:,[7,5,6]))./sqrt(size(sums,1)).*100;

fh=figure();
hold on
bh=bar(mm.*eye(3),'stacked');
errorbar(1:3,mm,sem,'k.')
bh(1).FaceColor=[0.5,0.5,0.5];
bh(2).FaceColor='k';
bh(3).FaceColor='w';
set(gca,'XTick',[]);
legend(bh,{'Within','Consistent','Inconsistent'},'Location','northoutside','Orientation','horizontal')

consistent_incon=ranksum(sums(:,5),sums(:,6));
% within_consis=ranksum(sums(:,5),sums(:,7));
% within_incon=ranksum(sums(:,6),sums(:,7));

title(sprintf('consis-incon%.4f',consistent_incon));
ylabel('Proportion of all chains (%)')
ylim([0,50])
if ~opt.skip_save
    % appendfig('tag','chain consistent inconsistent,chain_plots.m')
    if opt.non_mem
        savefig(fh,fullfile('binary','chains_consist_incon_nonmem.fig'))
    else
        savefig(fh,fullfile('binary','chains_consist_incon.fig'))
    end
end

if opt.shuf
    fwdcnt=nnz(chains_uf_all.cross_reg);
    wtncnt=nnz(~chains_uf_all.cross_reg);
    revcnt=nnz(chains_uf_rev_all.cross_reg);

    [shuf_fwdcnt,shuf_wtncnt,shuf_revcnt]=deal([]);
    load(fullfile('binary','bz_ring_shufs.mat'),'shufs');
    for shufidx=1:numel(shufs)
        chains_shuf_all=wave.COM_chain_reg(su_meta,sel_meta,reg_com_maps,'shuf',true,'shuf_data',shufs{shufidx},'cross_only',cross_only,'non_mem',opt.non_mem);
        chains_shuf_rev_all=wave.COM_chain_reg(su_meta,sel_meta,reg_com_maps,'reverse',true,'shuf',true,'shuf_data',shufs{shufidx},'cross_only',cross_only,'non_mem',opt.non_mem);
        shuf_fwdcnt=[shuf_fwdcnt;nnz(chains_shuf_all.cross_reg)];
        shuf_wtncnt=[shuf_wtncnt;nnz(~chains_shuf_all.cross_reg)];
        shuf_revcnt=[shuf_revcnt;nnz(chains_shuf_rev_all.cross_reg)];
    end
    
    boot_delta_mm=bootstrp(100,@(x) [wtncnt,fwdcnt,revcnt]-mean(x),[shuf_wtncnt,shuf_fwdcnt,shuf_revcnt]);
    shufstd=std([shuf_wtncnt,shuf_fwdcnt,shuf_revcnt]);
    zscores=boot_delta_mm./shufstd;
    mm=mean(zscores);
    sem=std(zscores)./sqrt(size(zscores,1));

    pp=ranksum(zscores(:,2),zscores(:,3));

    fh=figure();
    tiledlayout(1,2);
    nexttile();
    hold on
    bh=bar(mm(1));
    errorbar(1,mm(1),sem(1),'k.','CapSize',12)
    set(gca,'XTick',1,'XtickLabel',{'Within'})
    if opt.non_mem
        bh.BaseValue=0;
        bh.BaseLine.LineStyle = ':';
        bh.BaseLine.Color = 'k';
        bh.BaseLine.LineWidth = 1;
        ylim([-2,1]);
    else
        ylim([0,110])
    end
    ylabel('Normalized occurrence (Z-Score)')

    nexttile();
    hold on
    bh=bar(mm(2:3));
    errorbar(1:2,mm(2:3),sem(2:3),'k.','CapSize',12)
    if opt.non_mem
        bh.BaseValue=0;
        bh.BaseLine.LineStyle = ':';
        bh.BaseLine.Color = 'k';
        bh.BaseLine.LineWidth = 1;
        ylim([-8,1]);
    else
    ylim([0,50])
    end
    set(gca,'XTick',1:2,'XtickLabel',{'Consistent','Inconsistent'})
    ylabel('Normalized occurrence (Z-Score)')

    sgtitle(sprintf('consist-vs-incong %.4f',pp))

    if ~opt.skip_save
        if opt.non_mem
            savefig(fh,fullfile('binary','chain_consist_incong_z_score_nonmem.fig'));
        else
            savefig(fh,fullfile('binary','chain_consist_incong_z_score.fig'));
        end
    end
end
end

%% All reverse within is either fwd within or fwd partial
% nonovl=cell(0,3);
% wisel=find(~chains_uf_rev_all.cross_reg);
% for ii=reshape(wisel,1,[])
%     sess=chains_uf_rev_all.sess(ii);
%     rcid=chains_uf_rev_all.cids{ii};
%     fcids=chains_uf_all.cids(chains_uf_all.sess==sess);
%     if any(cellfun(@(x) all(ismember(rcid,x),'all'), fcids))
%         continue
%     else
%         nonovl=[nonovl;{ii,sess,rcid}];
%     end
% end
