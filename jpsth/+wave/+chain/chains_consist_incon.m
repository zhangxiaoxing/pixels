%%
function fh=chains_consist_incon(su_meta,sel_meta,opt)
arguments
    su_meta = []
    sel_meta = []
    opt.skip_save (1,1) logical = true
    opt.shuf (1,1) logical = false
    opt.non_mem (1,1) logical = false
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
    opt.poolsize (1,1) double = 2
end
global_init;
cross_only=false;
if isempty(su_meta)
    switch opt.criteria
        case 'WT'
            load(fullfile('binary','su_meta.mat'),'su_meta');
        case 'Learning'
            su_meta=ephys.util.load_meta("save_file",false,"adjust_white_matter",true,"criteria","Learning","load_file",false,"skip_stats",true);
        case 'any'
            keyboard()
    end
end
if isempty(sel_meta)
    switch opt.criteria
        case 'WT'
            fstr=load(fullfile('binary','wrs_mux_meta.mat'));
            sel_meta=fstr.wrs_mux_meta;
            clear fstr
        case 'Learning'
            sel_meta=ephys.get_wrs_mux_meta('load_file',false,'save_file',false,'criteria','Learning','extend6s',true);
        case 'any'
            keyboard()
    end
end
reg_com_maps=wave.get_reg_com_maps(sel_meta,'criteria',opt.criteria);

poolh=parpool(opt.poolsize);
Fconsist=parfeval(poolh,@wave.COM_chain_reg,1,su_meta,sel_meta,reg_com_maps,'cross_only',cross_only,'non_mem',opt.non_mem,'criteria',opt.criteria);
Fincong=parfeval(poolh,@wave.COM_chain_reg,1,su_meta,sel_meta,reg_com_maps,'reverse',true,'cross_only',cross_only,'non_mem',opt.non_mem,'criteria',opt.criteria);

chains_uf_all=fetchOutputs(Fconsist);
chains_uf_rev_all=fetchOutputs(Fincong);

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
        switch opt.criteria
            case 'WT'
                savefig(fh,fullfile('binary','chains_consist_incon_nonmem.fig'))
            case 'Learning'
                savefig(fh,fullfile('binary','LN_chains_consist_incon_nonmem.fig'))
            otherwise
                keyboard()
        end
    else
        switch opt.criteria
            case 'WT'
                savefig(fh,fullfile('binary','chains_consist_incon.fig'))
            case 'Learning'
                savefig(fh,fullfile('binary','LN_chains_consist_incon.fig'))
            otherwise
                keyboard()
        end
    end
end

if opt.shuf
    fwdcnt=nnz(chains_uf_all.cross_reg);
    wtncnt=nnz(~chains_uf_all.cross_reg);
    revcnt=nnz(chains_uf_rev_all.cross_reg);

    
    switch opt.criteria
        case 'WT'
            load(fullfile('binary','bz_ring_shufs.mat'),'shufs');
        case 'Learning'
            load(fullfile('binary','LN_bz_shufs.mat'),'shufs');
        otherwise
            keyboard()
    end

    Fconsist=parallel.FevalFuture.empty(0,numel(shufs));
    Fincong=parallel.FevalFuture.empty(0,numel(shufs));
    for shufidx=1:numel(shufs) 
        Fconsist(shufidx)=parfeval(poolh,@wave.COM_chain_reg,1,su_meta,sel_meta,reg_com_maps,'shuf',true,'shuf_data',shufs{shufidx},'cross_only',cross_only,'non_mem',opt.non_mem,'criteria',opt.criteria);
        Fincong(shufidx)=parfeval(poolh,@wave.COM_chain_reg,1,su_meta,sel_meta,reg_com_maps,'reverse',true,'shuf',true,'shuf_data',shufs{shufidx},'cross_only',cross_only,'non_mem',opt.non_mem,'criteria',opt.criteria);
    end
    
    [shuf_fwdcnt,shuf_wtncnt,shuf_revcnt]=deal([]);
    for shufidx=1:numel(shufs) 
        chains_shuf_all=fetchOutputs(Fconsist(shufidx));
        chains_shuf_rev_all=fetchOutputs(Fincong(shufidx));
        shuf_fwdcnt=[shuf_fwdcnt;nnz(chains_shuf_all.cross_reg)];
        shuf_wtncnt=[shuf_wtncnt;nnz(~chains_shuf_all.cross_reg)];
        shuf_revcnt=[shuf_revcnt;nnz(chains_shuf_rev_all.cross_reg)];
    end

    if false
        if opt.non_mem
            chain_consist=cell2struct({wtncnt;fwdcnt;revcnt;shuf_wtncnt;shuf_fwdcnt;shuf_revcnt},...
                {'Within_observed','Consistent_observed','Inconsistent_observed','Within_shuffled','Consistent_shuffled','Inconsistent_shuffled'});
            fid=fopen(fullfile('binary','upload','F2IR_chain_nonmemory_consistent_inconsistent.json'),'w');
            fprintf(fid,jsonencode(chain_consist));
            fclose(fid);
        else
            chain_consist=cell2struct({wtncnt;fwdcnt;revcnt;shuf_wtncnt;shuf_fwdcnt;shuf_revcnt},...
                {'Within_observed','Consistent_observed','Inconsistent_observed','Within_shuffled','Consistent_shuffled','Inconsistent_shuffled'});
            fid=fopen(fullfile('binary','upload','F2IL_chain_congruent_consistent_inconsistent.json'),'w');
            fprintf(fid,jsonencode(chain_consist));
            fclose(fid);
        end
    end


    boot_delta_mm=bootstrp(1000,@(x) [wtncnt,fwdcnt,revcnt]-mean(x),[shuf_wtncnt,shuf_fwdcnt,shuf_revcnt]);
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
    ylim([0,70])
    end
    set(gca,'XTick',1:2,'XtickLabel',{'Consistent','Inconsistent'})
    ylabel('Normalized occurrence (Z-Score)')

    sgtitle(sprintf('consist-vs-incong %.4f',pp))

    if ~opt.skip_save
        if opt.non_mem
            switch opt.criteria
                case 'WT'
                    savefig(fh,fullfile('binary','chain_consist_incong_z_score_nonmem.fig'));
                case 'Learning'
                    savefig(fh,fullfile('binary','LN_chain_consist_incong_z_score_nonmem.fig'));
                otherwise
                    keyboard()
            end
        else
            switch opt.criteria
                case 'WT'
                    savefig(fh,fullfile('binary','chain_consist_incong_z_score.fig'));
                case 'Learning'
                    savefig(fh,fullfile('binary','LN_chain_consist_incong_z_score.fig'));
                otherwise
                    keyboard();
            end
        end
    end
end
delete(poolh);
end

