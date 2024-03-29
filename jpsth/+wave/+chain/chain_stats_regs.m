function [fh3,fhf]=chain_stats_regs(su_meta,sel_meta,len_thresh,opt)
arguments
    su_meta = []
    sel_meta = []
    len_thresh (1,1) double =3
    opt.odor_only (1,1) logical = true
end

if isempty(su_meta)
    load(fullfile('binary','su_meta.mat'))
end

if isempty(sel_meta)
    fstr=load(fullfile('binary','wrs_mux_meta.mat'));
    sel_meta=fstr.wrs_mux_meta;
    clear fstr
end

reg_com_maps=wave.get_reg_com_maps(sel_meta);

shuf_sum=struct([]);
%% data collection
greys=ephys.getGreyRegs('range','grey');
% vs shuffle
% jpsth/+bz/+rings/shuffle_conn_bz_alt.m
global gather_config
load(fullfile('binary','bz_ring_shufs.mat'),'shufs');

for shufidx=[1:numel(shufs),0] % 0 has to be last to keep value for later pipeline
    if shufidx>0
        chains_uf=wave.COM_chain_reg(su_meta,sel_meta,reg_com_maps,'shuf',true,'shuf_data',shufs{shufidx},'cross_only',false,'non_mem',false);
    else
        chains_uf=wave.COM_chain_reg(su_meta,sel_meta,reg_com_maps,'cross_only',false,'non_mem',false);
    end

    % region tag
    chains_uf.reg=cell(size(chains_uf.cids));
    chains_uf.reg_sel=addcats(categorical(NaN(size(chains_uf.cids))),["within","cross"]);

    % Screen through corresponding regions.
    % TODO: separate processing of within region and cross region chains.
    for sess=reshape(unique(chains_uf.sess),1,[])
        sesscid=su_meta.allcid(su_meta.sess==sess);
        sessreg=su_meta.reg_tree(5,su_meta.sess==sess);
        for cnid=reshape(find(chains_uf.sess==sess),1,[])
            [~,supos]=ismember(uint16(chains_uf.cids{cnid}),sesscid);
            chains_uf.reg{cnid}=sessreg(supos);
            if all(ismember(sessreg(supos),greys),'all')
                if all(strcmp(sessreg(supos(2:end)),sessreg(supos(1))),'all')
                    % within region
                    chains_uf.reg_sel(cnid)="within";
                else
                    % cross region
                    chains_uf.reg_sel(cnid)="cross";
                end
            end
        end
    end %sess
    % within & cross count
%     disp([...
%         nnz(chains_uf.reg_sel=="within" & cellfun(@(x) numel(x)>=opt.len_thresh,chains_uf.cids)),...
%         nnz(chains_uf.reg_sel=="cross" & cellfun(@(x) numel(x)>=opt.len_thresh,chains_uf.cids))...
%         ]);

    % repeated occurrence
    chains_uf.uid=arrayfun(@(x) chains_uf.sess(x)*100000+int32(chains_uf.cids{x}), 1:numel(chains_uf.sess),'UniformOutput',false);
    len_sel=cellfun(@(x) numel(x),chains_uf.cids)>=len_thresh;

    for curr_tag=["within","cross"]

        if opt.odor_only
            olf_sel=ismember(chains_uf.wave,["olf_s1","olf_s2","s1d3","s2d3","s1d6","s2d6"]) & len_sel & (chains_uf.reg_sel==curr_tag);
            both_u_reg_cat=struct();
        else
            olf_sel=ismember(chains_uf.wave,["olf_s1","olf_s2"]) & len_sel & (chains_uf.reg_sel==curr_tag);
        end
        olf_uid=[chains_uf.uid{olf_sel}];
        [chain_olf_uid,usel]=unique(olf_uid);

        olf_reg=[chains_uf.reg{olf_sel}];
        olf_reg_cat.(curr_tag)=categorical(olf_reg);
        olf_u_reg_cat.(curr_tag)=categorical(olf_reg(usel));

        if ~opt.odor_only
            both_sel=ismember(chains_uf.wave,["s1d3","s2d3","s1d6","s2d6"]) & len_sel & chains_uf.reg_sel==curr_tag;
            both_uid=[chains_uf.uid{both_sel}];
            [chain_both_uid,usel]=unique(both_uid);
            both_reg=[chains_uf.reg{both_sel}];
            both_reg_cat.(curr_tag)=categorical(both_reg);
            both_u_reg_cat.(curr_tag)=categorical(both_reg(usel));
        end

        olf_ratio.(curr_tag)=[];
        both_ratio.(curr_tag)=[];
        if shufidx>0
            all_chain_regs.(curr_tag)=greys;
        else
            all_chain_regs.(curr_tag)=unique([struct2array(olf_u_reg_cat),struct2array(both_u_reg_cat)]);
        end
        for rr=all_chain_regs.(curr_tag)
            su_cnt=nnz(strcmp(su_meta.reg_tree(5,:),char(rr)));
            chains_occur_olf=nnz(olf_reg_cat.(curr_tag)==rr);
            olf_ratio.(curr_tag)=[olf_ratio.(curr_tag);chains_occur_olf./su_cnt];
            if ~opt.odor_only
                chains_occur_both=nnz(both_reg_cat.(curr_tag)==rr);
                both_ratio.(curr_tag)=[both_ratio.(curr_tag);chains_occur_both./su_cnt];
            end
        end
        if shufidx==0
            [olf_ratio.(curr_tag),olf_srt.(curr_tag)]=sort(olf_ratio.(curr_tag),'descend');
            if ~opt.odor_only
                [both_ratio.(curr_tag),both_srt.(curr_tag)]=sort(both_ratio.(curr_tag),'descend');
            end
        end
    end
    if shufidx>0
        shuf_sum=[shuf_sum;struct('all_chain_regs',all_chain_regs,'olf_ratio',olf_ratio,'both_ratio',both_ratio)];
    end
end

%% occurrence of chain neuron per region neuron
if false
    fh3=figure();
    if opt.odor_only
        tiledlayout(1,2)
    else
        tiledlayout(2,2)
    end
    for curr_tag=["within","cross"]
        nexttile()
        bary=zeros(1,3);
        barn=min(numel(olf_ratio.(curr_tag)),3);
        bary(1:barn)=olf_ratio.(curr_tag)(1:barn);
        bar(bary);
        set(gca(),'XTick',1:barn,'XTickLabel',all_chain_regs.(curr_tag)(olf_srt.(curr_tag)(1:barn)),'YScale','log','XTickLabelRotation',90);
        ylim([0.001,10]);
        ylabel('occurrence per neuron')
        title("olf "+curr_tag)
        if opt.odor_only
            continue
        end
        nexttile()
        bary=zeros(1,3);
        barn=min(numel(both_ratio.(curr_tag)),3);
        bary(1:barn)=both_ratio.(curr_tag)(1:barn);
        bar(bary);
        set(gca(),'XTick',1:barn,'XTickLabel',all_chain_regs.(curr_tag)(both_srt.(curr_tag)(1:barn)),'YScale','log','XTickLabelRotation',90);
        ylim([0.001,10]);
        ylabel('occurrence per neuron')
        title("both "+curr_tag)
    end
end

%% full list vs shuffle

fhf=figure();
if opt.odor_only
    tiledlayout(1,2)
else
    tiledlayout(2,2)
end
chain_reg=struct();
for curr_tag=["within","cross"]
    nexttile()
    hold on
    barsel=olf_ratio.(curr_tag)>0;
    regsel=all_chain_regs.(curr_tag)(olf_srt.(curr_tag)(barsel));
    shufmat=nan(numel(shuf_sum),numel(regsel));
    regidx=1;
    for onereg=reshape(regsel,1,[])
        for ii=1:numel(shuf_sum)
            shufmat(ii,regidx)=shuf_sum(ii).olf_ratio.(curr_tag)(strcmp(shuf_sum(ii).all_chain_regs.(curr_tag),char(onereg)));
        end
        regidx=regidx+1;
    end
    bh=bar(1:numel(regsel),[olf_ratio.(curr_tag)(barsel),mean(shufmat).'],'grouped');
    bh(1).FaceColor='k';
    if numel(bh)>1
        bh(2).FaceColor='w';
    end
    shufstd=std(shufmat);
    shufsem=shufstd./sqrt(size(shufmat,1));
    if numel(bh)>1
        errorbar(bh(2).XEndPoints,bh(2).YData,shufsem,'k.')
    end
    % olf_ratio.(curr_tag)(barsel)-mean(shufmat).'


    set(gca(),'XTick',1:numel(regsel),'XTickLabel',all_chain_regs.(curr_tag)(olf_srt.(curr_tag)(barsel)),'YScale','log','XTickLabelRotation',90);
    ylim([1e-4,10]);
    ylabel('occurrence per neuron')
    title("olf "+curr_tag)

    zscores=abs(olf_ratio.(curr_tag)(barsel)-mean(shufmat).')./std(shufmat).';
    chain_reg.(curr_tag).zscore=zscores;
    stars3=zscores>norminv(0.9995);
    stars2=zscores>norminv(0.995) & zscores<=norminv(0.9995);
    stars1=zscores>norminv(0.975) & zscores<=norminv(0.995);
    text(find(stars3),repmat(9,nnz(stars3),1),'***','HorizontalAlignment','center')
    text(find(stars2),repmat(3.16,nnz(stars2),1),'**','HorizontalAlignment','center')
    text(find(stars1),repmat(1,nnz(stars1),1),'*','HorizontalAlignment','center')
    
    if opt.odor_only
        chainreg.(curr_tag)=table(all_chain_regs.(curr_tag)(olf_srt.(curr_tag)(barsel)).',olf_ratio.(curr_tag)(barsel),shufmat.','VariableNames',{'Region','Observed','Shuffled'});
        continue
    end

    nexttile()
    hold on
    barsel=both_ratio.(curr_tag)>0;
    regsel=all_chain_regs.(curr_tag)(both_srt.(curr_tag)(barsel));
    shufmat=nan(numel(shuf_sum),numel(regsel));
    regidx=1;
    for onereg=reshape(regsel,1,[])
        for ii=1:numel(shuf_sum)
            shufmat(ii,regidx)=shuf_sum(ii).both_ratio.(curr_tag)(strcmp(shuf_sum(ii).all_chain_regs.(curr_tag),char(onereg)));
        end
        regidx=regidx+1;
    end
    bh=bar(1:numel(regsel),[both_ratio.(curr_tag)(barsel),mean(shufmat).'],'grouped');
    bh(1).FaceColor='k';
    bh(2).FaceColor='w';
    shufsem=std(shufmat)./sqrt(size(shufmat,1));
    errorbar(bh(2).XEndPoints,bh(2).YData,shufsem,'k.')
    set(gca(),'XTick',1:numel(regsel),'XTickLabel',all_chain_regs.(curr_tag)(both_srt.(curr_tag)),'YScale','log','XTickLabelRotation',90);
    ylim([1e-4,10]);
    ylabel('occurrence per neuron')
    title("both "+curr_tag)
    zscores=abs(both_ratio.(curr_tag)(barsel)-mean(shufmat).')./std(shufmat).';
    disp(curr_tag)
    disp(zscores)
    stars3=zscores>norminv(0.9995);
    stars2=zscores>norminv(0.995) & zscores<=norminv(0.9995);
    stars1=zscores>norminv(0.975) & zscores<=norminv(0.995);
    text(find(stars3),repmat(9,nnz(stars3),1),'***','HorizontalAlignment','center')
    text(find(stars2),repmat(3.16,nnz(stars2),1),'**','HorizontalAlignment','center')
    text(find(stars1),repmat(1,nnz(stars1),1),'*','HorizontalAlignment','center')
end
savefig(fhf,fullfile('binary','chains_occur_per_neuron_per_reg.fig'));
blame=vcs.blame();
save(fullfile('binary','chains_occur_per_neuron_per_reg.mat'),'chain_reg','blame');
if false
    fid=fopen(fullfile('binary','upload','F2P_SF6B_SF6C_Chain_occurrence_across_regions.json'),'w');
    fprintf(fid,jsonencode(chainreg));
    fclose(fid);
end
end