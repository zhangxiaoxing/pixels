shuf_sum=struct([]);

%% data collect
global_init;
wrs_mux_meta=ephys.get_wrs_mux_meta();
greys=ephys.getGreyRegs('range','grey');
su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);

% vs shuffle
% jpsth/+bz/+rings/shuffle_conn_bz_alt.m
% wave.COM_chain_shuf(wrs_mux_meta);
load('chains_shuf.mat','shuf_chains')
for shuf_idx=[1:100,0]
    if shuf_idx>0
        chains_uf=shuf_chains{shuf_idx};
    else
        load(fullfile('bzdata','chains_mix.mat'),'chains_uf');
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
    disp([...
        nnz(chains_uf.reg_sel=="within" & cellfun(@(x) numel(x)>4,chains_uf.cids)),...
        nnz(chains_uf.reg_sel=="cross" & cellfun(@(x) numel(x)>4,chains_uf.cids))...
        ]);

    % repeated occurrence
    chains_uf.uid=arrayfun(@(x) chains_uf.sess(x)*100000+int32(chains_uf.cids{x}), 1:numel(chains_uf.sess),'UniformOutput',false);
    len_sel=cellfun(@(x) numel(x),chains_uf.cids)>4;

    olf_sel=ismember(chains_uf.wave,["olf_s1","olf_s2"]) & len_sel; % & (chains_uf.reg_sel==curr_tag)
    olf_uid=[chains_uf.uid{olf_sel}];


    both_sel=ismember(chains_uf.wave,["s1d3","s2d3","s1d6","s2d6"]) & len_sel; % & chains_uf.reg_sel==curr_tag
    both_uid=[chains_uf.uid{both_sel}];

    dur_sel=ismember(chains_uf.wave,["dur_d3","dur_d6"]) & len_sel; % & chains_uf.reg_sel==curr_tag
    dur_uid=[chains_uf.uid{dur_sel}];


    for curr_tag=["within","cross"]
        olf_sel=ismember(chains_uf.wave,["olf_s1","olf_s2"]) & len_sel & (chains_uf.reg_sel==curr_tag);

        olf_uid=[chains_uf.uid{olf_sel}];
        [chain_olf_uid,usel]=unique(olf_uid);

        olf_reg=[chains_uf.reg{olf_sel}];
        olf_reg_cat.(curr_tag)=categorical(olf_reg);
        olf_u_reg_cat.(curr_tag)=categorical(olf_reg(usel));

        both_sel=ismember(chains_uf.wave,["s1d3","s2d3","s1d6","s2d6"]) & len_sel & chains_uf.reg_sel==curr_tag;
        both_uid=[chains_uf.uid{both_sel}];
        [chain_both_uid,usel]=unique(both_uid);
        both_reg=[chains_uf.reg{both_sel}];
        both_reg_cat.(curr_tag)=categorical(both_reg);
        both_u_reg_cat.(curr_tag)=categorical(both_reg(usel));

        % dur_sel=ismember(chains_uf.wave,["dur_d3","dur_d6"]) & len_sel; % & chains_uf.reg_sel==curr_tag
        % dur_uid=[chains_uf.uid{dur_sel}];
        % [chain_dur_uid,usel]=unique(dur_uid);
        % dur_reg=[chains_uf.reg{dur_sel}];
        % dur_reg_cat=categorical(dur_reg);
        % dur_u_reg_cat=categorical(dur_reg(usel));

        olf_ratio.(curr_tag)=[];
        both_ratio.(curr_tag)=[];
        if shuf_idx>0
            all_chain_regs.(curr_tag)=greys;
        else
            all_chain_regs.(curr_tag)=unique([struct2array(olf_u_reg_cat),struct2array(both_u_reg_cat)]);
        end
        for rr=all_chain_regs.(curr_tag)
            su_cnt=nnz(strcmp(su_meta.reg_tree(5,:),char(rr)));
            chains_occur_olf=nnz(olf_reg_cat.(curr_tag)==rr);
            chains_occur_both=nnz(both_reg_cat.(curr_tag)==rr);
            olf_ratio.(curr_tag)=[olf_ratio.(curr_tag);chains_occur_olf./su_cnt];
            both_ratio.(curr_tag)=[both_ratio.(curr_tag);chains_occur_both./su_cnt];
        end
        if shuf_idx==0
            [olf_ratio.(curr_tag),olf_srt.(curr_tag)]=sort(olf_ratio.(curr_tag),'descend');
            [both_ratio.(curr_tag),both_srt.(curr_tag)]=sort(both_ratio.(curr_tag),'descend');
        end
    end
    if shuf_idx>0
        shuf_sum=[shuf_sum;struct('all_chain_regs',all_chain_regs,'olf_ratio',olf_ratio,'both_ratio',both_ratio)];
    end
end

%% occurrence of chain neuron per region neuron
figure()
tiledlayout(2,2)
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
 

%% full list vs shuffle

figure()
tiledlayout(2,2)
for curr_tag=["within","cross"]
    nexttile()
    hold on
    barsel=olf_ratio.(curr_tag)>0;
    regsel=all_chain_regs.(curr_tag)(olf_srt.(curr_tag)(barsel));
    shufmat=nan(numel(shuf_chains),numel(regsel));
    regidx=1;
    for onereg=reshape(regsel,1,[])
        for ii=1:numel(shuf_sum)
            shufmat(ii,regidx)=shuf_sum(ii).olf_ratio.(curr_tag)(strcmp(shuf_sum(ii).all_chain_regs.(curr_tag),char(onereg)));
        end             
        regidx=regidx+1;
    end
    bh=bar(1:numel(regsel),[olf_ratio.(curr_tag)(barsel),mean(shufmat).'],'grouped');
    bh(1).FaceColor='w';
    bh(2).FaceColor='k';
    shufstd=std(shufmat);
    shufsem=shufstd./sqrt(size(shufmat,1));
    errorbar(bh(2).XEndPoints,bh(2).YData,shufsem,'k.')
    olf_ratio.(curr_tag)(barsel)-mean(shufmat).'


    set(gca(),'XTick',1:numel(regsel),'XTickLabel',all_chain_regs.(curr_tag)(olf_srt.(curr_tag)),'YScale','log','XTickLabelRotation',90);
    ylim([1e-5,10]);
    ylabel('occurrence per neuron')
    title("olf "+curr_tag)

    zscores=abs(olf_ratio.(curr_tag)(barsel)-mean(shufmat).')./std(shufmat).';
    disp(curr_tag)
    disp(zscores)
    text(find(zscores>norminv(0.9995)),repmat(9,nnz(zscores>norminv(0.9995)),1),'***','HorizontalAlignment','center')
    text(find(zscores>norminv(0.995)),repmat(3.16,nnz(zscores>norminv(0.995)),1),'**','HorizontalAlignment','center')
    text(find(zscores>norminv(0.975)),repmat(1,nnz(zscores>norminv(0.975)),1),'*','HorizontalAlignment','center')


    nexttile()
    hold on
    barsel=both_ratio.(curr_tag)>0;
    regsel=all_chain_regs.(curr_tag)(both_srt.(curr_tag)(barsel));
    shufmat=nan(numel(shuf_chains),numel(regsel));
    regidx=1;
    for onereg=reshape(regsel,1,[])
        for ii=1:numel(shuf_sum)
            shufmat(ii,regidx)=shuf_sum(ii).both_ratio.(curr_tag)(strcmp(shuf_sum(ii).all_chain_regs.(curr_tag),char(onereg)));
        end
        regidx=regidx+1;
    end
    bh=bar(1:numel(regsel),[both_ratio.(curr_tag)(barsel),mean(shufmat).'],'grouped');
    bh(1).FaceColor='w';
    bh(2).FaceColor='k';
    shufsem=std(shufmat)./sqrt(size(shufmat,1));
    errorbar(bh(2).XEndPoints,bh(2).YData,shufsem,'k.')
    set(gca(),'XTick',1:numel(regsel),'XTickLabel',all_chain_regs.(curr_tag)(both_srt.(curr_tag)),'YScale','log','XTickLabelRotation',90);
    ylim([1e-5,10]);
    ylabel('occurrence per neuron')
    title("both "+curr_tag)
    zscores=abs(both_ratio.(curr_tag)(barsel)-mean(shufmat).')./std(shufmat).';
    disp(curr_tag)
    disp(zscores)
    text(find(zscores>norminv(0.9995)),repmat(9,nnz(zscores>norminv(0.9995)),1),'***','HorizontalAlignment','center')
    text(find(zscores>norminv(0.995)),repmat(3.16,nnz(zscores>norminv(0.995)),1),'**','HorizontalAlignment','center')
    text(find(zscores>norminv(0.975)),repmat(1,nnz(zscores>norminv(0.975)),1),'*','HorizontalAlignment','center')
end


