function [fc_com_pvsst_stats,screen_stats]=fc_conn_screen(com_map,pct_meta,opt)
arguments
    com_map
    pct_meta
    opt.title_suffix (1,:) char = []
end
[sig,~]=bz.load_sig_sums_conn_file();
sig=bz.join_fc_waveid(sig,pct_meta.wave_id);
idmap=load(fullfile('..','align','reg_ccfid_map.mat'));

sink_ccfid=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/grey_targets');
src_ccfid=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/grey_srcs');
sink_src_mat=h5read(fullfile('..','allensdk','proj_mat.hdf5'),'/src_target_matrix');

screen_stats=nan(numel(src_ccfid),13);% src_ccfid, #olf_l2h, #olf_l2h_consist, #olf_h2l, #olf_h2l_consist, {DUR x 4},{MIX x 4}

for srcidx=1:numel(src_ccfid)
    connv=sink_src_mat(:,srcidx);
%     hiermap=containers.Map(cellfun(@(x) x,idmap.ccfid2reg.values(num2cell(sink_ccfid))),...
    hiermap=containers.Map(sink_ccfid,connv);
    hiermap(0)=nan;

    fc_com_pvsst_stats=[];
    usess=unique(sig.sess);
    
    for sii=reshape(usess,1,[]) %iterate through sessions
        sesssel=sig.sess==sii;
        
        suid=sig.suid(sesssel,:);
        waveid=sig.waveid(sesssel,:);
        regsess=squeeze(sig.reg(sesssel,5,:));
        hier_v=cell2mat(hiermap.values(num2cell(regsess)));
        com_sess_pct=nan(size(suid));
        for ff=["s1d3","s1d6","s2d3","s2d6","olf_s1","olf_s2","dur_d3","dur_d6"]
            sukeys=com_map.(['s',num2str(sii)]).(ff).com.keys(); % prefered SUid
            susel=ismember(suid,int32(cell2mat(sukeys)));% same dim as suid
            com_sess_pct(susel)=cell2mat(com_map.(['s',num2str(sii)]).(ff).com.values(num2cell(suid(susel)))); % out put is nx2 in dim
        end
        fc_com_pvsst_stats=[fc_com_pvsst_stats;double(sii).*ones(size(suid(:,1))),double(suid),com_sess_pct,hier_v,double(regsess),double(waveid)];
        %==================================================sess=====================suid=======COM_context1=hiermap_v====ccfid===========waveid======
    end
    % olf and dur
    finite_sel=all(isfinite(fc_com_pvsst_stats(:,4:7)),2);
    same_sel=fc_com_pvsst_stats(:,6)==fc_com_pvsst_stats(:,7);
    m2l_sel=fc_com_pvsst_stats(:,6)>fc_com_pvsst_stats(:,7);
    l2m_sel=fc_com_pvsst_stats(:,6)<fc_com_pvsst_stats(:,7);
    congrusel=pct.su_pairs.get_congru(fc_com_pvsst_stats(:,10:11));
    
    olf_m2l_sel=finite_sel & m2l_sel & congrusel & ~any(ismember(fc_com_pvsst_stats(:,10:11),7:8),2);  % no 7 or 8
    olf_l2m_sel=finite_sel & l2m_sel & congrusel & ~any(ismember(fc_com_pvsst_stats(:,10:11),7:8),2);  % no 7 or 8
    dur_m2l_sel=finite_sel & m2l_sel & congrusel & ~any(ismember(fc_com_pvsst_stats(:,10:11),5:6),2);  % no 5 or 6
    dur_l2m_sel=finite_sel & l2m_sel & congrusel & ~any(ismember(fc_com_pvsst_stats(:,10:11),5:6),2);  % no 5 or 6
    mix_m2l_sel=finite_sel & m2l_sel & congrusel & ~any(ismember(fc_com_pvsst_stats(:,10:11),5:8),2);  % no 7 or 8
    mix_l2m_sel=finite_sel & l2m_sel & congrusel & ~any(ismember(fc_com_pvsst_stats(:,10:11),5:8),2);  % no 7 or 8    
    
    if strcmp(idmap.ccfid2reg(src_ccfid(srcidx)),'COA')
        olf_sel=finite_sel & congrusel & ~any(ismember(fc_com_pvsst_stats(:,10:11),7:8),2);  % no 7 or 8
        olf_same_sel=finite_sel & same_sel & congrusel & ~any(ismember(fc_com_pvsst_stats(:,10:11),7:8),2);  % no 7 or 8
        histo_panel(fc_com_pvsst_stats(olf_sel,:),fc_com_pvsst_stats(olf_same_sel,:),...
            fc_com_pvsst_stats(olf_m2l_sel,:),fc_com_pvsst_stats(olf_l2m_sel,:));
        sgtitle('Olfactory congruent FC, COA projection gradient')
    elseif strcmp(idmap.ccfid2reg(src_ccfid(srcidx)),'MTN')
        mix_sel=finite_sel & congrusel & ~any(ismember(fc_com_pvsst_stats(:,10:11),5:8),2);  % no 7 or 8
        mix_same_sel=finite_sel & same_sel & congrusel & ~any(ismember(fc_com_pvsst_stats(:,10:11),5:8),2);  % no 7 or 8
        histo_panel(fc_com_pvsst_stats(mix_sel,:),fc_com_pvsst_stats(mix_same_sel,:),...
            fc_com_pvsst_stats(mix_m2l_sel,:),fc_com_pvsst_stats(mix_l2m_sel,:));
        sgtitle('Mixed congruent FC, MTN projection gradient')
    elseif strcmp(idmap.ccfid2reg(src_ccfid(srcidx)),'ILA')
        dur_sel=finite_sel & congrusel & ~any(ismember(fc_com_pvsst_stats(:,10:11),5:6),2);  % no 7 or 8
        dur_same_sel=finite_sel & same_sel & congrusel & ~any(ismember(fc_com_pvsst_stats(:,10:11),5:6),2);  % no 7 or 8
        histo_panel(fc_com_pvsst_stats(dur_sel,:),fc_com_pvsst_stats(dur_same_sel,:),...
            fc_com_pvsst_stats(dur_m2l_sel,:),fc_com_pvsst_stats(dur_l2m_sel,:));
        sgtitle('Duration congruent FC, ILA projection gradient')
    end



    % src_ccfid, #olf_total, #olf_consist, olf_com_latency,#dur_total, #dur_consist, dur_com_latency    
    screen_stats(srcidx,:)=[src_ccfid(srcidx),...
        nnz(olf_m2l_sel),...
        nnz(diff(fc_com_pvsst_stats(olf_m2l_sel,4:5),1,2)>0),...%early to late
        nnz(olf_l2m_sel),...
        nnz(diff(fc_com_pvsst_stats(olf_l2m_sel,4:5),1,2)>0),...
        nnz(dur_m2l_sel),...
        nnz(diff(fc_com_pvsst_stats(dur_m2l_sel,4:5),1,2)>0),...
        nnz(dur_l2m_sel),...
        nnz(diff(fc_com_pvsst_stats(dur_l2m_sel,4:5),1,2)>0),...
        nnz(mix_m2l_sel),...
        nnz(diff(fc_com_pvsst_stats(mix_m2l_sel,4:5),1,2)>0),...%early to late
        nnz(mix_l2m_sel),...
        nnz(diff(fc_com_pvsst_stats(mix_l2m_sel,4:5),1,2)>0)...
        ];

end

for shift=[0,2]
    olf_consist_ratio=screen_stats(:,3+shift)./screen_stats(:,2+shift);
    [~,sidx]=sort(olf_consist_ratio);
    figure('Color','w','Position',[4,4,1440,720])
    bar(olf_consist_ratio(sidx));
    regs=idmap.ccfid2reg.values(num2cell(screen_stats(sidx,1)));
    set(gca,'XTick',1:numel(sidx),'XTickLabel',[regs{:}])
    ylim([0.3,0.7])
    yline(0.5,'r--')
    if shift==0
        title("Olfactory FC early to late v.s more to less projection "+opt.title_suffix)
    else
        title("Olfactory FC early to late v.s less to more projection "+opt.title_suffix)
    end

    dur_consist_ratio=screen_stats(:,7+shift)./screen_stats(:,6+shift);
    [~,sidx]=sort(dur_consist_ratio);
    figure('Color','w','Position',[4,4,1440,720])
    bar(dur_consist_ratio(sidx));
    regs=idmap.ccfid2reg.values(num2cell(screen_stats(sidx,1)));
    set(gca,'XTick',1:numel(sidx),'XTickLabel',[regs{:}])
    ylim([0.3,0.7])
    yline(0.5,'r--')
    if shift==0
        title("Duration FC early to late v.s more to less projection "+opt.title_suffix)
    else
        title("Duration FC early to late v.s less to more projection "+opt.title_suffix)
    end

    mix_consist_ratio=screen_stats(:,11+shift)./screen_stats(:,10+shift);
    [~,sidx]=sort(mix_consist_ratio);
    figure('Color','w','Position',[4,4,1440,720])
    bar(mix_consist_ratio(sidx));
    regs=idmap.ccfid2reg.values(num2cell(screen_stats(sidx,1)));
    set(gca,'XTick',1:numel(sidx),'XTickLabel',[regs{:}])
    ylim([0.3,0.7])
    yline(0.5,'r--')
    if shift==0
        title("Mixed FC early to late v.s more to less projection "+opt.title_suffix)
    else
        title("Mixed FC early to late v.s less to more projection "+opt.title_suffix)
    end

end

end



function fh=histo_panel(fc_com_pvsst_stats_ovall,fc_com_pvsst_stats_same,fc_com_pvsst_stats_m2l,fc_com_pvsst_stats_l2m)
fh=figure('Color','w');
tiledlayout(2,2)
% total
fcs={fc_com_pvsst_stats_ovall,fc_com_pvsst_stats_same,fc_com_pvsst_stats_m2l,fc_com_pvsst_stats_l2m};
sub_title={'Overall','Same region','More to less innervated','Less to more innervated'}
for fii=1:4
    fc=fcs{fii};
    nexttile
    latency=diff(fc(:,4:5),1,2)./4;
    [h,p]=ttest(latency);
    histogram(latency,-1:0.1:1,'Normalization','probability')
    xline(0,'--k','LineWidth',1)
    xline(mean(latency),'--r','LineWidth',1)
    xlabel({'Lead-Follow TCOM latency (s)',sprintf('mean=%.3f, p=%.3f',mean(latency),p)})
    ylabel('Probability')
    title(sub_title{fii});
end
% % same
% nexttile
% latency=diff(fc_com_pvsst_stats_same(:,4:5),1,2)./4;
% histogram(latency,-0.5:0.1:0.5)
% xline(0,'--k')
% xline(mean(latency),'--r')
% % more to less
% nexttile
% latency=diff(fc_com_pvsst_stats_m2l(:,4:5),1,2)./4;
% histogram(latency,-0.5:0.1:0.5)
% xline(0,'--k')
% xline(mean(latency),'--r')
% % less to more
% nexttile
% latency=diff(fc_com_pvsst_stats_l2m(:,4:5),1,2)./4;
% histogram(latency,-0.5:0.1:0.5)
% xline(0,'--r')
% xline(mean(latency),'--k')
end