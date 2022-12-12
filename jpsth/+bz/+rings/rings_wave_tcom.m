function pstats=rings_wave_tcom(tcom3_maps,tcom6_maps)
persistent sums_all
if isempty(sums_all)
    load(fullfile('bzdata','sums_ring_stats_all.mat'));
end
pstats=struct();
[pstats.both3,pstats.olf3,pstats.dur3,pstats.both6,pstats.olf6,pstats.dur6]=deal(cell(0));

su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
wrs_mux_meta=ephys.get_wrs_mux_meta();

for rsize=3:5
    rstats=cell(0,10);
    one_rsize=sums_all{rsize-2};
    curr_sess=-1;
    for ridx=1:size(one_rsize,1)
        if curr_sess~=one_rsize{ridx,1}
            curr_sess=one_rsize{ridx,1};
            sesscid=su_meta.allcid(su_meta.sess==curr_sess);
            sesswaveid=wrs_mux_meta.wave_id(su_meta.sess==curr_sess);
            sess_part=su_meta.reg_tree(1,su_meta.sess==curr_sess);
            sess_reg=su_meta.reg_tree(5,su_meta.sess==curr_sess);

            sess_wave_map=containers.Map(num2cell(sesscid),num2cell(sesswaveid));
            sess_part_map=containers.Map(num2cell(sesscid),sess_part);
            sess_reg_map=containers.Map(num2cell(sesscid),sess_reg);
        end
        curr_waveid=cell2mat(sess_wave_map.values(num2cell(one_rsize{ridx,3})));
        curr_part=sess_part_map.values(num2cell(one_rsize{ridx,3}));
        curr_reg=sess_reg_map.values(num2cell(one_rsize{ridx,3}));
        if (~all(ismember(curr_part,{'CH','BS'}),'all')) || any(ismissing(curr_reg),'all')
            continue
        end
        [rwid,seltype]=bz.rings.ring_wave_type(curr_waveid);
        if strcmp(rwid,'congru')
            rstats=[rstats;one_rsize(ridx,:),curr_waveid,{curr_reg},seltype];
        end
    end

    rstats=rstats(cell2mat(rstats(:,6))>0.1 & cellfun(@(x) numel(unique(x)),rstats(:,3))==rsize,:);
    usess=unique(cell2mat(rstats(:,1)));
    
    for sessid=usess.'
        [spkID,spkTS,trials,~,~,FT_SPIKE]=ephys.getSPKID_TS(sessid,'keep_trial',true);
        rid=find(cell2mat(rstats(:,1))==sessid);
        for ri=reshape(rid,1,[])
            ts_id=[];
            cids=rstats{ri,3};
            disp({sessid,ri});

            per_cid_spk_cnt=cids;
            for in_ring_pos=1:numel(cids) % TODO, 1:rsize
                one_ring_sel=spkID==cids(in_ring_pos);
                per_cid_spk_cnt(in_ring_pos)=nnz(one_ring_sel);
                rawts=spkTS(one_ring_sel);

                ft_sel=strcmp(FT_SPIKE.label,num2str(cids(in_ring_pos)));
                ft_ts=FT_SPIKE.timestamp{ft_sel};
                ft_trl_time=FT_SPIKE.time{ft_sel};
                ft_trl=FT_SPIKE.trial{ft_sel};

                [~,tspos]=ismember(ft_ts,rawts);
                ext_time=repmat(-realmax,numel(rawts),1);
                ext_time(tspos)=ft_trl_time;

                ext_trl=repmat(-realmax,numel(rawts),1);
                ext_trl(tspos)=ft_trl;

                ts_id=cat(1,ts_id,[rawts,... % 1
                    repmat(cids(in_ring_pos),per_cid_spk_cnt(in_ring_pos),1),...  % 2
                    ones(per_cid_spk_cnt(in_ring_pos),1)*in_ring_pos,...  % 3
                    ext_time,...  % 4
                    ext_trl]); % 5
            end
            ts_id=sortrows(ts_id,1);
            ts_id=[ts_id,full(rstats{ri,5}.tags)]; % join TS, ring tag % 6


            [pref3,pref6]=bz.rings.preferred_trials_rings(rstats{ri,7},trials);

            [~,sel_type_idx]=ismember(rstats(ri,9),{'both','olf','dur'});
            if ~isempty(pref3)
                if all(ismember(rstats{ri,8},tcom3_maps{sel_type_idx}.keys()),'all')
                    reg_tcom=cell2mat(tcom3_maps{sel_type_idx}.values(rstats{ri,8}));
                    tssel=ismember(ts_id(:,5),pref3) & ts_id(:,4)>=1 & ts_id(:,4)<4 & ts_id(:,6)~=0;
                    rsums3=ts_id(tssel,:);
                    tag=[rstats{ri,9},'3'];
                    pstats.(tag)=[pstats.(tag);{rsums3(:,4)},rstats(ri,8),{reg_tcom}];
                end
            end
            if ~isempty(pref6)
                if all(ismember(rstats{ri,8},tcom6_maps{sel_type_idx}.keys()),'all')
                    reg_tcom=cell2mat(tcom6_maps{sel_type_idx}.values(rstats{ri,8}));
                    tssel=ismember(ts_id(:,5),pref6) & ts_id(:,4)>=1 & ts_id(:,4)<7 & ts_id(:,6)~=0;
                    rsums6=ts_id(tssel,:);
                    tag=[rstats{ri,9},'6'];
                    pstats.(tag)=[pstats.(tag);{rsums6(:,4)},rstats(ri,8),{reg_tcom}];
                end
            end
        end
    end
end
end

function post_analysis()
if false % overall descriptive
    sum([size(pstats.both3,1);size(pstats.olf3,1);size(pstats.both6,1);size(pstats.olf6,1)])

    nnz([cellfun(@(x) all(strcmp(x,'HIP'),'all'),pstats.both3(:,2));...
        cellfun(@(x) all(strcmp(x,'HIP'),'all'),pstats.olf3(:,2));...
        cellfun(@(x) all(strcmp(x,'HIP'),'all'),pstats.both6(:,2));...
        cellfun(@(x) all(strcmp(x,'HIP'),'all'),pstats.olf6(:,2))])

    nnz([cellfun(@(x) all(strcmp(x,'ORB'),'all'),pstats.both3(:,2));...
        cellfun(@(x) all(strcmp(x,'ORB'),'all'),pstats.olf3(:,2));...
        cellfun(@(x) all(strcmp(x,'ORB'),'all'),pstats.both6(:,2));...
        cellfun(@(x) all(strcmp(x,'ORB'),'all'),pstats.olf6(:,2))])
end

figure()
tiledlayout(2,3)

both6tcom=cellfun(@(x) sum(histcounts(x,1:0.25:7).*(1:24))./sum(histcounts(x,1:0.25:7))./4, pstats.both6(:,1));
both6_reg_wave=cellfun(@(x) mean(x), pstats.both6(:,3));
nexttile()
hold on
scatter(both6_reg_wave,both6tcom,16,'o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha',0.25);
c=[both6_reg_wave,ones(numel(both6tcom),1)]\both6tcom;
minmaxx=[min(both6_reg_wave),max(both6_reg_wave)];
plot(minmaxx,c(1).*minmaxx+c(2),'k--')
[r,p]=corr(both6_reg_wave,both6tcom);
xlabel('Averaged loops region wave timing (s)');
ylabel('Loops wave timing');
xlim([2.55,2.75])
text(max(xlim()),max(ylim()),sprintf('r=%.2f,p=%.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top')
title('both, 6s')

olf6tcom=cellfun(@(x) sum(histcounts(x,1:0.25:7).*(1:24))./sum(histcounts(x,1:0.25:7))./4, pstats.olf6(:,1));
olf6_reg_wave=cellfun(@(x) mean(x), pstats.olf6(:,3));
nexttile()
hold on
scatter(olf6_reg_wave,olf6tcom,16,'o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha',0.25);
c=[olf6_reg_wave,ones(numel(olf6tcom),1)]\olf6tcom;
minmaxx=[min(olf6_reg_wave),max(olf6_reg_wave)];
plot(minmaxx,c(1).*minmaxx+c(2),'k--')
[r,p]=corr(olf6_reg_wave,olf6tcom);
xlabel('Averaged loops region wave timing (s)');
ylabel('Loops wave timing');
% xlim([2.55,2.75])
text(max(xlim()),max(ylim()),sprintf('r=%.2f,p=%.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top')
title('Olf, 6s')

merge6tcom=[both6tcom;olf6tcom];
merge6_reg_wave=[both6_reg_wave;olf6_reg_wave];

nexttile()
hold on
scatter(merge6_reg_wave,merge6tcom,16,'o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha',0.25);
c=[merge6_reg_wave,ones(numel(merge6tcom),1)]\merge6tcom;
minmaxx=[min(merge6_reg_wave),max(merge6_reg_wave)];
plot(minmaxx,c(1).*minmaxx+c(2),'k--')
[r,p]=corr(merge6_reg_wave,merge6tcom);
xlabel('Averaged loops region wave timing (s)');
ylabel('Loops wave timing');
% xlim([2.55,2.75])
text(max(xlim()),max(ylim()),sprintf('r=%.2f,p=%.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top')
title('Merged both-olf, 6s')


both3tcom=cellfun(@(x) sum(histcounts(x,1:0.25:4).*(1:12))./sum(histcounts(x,1:0.25:4))./4, pstats.both3(:,1));
both3_reg_wave=cellfun(@(x) mean(x), pstats.both3(:,3));
nexttile()
hold on
scatter(both3_reg_wave,both3tcom,16,'o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha',0.25);
c=[both3_reg_wave,ones(numel(both3tcom),1)]\both3tcom;
minmaxx=[min(both3_reg_wave),max(both3_reg_wave)];
plot(minmaxx,c(1).*minmaxx+c(2),'k--')
[r,p]=corr(both3_reg_wave,both3tcom);
xlabel('Averaged loops region wave timing (s)');
ylabel('Loops wave timing');
% xlim([2.55,2.75])
text(max(xlim()),max(ylim()),sprintf('r=%.2f,p=%.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top')
title('both, 3s')


olf3tcom=cellfun(@(x) sum(histcounts(x,1:0.25:4).*(1:12))./sum(histcounts(x,1:0.25:4))./4, pstats.olf3(:,1));
olf3_reg_wave=cellfun(@(x) mean(x), pstats.olf3(:,3));
nexttile()
hold on
scatter(olf3_reg_wave,olf3tcom,16,'o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha',0.25);
c=[olf3_reg_wave,ones(numel(olf3tcom),1)]\olf3tcom;
minmaxx=[min(olf3_reg_wave),max(olf3_reg_wave)];
plot(minmaxx,c(1).*minmaxx+c(2),'k--')
[r,p]=corr(olf3_reg_wave,olf3tcom);
xlabel('Averaged loops region wave timing (s)');
ylabel('Loops wave timing');
% xlim([2.55,2.75])
text(max(xlim()),max(ylim()),sprintf('r=%.2f,p=%.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top')
title('Olf, 3s')

merge3tcom=[both3tcom;olf3tcom];
merge3_reg_wave=[both3_reg_wave;olf3_reg_wave];

nexttile()

hold on
scatter(merge3_reg_wave,merge3tcom,16,'o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha',0.25);
c=[merge3_reg_wave,ones(numel(merge3tcom),1)]\merge3tcom;
minmaxx=[min(merge3_reg_wave),max(merge3_reg_wave)];
plot(minmaxx,c(1).*minmaxx+c(2),'k--')
[r,p]=corr(merge3_reg_wave,merge3tcom);
xlabel('Averaged loops region wave timing (s)');
ylabel('Loops wave timing');
% xlim([2.55,2.75])
text(max(xlim()),max(ylim()),sprintf('r=%.2f,p=%.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top')
title('Merged both-olf, 3s')


end