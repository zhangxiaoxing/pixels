function pstats=rings_su_wave_tcom_corr(sums_all,opt)
arguments
    sums_all
    opt.load_file (1,1) logical = true
end
% function cong_dir=rings_su_tcom_order(sums_all) % modified from
% persistent sums_all
% if isempty(sums_all)
%     load(fullfile('bzdata','sums_ring_stats_all.mat'));
% end
if ~opt.load_file
    pstats=struct();
    [pstats.both3,pstats.olf3,pstats.dur3,pstats.both6,pstats.olf6,pstats.dur6]=deal(cell(0));

    su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
    wrs_mux_meta=ephys.get_wrs_mux_meta();
    com_map=wave.get_pct_com_map(wrs_mux_meta,'early_smooth',false);

    rstats=cell(0,11);
    for rsize=3:5
        one_rsize=sums_all{rsize-2};
        curr_sess=-1;
        for ridx=1:size(one_rsize,1)
            if curr_sess~=one_rsize{ridx,1}
                curr_sess=one_rsize{ridx,1};
                sesscid=su_meta.allcid(su_meta.sess==curr_sess);
                sesswaveid=wrs_mux_meta.wave_id(su_meta.sess==curr_sess);
                sess_wave_map=containers.Map(num2cell(sesscid),num2cell(sesswaveid));
                sess_su_com_map=com_map.(['s',num2str(curr_sess)]);
            end
            curr_waveid=cell2mat(sess_wave_map.values(num2cell(one_rsize{ridx,3})));
            %         if (~all(ismember(curr_part,{'CH','BS'}),'all'))
            %             continue
            %         end
            [rwid,seltype]=bz.rings.ring_wave_type(curr_waveid);

            if ~strcmp(rwid,'congru')
                continue
            end
            [su_tcom3,su_tcom6]=deal([]);
            if strcmp(seltype,'both')
                if ~any(ismember(curr_waveid,[2 4 8]),'all') % 3s
                    su_tcom3=nan(size(curr_waveid));
                    for ii=1:rsize
                        switch curr_waveid(ii)
                            case 1
                                su_tcom3(ii)=sess_su_com_map.s1d3.com3(one_rsize{ridx,3}(ii));
                            case 3
                                su_tcom3(ii)=sess_su_com_map.s2d3.com3(one_rsize{ridx,3}(ii));
                            case 5
                                su_tcom3(ii)=sess_su_com_map.olf_s1.com3(one_rsize{ridx,3}(ii));
                            case 6
                                su_tcom3(ii)=sess_su_com_map.olf_s2.com3(one_rsize{ridx,3}(ii));
                            case 7
                                su_tcom3(ii)=sess_su_com_map.dur_d3.com3(one_rsize{ridx,3}(ii));
                        end
                    end
                elseif ~any(ismember(curr_waveid,[1 3 7]),"all") % 6s
                    su_tcom6=nan(size(curr_waveid));
                    for ii=1:rsize
                        switch curr_waveid(ii)
                            case 2
                                su_tcom6(ii)=sess_su_com_map.s1d6.com6(one_rsize{ridx,3}(ii));
                            case 4
                                su_tcom6(ii)=sess_su_com_map.s2d6.com6(one_rsize{ridx,3}(ii));
                            case 5
                                su_tcom6(ii)=sess_su_com_map.olf_s1.com6(one_rsize{ridx,3}(ii));
                            case 6
                                su_tcom6(ii)=sess_su_com_map.olf_s2.com6(one_rsize{ridx,3}(ii));
                            case 8
                                su_tcom6(ii)=sess_su_com_map.dur_d6.com6(one_rsize{ridx,3}(ii));
                        end
                    end
                end

            elseif strcmp(seltype,'olf')
                if all(ismember(one_rsize{ridx,3},...
                        cell2mat(sess_su_com_map.olf_s1.com3.keys())),"all")
                    % 3s delay
                    su_tcom3=cell2mat(sess_su_com_map.olf_s1.com3.values(num2cell(one_rsize{ridx,3})));
                    % 6s delay
                    su_tcom6=cell2mat(sess_su_com_map.olf_s1.com6.values(num2cell(one_rsize{ridx,3})));
                elseif all(ismember(one_rsize{ridx,3},...
                        cell2mat(sess_su_com_map.olf_s2.com3.keys())),"all")
                    % 3s delay
                    su_tcom3=cell2mat(sess_su_com_map.olf_s2.com3.values(num2cell(one_rsize{ridx,3})));
                    % 6s delay]
                    su_tcom6=cell2mat(sess_su_com_map.olf_s2.com6.values(num2cell(one_rsize{ridx,3})));
                end
            elseif strcmp(seltype,'dur')
                if all(ismember(one_rsize{ridx,3},...
                        cell2mat(sess_su_com_map.dur_d3.com3.keys())),"all")
                    % 3s delay
                    su_tcom3=cell2mat(sess_su_com_map.dur_d3.com3.values(num2cell(one_rsize{ridx,3})));
                elseif all(ismember(one_rsize{ridx,3},...
                        cell2mat(sess_su_com_map.dur_d6.com3.keys())),"all")
                    % 6s delay
                    su_tcom6=cell2mat(sess_su_com_map.dur_d6.com6.values(num2cell(one_rsize{ridx,3})));
                end
            end
            rstats=[rstats;one_rsize(ridx,:),curr_waveid,seltype,{su_tcom3./4},{su_tcom6./4},rsize];
        end
    end

    rstats=rstats(cell2mat(rstats(:,6))>0.1 & cellfun(@(x) numel(unique(x)),rstats(:,3))==cell2mat(rstats(:,11)),:);
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

            if ~isempty(pref3)
                su_tcom3=rstats{ri,9};
                tssel=ismember(ts_id(:,5),pref3) & ts_id(:,4)>=1 & ts_id(:,4)<4 & ts_id(:,6)~=0;
                rsums3=ts_id(tssel,:);
                tag=[rstats{ri,8},'3'];
                pstats.(tag)=[pstats.(tag);{rsums3(:,4)},{su_tcom3}];
            end
            if ~isempty(pref6)
                su_tcom6=rstats{ri,10};
                tssel=ismember(ts_id(:,5),pref6) & ts_id(:,4)>=1 & ts_id(:,4)<7 & ts_id(:,6)~=0;
                rsums6=ts_id(tssel,:);
                tag=[rstats{ri,8},'6'];
                pstats.(tag)=[pstats.(tag);{rsums6(:,4)},{su_tcom6}];
            end
        end
    end
    % keyboard();
    blame=vcs.blame();
    save(fullfile('bzdata','loop_tcom_su_tcom_corr.mat'),'pstats','blame')
else
    load(fullfile('bzdata','loop_tcom_su_tcom_corr.mat'),'pstats')
end
sum_stats(pstats);
end

function sum_stats(pstats)
figure()
tiledlayout(2,3)
one_panel(pstats.both6,6,'both, 6s')
one_panel(pstats.olf6,6,'olf, 6s')
one_panel([pstats.both6;pstats.olf6;pstats.dur6],6,'merged, 6s')
% one_panel(pstats.olf6,6,'olf, 6s')
one_panel(pstats.both3,3,'both, 3s')
one_panel(pstats.olf3,3,'olf, 3s')
one_panel([pstats.both3;pstats.olf3;pstats.dur3],3,'merged, 3s')
% one_panel(pstats.olf3,3,'olf, 3s')
end

function one_panel(one_stat,delay,ttitle)
if delay==6
    loop_tcom=cellfun(@(x) sum(histcounts(x,1:0.25:7).*(1:24))./sum(histcounts(x,1:0.25:7))./4, one_stat(:,1));
elseif delay==3
    loop_tcom=cellfun(@(x) sum(histcounts(x,1:0.25:4).*(1:12))./sum(histcounts(x,1:0.25:4))./4, one_stat(:,1));
end
su_tcom=cellfun(@(x) mean(x,"all"), one_stat(:,2));
nexttile()
hold on
scatter(loop_tcom,su_tcom,16,'o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha',0.25);
c=[loop_tcom,ones(numel(loop_tcom),1)]\su_tcom;
minmaxx=[min(loop_tcom),max(loop_tcom)];
plot(minmaxx,c(1).*minmaxx+c(2),'k--')
[r,p]=corr(loop_tcom,su_tcom);
xlabel('Loops wave timing (s)');
ylabel('Mean neuron wave timing (s)');
text(max(xlim()),max(ylim()),sprintf('r=%.2f,p=%.3f',r,p),'HorizontalAlignment','right','VerticalAlignment','top')
title(ttitle)

end


% end

