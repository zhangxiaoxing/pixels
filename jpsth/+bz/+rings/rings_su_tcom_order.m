function cong_dir=rings_su_tcom_order()
persistent sums_all
if isempty(sums_all)
    load(fullfile('bzdata','sums_ring_stats_all.mat'));
end
% pstats=struct();
% [pstats.both3,pstats.olf3,pstats.dur3,pstats.both6,pstats.olf6,pstats.dur6]=deal(cell(0));

cong_dir=struct();
[cong_dir.both3,cong_dir.both6,cong_dir.olf3,cong_dir.olf6,cong_dir.dur3,cong_dir.dur6]=deal([]);

su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
wrs_mux_meta=ephys.get_wrs_mux_meta();
com_map=wave.get_pct_com_map(wrs_mux_meta,'curve',true,'early_smooth',false);

for rsize=3:5
%     rstats=cell(0,10);
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
                dir_congru3=[arrayfun(@(x) su_tcom3(x)>su_tcom3(x-1),2:rsize),su_tcom3(1)>su_tcom3(end)];
                cong_dir.both3=[cong_dir.both3;rsize,mean(dir_congru3,'all')];
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
                dir_congru6=[arrayfun(@(x) su_tcom6(x)>su_tcom6(x-1),2:rsize),su_tcom6(1)>su_tcom6(end)];
                cong_dir.both6=[cong_dir.both6;rsize,mean(dir_congru6,'all')];
            else
                keyboard();
            end

        elseif strcmp(seltype,'olf')
            if all(ismember(one_rsize{ridx,3},...
                    cell2mat(sess_su_com_map.olf_s1.com3.keys())),"all")
                % 3s delay
                su_tcom3=cell2mat(sess_su_com_map.olf_s1.com3.values(num2cell(one_rsize{ridx,3})));
                dir_congru3=[arrayfun(@(x) su_tcom3(x)>su_tcom3(x-1),2:rsize),su_tcom3(1)>su_tcom3(end)];
                cong_dir.olf3=[cong_dir.olf3;rsize,mean(dir_congru3,'all')];
                % 6s delay
                su_tcom6=cell2mat(sess_su_com_map.olf_s1.com6.values(num2cell(one_rsize{ridx,3})));
                dir_congru6=[arrayfun(@(x) su_tcom6(x)>su_tcom6(x-1),2:rsize),su_tcom6(1)>su_tcom6(end)];
                cong_dir.olf6=[cong_dir.olf6;rsize,mean(dir_congru6,'all')];
            elseif all(ismember(one_rsize{ridx,3},...
                    cell2mat(sess_su_com_map.olf_s2.com3.keys())),"all")
                % 3s delay
                su_tcom3=cell2mat(sess_su_com_map.olf_s2.com3.values(num2cell(one_rsize{ridx,3})));
                dir_congru3=[arrayfun(@(x) su_tcom3(x)>su_tcom3(x-1),2:rsize),su_tcom3(1)>su_tcom3(end)];
                cong_dir.olf3=[cong_dir.olf3;rsize,mean(dir_congru3,'all')];
                % 6s delay]
                su_tcom6=cell2mat(sess_su_com_map.olf_s2.com6.values(num2cell(one_rsize{ridx,3})));
                dir_congru6=[arrayfun(@(x) su_tcom6(x)>su_tcom6(x-1),2:rsize),su_tcom6(1)>su_tcom6(end)];
                cong_dir.olf6=[cong_dir.olf6;rsize,mean(dir_congru6,'all')];
            else
                keyboard();
            end
        elseif strcmp(seltype,'dur')
            if all(ismember(one_rsize{ridx,3},...
                    cell2mat(sess_su_com_map.dur_d3.com3.keys())),"all")
                % 3s delay
                su_tcom3=cell2mat(sess_su_com_map.dur_d3.com3.values(num2cell(one_rsize{ridx,3})));
                dir_congru3=[arrayfun(@(x) su_tcom3(x)>su_tcom3(x-1),2:rsize),su_tcom3(1)>su_tcom3(end)];
                cong_dir.dur3=[cong_dir.dur3;rsize,mean(dir_congru3,'all')];
            elseif all(ismember(one_rsize{ridx,3},...
                    cell2mat(sess_su_com_map.dur_d6.com3.keys())),"all")
                % 6s delay
                su_tcom6=cell2mat(sess_su_com_map.dur_d6.com6.values(num2cell(one_rsize{ridx,3})));
                dir_congru6=[arrayfun(@(x) su_tcom6(x)>su_tcom6(x-1),2:rsize),su_tcom6(1)>su_tcom6(end)];
                cong_dir.dur6=[cong_dir.dur6;rsize,mean(dir_congru6,'all')];
            else
                keyboard();
            end
        end
    end
end
end

function plot()
figure()
bar([mean(cong_dir.both3(:,2)),mean(cong_dir.both6(:,2)),mean(cong_dir.olf3(:,2)),mean(cong_dir.olf6(:,2))])
yline(0.5,'k--')
set(gca(),'XTick',1:4,'XTickLabel',{'Both-3s','Both-6s','Olf.-3s','Olf.-6s'},'YTick',0:0.25:0.75,'YTickLabel',0:25:75)
ylim([0,0.8])
ylabel('Loop-order follows neuron-TCOM (%)')
end

function shuffle()

rr=rand(10000,3);
shuf3=[rr(:,2)>rr(:,1),rr(:,3)>rr(:,2),rr(:,1)>rr(:,3)];
mean(shuf3,'all')

rr=rand(10000,4);
shuf4=[rr(:,2)>rr(:,1),rr(:,3)>rr(:,2),rr(:,4)>rr(:,3),rr(:,1)>rr(:,4)];
mean(shuf4,'all')

rr=rand(10000,5);
shuf5=[rr(:,2)>rr(:,1),rr(:,3)>rr(:,2),rr(:,4)>rr(:,3),rr(:,5)>rr(:,4),rr(:,1)>rr(:,5)];
mean(shuf5,'all')

end

function obsolete
    rstats=rstats(cell2mat(rstats(:,6))>0.1 & cellfun(@(x) numel(unique(x)),rstats(:,3))==rsize,:);
    usess=unique(cell2mat(rstats(:,1)));
    
    for sessid=usess.'
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
% end

