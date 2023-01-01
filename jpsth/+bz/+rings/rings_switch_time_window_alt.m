% TODO: loop-loop switch

function single_su_multi_ring=rings_switch_time_window_alt(pstats)
%data file from rings_time_constant.m
if false
    load(fullfile('bzdata','rings_spike_trial_tagged.mat'),'pstats')
    rstats.congru=pstats.congru;
    pstats=rstats;
    clear rstats
end

% for each session
single_su_multi_ring=struct();

all_rings=fieldnames(pstats.congru);
waveids={[1 5],[2 5],[3 6],[4 6],[1 7],[2 8],[3 7],[4 8]};
for wid=1:8
    waveid=waveids{wid};
    for sessid=33%[100,18,33]
        sess=str2double(string(regexp(all_rings,'(?<=s)\d+(?=r.*)','match')));
        sess_ring=all_rings(sess==sessid);
        trials=pstats.congru.(sess_ring{1}).trials;
        if ismember(3,waveid)
            trial_sel=find(trials(:,5)==8 & trials(:,8)==3 & all(trials(:,9:10)>0,2));
        elseif ismember(1,waveid)
            trial_sel=find(trials(:,5)==4 & trials(:,8)==3 & all(trials(:,9:10)>0,2));
        elseif ismember(2,waveid)
            trial_sel=find(trials(:,5)==4 & trials(:,8)==6 & all(trials(:,9:10)>0,2));
        elseif ismember(4,waveid)
            trial_sel=find(trials(:,5)==8 & trials(:,8)==6 & all(trials(:,9:10)>0,2));
        end

        for rr=reshape(sess_ring,1,[])
            if ~all(ismember(pstats.congru.(rr{1}).rstats{4},waveid),2)
                continue
            end
            ts_id=pstats.congru.(rr{1}).ts_id;
            rcids=pstats.congru.(rr{1}).rstats{3};
            for su=rcids
                sutag="s"+sessid+"w"+wid+"u"+num2str(su);
                if ~isfield(single_su_multi_ring,sutag)
                    single_su_multi_ring.(sutag)=ts_id(ts_id(:,2)==su & ismember(ts_id(:,5),trial_sel),4:5);
                end
                single_su_multi_ring.(sutag)(:,end+1)=ts_id(ts_id(:,2)==su & ismember(ts_id(:,5),trial_sel) ,6);
            end
        end
    end
end

end



function overlap_ring=overlap_map(pstats,sess_ring)
    overlap_ring=struct();
    for rr=reshape(sess_ring,1,[])
        overlap_ring.(rr{1})=cell(0);
        ring_cid=pstats.congru.(rr{1}).rstats{3}; % skip per-wave seperation, no intersect anyway
        for rrm=reshape(sess_ring,1,[])
            if isequal(rr,rrm)
                continue
            end
            ov_cid=pstats.congru.(rrm{1}).rstats{3}; % skip per-wave seperation, no intersect anyway
            if any(ismember(ring_cid,ov_cid),2)
                overlap_ring.(rr{1})=[overlap_ring.(rr{1}),rrm];
            end
        end
    end
end
 