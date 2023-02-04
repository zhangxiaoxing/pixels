% TODO: loop-loop switch

function rings_switch_time_window(sums_all)
%data file from rings_time_constant.m
if true
    load(fullfile('bzdata','rings_spike_trial_tagged.mat'),'pstats')
    rstats.congru=pstats.congru;
    pstats=rstats;
    clear rstats
end

% for each session
fn=fieldnames(pstats.congru);
sess=str2double(string(regexp(fn,'(?<=s)\d+(?=r.*)','match')));
for sessid=33%[100,18,33]
    sess_ring=fn(sess==sessid);
    %% build overlap-ring map
    overlap_ring=overlap_map(pstats,sess_ring);
    %%  for each ring
    trials=pstats.congru.(sess_ring{1}).trials;
    for rr=reshape(sess_ring,1,[])
    %       for each trial
    %           for each loop activity
    %               for time-window
    %                   self repeat count
    %                   for each over-lap-ring
    %                       repeat count
    end
end
% sum-up and histograph
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
 