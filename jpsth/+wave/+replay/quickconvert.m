function rings_replay_tbl=quickconvert(ring_replay)
rings_replay_tbl=[];
for dd=["d3","d6"]
    for ss=["olf_s1","olf_s2"]
        for fn=reshape(fieldnames(ring_replay.(dd).(ss)),1,[])
            oner=ring_replay.(dd).(ss).(fn{1});
            delay=str2double(replace(dd,"d",""));
            wid=replace(ss,"olf_","")+dd;
            cid=string(fn{1});
            c={oner.meta{1},delay,wid,cid,oner.ts,oner.meta,[],oner.trl_align,oner.freqstats};
            rings_replay_tbl=[rings_replay_tbl;cell2table(c,'VariableNames',{'session','delay','wave','ring_id','ts','meta','ts_id','trl_align','freqstats'})];
        end
    end
end
end

