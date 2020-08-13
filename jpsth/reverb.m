prefix='0810_selec';
gen_conn_chain=true;
if gen_conn_chain
    load(fullfile('reg_keep.mat'));
%    reg_set=join_reg;
    
    reg_set=reg_set(~strcmp(reg_set,'Unlabeled'));
    reg_set=reg_set(~strcmp(reg_set,'root'));
    
    for bin=1:6
%        load(sprintf('0626_selec_XCORR_stats_delay_6_%d_%d_2msbin.mat',bin,bin+1));
        stats=stats_bins{bin};
        for i=1:length(stats)
            s=stats{i};
            s.uid1=s.fileidx*100000+s.su1_clusterid;
            s.uid2=s.fileidx*100000+s.su2_clusterid;
            stats{i}=s;
        end

        pair_chain=zeros(0,2);
        pair_reg=zeros(0,2);

        conn_chain_S1=zeros(0,2);
        reg_chain_S1=zeros(0,2);
        pref_chain_S1=zeros(0,12);

        conn_chain_S2=zeros(0,2);
        reg_chain_S2=zeros(0,2);
        pref_chain_S2=zeros(0,12);

        conn_chain_both=zeros(0,2);
        reg_chain_both=zeros(0,2);
        pref_chain_both=zeros(0,12);
        tbin=bin;
        for pidx=1:length(stats)
            s=stats{pidx};

            su1reg_idx=find(strcmp(s.reg_su1,reg_set));
            su2reg_idx=find(strcmp(s.reg_su2,reg_set));

            if isempty(su1reg_idx) || isempty(su2reg_idx)
                continue;
            end

            if s.s1_trials<20 || s.s2_trials<20
                continue
            end

            pair_chain(end+1,:)=[s.uid1,s.uid2];
            pair_reg(end+1,:)=[su1reg_idx,su2reg_idx];

%            if isempty(su1reg_idx) || isempty(su2reg_idx)
%                fprintf('%s, %s\n', s.reg_su1, s.reg_su2);
%                continue
%            end
            if s.totalcount<250
                continue
            end
            if s.prefered_sample_su1(bin+1) && s.prefered_sample_su2(bin+1) && s.prefered_sample_su1(bin+1)==s.prefered_sample_su2(bin+1)  
                sel_flag=true;
            else
                sel_flag=false;
            end

            %%S1
            if s.s1_peak_significant
                if s.AIs1>0.4
                    conn_chain_S1(end+1,:)=[s.uid1,s.uid2];
                    reg_chain_S1(end+1,:)=[su1reg_idx,su2reg_idx];
                    pref_chain_S1(end+1,:)=[s.prefered_sample_su1(2:end),s.prefered_sample_su2(2:end)];
                    if sel_flag
                    end
                elseif s.AIs1<-0.4
                    conn_chain_S1(end+1,:)=[s.uid2,s.uid1];
                    reg_chain_S1(end+1,:)=[su2reg_idx,su1reg_idx];
                    pref_chain_S1(end+1,:)=[s.prefered_sample_su2(2:end),s.prefered_sample_su1(2:end)];
                    if sel_flag
                    end
                end
            end
            if s.s2_peak_significant
                if s.AIs2>0.4
                    conn_chain_S2(end+1,:)=[s.uid1,s.uid2];
                    reg_chain_S2(end+1,:)=[su1reg_idx,su2reg_idx];
                    pref_chain_S2(end+1,:)=[s.prefered_sample_su1(2:end),s.prefered_sample_su2(2:end)];
                    if sel_flag
                    end
                elseif s.AIs2<-0.4
                    conn_chain_S2(end+1,:)=[s.uid2,s.uid1];
                    reg_chain_S2(end+1,:)=[su2reg_idx,su1reg_idx];
                    pref_chain_S2(end+1,:)=[s.prefered_sample_su2(2:end),s.prefered_sample_su1(2:end)];
                    if sel_flag
                    end
                end
            end
            if s.s1_peak_significant && s.s2_peak_significant 
               if  s.AIs1>0.4 & s.AIs2>0.4
                    conn_chain_both(end+1,:)=[s.uid1,s.uid2];
                    reg_chain_both(end+1,:)=[su1reg_idx,su2reg_idx];
                    pref_chain_both(end+1,:)=[s.prefered_sample_su1(2:end),s.prefered_sample_su2(2:end)];
               end

               if  s.AIs1<-0.4 & s.AIs2<-0.4
                    conn_chain_both(end+1,:)=[s.uid2,s.uid1];
                    reg_chain_both(end+1,:)=[su2reg_idx,su1reg_idx];
                    pref_chain_both(end+1,:)=[s.prefered_sample_su2(2:end),s.prefered_sample_su1(2:end)];
               end
            end
        end
        disp('check file name')
        save(sprintf('%s_conn_chain_duo_6s_%d_%d.mat',prefix,bin,bin+1),'conn_chain_S1','reg_chain_S1','pref_chain_S1','conn_chain_S2','reg_chain_S2','pref_chain_S2','conn_chain_both','reg_chain_both','pref_chain_both','pair_chain','pair_reg')    
    end
return
end

return

count2order=false;
if count2order
    done=[];
    found_2nd=[];
    reg_2nd=[];
    for i=1:length(conn_chain)
        source=conn_chain(i,1);
        if ~ismember(source,done)
            %% do chain stuff
            sec_order=conn_chain(conn_chain(:,1)==source,2);
            matches=conn_chain(ismember(conn_chain(:,1),sec_order) & conn_chain(:,2)==source,1);
            match_reg=reg_chain(ismember(conn_chain(:,1),sec_order) & conn_chain(:,2)==source,1);
            if ~isempty(matches)
                for mi=1:size(matches,1)
                    found_2nd(end+1,:)=[source,matches(mi)];
                    reg_2nd(end+1,:)=[reg_chain(i,1),match_reg(mi)];
                end
            end
            done=[done,source];
        end
    end
end

count3order=false;
if count3order
    done=[];
    found_3rd=[];
    for i=1:length(conn_chain)
        source=conn_chain(i,1);
        if ~ismember(source,done)
            %% do chain stuff
            sec_order=conn_chain(conn_chain(:,1)==source,2);
            third_order=conn_chain(ismember(conn_chain(:,1),sec_order),2);
            matches=conn_chain(ismember(conn_chain(:,1),third_order) & conn_chain(:,2)==source,:);
            if ~isempty(matches)
                for i=1:size(matches,1)
                    found_3rd(end+1,:)=[source,matches(i,:)];
                end
            end
            done=[done,source];
        end
    end
end

uids=zeros(1,0)
reg_uids=cell(1,0)
for i=1:length(stats)
    s=stats{i};
    uids(end+1)=s.uid1;
    reg_uids{end+1}=s.reg_su1;
    uids(end+1)=s.uid2;
    reg_uids{end+1}=s.reg_su2;
end

[uids,IA,IC]=unique(uids);
reg_uids=reg_uids(IA);

reg_2nd_chain=cell(0,2);
for i=1:length(found)
    reg_2nd_chain(i,:)={reg_uids(uids==found(i,1)),reg_uids(uids==found(i,2))};
end
reg_3rd_chain=cell(0,3);
for i=1:length(found_3rd)
    reg_3rd_chain(i,:)={reg_uids(uids==found_3rd(i,1)),reg_uids(uids==found_3rd(i,2)),reg_uids(uids==found_3rd(i,3))};
end



c=0;
for i=1:length(stats)
    s=stats{i};
    if strcmp(s.reg_su1,'TTv') && strcmp(s.reg_su2,'TTv')
        c=c+1;
        disp(c)
    end
end

