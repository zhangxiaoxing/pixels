%for i=1:length(stats)
%    s=stats{i};
%    s.uid1=s.fileidx*1000+s.su1_label_idx;
%    s.uid2=s.fileidx*1000+s.su2_label_idx;
%    stats{i}=s;
%end
%
currbin=1;
gen_conn_chain=false;
if gen_conn_chain
    keyboard
%    conn_chain=zeros(0,2);
%    reg_chain=cell(0,2);
    if ~exist('join_reg_set','var')
        load(fullfile('..','join_reg_set.mat'));
        reg_set=join_reg_set;
    end
    
    reg_set=reg_set(~strcmp(reg_set,'Unlabeled'));
    reg_set=reg_set(~strcmp(reg_set,'root'));
    
    for bin=currbin
%         load(sprintf('XCORR_stats_delay_6_%d_%d_2msbin.mat',bin,bin+1));
        conn_chain=zeros(0,2);
        reg_chain=zeros(0,2);
        pref_chain=zeros(0,12);
        tbin=bin;
        for pidx=1:length(stats)
            s=stats{pidx};
            if s.totalcount<1000 || s.s1_trials<20 || s.s2_trials<20 || strcmp(s.reg_su1,'Unlabeled') || strcmp(s.reg_su2,'Unlabeled') || strcmp(s.reg_su1,'root') || strcmp(s.reg_su2,'root')
                continue
            end
            
            su1reg_idx=find(strcmp(s.reg_su1,reg_set));
            su2reg_idx=find(strcmp(s.reg_su2,reg_set));
            if isempty(su1reg_idx) || isempty(su2reg_idx)
                fprintf('%s, %s\n', s.reg_su1, s.reg_su2);
                %keyboard
                continue
            end
            
%             not helpful in selective or non-selective type of datasets
%             if ~any(s.prefered_sample_su1(2:end)) || ~any(s.prefered_sample_su2(2:end))
%                 continue
%             end
            
            if s.prefered_sample_su1(currbin+1) && s.prefered_sample_su2(currbin+1) && s.prefered_sample_su1(currbin+1)==s.prefered_sample_su2(currbin+1)  
                sel_flag=true;
            else
                sel_flag=false;
            end

            if s.s1_peak_significant && s.s2_peak_significant
                if (s.AIs1>0 && s.AIs2>=0) || (s.AIs1>=0 && s.AIs2>0) %2 to 1
                    conn_chain(end+1,:)=[s.uid1,s.uid2];
                    reg_chain(end+1,:)=[su1reg_idx,su2reg_idx];
                    pref_chain(end+1,:)=[s.prefered_sample_su1(2:end),s.prefered_sample_su2(2:end)];
                elseif (s.AIs1<0 && s.AIs2<=0) || (s.AIs1<=0 && s.AIs2<0)
                    conn_chain(end+1,:)=[s.uid2,s.uid1];
                    reg_chain(end+1,:)=[su2reg_idx,su1reg_idx];
                    pref_chain(end+1,:)=[s.prefered_sample_su2(2:end),s.prefered_sample_su1(2:end)];
                elseif (s.AIs1<0 && s.AIs2>0) || (s.AIs1>0 && s.AIs2<0)
                    conn_chain(end+1,:)=[s.uid1,s.uid2];
                    reg_chain(end+1,:)=[su1reg_idx,su2reg_idx];
                    pref_chain(end+1,:)=[s.prefered_sample_su1(2:end),s.prefered_sample_su2(2:end)];
                    conn_chain(end+1,:)=[s.uid2,s.uid1];
                    reg_chain(end+1,:)=[su2reg_idx,su1reg_idx];
                    pref_chain(end+1,:)=[s.prefered_sample_su2(2:end),s.prefered_sample_su1(2:end)];
                    if sel_flag
                    end
                end
            elseif s.s1_peak_significant
                if s.AIs1>0 % su2 to su1
                    conn_chain(end+1,:)=[s.uid1,s.uid2];
                    reg_chain(end+1,:)=[su1reg_idx,su2reg_idx];
                    pref_chain(end+1,:)=[s.prefered_sample_su1(2:end),s.prefered_sample_su2(2:end)];
                elseif s.AIs1<0  %=0 is possible
                    conn_chain(end+1,:)=[s.uid2,s.uid1];
                    reg_chain(end+1,:)=[su2reg_idx,su1reg_idx];
                    pref_chain(end+1,:)=[s.prefered_sample_su2(2:end),s.prefered_sample_su1(2:end)];
                end
                
            elseif s.s2_peak_significant
                if s.AIs2>0 % su2 to su1
                    conn_chain(end+1,:)=[s.uid1,s.uid2];
                    reg_chain(end+1,:)=[su1reg_idx,su2reg_idx];
                    pref_chain(end+1,:)=[s.prefered_sample_su1(2:end),s.prefered_sample_su2(2:end)];
                elseif s.AIs2<0  %=0 is possible
                    conn_chain(end+1,:)=[s.uid2,s.uid1];
                    reg_chain(end+1,:)=[su2reg_idx,su1reg_idx];
                    pref_chain(end+1,:)=[s.prefered_sample_su2(2:end),s.prefered_sample_su1(2:end)];
                end
            end
                
        end
    end
disp('check file name')
% keyboard
save(sprintf('%s_conn_chain_duo_6s_%d_%d.mat',prefix,bin_range(1),bin_range(2)),'conn_chain','reg_chain','pref_chain')    
return
end

gen_conn_mat=false;
if gen_conn_mat
    conn_mat_all=cell(0);
    if ~exist('join_reg_set','var')
        load(fullfile('..','join_reg_set.mat'));
        reg_set=join_reg_set;
    end
    
    reg_set=reg_set(~strcmp(reg_set,'Unlabeled'));
    reg_set=reg_set(~strcmp(reg_set,'root'));
    
    for bin=currbin
%         load(sprintf('XCORR_stats_delay_6_%d_%d_2msbin.mat',bin,bin+1));
        conn_mat=zeros(length(reg_set),length(reg_set));
        conn_sel_mat=zeros(length(reg_set),length(reg_set));
        tbin=bin;
        for pidx=1:length(stats)
            s=stats{pidx};
            if s.totalcount<1000 || s.s1_trials<20 || s.s2_trials<20 || strcmp(s.reg_su1,'Unlabeled') || strcmp(s.reg_su2,'Unlabeled') || strcmp(s.reg_su1,'root') || strcmp(s.reg_su2,'root')
                continue
            end
            
            su1reg_idx=find(strcmp(s.reg_su1,reg_set));
            su2reg_idx=find(strcmp(s.reg_su2,reg_set));
            if isempty(su1reg_idx) || isempty(su2reg_idx)
                fprintf('%s, %s\n', s.reg_su1, s.reg_su2);
                %keyboard
                continue
            end
            
%             not helpful in selective or non-selective type of datasets
%             if ~any(s.prefered_sample_su1(2:end)) || ~any(s.prefered_sample_su2(2:end))
%                 continue
%             end
            
            if s.prefered_sample_su1(currbin+1) && s.prefered_sample_su2(currbin+1) && s.prefered_sample_su1(currbin+1)==s.prefered_sample_su2(currbin+1)  
                sel_flag=true;
            else
                sel_flag=false;
            end

            if s.s1_peak_significant && s.s2_peak_significant
                if (s.AIs1>0 && s.AIs2>=0) || (s.AIs1>=0 && s.AIs2>0) %2 to 1
                    conn_mat(su1reg_idx,su2reg_idx)=conn_mat(su1reg_idx,su2reg_idx)+1;
                    if sel_flag
                        conn_sel_mat(su1reg_idx,su2reg_idx)=conn_sel_mat(su1reg_idx,su2reg_idx)+1;
                    end
                elseif (s.AIs1<0 && s.AIs2<=0) || (s.AIs1<=0 && s.AIs2<0)
                    conn_mat(su2reg_idx,su1reg_idx)=conn_mat(su2reg_idx,su1reg_idx)+1;
                    if sel_flag
                        conn_sel_mat(su2reg_idx,su1reg_idx)=conn_sel_mat(su2reg_idx,su1reg_idx)+1;
                    end                    
                elseif (s.AIs1<0 && s.AIs2>0) || (s.AIs1>0 && s.AIs2<0)
                    conn_mat(su2reg_idx,su1reg_idx)=conn_mat(su2reg_idx,su1reg_idx)+1;
                    conn_mat(su1reg_idx,su2reg_idx)=conn_mat(su1reg_idx,su2reg_idx)+1;
                    if sel_flag
                        conn_sel_mat(su1reg_idx,su2reg_idx)=conn_sel_mat(su1reg_idx,su2reg_idx)+1;
                        conn_sel_mat(su2reg_idx,su1reg_idx)=conn_sel_mat(su2reg_idx,su1reg_idx)+1;
                    end
                end
            elseif s.s1_peak_significant
                if s.AIs1>0 % su2 to su1
                    conn_mat(su1reg_idx,su2reg_idx)=conn_mat(su1reg_idx,su2reg_idx)+1;
                    if sel_flag
                        conn_sel_mat(su1reg_idx,su2reg_idx)=conn_sel_mat(su1reg_idx,su2reg_idx)+1;
                    end                    
                elseif s.AIs1<0  %=0 is possible
                    conn_mat(su2reg_idx,su1reg_idx)=conn_mat(su2reg_idx,su1reg_idx)+1;
                    if sel_flag
                        conn_sel_mat(su2reg_idx,su1reg_idx)=conn_sel_mat(su2reg_idx,su1reg_idx)+1;
                    end
                end
                
            elseif s.s2_peak_significant
                if s.AIs2>0 % su2 to su1
                    conn_mat(su1reg_idx,su2reg_idx)=conn_mat(su1reg_idx,su2reg_idx)+1;
                    if sel_flag
                        conn_sel_mat(su1reg_idx,su2reg_idx)=conn_sel_mat(su1reg_idx,su2reg_idx)+1;
                    end
                elseif s.AIs2<0  %=0 is possible
                    conn_mat(su2reg_idx,su1reg_idx)=conn_mat(su2reg_idx,su1reg_idx)+1;
                    if sel_flag
                        conn_sel_mat(su2reg_idx,su1reg_idx)=conn_sel_mat(su2reg_idx,su1reg_idx)+1;
                    end
                end
            end
                
        end
    end
disp('check file name')
% keyboard
save(sprintf('%s_conn_mat_duo_6s_%d_%d.mat',prefix,bin_range(1),bin_range(2)),'conn_mat','conn_sel_mat')    
% return
end

count2order=true;
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

keyboard
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

