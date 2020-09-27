sample=2;
for bin=1:6
    fstr{bin}=load(sprintf('0831_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
end
load hebb_pattern.mat
if sample==1
    hebbPattern=hebbPatternS1;
elseif sample==2
    hebbPattern=hebbPatternS2;
else
    keyboard
end
peak_stats=nan(size(hebbPattern,1),11);
count_stats=nan(size(hebbPattern,1),11);
for i=1:length(hebbPattern)
    if ~rem(i,1000)
        disp(i)
    end
    bins=hebbPattern{i,4};
    [cc,pp,tt]=sel_data(fstr,bins(1),sample);
    peak_stats(i,1:3)=[pp(cc(:,1)==hebbPattern{i,1}(1) & cc(:,2)==hebbPattern{i,1}(2)),...
        pp(cc(:,1)==hebbPattern{i,1}(2) & cc(:,2)==hebbPattern{i,1}(3)),...
        pp(cc(:,1)==hebbPattern{i,1}(3) & cc(:,2)==hebbPattern{i,1}(1))];
    count_stats(i,1:3)=[tt(cc(:,1)==hebbPattern{i,1}(1) & cc(:,2)==hebbPattern{i,1}(2)),...
        tt(cc(:,1)==hebbPattern{i,1}(2) & cc(:,2)==hebbPattern{i,1}(3)),...
        tt(cc(:,1)==hebbPattern{i,1}(3) & cc(:,2)==hebbPattern{i,1}(1))];
    
    [cc,pp,tt]=sel_data(fstr,bins(2),sample);
    peak_stats(i,4:7)=[pp(cc(:,1)==hebbPattern{i,2}(1) & cc(:,2)==hebbPattern{i,2}(2)),...
        pp(cc(:,1)==hebbPattern{i,2}(2) & cc(:,2)==hebbPattern{i,2}(3)),...
        pp(cc(:,1)==hebbPattern{i,2}(3) & cc(:,2)==hebbPattern{i,2}(4)),...
        pp(cc(:,1)==hebbPattern{i,2}(4) & cc(:,2)==hebbPattern{i,2}(1))];
    
    count_stats(i,4:7)=[tt(cc(:,1)==hebbPattern{i,2}(1) & cc(:,2)==hebbPattern{i,2}(2)),...
        tt(cc(:,1)==hebbPattern{i,2}(2) & cc(:,2)==hebbPattern{i,2}(3)),...
        tt(cc(:,1)==hebbPattern{i,2}(3) & cc(:,2)==hebbPattern{i,2}(4)),...
        tt(cc(:,1)==hebbPattern{i,2}(4) & cc(:,2)==hebbPattern{i,2}(1))];
    
    [cc,pp,tt]=sel_data(fstr,bins(3),sample);
    peak_stats(i,8:11)=[pp(cc(:,1)==hebbPattern{i,3}(1) & cc(:,2)==hebbPattern{i,3}(2)),...
        pp(cc(:,1)==hebbPattern{i,3}(2) & cc(:,2)==hebbPattern{i,3}(3)),...
        pp(cc(:,1)==hebbPattern{i,3}(3) & cc(:,2)==hebbPattern{i,3}(4)),...
        pp(cc(:,1)==hebbPattern{i,3}(4) & cc(:,2)==hebbPattern{i,3}(1))];
    
    count_stats(i,8:11)=[tt(cc(:,1)==hebbPattern{i,3}(1) & cc(:,2)==hebbPattern{i,3}(2)),...
        tt(cc(:,1)==hebbPattern{i,3}(2) & cc(:,2)==hebbPattern{i,3}(3)),...
        tt(cc(:,1)==hebbPattern{i,3}(3) & cc(:,2)==hebbPattern{i,3}(4)),...
        tt(cc(:,1)==hebbPattern{i,3}(4) & cc(:,2)==hebbPattern{i,3}(1))];
end

mtt=mean(count_stats,2);
mm=mean(peak_stats,2);
% miin=min(peak_stats,[],2);
bins=cell2mat(hebbPattern(:,4));
[mm,idces]=sort(mm,'descend');
if sample==1
    cmp_stats_s1=[idces,mtt(idces),mm,bins(idces,:)];
    % cmp_stats_sel=cmp_stats(cmp_stats(:,3)>3.5 & cmp_stats(:,2)>6,:);
    % crs_bin=cmp_stats(cmp_stats(:,4)~=cmp_stats(:,5) & cmp_stats(:,5)~=cmp_stats(:,6),:);
    save('hebb_showcase_stats.mat','cmp_stats_s1','-append')
elseif sample==2
    cmp_stats_s2=[idces,mtt(idces),mm,bins(idces,:)];
    save('hebb_showcase_stats.mat','cmp_stats_s2','-append')
end

function [cc,pp,tt]=sel_data(fstr,bin,sample)
    if sample==1
        cc=fstr{bin}.conn_chain_S1;
        pp=fstr{bin}.peaks1;
        tt=fstr{bin}.totalcount_S1;
    else
        cc=fstr{bin}.conn_chain_S2;
        pp=fstr{bin}.peaks2;
        tt=fstr{bin}.totalcount_S2;
    end
end