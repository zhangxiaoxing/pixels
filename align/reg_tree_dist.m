function dist = reg_tree_dist(reg)
%REG_TREE_DIST Summary of this function goes here
%   Detailed explanation goes here

arguments
    reg (:,6,2) int32
end

% idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));
dist=nan(size(reg,1),1);
for i=1:size(reg,1)
    if all(ismember(reg(i,1,:),[343,567]))
        mx_eql=find(reg(i,:,1)==reg(i,:,2) & ~any(reg(i,:,:)==0,3),1,'last');
        if isempty(mx_eql)
            dist(i)=6;
        else
            dist(i)=6-mx_eql;
        end
    else
        dist(i)=-1;
    end
end

end

