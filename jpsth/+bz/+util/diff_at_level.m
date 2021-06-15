function [is_diff,is_same]=diff_at_level(reg)
is_diff=false(size(reg,1),6);
is_same=false(size(reg,1),6);
graysel=all(reg(:,1,:)==567 | reg(:,1,:)==343,3);
for dep=1:6
    is_diff(:,dep)=graysel & reg(:,dep,1)~=reg(:,dep,2);
    is_same(:,dep)=graysel & reg(:,dep,1)==reg(:,dep,2);
end
end