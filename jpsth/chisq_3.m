function p=chisq_3(pos1,cnt1,pos2,cnt2,pos3,cnt3)
vec1=[ones(cnt1,1);...
    2*ones(cnt2,1);...
    3*ones(cnt3,1)];
vec2=[(1:cnt1).'>pos1;...
    (1:cnt2).'>pos2;...
    (1:cnt3).'>pos3];

[~,~,p]=crosstab(vec1,vec2);

end




