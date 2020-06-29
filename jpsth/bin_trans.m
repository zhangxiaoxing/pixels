pre_post_bin=[];
for i=1:length(pref_chain_S1)
%     pre=(pref_chain_S1(i,1:6)>0).*(1:6);  %TODO should do both pref >0 and ==1 test
%     post=(pref_chain_S1(i,7:12)>0).*(1:6);
      pre1st=find(pref_chain_S1(i,1:6)>0,1);
      post1st=find(pref_chain_S1(i,7:12)>0,1);
      pre_post_bin(end+1,:)=[pre1st,post1st];
end
% doesn't really work
% scatter(pre_post_bin(:,1),pre_post_bin(:,2),20,'o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.01);

dbin=diff(pre_post_bin,1,2);
histogram(dbin)

pre_post_bin=[];
for i=1:length(pref_chain_S1)
     pre=(pref_chain_S1(i,1:6)>0)*(1:6)';  %TODO should do both pref >0 and ==1 test
     post=(pref_chain_S1(i,7:12)>0)*(1:6)';
      pre_post_bin(end+1,:)=[pre,post];
end
% doesn't really work
% scatter(pre_post_bin(:,1),pre_post_bin(:,2),20,'o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.01);

dbin=diff(pre_post_bin,1,2);
histogram(dbin)




pre_post_bin=[];
for i=1:length(pref_chain_S2)
%     pre=(pref_chain_S1(i,1:6)>0).*(1:6);  %TODO should do both pref >0 and ==1 test
%     post=(pref_chain_S1(i,7:12)>0).*(1:6);
      pre1st=find(pref_chain_S2(i,1:6)>0,1);
      post1st=find(pref_chain_S2(i,7:12)>0,1);
      pre_post_bin(end+1,:)=[pre1st,post1st];
end
% doesn't really work
% scatter(pre_post_bin(:,1),pre_post_bin(:,2),20,'o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.01);

dbin=diff(pre_post_bin,1,2);
histogram(dbin)