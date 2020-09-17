
[geplist1,reglist1]=sortdata(1,fstr);
[geplist2,reglist2]=sortdata(2,fstr);

geplist=[geplist1;geplist2];
reglist=[reglist1;reglist2];
%% export
edgecell=num2cell(geplist);
% csvcell(1,:)={'Source','Target','Bin','Weight'};
edgecell(:,5)=arrayfun(@(x) sprintf('<[%d,%d]>',x,x+1),geplist(:,3),'UniformOutput',false);
edgecell=[{'Source','Target','Label','Peak','timeset'};edgecell];
% [~,idx]=sort(cell2mat(csvcell(:,3)),'descend');
% csvcell=csvcell(idx,:);
writecell(edgecell,'chain_plot_S1_edge.csv');

load('reg_keep.mat','reg_set')
% nodecell=unique(reglist,'rows');
sess=num2cell(idivide(int32(nodecell(:,1)),int32(100000)));
reg=arrayfun(@(x) reg_set{x}, nodecell(:,2),'UniformOutput',false);
nodecell=num2cell(nodecell);
nodecell=[nodecell(:,[1,1]), reg,sess];
nodecell=[{'Id','Label','Reg','Session'};nodecell];
writecell(nodecell,'chain_plot_S1_node.csv');


function [geplist,reglist]=sortdata(prefsamp,fstr)
if prefsamp==1
    conn_chain_SC='conn_chain_S1';
    reg_chain_SC='reg_chain_S1';
    pref_chain_SC='pref_chain_S1';
    peakSC='peaks1';
elseif prefsamp==2
    conn_chain_SC='conn_chain_S2';
    reg_chain_SC='reg_chain_S2';
    pref_chain_SC='pref_chain_S2';
    peakSC='peaks2';
end
geplist=[];
reglist=[];
bin1sel=(fstr{1}.(pref_chain_SC)(:,1)==prefsamp & all(fstr{1}.(pref_chain_SC)(:,7:8)==prefsamp,2));
bin1postu=unique(fstr{1}.(conn_chain_SC)(bin1sel,2));
geplist=[geplist;fstr{1}.(conn_chain_SC)(bin1sel,:),ones(nnz(bin1sel),1),fstr{1}.(peakSC)(bin1sel)'];
reglist=[reglist;reshape(fstr{1}.(conn_chain_SC)(bin1sel,:),[],1),reshape(fstr{1}.(reg_chain_SC)(bin1sel,:),[],1)];

bin2sel=ismember(fstr{2}.(conn_chain_SC)(:,1),bin1postu) & all(fstr{2}.(pref_chain_SC)(:,8:9)==prefsamp,2);
bin2postu=unique(fstr{2}.(conn_chain_SC)(bin2sel,2));
geplist=[geplist;fstr{2}.(conn_chain_SC)(bin2sel,:),ones(nnz(bin2sel),1)*2,fstr{2}.(peakSC)(bin2sel)'];
reglist=[reglist;reshape(fstr{2}.(conn_chain_SC)(bin2sel,:),[],1),reshape(fstr{2}.(reg_chain_SC)(bin2sel,:),[],1)];

bin3sel=ismember(fstr{3}.(conn_chain_SC)(:,1),bin2postu) & all(fstr{3}.(pref_chain_SC)(:,9:10)==prefsamp,2);
bin3postu=unique(fstr{3}.(conn_chain_SC)(bin3sel,2));
geplist=[geplist;fstr{3}.(conn_chain_SC)(bin3sel,:),ones(nnz(bin3sel),1)*3,fstr{3}.(peakSC)(bin3sel)'];
reglist=[reglist;reshape(fstr{3}.(conn_chain_SC)(bin3sel,:),[],1),reshape(fstr{3}.(reg_chain_SC)(bin3sel,:),[],1)];

bin4sel=ismember(fstr{4}.(conn_chain_SC)(:,1),bin3postu) & all(fstr{4}.(pref_chain_SC)(:,10:11)==prefsamp,2);
bin4postu=unique(fstr{4}.(conn_chain_SC)(bin4sel,2));
geplist=[geplist;fstr{4}.(conn_chain_SC)(bin4sel,:),ones(nnz(bin4sel),1)*4,fstr{4}.(peakSC)(bin4sel)'];
reglist=[reglist;reshape(fstr{4}.(conn_chain_SC)(bin4sel,:),[],1),reshape(fstr{4}.(reg_chain_SC)(bin4sel,:),[],1)];

bin5sel=ismember(fstr{5}.(conn_chain_SC)(:,1),bin4postu) & all(fstr{5}.(pref_chain_SC)(:,11:12)==prefsamp,2);
geplist=[geplist;fstr{5}.(conn_chain_SC)(bin5sel,:),ones(nnz(bin5sel),1)*5,fstr{5}.(peakSC)(bin5sel)'];
reglist=[reglist;reshape(fstr{5}.(conn_chain_SC)(bin5sel,:),[],1),reshape(fstr{5}.(reg_chain_SC)(bin5sel,:),[],1)];
gplist(:,5:6)=prefsamp;
reglist=unique(reglist,'rows');
end

function unchained()


return
end