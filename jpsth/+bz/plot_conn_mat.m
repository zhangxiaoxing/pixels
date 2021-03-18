function fh=plot_conn_mat(in,depth,type)
arguments
    in (1,1) struct
    depth (1,1) double
    type (1,:) char
end
idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));
uccfid=unique(in.netedge(:));
regjoin=@(x) strjoin(x(3:end),'-');
unodes=arrayfun(@(x) regjoin(idmap.reg2tree(char(idmap.ccfid2reg(x)))),uccfid,'UniformOutput',false);
leafnodes=arrayfun(@(x) char(idmap.ccfid2reg(x)),uccfid,'UniformOutput',false);
[unodes,idx]=sort(unodes);
leafnodes=leafnodes(idx);
uccfid=uccfid(idx)';
conn_mat=nan(numel(uccfid));
for i=1:numel(uccfid)
    for j=1:numel(uccfid)
        connsel=find(in.netedge(:,1)==uccfid(i) & in.netedge(:,2)==uccfid(j));
        if ~isempty(connsel)
            conn_mat(i,j)=in.net_count(connsel,3);
        end
    end
end
% set(groot, 'DefaultTextFontSize', 20);
fh=figure('Color','w','Position',[100,100,850,600]);
imagesc(conn_mat.*100,[0,3])
colormap('jet')
set(gca(),'YDir','normal',...
    'XTick',1:numel(unodes),'YTick',1:numel(unodes),...
    'XTickLabel',leafnodes,...
    'XTickLabelRotation',90,...
    'YTickLabel',unodes);
xlabel('Target')
ylabel('Source')
cb=colorbar();
cb.Label.String='Connected fraction (%)';
title(sprintf('%s net, level %d',type, depth));

end