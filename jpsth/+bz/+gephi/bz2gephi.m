function bz2gephi(in_m,in_n,in_a,depth,type,opt)
arguments
    in_m (1,1) struct
    in_n (1,1) struct
    in_a (1,1) struct
    depth (1,1) double
    type (1,:) char = 'diff'
    opt.min_conn (1,1) double = 20
end
gephisel=in_m.net_count(:,1)>=opt.min_conn & in_m.netedge(:,1)~=in_m.netedge(:,2); % [pair_count, sig_count, ratio]
csvcell={'Source','Target','Weight','Diff','Chisq_p'};
keepedge=in_m.netedge(gephisel,:);
keepnode=in_m.netnode(gephisel,:);
keepcount=in_m.net_count(gephisel,:);
for i=1:size(keepedge,1)
    nidx=find(in_n.netedge(:,1)==keepedge(i,1) & in_n.netedge(:,2)==keepedge(i,2));
    aidx=find(in_a.netedge(:,1)==keepedge(i,1) & in_a.netedge(:,2)==keepedge(i,2));
    %TODO significant memory assemblies
    if ~isempty(nidx) && in_n.net_count(nidx,1)>=opt.min_conn
        [~,~,csqp]=crosstab([zeros(1,keepcount(i,1)),ones(1,in_n.net_count(nidx,1))],...
            [(1:keepcount(i,1))>keepcount(i,2),(1:in_n.net_count(nidx,1))>in_n.net_count(nidx,2)]);%chisq test
    else, continue;
    end
    
    if ~isempty(aidx)&& in_a.net_count(aidx,3)>0
        csvcell(end+1,:)={char(keepnode{i,1}),...
            char(keepnode{i,2}),...
            in_a.net_count(aidx,3),...
            (keepcount(i,3)-in_n.net_count(nidx,3))/(keepcount(i,3)+in_n.net_count(nidx,3)),...
            csqp};
    end
end
writecell(csvcell,fullfile('bzdata',sprintf('conn4gephi_%s_%d.csv',type,depth)));
end