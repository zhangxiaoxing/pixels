sessid=100;
maxid=4;
[mmax,ttrial,meta]=wave.chain_loop_SC_spk(false,sessid,maxid,'skip_plot',true);
% predefined layout using gephi and Ordered-Graph-Layout
%% nodes
gephilayout=jsondecode(fileread(fullfile("+gephi","SS_chain_loop_230315.json")));
nodecell={'Id','Label','XX','YY'};
for nn=1:numel(gephilayout.nodes)
    nodecell=[nodecell;...
        {gephilayout.nodes(nn).attributes.label,...
        gephilayout.nodes(nn).attributes.label,...
        gephilayout.nodes(nn).attributes.x,...
        gephilayout.nodes(nn).attributes.y}];
end
writecell(nodecell,fullfile('bzdata',sprintf('SnapShotNode4gephi_s%dm_%d.csv',sessid,maxid)));

%% edges
spkmat=cell2mat(mmax(:,4));
mminTS=min(spkmat(:,2));
mmaxTS=max(spkmat(:,2));

% snapshot bins
% 10ms one-off test
binw=10;
binidx=1;
conmat=[];
concell={};
idxmap=struct();
for onset=mminTS:(binw*30):mmaxTS
    idxmap.("B"+binidx)=containers.Map('KeyType','char','ValueType','any');
    %>=onset, <onset+300
    for patIdx=1:size(mmax,1)
        tsmat=mmax{patIdx,4};
        for ii=2:size(tsmat,1)
            if tsmat(ii,2)>=onset && tsmat(ii,2)<onset+binw*30
                conmat=[conmat;tsmat(ii-1,1),tsmat(ii,1),ii-1];
                key=sprintf('%d_%d',tsmat(ii-1,1),tsmat(ii,1));
                ptag=regexp(mmax{patIdx,2},'(?<=s\d*)[cr]\d*','match','once');
                lbl=string(ptag{1})+"#"+num2str(ii-1);
                if idxmap.("B"+binidx).isKey(key)
                    idxmap.("B"+binidx)(key)=idxmap.("B"+binidx)(key)+", "+lbl;
                else
                    idxmap.("B"+binidx)(key)=lbl;
                end
            end
        end
    end
    binidx=binidx+1;
end
    
% write conn spreadsheet with maps

csvcell={'Source','Target','TimeSet','Label'};
for bbinidx=1:binidx-1
    binkeys=idxmap.("B"+bbinidx).keys;
    for kk=1:numel(binkeys)
        source=regexp(binkeys{kk},'\d*(?=_\d*)','match','once');
        target=regexp(binkeys{kk},'(?<=\d*_)\d*','match','once');
%         time=sprintf('<[%d,%d)>',[mminTS,mminTS+binw*30]+(bbinidx-1)*binw*30);
        time=sprintf('<[%d,%d)>',bbinidx,bbinidx+1);
        label=idxmap.("B"+bbinidx)(binkeys{kk});
        csvcell=[csvcell;...
            {source,target,time,label}];
    end
end
writecell(csvcell,fullfile('bzdata',sprintf('SnapShotConn4gephi_s%dm_%d.csv',sessid,maxid)));


