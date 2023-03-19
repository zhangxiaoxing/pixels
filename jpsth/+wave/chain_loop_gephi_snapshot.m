[mmax,ttrial]=wave.chain_loop_SC_spk(false,100,4,'skip_plot',false);
% predefined layout using gephi and Ordered-Graph-Layout
gephilayout=jsondecode(fileread(fullfile("+gephi","SS_chain_loop_230315.json")));

%TODO: node coord from previously exported json
nodecell=[{'Id','Label','XX','YY'};...
    num2cell([ucid,ucid,(numel(ucid):-1:1).',reg_prop(:,1),topopos(tcomidx)])];
writecell(nodecell,fullfile('bzdata',sprintf('SingleSpkNode4gephi_s%dm_%d.csv',sessid,maxid)));

% Source, Target, Interval, Activity-pattern-spike-Label
csvcell=[{'Source','Target','Bin','PatternLabel'};num2cell(allconns)];
writecell(csvcell,fullfile('bzdata',sprintf('SingleSpkConn4gephi_s%dm_%d.csv',sessid,maxid)));

spkmat=cell2mat(mmax(:,4));

mminTS=min(spkmat(:,2));
mmaxTS=max(spkmat(:,2));
% snapshot bins
% 10ms one-off test
binw=10;

binidx=1;
conmat=[];
concell={};
idxmap=containers.Map('KeyType','char','ValueType','double');
for onset=mminTS:(binw*30):mmaxTS
    %>=onset, <onset+300
    for patIdx=1:size(mmax,1)
        tsmat=mmax{patIdx,4};
        for ii=2:size(tsmat,1)
            if tsmat(ii,2)>=onset && tsmat(ii,2)<onset+binw*10
                conmat=[conmat;tsmat(ii-1,1),tsmat(ii,1),ii-1];
                key=sprintf('%d_%d',tsmat(ii-1,1),tsmat(ii,1));
                lbl=regexp(mmax{patIdx,2},'(?<=s\d*)[cr]\d*','match','once');
                
                if idxmap.isKey(key)
                    idxmap(key)=strjoin({idxmap(key),''})
                else
                    idxmap(key)=
                end
            end
        end
    end
    binidx=binidx+1;
end
    

                                                                                                                                                          