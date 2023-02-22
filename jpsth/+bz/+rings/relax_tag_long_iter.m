% Algorithm for replacing recursive relax_tag_long function with 
% iterations due to performance concerns. 23.01.30
% TODO: Merge single spike loops search

% in_={[10:450:920,2150:450:2700,([10:450:920,2150:450:2700])+30000].',[500:450:1500,(500:450:1500)+30000].',[1550:450:2200,(1550:450:2200)+30000].'}
% bz.rings.relax_tag_long_iter(in_);
function out=relax_tag_long_iter(in,opt)
arguments
    in
    opt.burstInterval (1,1) double = 600 % ticks, at 30k sps
    opt.burst (1,1) logical = true % option to merge single spike loop algorithm
end
out=cell(0); %{[depth,id]}
already=[0,0,0,0];

rsize=size(in,2);
[~,loopIdx]=min(cellfun(@(x) x(1),in));
recDepth=1;
loopCnt=1;
perSU=ones(1,rsize);
minIdx=ones(1,rsize);

%% iterative
istack={loopIdx,recDepth,loopCnt,perSU,out,'burst'};

% TODO: where are the results?
while ~isempty(istack)
    % get parameters from stack
    loopIdx=istack{1,1};
    recDepth=istack{1,2};
    loopCnt=istack{1,3};
    perSU=istack{1,4};
    % ??
    nextSect=istack{1,6};
    istack(1,:)=[]; % pop
    switch nextSect
        case 'burst'
            link=[loopIdx,perSU(loopIdx),in{loopIdx}(perSU(loopIdx))];
            % push return
            istack=[{loopIdx,recDepth,loopCnt,perSU,link,'branch'};...
                istack];
            if ~opt.burst
                continue
            end
            if perSU(loopIdx)<=numel(in{loopIdx})-1 && in{loopIdx}(perSU(loopIdx)+1)-in{loopIdx}(perSU(loopIdx))<opt.burstInterval
                % push recursive
                perSU(loopIdx)=perSU(loopIdx)+1;
                istack=[{loopIdx,recDepth+1,loopCnt,perSU,link,'burst'};...
                    istack];
            else
                if loopCnt>rsize
                    cstack=flip(cell2mat(istack(:,5)));
                    currconn=[cstack(1:end-1,1:2),cstack(2:end,1:2)];
                    if ~all(ismember(currconn,already,'rows'))
                        out(end+1)={cstack};
                    end
                    already=unique([already;currconn],'rows');
                end
            end
        case 'branch'
            if loopIdx==rsize
                nxtd=1;
            else
                nxtd=loopIdx+1;
            end
            
            link=[loopIdx,perSU(loopIdx),in{loopIdx}(perSU(loopIdx))];
            % push return
            istack=[{loopIdx,recDepth,loopCnt,perSU,link,'housekeeping'};...
                istack];
            perSU(nxtd)=findFirst(in{nxtd},in{loopIdx}(perSU(loopIdx))+24);
            if perSU(nxtd)<=numel(in{nxtd}) && in{nxtd}(perSU(nxtd))<in{loopIdx}(perSU(loopIdx))+300 % in window
                % push recursive
                istack=[{nxtd,recDepth+1,loopCnt+1,perSU,link,'burst'};...
                    istack];
            end
        case 'housekeeping'
            %% Unencumbered by any FC
            % TODO: push?
            if recDepth==1
                nxtIdx=ones(1,rsize);
                if size(already,1)<rsize+1
                    nxtIdx(loopIdx)=perSU(loopIdx)+1;
                else
                    nxtIdx(loopIdx)=max(already(already(:,3)==loopIdx,4)+1);
                end
                for jj=setdiff(1:rsize,loopIdx)
                    remainset=setdiff(1:numel(in{jj})+1,already(already(:,3)==jj,4));
                    nxtIdx(jj)=min(remainset);
                end
                nxtIdx=max([nxtIdx;minIdx]);
                if ~all(nxtIdx<cellfun(@(x) numel(x),in),'all')
                    break % valid break of control loop
                else
                    nextTS=arrayfun(@(x) in{x}(nxtIdx(x)),1:rsize);
                    [~,loopIdx]=min(nextTS);
                    perSU=nxtIdx;
                    minIdx=nxtIdx;                    
                    minIdx(loopIdx)=minIdx(loopIdx)+1;
                    istack={loopIdx,1,1,perSU,[],'burst'};
                end
            end
    end
end
end


function idx = findFirst(nxt, curr) %nxt should be first > curr
l = 1;
r = length(nxt);
idx = 1;
while l < r
    idx = 1 + floor((l + r - 1) / 2);
    if nxt(idx) > curr
        r = idx - 1;
    elseif nxt(idx) <= curr
        l = idx;
    end
end
if l == r
    idx = r;
end
if nxt(idx) <= curr
    idx = idx+1;
end
end