%% Repurosed for burst-loop algorithm as of 23.01.30
function out=relax_tag_long(in,loopIdx,recDepth,loopCnt,perSU,opt)
arguments
    in
    loopIdx % current neuron idx
    recDepth
    loopCnt % accumated loop neuron count
    perSU
    opt.burstInterval (1,1) double = 600 % ticks, at 30k sps
    opt.burst (1,1) logical = true % option to merge single spike loop algorithm
end
out=cell(0); %{[depth,id]}
rsize=size(in,2);
if ~isempty(perSU) && rem(perSU(1),1000)==0
    disp([perSU(1),numel(in{1})]);
end

if isempty(loopIdx)
    [~,loopIdx]=min(cellfun(@(x) x(1),in));
    recDepth=1;
    loopCnt=1;
    perSU=ones(1,rsize);
end

% disp([loopCnt,loopIdx,perSU(loopIdx),recDepth]);


while true
    %% same neuron
    rec_same=cell(0);
    if perSU(loopIdx)<=numel(in{loopIdx})-1 && in{loopIdx}(perSU(loopIdx)+1)-in{loopIdx}(perSU(loopIdx))<opt.burstInterval
        link=[loopIdx,perSU(loopIdx),in{loopIdx}(perSU(loopIdx))];
        delta=zeros(1,numel(perSU));
        delta(loopIdx)=1;
        rec_same=bz.rings.relax_tag_long(in,loopIdx,recDepth+1,loopCnt,perSU+delta,'burstInterval',opt.burstInterval);
        for ii=1:numel(rec_same)
            out(end+1)={[link;rec_same{ii}]};
        end
    else % no branch last ts
        if loopCnt>rsize
            out={[loopIdx,perSU(loopIdx),in{loopIdx}(perSU(loopIdx))]};
        end
    end
    %% branch
    rec_branch=cell(0);
    if perSU(loopIdx)>numel(in{loopIdx})
        break
    end
    nxtd=loopIdx+1;
    if nxtd>rsize
        nxtd=1;
    end
    perSU(nxtd)=1;
%     while perSU(nxtd)<=numel(in{nxtd}) && in{nxtd}(perSU(nxtd))<in{loopIdx}(perSU(loopIdx))+24
%         perSU(nxtd)=perSU(nxtd)+1; % TODO: improve performance with binary tree search
%     end
    perSU(nxtd)=findFirst(in{nxtd},in{loopIdx}(perSU(loopIdx))+24);
    if perSU(nxtd)<=numel(in{nxtd}) && in{nxtd}(perSU(nxtd))<in{loopIdx}(perSU(loopIdx))+300 % in window
        link=[loopIdx,perSU(loopIdx),in{loopIdx}(perSU(loopIdx))];
        rec_branch=bz.rings.relax_tag_long(in,nxtd,recDepth+1,loopCnt+1,perSU,'burstInterval',opt.burstInterval);
        for ii=1:numel(rec_branch)
            out(end+1)={[link;rec_branch{ii}]};
        end
    end

    % TODO: consider a to-be-skipped list?
    if recDepth==1
        if ~isempty(rec_branch) || ~isempty(rec_same)
            for ss=1:numel(perSU)
                perSU(ss)=max(cellfun(@(x) max(x(x(:,1)==ss,2)),[rec_same,rec_branch]),[],"all");
            end
        end
        if ~all(perSU+1<cellfun(@(x) numel(x),in),'all')
            break
        else
            nextTS=arrayfun(@(x) in{x}(perSU(x)+1),1:rsize);
            [~,loopIdx]=min(nextTS);
            perSU(loopIdx)=perSU(loopIdx)+1;
        end
    else
        break
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