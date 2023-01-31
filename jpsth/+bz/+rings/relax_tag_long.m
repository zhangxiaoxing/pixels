%% Repurosed for burst-loop algorithm as of 23.01.30

function [out,perSU]=relax_tag_long(in,loopIdx,recDepth,loopCnt,perSU,opt)
arguments
    in
    loopIdx % current neuron idx
    recDepth
    loopCnt % accumated loop neuron count
    perSU
    opt.burstInterval (1,1) double = 300 % ticks, at 30k sps
end

out=cell(0); %{[depth,id]}
rsize=size(in,2);

if isempty(loopIdx)
    [~,loopIdx]=min(cellfun(@(x) x(1),in));
    recDepth=1;
    loopCnt=1;
    perSU=ones(1,rsize);
end

% disp([loopCnt,loopIdx,perSU(loopIdx),recDepth]);



while true
    %% same neuron
    if perSU(loopIdx)<=numel(in{loopIdx})-1 && in{loopIdx}(perSU(loopIdx)+1)-in{loopIdx}(perSU(loopIdx))<opt.burstInterval
        link=[loopIdx,perSU(loopIdx),in{loopIdx}(perSU(loopIdx))];
        perSU(loopIdx)=perSU(loopIdx)+1;
        [rec,perSU]=bz.rings.relax_tag_long(in,loopIdx,recDepth+1,loopCnt,perSU,'burstInterval',opt.burstInterval);
        for ii=1:numel(rec)
            out(end+1)={[link;rec{ii}]};
        end
    else % no branch last ts
        if loopCnt>rsize
            out={[loopIdx,perSU(loopIdx),in{loopIdx}(perSU(loopIdx))]};
        end
    end
    %% branch
    %     if loopDepth<rsize
    % TODO varify loop index WIP
    nxtd=loopIdx+1;
    if nxtd>rsize
        nxtd=1;
    end
    perSU(nxtd)=1;
    while perSU(nxtd)<=numel(in{nxtd}) && in{nxtd}(perSU(nxtd))<in{loopIdx}(perSU(loopIdx))+24
        perSU(nxtd)=perSU(nxtd)+1; % TODO: improve performance with binary tree search
    end
    if perSU(nxtd)<=numel(in{nxtd}) && in{nxtd}(perSU(nxtd))<in{loopIdx}(perSU(loopIdx))+300 % in window
        link=[loopIdx,perSU(loopIdx),in{loopIdx}(perSU(loopIdx))];
        [rec,perSU]=bz.rings.relax_tag_long(in,nxtd,recDepth+1,loopCnt+1,perSU,'burstInterval',opt.burstInterval);
        for ii=1:numel(rec)
            out(end+1)={[link;rec{ii}]};
        end
    end

    % TODO: consider a to-be-skipped list?
    if recDepth==1
        if all(perSU<cellfun(@(x) numel(x),in),'all')
            nextTS=arrayfun(@(x) in{x}(perSU(x)+1),1:rsize);
            [~,loopIdx]=min(nextTS);
            perSU(loopIdx)=perSU(loopIdx)+1;
        else
            break
        end
    else
        break
    end
end
end