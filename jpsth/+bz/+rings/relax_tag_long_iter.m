%% WIP algorithm for replacing recursive relax_tag_long function with 
%% iterations due to performance concerns. 23.01.30


% in_={[10:450:920,2150:450:2700].',[950:450:1500].',[1550:450:2200].'}
% bz.rings.relax_tag_long_iter(in_);
function out=relax_tag_long_iter(in,opt)
arguments
    in
    opt.burstInterval (1,1) double = 600 % ticks, at 30k sps
    opt.burst (1,1) logical = true % option to merge single spike loop algorithm
end
out=cell(0); %{[depth,id]}
rsize=size(in,2);

[~,loopIdx]=min(cellfun(@(x) x(1),in));
recDepth=1;
loopCnt=1;
perSU=ones(1,rsize);

% disp([loopCnt,loopIdx,perSU(loopIdx),recDepth]);

%% iterative
istack={{loopIdx,recDepth,loopCnt,perSU,[]}};

% TODO: where is result?
while ~isempty(istack)
    % get parameters from stack
    currFrame= istack{1};
    loopIdx=currFrame{1};
    recDepth=currFrame{2};
    loopCnt=currFrame{3};
    perSU=currFrame{4};
    prevlinks=currFrame{5};
    if numel(istack)>1
        istack=istack(2:end);
    else
        istack=cell(0);
    end

    if perSU(loopIdx)<=numel(in{loopIdx})-1 && in{loopIdx}(perSU(loopIdx)+1)-in{loopIdx}(perSU(loopIdx))<opt.burstInterval
        link=[loopIdx,perSU(loopIdx),in{loopIdx}(perSU(loopIdx))];
        perSU(loopIdx)=perSU(loopIdx)+1;
        istack=[{{loopIdx,recDepth+1,loopCnt,perSU,[prevlinks;link]}};...
            istack];
        continue; % TODO varify jump point.
    else
        if loopCnt>rsize % proper dealing with end-of-line-return
            % TODO unfinished update of result
            % TODO check for repetition before udpate result
            out(end+1)=[prevlinks;loopIdx,perSU(loopIdx),in{loopIdx}(perSU(loopIdx))];
        end

        % push branch context into stack. i.e. depth first
        %                 rec_branch=cell(0);
        if perSU(loopIdx)>numel(in{loopIdx})
            break % TODO break from control loop or continue?
        end

        if loopIdx==rsize
            nxtd=1;
        else
            nxtd=loopIdx+1;
        end
        perSU(nxtd)=findFirst(in{nxtd},in{loopIdx}(perSU(loopIdx))+24);

        if perSU(nxtd)<=numel(in{nxtd}) && in{nxtd}(perSU(nxtd))<in{loopIdx}(perSU(loopIdx))+300 % in window
            link=[loopIdx,perSU(loopIdx),in{loopIdx}(perSU(loopIdx))];
            istack=[{{nxtd,recDepth+1,loopCnt+1,perSU,[prevlinks;link]}};...
                istack];
            continue;
        end
        % Unencumbered by any FC
        if recDepth==1
            % TODO: Skip processed spikes
            if ~isempty(rec_branch) || ~isempty(rec_same)
                for ss=1:numel(perSU)
                    perSU(ss)=max(cellfun(@(x) max(x(x(:,1)==ss,2)),[rec_same,rec_branch]),[],"all");
                end
            end
            if ~all(perSU+1<cellfun(@(x) numel(x),in),'all')
                break % valid break of control loop
            else
                nextTS=arrayfun(@(x) in{x}(perSU(x)+1),1:rsize);
                [~,loopIdx]=min(nextTS);
                perSU(loopIdx)=perSU(loopIdx)+1;

                % TODO should clear stack here?
                istack={{loopIdx,1,1,perSU,[]}};...
                    continue
            end
        else
            % Should modify stack
            continue
        end
    end
end
end


%             if recDepth==1
%                 % TODO: Skip processed spikes
%                 if ~isempty(rec_branch) || ~isempty(rec_same)
%                     for ss=1:numel(perSU)
%                         perSU(ss)=max(cellfun(@(x) max(x(x(:,1)==ss,2)),[rec_same,rec_branch]),[],"all");
%                     end
%                 end
%                 if ~all(perSU+1<cellfun(@(x) numel(x),in),'all')
%                     break % valid break of control loop
%                 else
%                     nextTS=arrayfun(@(x) in{x}(perSU(x)+1),1:rsize);
%                     [~,loopIdx]=min(nextTS);
%                     perSU(loopIdx)=perSU(loopIdx)+1;
% 
%                     % TODO should clear stack here?
%                     istack={{loopIdx,1,1,perSU,[]}};...
%                     continue
% 
%                 end
%             else
%                 break
%             end




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