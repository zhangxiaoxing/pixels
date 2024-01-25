function silly_sync(dpath)
display(dpath)
cd(dpath)
if ispc
    addpath(fullfile("k:","Lib","npy-matlab","npy-matlab"))
elseif isunix
    addpath('/home/zhangxx/code/npy-matlab/npy-matlab')
end
sr=readNPY('sync_raw.npy');

% 158us=4.74sp or 153us=4.59sp per signal bit, 9.18 or 9.33 sp for
% consecutive bits

delta=diff([0;bitand(sr,64);0]);
riseedge=find(delta>0); % sss
falledge=find(delta<0); % sss
edge_edge_dist=diff(riseedge);
word_hist=histcounts(edge_edge_dist,0.5:60.5);
[~,wordw]=max(word_hist);
highlevelw=falledge-riseedge; % sss

bith=histcounts(highlevelw,0.5:10.5);
[~,bitw]=max(bith);

ii=find(sr(1:100)==0,1,"first");
jj=find(sr(ii+1:ii+100)>0,1,"first");
ii=ii+jj;

evts=java.util.LinkedList;
% multievts=[];
previous=0;



while true
    if ii-previous>1000000
        disp(ii)
        previous=ii;
        % writeNPY(evts,'sync_events.npy');
    end
    % TODO: break if no more data
    if (ii+wordw+1)>size(sr,1)
        break
    end
    nxtrise=find(sr(ii:ii+wordw+1)==0,1,"last");
    if nxtrise>wordw-2 && nxtrise<wordw+2
        blocksample=sr(ii:ii+nxtrise-1);
        hicnt=nnz(blocksample>0);
        if ~(hicnt>bitw/2 && hicnt<bitw*1.5) % skip null events to improve performance 
            % currwordw=size(blocksample,1);
            evts.add([ii,...
                median(blocksample(7:9))>0,...
                median(blocksample(10:12))>0,...
                median(blocksample(14:16))>0,...
                median(blocksample(18:21))>0,...
                median(blocksample(23:25))>0]);
            % if sum(evts(end,2:end))>1
            %    multievts=[multievts;transpose(resize(blocksample,30))];
            %    if rem(size(multievts,1),1000)==0
            %        keyboard()
            %    end
            % end
        end
        ii=ii+nxtrise;
    else
        disp("non empty stop bit")
        ii=ii+nxtrise;
    end
end
evta=transpose(cell2mat(transpose(cell(evts.toArray()))));
writeNPY(evta,'sync_events.npy');
end

