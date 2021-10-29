%assume ran bz.load_sig_pair
load(fullfile('bzdata','rings_bz.mat'),'rings');
hebbPattern=cell(0,4);
[sig,~]=bz.load_sig_pair();
if ~exist('pstats3','var') || ~exist('pstats4','var')
    load('loops_proportion_stats_all.mat','pstats3','pstats4','mstats3','mstats4');
end

for sess=1:size(rings,1)
    if any(cellfun(@(x) isempty(x),rings(sess,1:2)),'all'), continue; end
    disp(sess);
    r3all=rings{sess,1};
    r4all=rings{sess,2};
    allconn=sig.suid(sig.sess==sess,:);

    for r4idx=1:size(r4all,1)
        for r3idx=1:size(r3all,1)
            if ~nnz(ismember(r4all(r4idx,:),r3all(r3idx,:)))==2, continue; end
            r4place=find(ismember(r4all(r4idx,:),r3all(r3idx,:)));
            if ~(isequal(r4place,[1 3]) || isequal(r4place,[2 4])), continue; end
            r3place=find(ismember(r3all(r3idx,:),r4all(r4idx,:)));
            switch sum(r3place)
                case 4 %1 3 or 3 1
                    final_ring124=r3all(r3idx,1:3);
                    %                                 
                case 3 % 1 2
                    final_ring124=r3all(r3idx,[2 3 1]);
                case 5 % 2 3
                    final_ring124=r3all(r3idx,[3 1 2]);
            end
            for b4pidx=1:length(r4all)
                if nnz(ismember(r4all(b4pidx,:),final_ring124))==3 ...
                        && isequal(...
                        bz.rings.flexsort(r4all(b4pidx,ismember(r4all(b4pidx,:),final_ring124))),...
                        bz.rings.flexsort(final_ring124))
                    pre14Idx=find(r4all(b4pidx,:)==final_ring124(3));
                    if (pre14Idx==1 && ~ismember(r4all(b4pidx,4),final_ring124)) ...
                            || (pre14Idx~=1 && ~ismember(r4all(b4pidx,pre14Idx-1),final_ring124))
                        post15idx=find(r4all(r4idx,:)==final_ring124(1))+1;
                        if post15idx>4
                            post15idx=post15idx-4;
                        end
                        if any(allconn(:,1)==final_ring124(1) & allconn(:,2)==r4all(r4idx,post15idx))
                            %                                             disp(r3all(r3idx,:))
                            %                                             disp(r4all(r4idx,:))
                            %                                             disp(r4all(b4pidx,:))
                            hebbPattern(end+1,:)={sess,r3all(r3idx,:),r4all(r4idx,:),r4all(b4pidx,:)};
                            %                                             keyboard
                        end
                    end
                end
            end
        end
    end
end

save('hebb_pattern_bz.mat','hebbPattern');
return
