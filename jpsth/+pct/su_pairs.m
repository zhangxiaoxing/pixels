classdef su_pairs
    methods (Static)
        function nonmem=get_nonmem(waveid)
            nonmem=any(waveid==0,2);
        end

        function congru=get_congru(waveid)
            congru=(waveid(:,1)==waveid(:,2) & waveid(:,1)>0)...
                | (ismember(waveid(:,1),1:2) & waveid(:,2)==5) | (ismember(waveid(:,2),1:2) & waveid(:,1)==5)...
                | (ismember(waveid(:,1),3:4) & waveid(:,2)==6) | (ismember(waveid(:,2),3:4) & waveid(:,1)==6)...
                | (ismember(waveid(:,1),[1,3]) & waveid(:,2)==7) | (ismember(waveid(:,2),[1,3]) & waveid(:,1)==7)...
                | (ismember(waveid(:,1),[2,4]) & waveid(:,2)==8) | (ismember(waveid(:,2),[2,4]) & waveid(:,1)==8);
        end

        function incong=get_incongru(waveid)
            incong=(waveid(:,1)~=waveid(:,2) & all(waveid>0,2) & all(waveid<5,2))...% 1 2 3 4
                | (waveid(:,1)~=waveid(:,2) & all(waveid>4,2) & all(waveid<7,2))...% 5 6
                | (waveid(:,1)~=waveid(:,2) & all(waveid>6,2))...% 7 8
                | (ismember(waveid(:,1),1:2) & waveid(:,2)==6) | (ismember(waveid(:,2),1:2) & waveid(:,1)==6)...
                | (ismember(waveid(:,1),3:4) & waveid(:,2)==5) | (ismember(waveid(:,2),3:4) & waveid(:,1)==5)...
                | (ismember(waveid(:,1),[1,3]) & waveid(:,2)==8) | (ismember(waveid(:,2),[1,3]) & waveid(:,1)==8)...
                | (ismember(waveid(:,1),[2,4]) & waveid(:,2)==7) | (ismember(waveid(:,2),[2,4]) & waveid(:,1)==7);
        end
    end
end