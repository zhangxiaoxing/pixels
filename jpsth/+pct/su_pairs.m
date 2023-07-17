classdef su_pairs
    methods (Static)
        function nonmem=get_nonmem(waveid)
            nonmem=any(waveid==0,2);
        end

        function congru=get_congru(waveid,opt)
            arguments
                waveid (:,2) double
                opt.asym_congru (1,1) logical = false
                opt.strict_mix (1,1) logical = false
                opt.odor_only (1,1) logical = false
            end

            congrumix=(waveid(:,1)==waveid(:,2) & waveid(:,1)>0);
            if opt.strict_mix
                congru=congrumix;
            elseif opt.asym_congru
                congru_single_mix=(ismember(waveid(:,1),1:2) & waveid(:,2)==5)...
                    | (ismember(waveid(:,1),3:4) & waveid(:,2)==6)...
                    | (ismember(waveid(:,1),[1,3]) & waveid(:,2)==7)...
                    | (ismember(waveid(:,1),[2,4]) & waveid(:,2)==8);
                congru=congrumix | congru_single_mix;
            else
                if opt.odor_only
                    congru=((waveid(:,1)==waveid(:,2) & ismember(waveid(:,1),1:6)))...
                        | (ismember(waveid(:,1),1:2) & waveid(:,2)==5) | (ismember(waveid(:,2),1:2) & waveid(:,1)==5)...
                        | (ismember(waveid(:,1),3:4) & waveid(:,2)==6) | (ismember(waveid(:,2),3:4) & waveid(:,1)==6);
                else
                    congru_single_mix=(ismember(waveid(:,1),1:2) & waveid(:,2)==5) | (ismember(waveid(:,2),1:2) & waveid(:,1)==5)...
                        | (ismember(waveid(:,1),3:4) & waveid(:,2)==6) | (ismember(waveid(:,2),3:4) & waveid(:,1)==6)...
                        | (ismember(waveid(:,1),[1,3]) & waveid(:,2)==7) | (ismember(waveid(:,2),[1,3]) & waveid(:,1)==7)...
                        | (ismember(waveid(:,1),[2,4]) & waveid(:,2)==8) | (ismember(waveid(:,2),[2,4]) & waveid(:,1)==8);
                    congru=congrumix | congru_single_mix;
                end
            end

        end

        function incong=get_incongru(waveid,opt)
            arguments
                waveid
                opt.odor_only (1,1) logical = false
            end
            if opt.odor_only
                incong=(any(ismember(waveid,1:2),2) & any(ismember(waveid,3:4),2))...
                    | (waveid(:,1)~=waveid(:,2) & all(waveid>4,2) & all(waveid<7,2))...% 5 6
                    | (ismember(waveid(:,1),1:2) & waveid(:,2)==6) | (ismember(waveid(:,2),1:2) & waveid(:,1)==6)...
                    | (ismember(waveid(:,1),3:4) & waveid(:,2)==5) | (ismember(waveid(:,2),3:4) & waveid(:,1)==5);
            else
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
end