function out=find_rings_bz(in,msize)
arguments
    in (:,2) double {mustBeNonempty(in)} %sig or pair?
    msize (1,1) double {mustBeMember(msize,3:5)} % size of ring
end
out=[];
pre_unit_set=unique(in(:,1));
for a_first=pre_unit_set(:)'
    seconds=in(in(:,1)==a_first,2);
    for a_second=seconds(:)'
        thirds=in(in(:,1)==a_second,2);
        thirds(thirds==a_first)=[];
        if msize==3
            ring=in(ismember(in(:,1),thirds) & in(:,2)==a_first,1);
            if ~isempty(ring)
                subgrp=[ring,ring,ring];
                subgrp(:,1:2)=repmat([a_first a_second],size(subgrp,1),1);
                out=[out;subgrp];
            end
        else
            for a_third=thirds(:)'
                fourths=in(in(:,1)==a_third,2);
                fourths(fourths==a_second)=[];
                if msize==4
                    ring=in(ismember(in(:,1),fourths) & in(:,2)==a_first,1);
                    if ~isempty(ring)
                        subgrp=[ring,ring,ring,ring];
                        subgrp(:,1:3)=repmat([a_first a_second a_third],size(subgrp,1),1);
                        out=[out;subgrp];
                    end
                elseif msize==5
                    fourths(fourths==a_first)=[];
                    for a_fourth=fourths(:)'
                        fifths=in(in(:,1)==a_fourth, 2);
                        fifths(ismember(fifths,[a_first a_second]))=[];
                        ring=in(ismember(in(:,1),fifths) & in(:,2)==a_first,1);
                        if ~isempty(ring)
                            subgrp=[ring,ring,ring,ring,ring];
                            subgrp(:,1:4)=repmat([a_first a_second a_third a_fourth],size(subgrp,1),1);
                            out=[out;subgrp];
                        end
                    end
                end
            end
        end
    end
end