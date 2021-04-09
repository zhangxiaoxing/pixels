function out=count_motif_congru(in,reg,pref,bin,msize,sample)
out=[];
pre_unit_set=unique(in(:,1));
for i=pre_unit_set(:)'
    mono_post=in(in(:,1)==i & reg(:,1)~=reg(:,2) & pref(:,bin)==sample & pref(:,bin+6)==sample & reg(:,1)<116 & reg(:,2)<116,2);
    for j=mono_post(:)'
        sec_post=in(in(:,1)==j & reg(:,1)~=reg(:,2) & pref(:,bin)==sample & pref(:,bin+6)==sample & reg(:,1)<116 & reg(:,2)<116,2);
        if msize==3
            ring=in(ismember(in(:,1),sec_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & pref(:,bin)==sample & pref(:,bin+6)==sample & reg(:,1)<116 & reg(:,2)<116,1);
            if ~isempty(ring)
                subgrp=[ring,ring,ring];
                subgrp(:,1:2)=repmat([i j],size(subgrp,1),1);
                out=[out;subgrp];
            end
        elseif msize==4
            for k=sec_post(:)'
                third_post=in(in(:,1)==k & reg(:,1)~=reg(:,2) & pref(:,bin)==sample & pref(:,bin+6)==sample & reg(:,1)<116 & reg(:,2)<116,2);
                ring=in(ismember(in(:,1),third_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & pref(:,bin)==sample & pref(:,bin+6)==sample & reg(:,1)<116 & reg(:,2)<116,1);
                if ~isempty(ring)
                    subgrp=[ring,ring,ring,ring];
                    subgrp(:,1:3)=repmat([i j k],size(subgrp,1),1);
                    out=[out;subgrp];
                end
            end
        elseif msize==5
            for k=sec_post(:)'
                third_post=in(in(:,1)==k & reg(:,1)~=reg(:,2) & pref(:,bin)==sample & pref(:,bin+6)==sample & reg(:,1)<116 & reg(:,2)<116,2);
                third_post(third_post==i)=[];
                for l=third_post(:)'
                    fourth_post=in(in(:,1)==l & reg(:,1)~=reg(:,2) & pref(:,bin)==sample & pref(:,bin+6)==sample & reg(:,1)<116 & reg(:,2)<116,2);
                    fourth_post(ismember(fourth_post,[i j]))=[];
                    ring=in(ismember(in(:,1),fourth_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & pref(:,bin)==sample & pref(:,bin+6)==sample & reg(:,1)<116 & reg(:,2)<116,1);
                    if ~isempty(ring)
                        subgrp=[ring,ring,ring,ring,ring];
                        subgrp(:,1:4)=repmat([i j k l],size(subgrp,1),1);
                        out=[out;subgrp];
                    end
                end
            end
        end
    end
end


end