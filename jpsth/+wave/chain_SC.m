for dur=reshape(fieldnames(out),1,[])
    for wv=reshape(fieldnames(out.(dur{1})),1,[])
        for lp=reshape(fieldnames(out.(dur{1}).(wv{1})),1,[])
            if size(out.(dur{1}).(wv{1}).(lp{1}).ts,2)<8
                continue
            end

            ccg_qual=out.(dur{1}).(wv{1}).(lp{1}).ccg_qual;
            % 1:Polarity 2:Time of peak 3:Noise peaks 4:FWHM 5:rising
            % edge 6:falling edge
            strict_sel=ccg_qual(:,2)>=252 & ccg_qual(:,4)>=2 & ccg_qual(:,4)<=40 & ccg_qual(:,5)>248;
            if all(strict_sel) % plot
                keyboard()
            else
                keyboard()
            end
        end
    end
end


