function [out,spk_peak_ai]=get_sig(xcorr,xshuf,opt)
    arguments
        xcorr (1,1) struct
        xshuf (1,1) struct
        opt.ntrial (1,1) double  {mustBePositive,mustBeInteger} = 20
        opt.nspk_thres (1,1) double {mustBePositive,mustBeInteger} = 500
        opt.thresh (1,1) double {mustBePositive,mustBeReal} = 2.8070
        opt.AI_thresh (1,1) double = 0.4
    end
    out=[];
    spk_peak_ai=[];
    for si=1:(size(xcorr.xcorr,1)-1)
        su1id=str2double(xcorr.label{si});
        for sj=(si+1):size(xcorr.xcorr,2)
            su2id=str2double(xcorr.label{sj});
            totalCount=nansum(squeeze(xcorr.xcorr(si,sj,:)));
            if numel(xcorr.cfg.trials)<opt.ntrial || totalCount<opt.nspk_thres
                continue
            end
            
            hist=squeeze(xcorr.xcorr(si,sj,:));
            shuf=squeeze(xshuf.shiftpredictor(si,sj,:));
            diff=hist-smooth(shuf);
            %diffs1(50:51)=0;
            stdv=std(shuf);
            score=diff(46:55)./stdv;
            %any score > thresh
            if any(score([1:4,6:10])>opt.thresh)
                peak=nanmax(score([1:4,6:10]));
                bincount=diff(46:55);
                bincount(score<=opt.thresh)=0;
                % A peak at a negative lag (I.E. AI>0) for stat.xcorr(chan1,chan2,:) means that chan1 is leading
                % chan2.
                sumdiff=(sum(bincount(1:4))-sum(bincount(7:10)));
                if sumdiff~=0
                    AsymIdx=sumdiff/(sum(bincount(1:4))+sum(bincount(7:end)));
                    if abs(AsymIdx)>=opt.AI_thresh
                        out=[out;su1id,su2id];
                        spk_peak_ai=[spk_peak_ai;totalCount,peak,AsymIdx];
                    end
                end
            end
        end
    end
end