function STF_SC(opt)
arguments
    opt.gen_STF (1,1) logical = false;
    opt.plot_STF (1,1) logical = true;
    opt.gen_COA (1,1) logical = false;
    opt.plot_COA (1,1) logical = false;
    opt.plotskip (1,1) double {mustBeInteger,mustBeNonnegative}=0
    opt.sessid (1,1) double = 18;
end
[spkID,spkTS,~,~,~]=ephys.getSPKID_TS(opt.sessid,'criteria','WT');
if opt.gen_STF || opt.gen_COA
    [correct_stats,correct_type]=bz.hist.get_stp_stats(6000,'BZWT','trialtype','correct');
    %% criteria in line with +bz\+hist\get_stats_by_mem_type.m
    sess_sel=correct_stats.correct.sess==opt.sessid & ~correct_stats.correct.skip & correct_stats.correct.rsq>0 & all(correct_stats.correct.mem_type>0,2);
    %%
    suids=correct_stats.correct.sess_suids(sess_sel,:);
    memtype=correct_stats.correct.mem_type(sess_sel,:);
    % fliplr(cmp_data.nonmem(:,2:end))
    %     coeff=fliplr(correct_stats.correct.coeff(sess_sel,2:end));
    
    %% traverse congruent A->C and incongruent B->C triplets
    sums_STF=cell(0);
    sums_COA=cell(0);
    processed=[];
    for fci=1:nnz(sess_sel)
        post_uid=suids(fci,2);
        if ismember(memtype(fci,2),[1,2])
            cmem=1:2;imem=3:4;
        else
            cmem=3:4;imem=1:2;
        end
        
        if ismember(post_uid,processed)
            continue
        end
        processed=[processed,post_uid];
        pre_sel=suids(:,2)==post_uid;
        if opt.gen_STF
            pre_c=find(ismember(memtype(:,1),cmem) & pre_sel);
            pre_i=find(ismember(memtype(:,1),imem) & pre_sel);
            if isempty(pre_c) || isempty(pre_i)
                continue
            end
            for cidx=1:numel(pre_c)
                for iidx=1:numel(pre_i)
                    cid1=suids(pre_c(cidx));
                    cid2=suids(pre_i(iidx));
                    cts1=spkTS(spkID==cid1);
                    cts2=spkTS(spkID==cid2);
                    pts=spkTS(spkID==post_uid);
                    disp([size(sums_STF,1),opt.sessid,cid1,cid2,post_uid]);
                    ts=findOneSTF(cts1,cts2,pts);
                    if ~isempty(ts)
                        sums_STF=[sums_STF;{[opt.sessid,cid1,cid2,post_uid],ts}];
                        save('STF_SC.mat','sums_STF');
                    end
                end
            end
        end
        
        if opt.gen_COA
            pre_c=find(ismember(memtype(:,1),cmem) & pre_sel);
            if numel(pre_c)<2
                continue
            end
            ccomb=nchoosek(pre_c,2);
            for cidx=1:size(ccomb,1)
                cid1=suids(ccomb(cidx,1));
                cts1=spkTS(spkID==cid1);
                cid2=suids(ccomb(cidx,2));
                cts2=spkTS(spkID==cid2);
                pts=spkTS(spkID==post_uid);
                disp([size(sums_COA,1),opt.sessid,cid1,cid2,post_uid]);
                ts=findOneCOA(cts1,cts2,pts);
                if ~isempty(ts)
                    sums_COA=[sums_COA;{[opt.sessid,cid1,cid2,post_uid],ts}];
                    save('COA_SC.mat','sums_COA');
                end
            end
        end
    end
end
%% plot all
if opt.plot_STF
    load('STF_SC.mat','sums_STF');
    for sidx=1:size(sums_STF,1)
        cid1=sums_STF{sidx,1}(2);
        cid2=sums_STF{sidx,1}(3);
        post_uid=sums_STF{sidx,1}(4);
        cts1=spkTS(spkID==cid1);
        cts2=spkTS(spkID==cid2);
        pts=spkTS(spkID==post_uid);
        tts=sums_STF{sidx,2};
        for ti=1:numel(tts)
            plotOneSTF(cts1,cts2,pts,opt.sessid,cid1,cid2,post_uid,tts(ti))
            keyboard()
        end
    end
end

if opt.plot_COA
    load('COA_SC.mat','sums_COA');
    cidx=0;
    for sidx=1:size(sums_COA,1)
        cid1=sums_COA{sidx,1}(2);
        cid2=sums_COA{sidx,1}(3);
        post_uid=sums_COA{sidx,1}(4);
        cts1=spkTS(spkID==cid1);
        cts2=spkTS(spkID==cid2);
        pts=spkTS(spkID==post_uid);
        tts=sums_COA{sidx,2};
        for ti=1:numel(tts)
            cidx=cidx+1;
            if opt.plotskip>0 && cidx<opt.plotskip
                continue
            end
            plotOneCOA(cts1,cts2,pts,opt.sessid,cid1,cid2,post_uid,tts(ti))
            keyboard()
        end
    end
end

end


function ts=findOneSTF(cts,its,pts)
sps=30000;
ts=[];
for t=0:0.25:(max([cts;its;pts])/sps-3)
    cc=cts(cts>=(t+0.5).*sps & cts<(t+3).*sps);
    ii=its(its>=(t+0.5).*sps & its<(t+3).*sps);
    pp=pts(pts>=(t+0.5).*sps & pts<(t+3).*sps);
    
    if nnz(cc>(t+1)*sps)>0 && nnz(cc<(t+1)*sps & cc>(t+0.5).*sps)==0 && nnz(cc>(t+1.5)*sps)==0 ...
            && nnz(ii>(t+2)*sps)>0 && nnz(ii<(t+2)*sps & ii>(t+0.5).*sps)==0 && nnz(ii>(t+2.5)*sps)==0
        tc=cc(1);
        ti=ii(1);
        avgbase=nnz(pp<tc)./(tc-(t+0.5).*sps);
        avgcongru=nnz(pp>tc & pp<=ti)./(ti-tc);
        avgincon=nnz(pp>ti)./((t+3).*sps-ti);
        
        if avgcongru>avgbase && avgincon<avgbase
            ts=[ts;t];
        end
    end
end
end

function ts=findOneCOA(cts1,cts2,pts)
sps=30000;
ts=[];
for t=0:0.02:(max([cts1;cts2;pts])/sps-0.5)
    cc1=cts1(cts1>=(t+0.1).*sps & cts1<(t+0.3).*sps);
    cc2=cts2(cts2>=(t+0.1).*sps & cts2<(t+0.3).*sps);
    pp=pts(pts>=(t+0.1).*sps & pts<(t+0.3).*sps);
    ppbase=nnz(pts>=t.*sps & pts<(t+0.1).*sps);
    if numel(cc1)>=10 || numel(cc2)>=10 || numel(pp)>=10 || numel(cc1)<2 || numel(cc2)<2 || numel(pp)<2 || ppbase > 1
        continue
    end
    comat=cc1-cc2.';
    [co1,co2]=find(comat<300 & comat>-300);
    if ~isempty(co1) && ~isempty(co2)
        ppa1=cc1(co1)-pp.';
        ppa2=cc2(co2)-pp.';
        if nnz(ppa1<0 & ppa1>-300 & ppa2<0 &ppa2>-300)>1
            ts=[ts;t];
        end
    end
end
end


function fh=plotOneSTF(cts,its,pts,sessid,cid,iid,pid,t)
sps=30000;
cc=cts(cts>=(t+0.5).*sps & cts<(t+3).*sps);
ii=its(its>=(t+0.5).*sps & its<(t+3).*sps);
pp=pts(pts>=(t+0.5).*sps & pts<(t+3).*sps);
fh=figure('Color','w','Position',[100,100,360,200]);
hold on;
plot([cc,cc].',repmat([2.1;2.5],1,numel(cc)),'-r');
plot([ii,ii].',repmat([1.6;2.0],1,numel(ii)),'-b');
plot([pp,pp].',repmat([1.1;1.5],1,numel(pp)),'-k');
ppext=pts(pts>=(t-1).*sps & pts<(t+2.5).*sps);
[pphc,ppedge]=histcounts(pp,((t+0.5):0.025:(t+3)).*sps);
pmvmm=movmean(pphc,10).*4;
pmvmm=pmvmm-min(pmvmm);
pmvmm=pmvmm./max(pmvmm);

cchc=histcounts(cc,((t+0.5):0.025:(t+3)).*sps);
cmvmm=movmean(cchc,10).*4;
cmvmm=cmvmm-min(cmvmm);
cmvmm=cmvmm./max(cmvmm);

iihc=histcounts(ii,((t+0.5):0.025:(t+3)).*sps);
imvmm=movmean(iihc,10).*4;
imvmm=imvmm-min(imvmm);
imvmm=imvmm./max(imvmm);

plot(ppedge(1:end-1)+0.0125*sps,cmvmm,'-r');
plot(ppedge(1:end-1)+0.0125*sps,imvmm,'-b');
plot(ppedge(1:end-1)+0.0125*sps,pmvmm,'-k');

xlim(([t+0.5,t+3]).*sps);
set(gca(),'YTick',1.3:0.5:2.3,'YTickLabel',{'Post','Diff.-Pre','Same-Pre'},'XTick',(t+0.5:t+2.5).*sps,'XTickLabel',0:2);
xlabel('Time (s)')
title(sprintf('%d-%d,%d,%d,%.4f',sessid,cid,iid,pid,t));
end


function fh=plotOneCOA(cts,its,pts,sessid,cid,iid,pid,t)
sps=30000;
cc=cts(cts>=t.*sps & cts<(t+0.3).*sps);
ii=its(its>=t.*sps & its<(t+0.3).*sps);
pp=pts(pts>=t.*sps & pts<(t+0.3).*sps);
fh=figure('Color','w','Position',[100,100,360,200]);
hold on;
plot([cc,cc].',repmat([2.1;2.5],1,numel(cc)),'-r');
plot([ii,ii].',repmat([1.6;2.0],1,numel(ii)),'-m');
plot([pp,pp].',repmat([1.1;1.5],1,numel(pp)),'-k');
% ppext=pts(pts>=t.*sps & pts<(t+0.3).*sps);
% [pphc,ppedge]=histcounts(pp,(t:0.005:(t+0.3)).*sps);
% pmvmm=movmean(pphc,10).*4;
% pmvmm=pmvmm-min(pmvmm);
% pmvmm=pmvmm./max(pmvmm);
% 
% cchc=histcounts(cc,(t:0.005:(t+0.3)).*sps);
% cmvmm=movmean(cchc,10).*4;
% cmvmm=cmvmm-min(cmvmm);
% cmvmm=cmvmm./max(cmvmm);
% 
% iihc=histcounts(ii,(t:0.005:(t+0.3)).*sps);
% imvmm=movmean(iihc,10).*4;
% imvmm=imvmm-min(imvmm);
% imvmm=imvmm./max(imvmm);
% 
% plot(ppedge(1:end-1)+0.0025*sps,smooth(cmvmm,5),'-r');
% plot(ppedge(1:end-1)+0.0025*sps,smooth(imvmm,5),'-m');
% plot(ppedge(1:end-1)+0.0025*sps,smooth(pmvmm,5),'-k');

xlim([t,t+0.3].*sps);
set(gca(),'YTick',1.3:0.5:2.3,'YTickLabel',{'Post','Same-Pre A','Same-Pre B'},'XTick',(t:0.1:t+0.3).*sps,'XTickLabel',0:0.1:0.3);
xlabel('Time (s)')
title(sprintf('%d-%d,%d,%d,%.4f',sessid,cid,iid,pid,t));
end

