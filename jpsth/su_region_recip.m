update_motif=false;
rpt=1000;

for bin=1:6
    fbin(bin)=load(sprintf('0831_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
end
fbase=load('0906_conn_chain_duo_6s_-2_-1.mat');
%% non-memory neurons temporarily skipped
%for bin=1:6
%    fstr=load(sprintf('0813_nonsel_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
%    fnbin(bin)=fstr;
%end
%    for bin=5
%        if ~isequal(fbin(bin).pair_chain,fbin(bin+1).pair_chain)
%            keyboard
%        end
%    end
load('reg_keep.mat');
% load('io_sel.mat');
greymatter=find(cellfun(@(x) ~isempty(regexp(x,'[A-Z]','match','once')), reg_set));
% any(cell2mat(cellfun(@(x) x(:,1)<100,ioselstats,'UniformOutput',false)))
% nansel=io_entire_delay(:,1)<100;
% sel=find(~nansel & greymatter);
su_set=unique(fbin(1).conn_chain_S1(:,1));
input_to_count=nan(length(reg_set),6);
output_from_count=nan(length(reg_set),6);
recip_count=nan(length(reg_set),6);
motif_count_actinact=nan(3,6);
motif_count_incong=nan(3,6);
motif_baseline_count=nan(3,1);
motif_baseline_count_incong=nan(3,1);
%triplet_nonsel_count=nan(1,6);

input_to_count_shuf=nan(length(reg_set),6,rpt);
output_from_count_shuf=nan(length(reg_set),6,rpt);
recip_count_shuf=nan(length(reg_set),6,rpt);
motif_count_shuf=nan(3,6,rpt);
motif_count_shuf_incong=nan(3,6,rpt);
motif_baseline_count_shuf_incong=nan(1,rpt);
motif_baseline_count_shuf=nan(3,rpt);
%% baseline
if exist('baseline_stats','var') && baseline_stats
    for midx=1:3
        disp(midx)
        msize=midx+2;
        motif_baseline_count(midx)=count_motif_congru_actinact(fbase.conn_chain_S1,fbase.reg_chain_S1,fbase.pref_chain_S1,bin,msize);
        %motif_baseline_count_incong=count_triplet_incong_actinact(fbase.conn_chain_S1,fbase.reg_chain_S1,fbase.pref_chain_S1,bin);
        parfor i=1:rpt
            [shuf,shufreg]=shuffle_conn_chain(fbase.conn_chain_S1,fbase.pair_chain,fbase.pair_reg);
            motif_baseline_count_shuf(midx,i)=count_motif_congru_actinact(shuf,shufreg,fbase.pref_chain_S1,bin,msize);
            %triplet_baseline_count_shuf_incong(i)=count_triplet_incong_actinact(shuf,shufreg,fbase.pref_chain_S1,bin);
        end
    end
    keyboard
    
end
%% candidate pairs
if exist('candidate_stats','var') && candidate_stats
    fut=parallel.FevalFuture.empty;
    for caIdx=1:4096 %7bit session, 3bit bin, 2bit msize
        if exist(sprintf('candidate_count_%d.mat',caIdx),'file')
            continue
        end
        session=bitshift(bitand(caIdx,bin2dec('111111100000')),-5);
        bin=bitshift(bitand(caIdx,bin2dec('000000011100')),-2);
        msize=bitand(caIdx,bin2dec('000000000011'))+2;
        if (~ismember(msize,3:5)) || (~ismember(bin,1:6))
            continue
        end
        mirror_chain=@(x) [x;x(:,[2 1])];
        cc=mirror_chain(fbin(bin).pair_chain);
        sel=cc(:,1)>session*100000 & cc(:,1)<(session+1)*100000;
        if nnz(sel)==0
            continue
        end
        rr=mirror_chain(fbin(bin).pair_reg);
        pp=[fbin(bin).pref_pair;fbin(bin).pref_pair(:,[7:12,1:6])];
        fut(end+1)=parfeval(@count_motif_congru,2,cc(sel,:),rr(sel,:),pp(sel,:),bin,msize,caIdx);
    end
    for futIdx=1:numel(fut)
        [cmpltIdx,motif_candidate_count,caIdx]=fetchNext(fut);
        savefun(sprintf('candidate_count_%d.mat',caIdx),session,bin,msize,motif_candidate_count);
        fprintf('%d, %d of %d\n',caIdx,futIdx,numel(fut));
    end
    keyboard
end
%% candidate baseline pairs
if exist('candidate_base_stats','var') && candidate_base_stats
    fut=parallel.FevalFuture.empty;
    bin=-2;
    for caIdx=1:512 %7bit session, 3bit bin, 2bit msize
        if exist(sprintf('candidate_base_count_%d.mat',caIdx),'file')
            continue
        end
        session=bitshift(bitand(caIdx,bin2dec('111111100')),-2);
        msize=bitand(caIdx,bin2dec('000000011'))+2;
        if ~ismember(msize,3:5)
            continue
        end
        mirror_chain=@(x) [x;x(:,[2 1])];
        cc=mirror_chain(fbase.pair_chain);
        sel=cc(:,1)>session*100000 & cc(:,1)<(session+1)*100000;
        if nnz(sel)==0
            continue
        end
        rr=mirror_chain(fbase.pair_reg);
        pp=[fbase.pref_pair;fbase.pref_pair(:,[7:12,1:6])];
        %% future
        fut(end+1)=parfeval(@count_motif_congru_actinact,2,cc(sel,:),rr(sel,:),pp(sel,:),bin,msize,caIdx);

    end
    %% get result
    for futIdx=1:numel(fut)
        [cmpltIdx,motif_candidate_count,caIdx]=fetchNext(fut);
        savefun(sprintf('candidate_base_count_%d.mat',caIdx),session,bin,msize,motif_candidate_count);
        fprintf('%d, %d of %d\n',caIdx,futIdx,numel(fut));
    end
    keyboard
end
%% memory delay
if true
    disp('memory delay')
    keyboard
    motif_count=nan(3,6);
    parfor bin=1:6
        disp(bin);
        for midx=1:3
            msize=midx+2;
            motif_count(midx,bin)=count_motif_congru(fbin(bin).conn_chain_S1,fbin(bin).reg_chain_S1,fbin(bin).pref_chain_S1,bin,msize);
            motif_count_incong(midx,bin)=count_motif_incong(fbin(bin).conn_chain_S1,fbin(bin).reg_chain_S1,fbin(bin).pref_chain_S1,bin,msize);
        end
        %triplet_count_actinact(bin)=count_triplet_congru_actinact(fbin(bin).conn_chain_S1,fbin(bin).reg_chain_S1,fbin(bin).pref_chain_S1,bin);
        %triplet_count_incong(bin)=count_triplet_incong_actinact(fbin(bin).conn_chain_S1,fbin(bin).reg_chain_S1,fbin(bin).pref_chain_S1,bin);
        %quadruplet_count_incong(bin)=count_triplet_incong_actinact(fbin(bin).conn_chain_S1,fbin(bin).reg_chain_S1,fbin(bin).pref_chain_S1,bin);
    end
    parfor i=1:rpt
        for midx=1:3
            msize=midx+2
            for bin=1:6
                [shuf,shufreg]=shuffle_conn_chain(fbin(bin).conn_chain_S1,fbin(bin).pair_chain,fbin(bin).pair_reg);
                motif_count_shuf(midx,bin,i)=count_motif_congru(shuf,shufreg,fbin(bin).pref_chain_S1,bin,msize);
                motif_count_shuf_incong(midx,bin,i)=count_motif_incong(shuf,shufreg,fbin(bin).pref_chain_S1,bin,msize);
                %motif_count_shuf(bin,i)=count_triplet_incong_actinact(shuf,shufreg,fbin(bin).pref_chain_S1,bin);
            end
        end
    end

    keyboard
    save('motif_count.mat','motif_count','motif_count_shuf','motif_count_incong','motif_count_shuf_incong')
end
%    parfor bin=1:6
%        disp(bin)
%        triplet_nonsel_count(bin)=count_triplet_congru(fnbin(bin).conn_chain_S1,fbin(bin).pair_reg);
%    end
%    save('triplet_count.mat','triplet_nonsel_count','-append')

parfor reg_target=1:length(reg_set)
    if ismember(reg_target,greymatter)
        disp(reg_target)
        for bin=1:6
            corrlist=[];
            for su=reshape(su_set,1,[])
                su2reg=nnz(fbin(bin).conn_chain_S1(:,1)==su & fbin(bin).reg_chain_S1(:,2)==reg_target);
                reg2su=nnz(fbin(bin).conn_chain_S1(:,2)==su & fbin(bin).reg_chain_S1(:,1)==reg_target);
                corrlist=[corrlist;su2reg,reg2su];
            end
            recip_count(reg_target,bin)=nnz(corrlist(:,1)>0 & corrlist(:,2)>0);
            output_from_count(reg_target,bin)=nnz(corrlist(:,2)>0);
            input_to_count(reg_target,bin)=nnz(corrlist(:,1)>0);
            
            for rr=1:rpt
                [shuf,shufreg]=shuffle_conn_chain(fbin(bin).conn_chain_S1,fbin(bin).pair_chain,fbin(bin).pair_reg);
                corrlist=[];
                for su=reshape(su_set,1,[])
                    su2reg=nnz(shuf(:,1)==su & shufreg(:,2)==reg_target);
                    reg2su=nnz(shuf(:,2)==su & shufreg(:,1)==reg_target);
                    corrlist=[corrlist;su2reg,reg2su];
                end
                recip_count_shuf(reg_target,bin,rr)=nnz(corrlist(:,1)>0 & corrlist(:,2)>0);
                output_from_count_shuf(reg_target,bin,rr)=nnz(corrlist(:,2)>0);
                input_to_count_shuf(reg_target,bin,rr)=nnz(corrlist(:,1)>0);
            end

        end
        %keyboard
    end
end
%keyboard
%save('su_region_recip.mat','per_reg_recip_idx') 
save('su_region_recip.mat','recip_count','input_to_count','output_from_count','recip_count_shuf','input_to_count_shuf','output_from_count_shuf')
keyboard
load('su_region_recip.mat')
wing=diff(ioselstats{1}(:,[3 7]),1,2);
per_reg_recip_idx=[per_reg_recip_idx,(1:140)',wing];
per_reg_recip_idx=per_reg_recip_idx(~isnan(per_reg_recip_idx(:,1)),:);

[~,idces]=sort(per_reg_recip_idx(:,1));
fh=figure('Color','w','Position',[0,100,1600,400]);
bar(1:length(idces),per_reg_recip_idx(idces,1))
set(gca,'XTick',1:length(idces),'XTickLabel',reg_set(per_reg_recip_idx(idces,2)),'XTickLabelRotation',90)
ylim([0,1])
ylabel('reciprocal fraction');
if ~verLessThan('matlab','9.8')
    exportgraphics(fh,'su_region_recip.pdf')
end

figure('Color','w','Position',[100,100,400,300]);
scatter(per_reg_recip_idx(:,1),per_reg_recip_idx(:,3),'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
xlabel('reciprocal fraction');
ylabel('WING');
[r,p]=corr(per_reg_recip_idx(:,1),per_reg_recip_idx(:,3));
legend(sprintf('r=%.3f,p=%.3f',r,p));
if ~verLessThan('matlab','9.8')
    exportgraphics(fh,'su_region_recip_WING.pdf')
end

function [out,caIdx]=count_motif_congru(in,reg,pref,bin,msize,caIdx)
out=0;
pre_unit_set=unique(in(:,1));
for i=reshape(pre_unit_set,1,[])
    mono_post=in(in(:,1)==i & reg(:,1)~=reg(:,2) & pref(:,bin)==pref(:,bin+6) & pref(:,bin)>0,2);
    for j=reshape(mono_post,1,[])
        sec_post=in(in(:,1)==j & reg(:,1)~=reg(:,2) & pref(:,bin)==pref(:,bin+6) & pref(:,bin)>0,2);
        if msize==3
            out=out+nnz(ismember(in(:,1),sec_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & pref(:,bin)==pref(:,bin+6) & pref(:,bin)>0);
        elseif msize==4
            for k=reshape(sec_post,1,[])
                third_post=in(in(:,1)==k & reg(:,1)~=reg(:,2) & pref(:,bin)==pref(:,bin+6) & pref(:,bin)>0,2);
                out=out+nnz(ismember(in(:,1),third_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & pref(:,bin)==pref(:,bin+6) & pref(:,bin)>0);
            end
        elseif msize==5
            for k=reshape(sec_post,1,[])
                third_post=in(in(:,1)==k & reg(:,1)~=reg(:,2) & pref(:,bin)==pref(:,bin+6) & pref(:,bin)>0,2);
                for l=reshape(third_post,1,[])
                    if l==i || l==j
                        continue;
                    end
                    fourth_post=in(in(:,1)==l & reg(:,1)~=reg(:,2) & pref(:,bin)==pref(:,bin+6) & pref(:,bin)>0,2);
                    fourth_post(ismember(fourth_post,[i j k l]))=[];
                    out=out+nnz(ismember(in(:,1),fourth_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & pref(:,bin)==pref(:,bin+6) & pref(:,bin)>0);
                end
            end
        end
    end
end
end

function out=count_motif_incong(in,reg,pref,bin,msize)
out=0;
pre_unit_set=unique(in(:,1));
for i=reshape(pre_unit_set,1,[])
    mono_post=in(in(:,1)==i & reg(:,1)~=reg(:,2) & sum(pref(:,[bin,bin+6]),2)==3,2);
    for j=reshape(mono_post,1,[])
        sec_post=in(in(:,1)==j & reg(:,1)~=reg(:,2) & sum(pref(:,[bin,bin+6]),2)==3,2);
        if msize==3
            out=out+nnz(ismember(in(:,1),sec_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & sum(pref(:,[bin,bin+6]),2)==3);
        elseif msize==4
            for k=reshape(sec_post,1,[])
                third_post=in(in(:,1)==k & reg(:,1)~=reg(:,2) & sum(pref(:,[bin,bin+6]),2)==3,2);
                out=out+nnz(ismember(in(:,1),third_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & sum(pref(:,[bin,bin+6]),2)==3);
            end
        elseif msize==5
            for k=reshape(sec_post,1,[])
                third_post=in(in(:,1)==k & reg(:,1)~=reg(:,2) & sum(pref(:,[bin,bin+6]),2)==3 ,2);
                for l=reshape(third_post,1,[])
                    if l==i || l==j
                        continue;
                    end
                    fourth_post=in(in(:,1)==l & reg(:,1)~=reg(:,2) & sum(pref(:,[bin,bin+6]),2)==3,2);
                    fourth_post(ismember(fourth_post,[i j k l]))=[];
                    out=out+nnz(ismember(in(:,1),fourth_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & sum(pref(:,[bin,bin+6]),2)==3);
                end
            end
        end
    end
end
end


function [out,caIdx]=count_motif_congru_actinact(in,reg,pref,bin,msize,caIdx)
out=0;
pre_unit_set=unique(in(:,1));
for i=reshape(pre_unit_set,1,[])
    %keyboard
    mono_post=in(in(:,1)==i & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)==max(pref(:,7:12),[],2) & max(pref(:,1:6),[],2)>0,2);
    for j=reshape(mono_post,1,[])
        sec_post=in(in(:,1)==j & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)==max(pref(:,7:12),[],2) & max(pref(:,1:6),[],2)>0,2);
        if msize==3
            out=out+nnz(ismember(in(:,1),sec_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)==max(pref(:,7:12),[],2) & max(pref(:,1:6),[],2)>0);
        elseif msize==4
            for k=reshape(sec_post,1,[])
                third_post=in(in(:,1)==k & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)==max(pref(:,7:12),[],2) & max(pref(:,1:6),[],2)>0,2);
                out=out+nnz(ismember(in(:,1),third_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)==max(pref(:,7:12),[],2) & max(pref(:,1:6),[],2)>0);
            end
        elseif msize==5
            for k=reshape(sec_post,1,[])
                third_post=in(in(:,1)==k & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)==max(pref(:,7:12),[],2) & max(pref(:,1:6),[],2)>0,2);
                for l=reshape(third_post,1,[])
                    if l==i || l==j
                        continue;
                    end
                    fourth_post=in(in(:,1)==l & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)==max(pref(:,7:12),[],2) & max(pref(:,1:6),[],2)>0,2);
                    fourth_post(ismember(fourth_post,[i j k l]))=[];
                    out=out+nnz(ismember(in(:,1),fourth_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)==max(pref(:,7:12),[],2) & max(pref(:,1:6),[],2)>0);
                end
            end
        end
    end
end
end

function out=count_triplet_incong_actinact(in,reg,pref,bin)
out=0;
pre_unit_set=unique(in(:,1));
for i=reshape(pre_unit_set,1,[])
    %keyboard
    mono_post=in(in(:,1)==i & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)~=max(pref(:,7:12),[],2) & max(pref(:,1:6),[],2)>0 & max(pref(:,7:12),[],2)>0,2);
    for j=reshape(mono_post,1,[])
        sec_post=in(in(:,1)==j & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)~=max(pref(:,7:12),[],2) & max(pref(:,1:6),[],2)>0 & max(pref(:,7:12),[],2)>0,2);
        out=out+nnz(ismember(in(:,1),sec_post) & in(:,2)==i & reg(:,1)~=reg(:,2) & max(pref(:,1:6),[],2)~=max(pref(:,7:12),[],2) & max(pref(:,1:6),[],2)>0 & max(pref(:,7:12),[],2)>0);
    end
end
end

function [outconn,outreg]=shuffle_conn_chain(in,inpair,inreg)
outconn=nan(size(in));
outreg=nan(size(in));
for i=1:114
    lbound=i*100000;
    ubound=(i+1)*100000;
    sel=in(:,1)>lbound & in(:,1)<ubound;
    if nnz(sel)>0
        selpair=find(inpair(:,1)>lbound & inpair(:,1)<ubound);
        shufsel=randperm(nnz(selpair));
%        keyboard
        shufdata=inpair(selpair(shufsel(1:nnz(sel))),:);
        flipsel=randi(2,size(shufdata,1),1)>1;
        shufdata(flipsel,:)=shufdata(flipsel,[2,1]);
        outconn(sel,:)=shufdata;        
        if exist('inreg','var') 
            outreg(sel,:)=inreg(selpair(shufsel(1:nnz(sel))),:);
        end
    end
end
%    outconn=in;
%    outreg=inreg;
%    for i=1:114
%        lbound=i*100000;
%        ubound=(i+1)*100000;
%        sel=in(:,1)>lbound & in(:,1)<ubound;
%        if nnz(sel)>0
%            seldata=in(sel,:);
%            shufsel=randperm(nnz(sel));
%            shufdata=seldata(shufsel,2);
%            out(sel,2)=shufdata;        
%
%            selreg=inreg(sel,:);
%            shufreg=selreg(shufsel,2);
%            outreg(sel,2)=shufreg;
%        end
%    end
end
function savefun(fn,session,bin,msize,motif_candidate_count)
    save(fn,'session','bin','msize','motif_candidate_count');
end
