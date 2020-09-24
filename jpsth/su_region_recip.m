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

greymatter=find(cellfun(@(x) ~isempty(regexp(x,'[A-Z]','match','once')), reg_set));

su_set=unique(fbin(1).conn_chain_S1(:,1));
input_to_count=nan(length(reg_set),6);
output_from_count=nan(length(reg_set),6);
recip_count=nan(length(reg_set),6);


input_to_count_shuf=nan(length(reg_set),6,rpt);
output_from_count_shuf=nan(length(reg_set),6,rpt);
recip_count_shuf=nan(length(reg_set),6,rpt);


%% candidate pairs
if exist('candidate_stats_easy','var') && candidate_stats_easy
%     candidate_delay=nan(2,6,3);
    for midx=1:3
        subtotal1=0;
        subtotal2=0;
        msize=midx+2;
        for bin=1:6
            for session=1:114 %7bit session, 3bit bin, 2bit msize
                disp(session)
                cc=fbin(bin).pair_chain;
                sel=cc(:,1)>=session*100000 & cc(:,1)<(session+1)*100000;
                if nnz(sel)==0
                    continue
                end
                connsel=cc(sel,:);
                regset=unique(fbin(bin).pair_reg(sel,:));
                reg_ring_set=nchoosek(regset,msize);
                
                for setIdx=1:size(reg_ring_set,1)
                    currcount1=1;
                    currcount2=1;
                    for mcounter=1:msize
                        if currcount1>0
                            sel1=fbin(bin).pair_reg(sel,:)==reg_ring_set(setIdx,mcounter) &...
                                [fbin(bin).pref_pair(sel,bin)==1, fbin(bin).pref_pair(sel,bin+6)==1];
                            count1=numel(unique(connsel(sel1)));
                            currcount1=currcount1*count1;
                        end
                        if currcount2>0
                            sel2=fbin(bin).pair_reg(sel,:)==reg_ring_set(setIdx,mcounter) &...
                                [fbin(bin).pref_pair(sel,bin)==2, fbin(bin).pref_pair(sel,bin+6)==2];
                            count2=numel(unique(connsel(sel2)));
                            currcount2=currcount2*count2;
                        end
                    end
                    
                    subtotal1=subtotal1+currcount1;
                    subtotal2=subtotal2+currcount2;
                end
            end
            candidate_delay(1,bin,midx)=subtotal1;
            candidate_delay(2,bin,midx)=subtotal2;
            save('candidate_count.mat','candidate_delay','-append')
        end
    end
    return
end





%% candidate baseline pairs easy
if exist('candidate_base_stats_easy','var') && candidate_base_stats_easy
    count_all=nan(2,3);
    for midx=1:3
        subtotal1=0;
        subtotal2=0;
        msize=midx+2;
        for session=1:114 %7bit session, 3bit bin, 2bit msize
            disp(session)
            cc=fbase.pair_chain;
            sel=cc(:,1)>=session*100000 & cc(:,1)<(session+1)*100000;
            if nnz(sel)==0
                continue
            end
            connsel=cc(sel,:);
            regset=unique(fbase.pair_reg(sel,:));
            reg_ring_set=nchoosek(regset,msize);
            
            for setIdx=1:size(reg_ring_set,1)
                currcount1=1;
                currcount2=1;
                for mcounter=1:msize
                    sel1=fbase.pair_reg(sel,:)==reg_ring_set(setIdx,mcounter) &...
                        [max(fbase.pref_pair(sel,1:6),[],2)==1, max(fbase.pref_pair(sel,7:12),[],2)==1];
                    count1=numel(unique(connsel(sel1)));
                    currcount1=currcount1*count1;
                    
                    sel2=fbase.pair_reg(sel,:)==reg_ring_set(setIdx,mcounter) &...
                        [max(fbase.pref_pair(sel,1:6),[],2)==2, max(fbase.pref_pair(sel,7:12),[],2)==2];
                    count2=numel(unique(connsel(sel2)));
                    currcount2=currcount2*count2;
                end
                
                subtotal1=subtotal1+currcount1;
                subtotal2=subtotal2+currcount2;
            end
        end
        count_all(1,midx)=subtotal1;
        count_all(2,midx)=subtotal2;
    end
return
end


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
