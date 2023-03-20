% SPADE related statistics were unsuitable for higher throughput recordings
% due to memory performance issues.
% Pending removal as of 2023.03.20


stpts=h5read('..\SPADE\stp_ts_export.hdf5','/spk_ts')'; %session,cid,sample,ts
stp_cnt=h5read('..\SPADE\stp_ts_export.hdf5','/spk_cnt')'; %session,cid,stp spk count, total spk count

fstr=cell(1,6);
for bin=1:6
    fstr{bin}=load(sprintf('0831_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
end

sess_ids=unique(stpts(:,1));
all_tagged=[];
for sessIdx=sess_ids'
    disp(sessIdx)
    patt_suids=unique(stpts(stpts(:,1)==sessIdx,2));
    fc_list=get_func_conn(sessIdx,patt_suids,fstr);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%
    for ssidx=1:size(fc_list,1)
        sample=fc_list(ssidx,3);
        bin=fc_list(ssidx,4);
        transIds=int32(rem(fc_list(ssidx,1:2),100000))';
        ts_id=[];
        for seqid=1:2
            ts=stpts(stpts(:,1)==sessIdx & stpts(:,2)==transIds(seqid) & stpts(:,3)==sample*4,4);
            ts(:,2)=seqid;
            ts_id=[ts_id;ts];
        end
        [~,s]=sort(ts_id(:,1));
        ts_id=ts_id(s,:);
        ts_id_tagged=fc_tag(ts_id);
        ts_id_tagged(:,4:7)=repmat([transIds',sample,sessIdx],size(ts_id_tagged,1),1);
        all_tagged=[all_tagged;ts_id_tagged];
    end
    
    %%TS,FC_SU_ID,FC_TAG,SU1_CID,SU2_CID,FC_SAMPLE,SESS_IDX
end
stp_cnt(:,5)=-1;
for i=1:size(stp_cnt,1)
    disp(i)
    sess=stp_cnt(i,1);
    cid=stp_cnt(i,2);
    stp_fc_spk=numel(unique(all_tagged(...
        all_tagged(:,7)==sess & ...
        all_tagged(:,3)==1  & ...
        ((all_tagged(:,2)==1 & all_tagged(:,4)==cid) | (all_tagged(:,2)==2 & all_tagged(:,5)==cid)),1)));
    stp_cnt(i,5)=stp_fc_spk;
    disp([stp_cnt(i,3:4),stp_fc_spk])
end
% save('all_tagged.mat','all_tagged','stp_cnt','-v7.3')
frac=stp_cnt(:,5)./stp_cnt(:,3).*100;
close all
fh=figure('Color','w','Position',[100,100,150,195]);
hold on;
plot(ones(size(frac))+rand(size(frac))*0.2-0.1,frac,'ko','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','none')
errorbar(mean(frac),std(frac)./sqrt(numel(frac)),'ro','MarkerSize',6)
set(gca,'YTick',0:25:100,'XTick',[])
xlim([0,2])
ylabel('Fraction of func. conn. spikes (%)')
exportgraphics(fh,'func_conn_frac_in_stp.pdf')



return






function out=get_func_conn(sessId,suids,fstr)
patt_su=suids+double(sessId)*100000;
out=[];
I=sessId;
lbound=100000*I;
ubound=100000*(I+1);
for bin=1:6
    sel11=fstr{bin}.conn_chain_S1(:,1)>=lbound & fstr{bin}.conn_chain_S1(:,1)<ubound & diff(fstr{bin}.reg_chain_S1,1,2);
    sel12=fstr{bin}.conn_chain_S2(:,1)>=lbound & fstr{bin}.conn_chain_S2(:,1)<ubound & diff(fstr{bin}.reg_chain_S2,1,2);
    out1=fstr{bin}.conn_chain_S1(sel11,:);
    out2=fstr{bin}.conn_chain_S2(sel12,:);
    out1(:,3)=1;
    out2(:,3)=2;
    out1(:,4)=bin;
    out2(:,4)=bin;
    out=[out;out1;out2];
    
end
out=unique(out,'rows');
patt_sel=ismember(out(:,1),patt_su) & ismember(out(:,2),patt_su);
out=out(patt_sel,:);

end


function out=fc_tag(in)
out=in;
out(:,3)=0;

for i=1:size(in,1)-1
    if in(i,2)==1
        j=i+1;
        while in(j,2)==2 && j<size(in,1)
            if in(j,1)<=in(i,1)+10 && in(j,1)>in(i,1)+2
                out(i,3)=1;
                out(j,3)=1;
            end
            j=j+1;
        end
    end
   
%             & in(rows,1)<tsseq(end,2)+0.01 ...
%             & in(rows,1)>tsseq(end,2)+0.0005 ... %matching time window, assuming 1kHz
end
end