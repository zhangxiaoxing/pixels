
function [overlapS1,overlapS2]=overlap_calc(mono_rez_S1,mono_rez_S2,sess,bin_range)
    if ~exist(sprintf('0831_selec_conn_chain_duo_6s_%d_%d.mat',bin_range(1),bin_range(2)),'file')
        overlapS1=[];
        overlapS2=[];
        return
    end
    fcon=load(sprintf('0831_selec_conn_chain_duo_6s_%d_%d.mat',bin_range(1),bin_range(2)));
    overlapS1=[0,0,size(mono_rez_S1.sig_con,1)];
    for i=1:size(mono_rez_S1.sig_con,1)
        t=mono_rez_S1.LUT(mono_rez_S1.sig_con(i,:))+sess*100000;
%         if ~any((fcon.pair_chain(:,1)==t(1) & fcon.pair_chain(:,2)==t(2)) |...
%             (fcon.pair_chain(:,1)==t(2) & fcon.pair_chain(:,2)==t(1)))
%             keyboard
%         end
%         disp(t);
        overlapS1(1)=overlapS1(1)+any(fcon.conn_chain_S1(:,1)==t(1) & fcon.conn_chain_S1(:,2)==t(2));
        overlapS1(2)=overlapS1(2)+any((fcon.pair_chain(:,1)==t(1) & fcon.pair_chain(:,2)==t(2))...
            |(fcon.pair_chain(:,1)==t(2) & fcon.pair_chain(:,2)==t(1)));
    end
    
    overlapS2=[0,0,size(mono_rez_S2.sig_con,1)];
    for i=1:size(mono_rez_S2.sig_con,1)
        t=mono_rez_S2.LUT(mono_rez_S2.sig_con(i,:))+sess*100000;
%         disp(t);
        overlapS2(1)=overlapS2(1)+any(fcon.conn_chain_S2(:,1)==t(1) & fcon.conn_chain_S2(:,2)==t(2));
        overlapS2(2)=overlapS2(2)+any((fcon.pair_chain(:,1)==t(1) & fcon.pair_chain(:,2)==t(2))...
        | (fcon.pair_chain(:,1)==t(2) & fcon.pair_chain(:,2)==t(1)));
    end
end