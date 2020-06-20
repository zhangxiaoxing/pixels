%assume loaded xcorr sums file here
%assume loaded sumsUniq file here

%needs plot xcorr per pair, psth per neuron
close all
bin1=unique(chainSumsUniq(:,1:2),'rows');
for i=1:length(bin1)
    %xcorr
    su1=bin1(i,1);
    su2=bin1(i,2);
    fidx=idivide(su1,int32(100000));
    xcorr=sums{fidx,5};
    su1_lbl=rem(su1,100000);
    su1_lbl_idx=find(strcmp(xcorr.label(:,1),num2str(su1_lbl)));
    su2_lbl=rem(su2,100000);
    su2_lbl_idx=find(strcmp(xcorr.label(:,1),num2str(su2_lbl)));
    
    figure();
    plot(squeeze(xcorr.xcorr(su1_lbl_idx,su2_lbl_idx,:)))
    keyboard
    %psth
end
