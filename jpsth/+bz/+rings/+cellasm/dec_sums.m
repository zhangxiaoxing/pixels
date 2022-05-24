function [phats,pcis,ucnts]=dec_sums(dec_data)
[phats,pcis,ucnts]=deal([]);
for c=dec_data.'
%     keyboard
    [phat,pci]=binofit(nnz(c{1}.result),numel(c{1}.result));
    ucnt=mean(c{1}.su_count);
    phats=[phats;phat];
    pcis=[pcis;pci];
    ucnts=[ucnts;ucnt];
end
end