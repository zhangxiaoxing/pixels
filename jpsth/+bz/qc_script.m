qcmat=cell2mat({sums_conn_str.qc}.');
figure()
subplot(1,3,1)
histogram(qcmat(qcmat(:,1)>0,4),0:1:51)
xlim([0,50])
title('FWHM')
subplot(1,3,2)
histogram(qcmat(qcmat(:,1)>0,3),0:1:51)
xlim([0,50])
title('Noise')

formsel=(qcmat(:,1)>0 & qcmat(:,3)<10 &qcmat(:,4)<5);
subplot(1,3,3)
histogram(qcmat(formsel,2),240:270)
xline(251)

title('Tpeak')

qcsel=(qcmat(:,1)>0 & qcmat(:,3)<10 &qcmat(:,4)<5 & qcmat(:,2)>253 & qcmat(:,2)<257);


%%%%%%%% following part is for visual inspection %%%%%%%%

scmat=cell2mat({sums_conn_str.ccg_sc}.');
ccgmat=scmat(:,7:end);
for ii=1:100:size(ccgmat,1)
    fh=figure('Color','w','Position',[100,100,800,600]);
    for jj=0:99
        subplot(3,4,jj+1)
        plot(-100:0.4:100,ccgmat(ii+jj,:),'-r')
        arrayfun(@(x) xline(x,'--k'),[-10,0,10]);
        xlim([-40,40]);
    end
    waitfor(fh);
end

