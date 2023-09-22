timedur=3600*1000;
FR1=10;
FR2=10;
FR3=10;

ts1=randi(timedur,FR1*3600,1)./1000+randi(10,FR1*3600,1)./10000;
ts2=[randsample(ts1,numel(ts1)./100)+0.002;randi(timedur,FR1*3600*0.99,1)./1000]+randi(10,FR1*3600,1)./10000;
ts3=randi(timedur,FR1*3600,1)./1000+randi(10,FR1*3600,1)./10000;
ephys.util.dependency(ft=false);
mono_res=bz.sortSpikeIDz([ts1;ts2;ts3],[ones(size(ts1));2.*ones(size(ts2));3.*ones(size(ts3))],'wallclocktime',true);

ccgfreq=mono_res.ccgR(:,1,2)./numel(ts1)./0.0004;
ccgbounds=squeeze(mono_res.Bounds(:,1,2,:)./numel(ts1)./0.0004);
figure()
hold on
plot(ccgfreq,'r-','LineWidth',1)
plot(ccgbounds,'m:','LineWidth',1)
xline([251-(10/0.4),251,251+(10/0.4)],'k--')
xlim([251-(12/0.4),251+(25/0.4)])
set(gca(),'XTick',226:25:310,'XTickLabel',-10:10:20)
xlabel('Time (ms)')
ylabel('Instant frequency (Hz)')
title('Butzaki')

ts1tsers=histcounts(ts1,0:0.0004:3600);
ts2tsers=histcounts(ts2,0:0.0004:3600);

gk=fspecial('gaussian',[1,201],25);
gc1=conv(ts1tsers,gk,"same");
gc2=conv(ts2tsers,gk,"same");

[r,lags]=xcorr(ts2tsers,ts1tsers,250);
[rc,lagsc]=xcorr(gc2,gc1,250);

figure();
hold on
plot(lags,r,'b-','LineWidth',1);
plot(lagsc,rc,'c--','LineWidth',1);
xlim([-12/0.4,25/0.4])
xline([-(10/0.4),0,(10/0.4)],'k--')
set(gca(),'XTick',-10/0.4:10/0.4:20/0.4,'XTickLabel',-10:10:20)
xlabel('Time (ms)')
ylabel('Dot product?')
title('Miyashita')
