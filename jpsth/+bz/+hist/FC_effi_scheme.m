t=0:1000;
b0=exp(-t/100);
b300=[zeros(1,300),b0(1:700)];
b330=[zeros(1,330),b0(1:670)];
b630=[zeros(1,630),b0(1:370)];
fh=figure('Color','w');
hold on

plot(b300,':r');
plot(b330,':r');
plot(b630,':r');
plot(b300+b330+b630,'-k');
set(gca(),'XTick',[],'YTick',[])
exportgraphics(fh,'FC_effi_schematics.pdf')