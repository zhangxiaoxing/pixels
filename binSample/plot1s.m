fo = fopen('sample1s.ap.bin');
mat = fread(fh,'int16');
fclose(fo);

per_ch=reshape(mat,385,[]);
for jj=1:3000:30000
    for ii=1:11:375
        fh=figure();
        ph=plot(per_ch(ii:ii+10,jj:jj+2999).'+(50:50:550));
        waitfor(fh)
    end
end