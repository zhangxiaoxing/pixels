function anovameta=get_dur_waveid()
%1:Samp,2:Dur,3:Bin,4:Samp*Dur,5:Samp*Bin,6:Dur*Bin,7:Samp*Dur*Bin
anovameta=ephys.selectivity_anova();
senssel=any(anovameta.anovap(:,[1 4 5 7])<0.05,2);
dursel=any(anovameta.anovap(:,[2 6])<0.05,2);
w=nan(size(anovameta.sess));
w(~dursel)=0;
w(dursel & ~senssel & anovameta.dur_selidx>=0)=3;
w(dursel & ~senssel & anovameta.dur_selidx<=0)=6;
anovameta.dur_waveid=w;

end