function anovameta=get_dur_waveid()
%1:Samp,2:Dur,3:Bin,4:Samp*Dur,5:Samp*Bin,6:Dur*Bin,7:Samp*Dur*Bin
anovameta=ephys.selectivity_anova();
dursel=any(anovameta.anovap(:,[2 6])<0.05,2);

anovameta.dur_waveid=zeros(size(anovameta.sess));
anovameta.dur_waveid(dursel & anovameta.dur_selidx>0)=3;
anovameta.dur_waveid(dursel & anovameta.dur_selidx<0)=6;

end