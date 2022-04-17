function anovameta=get_dur_waveid()
%1:Samp,2:Dur,3:Bin,4:Samp*Dur,5:Samp*Bin,6:Dur*Bin,7:Samp*Dur*Bin
anovameta=ephys.selectivity_anova();
%1:Samp,2:Dur,3:Bin,4:Samp*Dur,5:Samp*Bin,6:Dur*Bin,7:Samp*Dur*Bin
dur_sel_mix=any(anovameta.anovap(:,[2 4 6 7])<0.05,2);
sens_sel_mix=any(anovameta.anovap(:,[1 4 5 7])<0.05,2);
dur_sel_exclu=dur_sel_mix & ~sens_sel_mix;

anovameta.dur_waveid=zeros(size(anovameta.sess));
anovameta.dur_waveid(dur_sel_exclu & anovameta.dur_selidx>0)=3;
anovameta.dur_waveid(dur_sel_exclu & anovameta.dur_selidx<0)=6;
end