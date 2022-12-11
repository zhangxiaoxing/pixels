function [pref_trl3,pref_trl6]=preferred_trials_rings(curr_waveid,trials)
if any(curr_waveid==1,'all')
    pref_trl3=find(trials(:,5)==4 & trials(:,8)==3 & all(trials(:,9:10)>0,2));
    pref_trl6=[];
elseif any(curr_waveid==2,'all')
    pref_trl6=find(trials(:,5)==4 & trials(:,8)==6 & all(trials(:,9:10)>0,2));
    pref_trl3=[];
elseif any(curr_waveid==3,'all')
    pref_trl3=find(trials(:,5)==8 & trials(:,8)==3 & all(trials(:,9:10)>0,2));
    pref_trl6=[];
elseif any(curr_waveid==4,'all')
    pref_trl6=find(trials(:,5)==8 & trials(:,8)==6 & all(trials(:,9:10)>0,2));
    pref_trl3=[];
elseif any(curr_waveid==5,'all')
    pref_trl3=find(trials(:,5)==4 & trials(:,8)==3 & all(trials(:,9:10)>0,2));
    pref_trl6=find(trials(:,5)==4 & trials(:,8)==6 & all(trials(:,9:10)>0,2));
elseif any(curr_waveid==6,'all')
    pref_trl3=find(trials(:,5)==8 & trials(:,8)==3 & all(trials(:,9:10)>0,2));
    pref_trl6=find(trials(:,5)==8 & trials(:,8)==6 & all(trials(:,9:10)>0,2));
elseif any(curr_waveid==7,'all')
    pref_trl3=find(trials(:,8)==3 & all(trials(:,9:10)>0,2));
    pref_trl6=[];
elseif any(curr_waveid==8,'all')
    pref_trl3=[];
    pref_trl6=find(trials(:,8)==6 & all(trials(:,9:10)>0,2));
else
    pref_trl3=[];
    pref_trl6=[];
    keyboard();
end

end