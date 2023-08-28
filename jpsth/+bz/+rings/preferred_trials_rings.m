function [pref_trl3,pref_trl6]=preferred_trials_rings(curr_waveid,trials)
if isnumeric(curr_waveid)
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
elseif isstring(curr_waveid)
    switch curr_waveid
        case "s1d3"
            pref_trl3=find(trials(:,5)==4 & trials(:,8)==3 & all(trials(:,9:10)>0,2));
            pref_trl6=[];
        case "s1d6"
            pref_trl6=find(trials(:,5)==4 & trials(:,8)==6 & all(trials(:,9:10)>0,2));
            pref_trl3=[];
        case "s2d3"
            pref_trl3=find(trials(:,5)==8 & trials(:,8)==3 & all(trials(:,9:10)>0,2));
            pref_trl6=[];
        case "s2d6"
            pref_trl6=find(trials(:,5)==8 & trials(:,8)==6 & all(trials(:,9:10)>0,2));
            pref_trl3=[];
        otherwise
            pref_trl3=[];
            pref_trl6=[];
            keyboard();
    end
else
    error("Not implemented")
    keyboard()
end
end