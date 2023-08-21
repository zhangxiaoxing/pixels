function exception_check()
[spkID,spkTS,trials,SU_id,folder,FT_SPIKE]=ephys.getSPKID_TS(106,'keep_trial',true);
errsel=trials(:,10)==0;
sampsel=trials(:,5)==8;
dur3sel=trials(:,8)==3;
dur6sel=trials(:,8)==6;
susel=strcmp(FT_SPIKE.label,'10202');

corrct3=nnz(... % trial
    ismember(FT_SPIKE.trial{susel},find(~errsel & sampsel & dur3sel))...
    ... % time
    & FT_SPIKE.time{susel}>=1 & FT_SPIKE.time{susel}<4)...
    ... % total delay time
    ./(nnz(~errsel & sampsel & dur3sel).*3)

error3=nnz(... % trial
    ismember(FT_SPIKE.trial{susel},find(errsel & sampsel & dur3sel))...
    ... % time
    & FT_SPIKE.time{susel}>=1 & FT_SPIKE.time{susel}<4)...
    ... % total delay time
    ./(nnz(errsel & sampsel & dur3sel).*3)


corrct6=nnz(... % trial
    ismember(FT_SPIKE.trial{susel},find(~errsel & sampsel & dur3sel))...
    ... % time
    & FT_SPIKE.time{susel}>=1 & FT_SPIKE.time{susel}<7)...
    ... % total delay time
    ./(nnz(~errsel & sampsel & dur6sel).*6)

error6=nnz(... % trial
    ismember(FT_SPIKE.trial{susel},find(errsel & sampsel & dur3sel))...
    ... % time
    & FT_SPIKE.time{susel}>=1 & FT_SPIKE.time{susel}<7)...
    ... % total delay time
    ./(nnz(errsel & sampsel & dur6sel).*6)


for ii=1:size(trials,1)
    if rem(ii,40)==1
        figure()
        hold on
    end
    if all(trials(ii,9:10)==1,2)
        plot(FT_SPIKE.time{susel}(FT_SPIKE.trial{susel}==ii),ii,'k|')
    elseif trials(ii,10)==0
        plot(FT_SPIKE.time{susel}(FT_SPIKE.trial{susel}==ii),ii,'r|')
    else
        plot(FT_SPIKE.time{susel}(FT_SPIKE.trial{susel}==ii),ii,'|','Color',[0.5,0.5,0.5])
    end
end
end

