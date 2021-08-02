function onesess=fc_coding_one_sess(fidx,suids,opt)
arguments
    fidx (1,1) double {mustBeInteger,mustBePositive}
    suids (:,2) int32
    opt.circular_rpt (1,1) double {mustBeInteger,mustBePositive} = 100
end

%Fieldtrip routine
[~,~,trials,~,~,FT_SPIKE]=ephys.getSPKID_TS(fidx,'keep_trial',true,'suids',unique(suids),'only_delay',true);
%TOD: apply trial number criteria
onesess=struct();
onesess.trials=trials;
onesess.fidx=fidx;
onesess.fc=[];

for fci=1:size(suids,1)
    if rem(fci,100)==0,disp(fci);end
    %per fc stats
    [onefc,onefc_b]=deal(zeros(size(trials,1),1));
    [cshift,bshift]=deal(zeros(size(trials,1),opt.circular_rpt));
    for ti=1:size(trials,1)
        susel1=strcmp(FT_SPIKE.label,num2str(suids(fci,1)));
        susel2=strcmp(FT_SPIKE.label,num2str(suids(fci,2)));
        spk1=FT_SPIKE.time{susel1}(FT_SPIKE.trial{susel1}==ti & FT_SPIKE.time{susel1}<trials(ti,8));
        spk2=FT_SPIKE.time{susel2}(FT_SPIKE.trial{susel2}==ti & FT_SPIKE.time{susel2}<trials(ti,8));
        deltaT=spk2-spk1.';
        evts=deltaT>=0.0008 & deltaT<0.01;
        bevts=deltaT>=-0.01 & deltaT<=-0.0008;
        onefc(ti)=nnz(any(evts'));
        onefc_b(ti)=nnz(any(bevts'));
        for rpt=1:opt.circular_rpt
            spk2=rem(spk2+rand()*trials(ti,8),trials(ti,8));
            deltaT=spk2-spk1.';
            evts=deltaT>=0.0008 & deltaT<0.01;
            bevts=deltaT>=-0.01 & deltaT<=-0.0008;
            cshift(ti,rpt)=nnz(any(evts'));
            bshift(ti,rpt)=nnz(any(bevts'));
        end
    end
    onesess.fc=[onesess.fc;{suids(fci,:),onefc,onefc-mean(cshift,2),cshift,onefc_b,onefc_b-mean(bshift,2),bshift}];
end
blame=vcs.blame();
save(sprintf('fc_coding_%d.mat',fidx),'onesess','blame')
end

function plotonesess(onesess)
s1sel=onesess.trials(:,5)==4 & onesess.trials(:,8)==6 & all(onesess.trials(:,9:10),2);
s2sel=onesess.trials(:,5)==8 & onesess.trials(:,8)==6 & all(onesess.trials(:,9:10),2);
figure()
hold on
for i=1:size(onesess.fc)
    scatter(mean(onesess.fc{i,3}(s1sel)),mean(onesess.fc{i,3}(s2sel)),4,'r','Marker','.');
end
end