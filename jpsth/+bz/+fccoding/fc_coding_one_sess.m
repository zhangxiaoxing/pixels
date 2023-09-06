function onesess=fc_coding_one_sess(fidx,suids,opt)
arguments
    fidx (1,1) double {mustBeInteger,mustBePositive}
    suids (:,2) int32
    opt.circular_rpt (1,1) double {mustBeInteger,mustBeNonnegative} = 0%100
%     opt.window (1,2) double =[0.0008,0.01] % todo implement variable
%     window
    opt.per_trial (1,1) logical = false
end

ephys.util.dependency('buz',false);
%Fieldtrip routine
[~,~,trials,~,~,FT_SPIKE]=ephys.getSPKID_TS(fidx,'keep_trial',true,'suids',unique(suids),'only_delay',true);
%TODO: possible trial number criteria
onesess=struct();
onesess.trials=trials;
onesess.fidx=fidx;
onesess.fc=[];

for fci=1:size(suids,1)
    if rem(fci,100)==0,disp(fci);end
    %per fc stats
%     accuDelta=[];
    [fwd_fc,rev_fc,lead_evt,follow_evt]=deal(zeros(size(trials,1),1));
    [fwd_shift,rev_shift]=deal(zeros(size(trials,1),opt.circular_rpt));

    for ti=1:size(trials,1)
        susel1=strcmp(FT_SPIKE.label,num2str(suids(fci,1)));
        susel2=strcmp(FT_SPIKE.label,num2str(suids(fci,2)));
        %TODO per trial, per bin >>>>>>>>>>>>>>
        %% result-critical window width--------------------------------------------vvvv-----------------------vvvv
        spk1=FT_SPIKE.time{susel1}(FT_SPIKE.trial{susel1}==ti & FT_SPIKE.time{susel1}<3 & FT_SPIKE.time{susel1}>=0);
        spk2=FT_SPIKE.time{susel2}(FT_SPIKE.trial{susel2}==ti & FT_SPIKE.time{susel2}<3 & FT_SPIKE.time{susel2}>=0);
        
        deltaT=spk2-spk1.';
%         accuDelta=[accuDelta;deltaT(:)];
        fwd_evts=deltaT>=0.0008 & deltaT<0.01; 
        rev_evts=deltaT>=-0.01 & deltaT<=-0.0008;
        fwd_fc(ti)=nnz(any(fwd_evts)); % total number of evoked follow-spk
        rev_fc(ti)=nnz(any(rev_evts,2)); % total number of 'evoked' lead-spk
        lead_evt(ti)=numel(spk1);
        follow_evt(ti)=numel(spk2);
        if opt.circular_rpt>0
            for rpt=1:opt.circular_rpt
                spk2=rem(spk2+rand()*trials(ti,8),trials(ti,8));
                deltaT=spk2-spk1.';
                fwd_evts=deltaT>=0.0008 & deltaT<0.01;
                rev_evts=deltaT>=-0.01 & deltaT<=-0.0008;
                fwd_shift(ti,rpt)=nnz(any(fwd_evts));
                rev_shift(ti,rpt)=nnz(any(rev_evts,2));
            end
        end
    end
    onesess.fc=[onesess.fc;{suids(fci,:),fwd_fc,fwd_fc-mean(fwd_shift,2),fwd_shift,rev_fc,rev_fc-mean(rev_shift,2),rev_shift,lead_evt,follow_evt}];
end
blame=vcs.blame();
save(fullfile('binary','fccoding',sprintf('fc_coding_%d.mat',fidx)),'onesess','blame')
end

% function plotonesess(onesess)
% s1sel=onesess.trials(:,5)==4 & onesess.trials(:,8)==6 & all(onesess.trials(:,9:10),2);
% s2sel=onesess.trials(:,5)==8 & onesess.trials(:,8)==6 & all(onesess.trials(:,9:10),2);
% figure()
% hold on
% for i=1:size(onesess.fc)
%     scatter(mean(onesess.fc{i,3}(s1sel)),mean(onesess.fc{i,3}(s2sel)),4,'r','Marker','.');
% end
% end
