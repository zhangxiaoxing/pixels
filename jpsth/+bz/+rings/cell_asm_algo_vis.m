% TODO cell assembly as algorithm
% Visualization
steps=[10,20,50,100,200,500];

%% Rate selectivity index
% Selective
% SU

% FCSP

% Ring loops



%% Rate decoding
% Units vs decoding accuracy
% Selective congruent
% SU
su_data=bz.rings.cellasm.get_SU_data();

suprefs=cellfun(@(x) mean(x('pref')),su_data);
sunonps=cellfun(@(x) mean(x('nonp')),su_data);
sufreqsel=suprefs>6;
%decoding
su_dec=bz.rings.cellasm.decode_one(su_data(sufreqsel),'unit_count',steps,'pre_pca',true);


% FCSP
fcsp_data=bz.rings.cellasm.get_fcsp_data();

fcprefs=cellfun(@(x) mean(x('pref')),fcsp_data);
fcnonps=cellfun(@(x) mean(x('nonp')),fcsp_data);
fcfreqsel=fcprefs>6;
fcspdec=bz.rings.cellasm.decode_one(fcsp_data(fcfreqsel),'unit_count',steps,'pre_pca',true);


% Ring loops
rings_data=bz.rings.cellasm.get_rings_data();

rnprefs=cellfun(@(x) mean(x('pref')),rings_data);
rnnonps=cellfun(@(x) mean(x('nonp')),rings_data);
rnfreqsel=rnprefs>1;
%decoding
ringdec=bz.rings.cellasm.decode_one(rings_data(rnfreqsel),'unit_count',steps,'pre_pca',true);


%selectivity index
selidx=(prefs(freqsel)-nonps(freqsel))./(prefs(freqsel)+nonps(freqsel));
%total number of su in rings
meta=cellfun(@(x) x('meta'),rings_data,'UniformOutput',false);
numel(unique(cell2mat(cellfun(@(x) x{4}*100000+x{5}, (meta(freqsel)).','UniformOutput',false))))



% Random-ML

[suhat,suci,sucnt]=bz.rings.cellasm.dec_sums(su_dec);
[fchat,fcci,fccnt]=bz.rings.cellasm.dec_sums(fcspdec);
[rnhat,rnci,rncnt]=bz.rings.cellasm.dec_sums(ringdec);

% same number of total SU
fh1=figure('Color','w','Position',[32,32,600,300]);
hold on
fill([sucnt;flip(sucnt)],[suci(:,1);flip(suci(:,2))].*100,'k','EdgeColor','none','FaceAlpha',0.1);
fill([fccnt;flip(fccnt)],[fcci(:,1);flip(fcci(:,2))].*100,'b','EdgeColor','none','FaceAlpha',0.1);
fill([rncnt;flip(rncnt)],[rnci(:,1);flip(rnci(:,2))].*100,'r','EdgeColor','none','FaceAlpha',0.1);
psh=plot(sucnt,suhat*100,'-k');
pfh=plot(fccnt,fchat*100,'-b');
prh=plot(rncnt,rnhat*100,'-r');
ylabel('Classifcation accuracy (%)')
xlabel('Number of neurons in recording');
xlim([0,200])
legend([psh,pfh,prh],{'Neuron FR','Function-coupling event rate','Ring-asm loops rate'},'Location','eastoutside','Orientation','vertical');
title('Ring-assembly reduce front-end channels')
exportgraphics(fh1,'cellasm_algo_same_su.pdf','ContentType','vector')
% same number of total data-channel
fh2=figure('Color','w','Position',[32,32,360,300])
hold on
fill([sucnt;flip(sucnt)],[suci(:,1);flip(suci(:,2))].*100,'k','EdgeColor','none','FaceAlpha',0.1);
fill([sucnt;flip(sucnt)],[fcci(:,1);flip(fcci(:,2))].*100,'b','EdgeColor','none','FaceAlpha',0.1);
fill([sucnt;flip(sucnt)],[rnci(:,1);flip(rnci(:,2))].*100,'r','EdgeColor','none','FaceAlpha',0.1);
psh=plot(sucnt,suhat*100,'-k');
pfh=plot(sucnt,fchat*100,'-b');
prh=plot(sucnt,rnhat*100,'-r');
ylabel('Classifcation accuracy (%)')
xlabel('Number of feature entities for decoding');
xlim([0,200])
% legend([psh,pfh,prh],{'Neuron FR','Function-coupling event rate','Ring-asm loops rate'},'Location','northoutside','Orientation','horizontal');
title('FCSP reduce decoding feature input')
exportgraphics(fh2,'cellasm_algo_same_feat.pdf','ContentType','vector')



fh3=figure('Color','w','Position',[32,32,200,150]);
hold on
rn_su_sem=cellfun(@(x) std(x.su_count)./sqrt(numel(x.su_count)),ringdec);
fc_su_sem=cellfun(@(x) std(x.su_count)./sqrt(numel(x.su_count)),fcspdec);
fill([rncnt+rn_su_sem;flip(rncnt-rn_su_sem)],[sucnt;flip(sucnt)],'r','EdgeColor','none','FaceAlpha',0.2);
fill([fccnt+fc_su_sem;flip(fccnt-fc_su_sem)],[sucnt;flip(sucnt)],'b','EdgeColor','none','FaceAlpha',0.2);
prh=plot(rncnt,sucnt,'r');
pfh=plot(fccnt,sucnt,'b');
psh=plot([0,500],[0,500],'-k');
xlabel('Associated neuron number')
ylabel('Measured entity number')
% legend([psh,pfh,prh],{'Neuron','Functional-coupling','Ring-assembly'},'Location','northoutside','Orientation','horizontal');
exportgraphics(fh3,'cellasm_algo_su_vs_entity.pdf','ContentType','vector')

%% get data set, preferred vs non-preferred
% SU|feature:cell * label:map * trial:vector


