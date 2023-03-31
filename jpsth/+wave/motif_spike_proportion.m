burstinterval=600;
fl=dir(fullfile("bzdata","ChainedLoop"+burstinterval+"S*.mat"));
su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
wrs_mux_meta=ephys.get_wrs_mux_meta();
sums=[];
for fi=1:numel(fl)
    sessid=str2double(regexp(fl(fi).name,'(?<=S)\d{1,3}(?=\.mat)','match','once'));
    sess_cid=su_meta.allcid(su_meta.sess==sessid);
    sess_wave_id=wrs_mux_meta.wave_id(su_meta.sess==sessid);
    
    fstr=load(fullfile(fl(fi).folder,fl(fi).name));
    trl=fstr.FT_SPIKE.trialinfo;
    
    for suidx=1:numel(sess_cid)
        ftidx=strcmp(num2str(sess_cid(suidx)),fstr.FT_SPIKE.label);
        if ~(any(fstr.FT_SPIKE.lc_tag{ftidx},'all'))
            continue;
        end
         
        switch sess_wave_id(suidx)
            case 0
                continue
            case 1
                trlsel=find(trl(:,5)==4 & trl(:,8)==3 & trl(:,9)>0 & trl(:,10)>0);
                spksel=ismember(fstr.FT_SPIKE.trial{ftidx},trlsel) ...
                    & fstr.FT_SPIKE.time{ftidx}>=1 & fstr.FT_SPIKE.time{ftidx}<4;
            case 2
                trlsel=find(trl(:,5)==4 & trl(:,8)==6 & trl(:,9)>0 & trl(:,10)>0);
                spksel=ismember(fstr.FT_SPIKE.trial{ftidx},trlsel) ...
                    & fstr.FT_SPIKE.time{ftidx}>=1 & fstr.FT_SPIKE.time{ftidx}<7;                
            case 3
                trlsel=find(trl(:,5)==8 & trl(:,8)==3 & trl(:,9)>0 & trl(:,10)>0);
                spksel=ismember(fstr.FT_SPIKE.trial{ftidx},trlsel) ...
                    & fstr.FT_SPIKE.time{ftidx}>=1 & fstr.FT_SPIKE.time{ftidx}<4;
            case 4
                trlsel=find(trl(:,5)==8 & trl(:,8)==6 & trl(:,9)>0 & trl(:,10)>0);
                spksel=ismember(fstr.FT_SPIKE.trial{ftidx},trlsel) ...
                    & fstr.FT_SPIKE.time{ftidx}>=1 & fstr.FT_SPIKE.time{ftidx}<7;
            case 5
                trl3sel=find(trl(:,5)==4 & trl(:,8)==3 & trl(:,9)>0 & trl(:,10)>0);
                trl6sel=find(trl(:,5)==4 & trl(:,8)==6 & trl(:,9)>0 & trl(:,10)>0);
                spksel=(ismember(fstr.FT_SPIKE.trial{ftidx},trl3sel) ...
                    & fstr.FT_SPIKE.time{ftidx}>=1 & fstr.FT_SPIKE.time{ftidx}<4) ...
                    |(ismember(fstr.FT_SPIKE.trial{ftidx},trl6sel) ...
                    & fstr.FT_SPIKE.time{ftidx}>=1 & fstr.FT_SPIKE.time{ftidx}<7);
            case 6
                trl3sel=find(trl(:,5)==8 & trl(:,8)==3 & trl(:,9)>0 & trl(:,10)>0);
                trl6sel=find(trl(:,5)==8 & trl(:,8)==6 & trl(:,9)>0 & trl(:,10)>0);
                spksel=(ismember(fstr.FT_SPIKE.trial{ftidx},trl3sel) ...
                    & fstr.FT_SPIKE.time{ftidx}>=1 & fstr.FT_SPIKE.time{ftidx}<4) ...
                    |(ismember(fstr.FT_SPIKE.trial{ftidx},trl6sel) ...
                    & fstr.FT_SPIKE.time{ftidx}>=1 & fstr.FT_SPIKE.time{ftidx}<7);
            case 7
                trlsel=find(ismember(trl(:,5),[4,8]) & trl(:,8)==3 & trl(:,9)>0 & trl(:,10)>0);
                spksel=ismember(fstr.FT_SPIKE.trial{ftidx},trlsel) ...
                    & fstr.FT_SPIKE.time{ftidx}>=1 & fstr.FT_SPIKE.time{ftidx}<4;
            case 8
                trlsel=find(ismember(trl(:,5),[4,8]) & trl(:,8)==6 & trl(:,9)>0 & trl(:,10)>0);
                spksel=ismember(fstr.FT_SPIKE.trial{ftidx},trlsel) ...
                    & fstr.FT_SPIKE.time{ftidx}>=1 & fstr.FT_SPIKE.time{ftidx}<7;
        end
        % ratio: chain:1, loop:4, SSCL1|4, burst>0
        sums=[sums;...
            sessid,double(sess_cid(suidx)),sess_wave_id(suidx),...
            nnz(bitand(fstr.FT_SPIKE.lc_tag{ftidx}(spksel),1)),...
            nnz(bitand(fstr.FT_SPIKE.lc_tag{ftidx}(spksel),4)),...
            nnz(bitand(fstr.FT_SPIKE.lc_tag{ftidx}(spksel),5)),...
            nnz(fstr.FT_SPIKE.lc_tag{ftidx}(spksel)),...
            nnz(spksel)];
    end
end

% two panels olf|both, 4 box each, chain, loop, SSChL, BChL

olfsel=ismember(sums(:,3),5:6);
bothsel=ismember(sums(:,3),1:4);

olfboxdata=[reshape(sums(olfsel,4:6)./sums(olfsel,8),[],1),...
    reshape(repmat([1 2 3],nnz(olfsel),1),[],1)];
bothboxdata=[reshape(sums(bothsel,4:6)./sums(bothsel,8),[],1),...
    reshape(repmat([1 2 3],nnz(bothsel),1),[],1)];


figure()
tiledlayout(1,2)
nexttile()
hold on
boxplot(olfboxdata(:,1),olfboxdata(:,2),'Colors','k','Whisker',realmax)
set(gca(),'XTick',1:3,'XTickLabel',{'Ch','Lp','ChL'},'YTick',0:0.1:0.5,'YTickLabel',0:10:50)
ylim([0,0.5])
ylabel('Percentage of spikes during delay')
title('Odor only')
nexttile()
hold on
boxplot(bothboxdata(:,1),bothboxdata(:,2),'Colors','k','Whisker',realmax)
set(gca(),'XTick',1:3,'XTickLabel',{'Ch','Lp','ChL'},'YTick',0:0.1:0.5,'YTickLabel',0:10:50)
ylim([0,0.5])
ylabel('Percentage of spikes during delay')
title('Encode both')



olfboxdatab=[reshape(sums(olfsel,7)./sums(olfsel,8),[],1),...
    reshape(repmat(1,nnz(olfsel),1),[],1)];
bothboxdatab=[reshape(sums(bothsel,7)./sums(bothsel,8),[],1),...
    reshape(repmat(2,nnz(bothsel),1),[],1)];

bdata=[olfboxdatab;bothboxdatab];


figure()
hold on
boxplot(bdata(:,1),bdata(:,2),'Colors','k','Whisker',realmax)
set(gca(),'XTick',1:2,'XTickLabel',{'Odor','Both'},'YTick',0:0.2:1,'YTickLabel',0:20:100)
ylabel('Percentage of spikes during delay')
xlim([0.5,2.5])
title('burst motif spike proportion')

prctile(bdata(bdata(:,2)==1,1),[25 50 75])
prctile(bdata(bdata(:,2)==2,1),[25 50 75])
return

