burstinterval=600;
decmod=categorical({'olf'},{'olf','dur'});

trlN=20;
neuN=50;

poswithmat=[];
posw_o_mat=[];
negmat=[];


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
                if  decmod=='olf'
                    antitrlsel=find(trl(:,5)==8 & trl(:,8)==3 & trl(:,9)>0 & trl(:,10)>0);
                else
                    antitrlsel=find(trl(:,5)==4 & trl(:,8)==6 & trl(:,9)>0 & trl(:,10)>0);
                end
                
            case 2
                trlsel=find(trl(:,5)==4 & trl(:,8)==6 & trl(:,9)>0 & trl(:,10)>0);
                if  decmod=='olf'
                    antitrlsel=find(trl(:,5)==8 & trl(:,8)==6 & trl(:,9)>0 & trl(:,10)>0);
                else
                    antitrlsel=find(trl(:,5)==4 & trl(:,8)==3 & trl(:,9)>0 & trl(:,10)>0);
                end

            case 3
                trlsel=find(trl(:,5)==8 & trl(:,8)==3 & trl(:,9)>0 & trl(:,10)>0);
                if  decmod=='olf'
                    antitrlsel=find(trl(:,5)==4 & trl(:,8)==3 & trl(:,9)>0 & trl(:,10)>0);
                else
                    antitrlsel=find(trl(:,5)==8 & trl(:,8)==6 & trl(:,9)>0 & trl(:,10)>0);
                end

            case 4
                trlsel=find(trl(:,5)==8 & trl(:,8)==6 & trl(:,9)>0 & trl(:,10)>0);
                if  decmod=='olf'
                    antitrlsel=find(trl(:,5)==4 & trl(:,8)==6 & trl(:,9)>0 & trl(:,10)>0);
                else
                    antitrlsel=find(trl(:,5)==8 & trl(:,8)==3 & trl(:,9)>0 & trl(:,10)>0);
                end
            case 5
                trl3sel=find(trl(:,5)==4 & trl(:,8)==3 & trl(:,9)>0 & trl(:,10)>0);
                trl6sel=find(trl(:,5)==4 & trl(:,8)==6 & trl(:,9)>0 & trl(:,10)>0);
                antitrl3sel=find(trl(:,5)==8 & trl(:,8)==3 & trl(:,9)>0 & trl(:,10)>0);
                antitrl6sel=find(trl(:,5)==8 & trl(:,8)==6 & trl(:,9)>0 & trl(:,10)>0);
            case 6
                trl3sel=find(trl(:,5)==8 & trl(:,8)==3 & trl(:,9)>0 & trl(:,10)>0);
                trl6sel=find(trl(:,5)==8 & trl(:,8)==6 & trl(:,9)>0 & trl(:,10)>0);
                antitrl3sel=find(trl(:,5)==8 & trl(:,8)==3 & trl(:,9)>0 & trl(:,10)>0);
                antitrl6sel=find(trl(:,5)==8 & trl(:,8)==6 & trl(:,9)>0 & trl(:,10)>0);
            case 7
                trlsel=find(ismember(trl(:,5),[4,8]) & trl(:,8)==3 & trl(:,9)>0 & trl(:,10)>0);
                antitrlsel=find(ismember(trl(:,5),[4,8]) & trl(:,8)==6 & trl(:,9)>0 & trl(:,10)>0);
            case 8
                trlsel=find(ismember(trl(:,5),[4,8]) & trl(:,8)==6 & trl(:,9)>0 & trl(:,10)>0);
                antitrlsel=find(ismember(trl(:,5),[4,8]) & trl(:,8)==3 & trl(:,9)>0 & trl(:,10)>0);
        end
        
        if all(decmod=='olf')
            if ismember(sess_wave_id(suidx),1:4)
                if ~(numel(trlsel)>trlN && numel(antitrlsel)>trlN)
                    continue
                end
                postrl=randsample(trlsel,trlN);

                spksel=ismember(fstr.FT_SPIKE.trial{ftidx},trlsel) ...
                    & fstr.FT_SPIKE.time{ftidx}>=1 & fstr.FT_SPIKE.time{ftidx}<4;

                [gc,gr]=groupcounts(fstr.FT_SPIKE.trial{ftidx}(spksel).');


                poswithmat=[poswithmat,./trl(postrl,8)];
                posw_o_mat=[posw_o_mat;]; 
                negmat=[negmat;];
            end
        end
                % per preferred-trial spikes, before dealing with chain-loop tags

%                 spksel=(ismember(fstr.FT_SPIKE.trial{ftidx},trl3sel) ...
%                     & fstr.FT_SPIKE.time{ftidx}>=1 & fstr.FT_SPIKE.time{ftidx}<4) ...
%                     |(ismember(fstr.FT_SPIKE.trial{ftidx},trl6sel) ...
%                     & fstr.FT_SPIKE.time{ftidx}>=1 & fstr.FT_SPIKE.time{ftidx}<7);

% 
%                 antispksel=ismember(fstr.FT_SPIKE.trial{ftidx},antitrlsel) ...
%                     & fstr.FT_SPIKE.time{ftidx}>=1 & fstr.FT_SPIKE.time{ftidx}<4;
%                 spksel=ismember(fstr.FT_SPIKE.trial{ftidx},trlsel) ...
%                     & fstr.FT_SPIKE.time{ftidx}>=1 & fstr.FT_SPIKE.time{ftidx}<4;

%                 [gc,gr]=groupcounts(fstr.FT_SPIKE.trial{ftidx}(spksel).');

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



return

