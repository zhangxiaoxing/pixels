%deprecated as of 2021-may, only for historical reference


% %  A peak at a negative lag for stat.xcorr(chan1,chan2,:) means that chan1 is leading
% %  chan2. Thus, a negative lag represents a spike in the second dimension of
% %  stat.xcorr before the channel in the third dimension of stat.stat.
% if ~exist('bin_range','var')
%     disp('missing bin range');
%     return
% end
% % cd('~/pixels/jpsth')
% homedir='/home/zx/neupix/wyt';
% %currmodel='selec';
% %prefix='0831';
% delay=6;
% %bin_range=[4,5];
% % addpath(fullfile('npy-matlab-master','npy-matlab'))
% % addpath('fieldtrip-20200320')
% % ft_defaults
% sus_trans=h5read('../transient_6.hdf5','/sus_trans'); %export_arr = np.vstack((sust, transient, switched, unclassified, early_in_6s, late_in_6s, prefer_s))
% reg_list=h5read('../transient_6.hdf5','/reg');
% cid_list=h5read('../transient_6.hdf5','/cluster_id');
% path_list=h5read('../transient_6.hdf5','/path');
% 
% load('reg_keep.mat','reg_set');
% reg_keep=reg_set([1:112,114,115]);
% reg_sel=ismember(deblank(reg_list),reg_keep);
% if startsWith(currmodel,'memory')
%     sust=find(sus_trans(:,1) & reg_sel);
%     trans=find(any(sus_trans(:,[2,4]),2) & reg_sel);
%     supool=[sust;trans]';
% elseif startsWith(currmodel,'selec')
%     sust=find(sus_trans(:,1) & reg_sel);
%     trans=find(sus_trans(:,2) & reg_sel);
%     supool=[sust;trans]';
% elseif startsWith(currmodel,'nonsel')
%     sust=[];
%     trans=[];
%     nonsel_logic=~any(sus_trans(:,1:4),2);
%     supool=find(nonsel_logic & reg_sel);
%     supool=supool';
% end
% 
% 
% counter=[];
% done=[];
% 
% error_list=cell(0);
% 
% for i=1:length(supool)
%     disp(i)
%     skip=false;
%     avail=false;
%     if exist(sprintf('%s_%s_XCORR_duo_f%d_delay_%d_%d_%d_2msbin.mat',prefix,currmodel,i,delay,bin_range(1),bin_range(2)),'file')
%         disp('file exist');
%         fprintf('%s_%s_XCORR_duo_f%d_delay_%d_%d_%d_2msbin.mat\n',prefix,currmodel,i,delay,bin_range(1),bin_range(2));
%         skip=true;
%     end
% 
%     if ismember(supool(i),done)
%         continue
%     end
%     if exist('xcorrpause','file')
%         disp('paused by file')
%         disp(i)
%         keyboard
%     end
% 
%     folder=regexp(path_list{supool(i)},'(\w|\\|-)*','match','once');
%     [folderType,file,spkFolder,metaFolder,error_list]=jointFolder(folder,error_list);
%     if folderType<0
%         continue
%     end
%     wffile=fullfile(metaFolder,'wf_stats.hdf5');% posix
%     if isfile(wffile)
%         if startsWith(currmodel,'selec')
%             if folderType==1
%                 sustIds=cid_list(startsWith(path_list,folder) & sus_trans(:,1) & reg_sel);
%                 transIds=cid_list(startsWith(path_list,folder) & sus_trans(:,2) & reg_sel);
%                 sameFolder=find(startsWith(path_list,folder) & any(sus_trans(:,[1,2]),2) & reg_sel);
%                 done=[done;sameFolder];
%                 sustCount=numel(sustIds);
%                 transCount=numel(transIds);
%             elseif folderType==2
%                 fimec0=replace(folder,'imec1','imec0');
%                 fimec1=replace(folder,'imec0','imec1');
% 
%                 sustIds0=cid_list(startsWith(path_list,fimec0) & sus_trans(:,1) & reg_sel);
%                 transIds0=cid_list(startsWith(path_list,fimec0) & sus_trans(:,2) & reg_sel);
% 
%                 sustIds1=cid_list(startsWith(path_list,fimec1) & sus_trans(:,1) & reg_sel)+10000;
%                 transIds1=cid_list(startsWith(path_list,fimec1) & sus_trans(:,2) & reg_sel)+10000;
% 
%                 sameFolder=find((startsWith(path_list,fimec0)|startsWith(path_list,fimec1)) & any(sus_trans(:,[1,2]),2) & reg_sel);
%                 done=[done;sameFolder];
% 
%                 sustIds=[sustIds0(:);sustIds1(:)];
%                 transIds=[transIds0(:);transIds1(:)];
% 
%                 sustCount=numel(sustIds);
%                 transCount=numel(transIds);            
% 
%             end
% 
%             if transCount<1
%                 continue
%             end
%             nonselIds=[];
%             if ~skip
%                 [avail,spktrial]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,currmodel); % posix
%             end
%         elseif startsWith(currmodel,'memory')
%             if folderType==1
%                 sustIds=cid_list(startsWith(path_list,folder) & sus_trans(:,1) & reg_sel);
%                 transIds=cid_list(startsWith(path_list,folder) & any(sus_trans(:,[2,4]),2) & reg_sel);
%                 sameFolder=find(startsWith(path_list,folder) & any(sus_trans(:,[1,2,4]),2) & reg_sel);
%                 done=[done;sameFolder];
%                 sustCount=numel(sustIds);
%                 transCount=numel(transIds);
%             elseif folderType==2
%                 fimec0=replace(folder,'imec1','imec0');
%                 fimec1=replace(folder,'imec0','imec1');
% 
%                 sustIds0=cid_list(startsWith(path_list,fimec0) & sus_trans(:,1) & reg_sel);
%                 transIds0=cid_list(startsWith(path_list,fimec0) & any(sus_trans(:,[2,4]),2) & reg_sel);
% 
%                 sustIds1=cid_list(startsWith(path_list,fimec1) & sus_trans(:,1) & reg_sel)+10000;
%                 transIds1=cid_list(startsWith(path_list,fimec1) & any(sus_trans(:,[2,4]),2) & reg_sel)+10000;
% 
%                 sameFolder=find((startsWith(path_list,fimec0)|startsWith(path_list,fimec1)) & any(sus_trans(:,[1,2,4]),2) & reg_sel);
%                 done=[done;sameFolder];
% 
%                 sustIds=[sustIds0(:);sustIds1(:)];
%                 transIds=[transIds0(:);transIds1(:)];
% 
%                 sustCount=numel(sustIds);
%                 transCount=numel(transIds);            
% 
%             end
% 
%             if sustCount+transCount<2
%                 continue
%             end
%             nonselIds=[];
%             if ~skip
%                 [avail,spktrial]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,currmodel); % posix
%             end
%         elseif startsWith(currmodel,'nonsel')
%             if folderType==1
%                 sameFolder=find(startsWith(path_list,folder) & nonsel_logic & reg_sel);
%                 nonselIds=cid_list(startsWith(path_list,folder) & nonsel_logic & reg_sel);
%                 done=[done;sameFolder];
%             elseif folderType==2
%                 fimec0=replace(folder,'imec1','imec0');
%                 fimec1=replace(folder,'imec0','imec1');
% 
%                 nonselIds0=cid_list(startsWith(path_list,fimec0) & nonsel_logic & reg_sel);
% 
%                 nonselIds1=cid_list(startsWith(path_list,fimec1) & nonsel_logic & reg_sel)+10000;
% 
% 
%                 sameFolder=find((startsWith(path_list,fimec0)|startsWith(path_list,fimec1)) & nonsel_logic & reg_sel);
%                 done=[done;sameFolder];
% 
%                 nonselIds=[nonselIds0(:);nonselIds1(:)];
%                 sustIds=[];
%                 transIds=[];
% 
%             end
% 
%             if numel(nonselIds)<1
%                 continue
%             end
%             if ~skip
%                 [avail,spktrial]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,currmodel); % posix
%             end
%         
%         end
%         if avail
%             [xc_s1,xcshuf_s1,xc_s2,xcshuf_x2]=plotxcorr(spktrial,delay,bin_range)
%         else
%             continue
%         end
% 
%     else
%         continue
%     end
%     sums={i,folder,sustIds,transIds,xc_s1,xcshuf_s1,xc_s2,xcshuf_x2}; %per folder save
%     %return
%     save(sprintf('%s_%s_XCORR_duo_f%d_delay_%d_%d_%d_2msbin.mat',prefix,currmodel,i,delay,bin_range(1),bin_range(2)),'sums','-v7.3','sust','trans','supool','counter','done') %prefix
%  	fprintf('%d of %d\n',i,length(supool))
% end
% 
% quit(0) 
% 
% function [folderType,file,spkFolder,metaFolder,error_list]=jointFolder(folder,error_list)
%     metaFolder=replace(folder,'\','/');
%     metaFolder=fullfile('/home/zx/neupix/wyt/DataSum',metaFolder);
%     if isfolder(metaFolder)
%         spkFolder=replace(metaFolder,'imec1','imec0');
%         file=dir(fullfile(spkFolder,'spike_info.mat'));
%         if isempty(file)
%             folderType=-1;
%             file=[];
%             spkFolder=[];
%             disp('Error processing file 2-tracks');
%             disp(metaFolder);
%             error_list(end+1,:)={folderType,metaFolder};
% %             pause;
%             return
%         end
%         folderType=2;
%     else
%         metaFolder=replace(metaFolder,'DataSum','DataSum/singleProbe');
%         spkFolder=metaFolder;
%         file=dir(fullfile(spkFolder,'spike_times.npy'));
%         if isempty(file)
%             folderType=-1;
%             file=[];
%             spkFolder=[];
%             disp('Error processing file 1-track');
%             disp(metaFolder);
%             error_list(end+1,:)={folderType,metaFolder};
%             return
%         end
%         folderType=1;        
%     end
% end
% 
% 
% function [avail,out]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,model)
% sps=30000;
% trials=clearBadPerf(h5read(fullfile(metaFolder,'events.hdf5'),'/trials')',model);
% if isempty(trials)
%     avail=false;
%     out=[];
%     return
% end
% 
% %     trials=double(trials);
% %     info=[trials(:,1)/s1s,trials(:,2)/s1s,trials(:,5),trials(:,6),trials(:,7),trials(:,8)];
% if strcmp(model, 'full')
%     s1s=30000;
%     FR_Th=1.0;
% 
%     metaf=strtrim(ls(fullfile(metaFolder,'*.meta')));
%     fh=fopen(metaf);
%     ts=textscan(fh,'%s','Delimiter',{'\n'});
%     nSample=str2double(replace(ts{1}{startsWith(ts{1},'fileSizeBytes')},'fileSizeBytes=',''));
%     spkNThresh=nSample/385/s1s/2*FR_Th;
%     clusterInfo = readtable(fullfile(metaFolder,'cluster_info.tsv'),'FileType','text','Delimiter','tab');
%     waveformGood=strcmp(clusterInfo{:,4},'good');
%     freqGood=clusterInfo{:,10}>spkNThresh;
%     cluster_ids = table2array(clusterInfo(waveformGood & freqGood,1));
% elseif startsWith(model,'selec') || startsWith(model,'memory')
%     cluster_ids=[sustIds;transIds];
% elseif startsWith(model,'nonsel')
%     cluster_ids=nonselIds;
% else
%     keyboard
% end
% 
% %  single-unit candidate
% 
% if folderType==1
%     spkTS=readNPY(fullfile(spkFolder,'spike_times.npy'));
%     spkId=readNPY(fullfile(spkFolder,'spike_clusters.npy'));
% elseif folderType==2
%     fstr=load(fullfile(spkFolder,'spike_info.mat'));
%     spkId=double([fstr.spike_info{1}{1};fstr.spike_info{1}{2}]);
%     spkTS=double([fstr.spike_info{2}{1};fstr.spike_info{2}{2}]);
% end
% FT_SPIKE=struct();
% 
% FT_SPIKE.label=strtrim(cellstr(num2str(cluster_ids)));
% FT_SPIKE.timestamp=cell(1,numel(cluster_ids));
% for i=1:numel(cluster_ids)
%     FT_SPIKE.timestamp{i}=spkTS(spkId==cluster_ids(i))';
% end
% %  continuous format F T struct file
% cfg=struct();
% cfg.trl=[trials(:,1)-3*sps,trials(:,1)+11*sps,zeros(size(trials,1),1)-3*sps,trials];
% cfg.trlunit='timestamps';
% cfg.timestampspersecond=sps;
% 
% FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);
% 
% out=FT_SPIKE;
% avail=true;
% end
% 
% 
% 
% function out=clearBadPerf(facSeq, model)
% if strcmp(model, 'error')
%     if length(facSeq)>=40
%         errorsel=~xor(facSeq(:,5)==facSeq(:,6) , facSeq(:,7)>0);
%         out=facSeq(errorsel,:);
%     else
%         out=[];
%     end
% else
%     if length(facSeq)>=40
%         facSeq(:,9)=0;
%         i=40;
%         while i<=length(facSeq)
%             goodOff=nnz(xor(facSeq(i-39:i,5)==facSeq(i-39:i,6) , facSeq(i-39:i,7)>0));
%             if goodOff>=30 %.75 correct rate
%                 facSeq(i-39:i,9)=1;
%             end
%             i=i+1;
%         end
%         out=facSeq(facSeq(:,9)==1,:);
%     else
%         out=[];
%     end
% end
% end
% 
% 
% function [Xc_S1,Xshuff_S1,Xc_S2,Xshuff_S2]=plotxcorr(spikeTrials,delay,bin_range)
% % https://www.nature.com/articles/nn799
% % A role for inhibition in shaping the temporal flow of information in prefrontal cortex
% % Christos Constantinidis, Graham V. Williams & Patricia S. Goldman-Rakic 
% % Nature Neuroscience volume 5, pages175-180(2002)
% % 
% % Neuron, Volume 76
% % Functional Microcircuit Recruited during Retrieval of Object Association Memory in Monkey Perirhinal Cortex
% % Toshiyuki Hirabayashi, Daigo Takeuchi, Keita Tamura, and Yasushi Miyashita
% 
% 
% cfg             = [];
% cfg.maxlag      = 0.1; % maximum 100 ms
% cfg.binsize     = 0.002; % bins of 2 ms
% cfg.outputunit  = 'raw'; % make unit area
% cfg.latency     = bin_range; % time bin based on sample onset
% cfg.vartriallen = 'no'; % allow variable trial lengths
% cfg.debias      = 'no';
% 
% cfg.trials      = find(spikeTrials.trialinfo(:,5)==4 & spikeTrials.trialinfo(:,8)==delay);
% if numel(cfg.trials)<2
%     Xc_S1=[];
%     Xshuff_S1=[];
% else
%     cfg.method      = 'xcorr'; % compute the normal cross-correlogram
%     Xc_S1 = ft_spike_xcorr(cfg,spikeTrials);
%     cfg.method      = 'shiftpredictor'; % compute the shift predictor
%     Xshuff_S1 = ft_spike_xcorr(cfg,spikeTrials);
% end
% 
% cfg.trials      = find(spikeTrials.trialinfo(:,5)==8 & spikeTrials.trialinfo(:,8)==delay);
% if numel(cfg.trials)<2
%     Xc_S2=[];
%     Xshuff_S2=[];
% else
%     cfg.method      = 'xcorr'; % compute the normal cross-correlogram
%     Xc_S2 = ft_spike_xcorr(cfg,spikeTrials);
%     cfg.method      = 'shiftpredictor'; % compute the shift predictor
%     Xshuff_S2 = ft_spike_xcorr(cfg,spikeTrials);
% end
% % compute the shuffled correlogram
% end
% 
% 
% 
