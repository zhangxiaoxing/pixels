%TODO Waveform based SU filter

function su_meta_=load_meta(opt)
arguments
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO','MY'})}='neupix'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
    % opt.n_bin (1,1) double {mustBeInteger,mustBePositive} = 3
    opt.skip_stats (1,1) logical = true
    opt.adjust_white_matter (1,1) logical = true
    opt.save_file (1,1) logical = false
    opt.load_file (1,1) logical = false
    opt.filename (1,:) char = 'su_meta.mat'
    opt.more_folder (1,1) logical = false
end
assert(opt.skip_stats,"stats have been relocated")
persistent su_meta opt_
if isempty(su_meta) || ~isequaln(opt,opt_)
    global gather_config
    if ~isempty(gather_config)
        opt.adjust_white_matter=gather_config.adjust_white_matter;
    end
    if opt.load_file
        load(fullfile('binary',opt.filename),'su_meta');
    else
        [homedir,homedir_all]=ephys.util.getHomedir();
        if opt.more_folder
            fl=[];
            for ii=1:numel(homedir_all)
                fl=[fl;struct2table(dir(fullfile(homedir_all{ii},'**','FR_All*.hdf5')))];
            end
            fl=fl(~contains(fl.name,'250'),:);
        else
            fl=dir(fullfile(homedir,'**','FR_All_1000.hdf5'));
        end
        % [~,sidx]=sort({fl.folder}); sorted by default
        su_meta=cell2struct({cell(0);[];{};[]},{'allpath','allcid','reg_tree','sess'});
        wtsessidx=1;
        for fi=1:numel(fl)
            suids=h5read(fullfile(fl(fi).folder,fl(fi).name),'/SU_id');
            if isempty(suids)
                keyboard
            end
            trials=h5read(fullfile(fl(fi).folder,fl(fi).name),'/Trials');
            if strcmp(opt.criteria,'Learning')
                ltrials=behav.procPerf(trials,"criteria","Learning");
            end
            if (strcmp(opt.criteria,'WT') && sum(trials(:,9))<40)...
                    || (strcmp(opt.criteria,'Learning') && (sum(trials(:,9))>=40 || sum(ltrials(:,9))<40))
                continue
            end
            if exist(fullfile(fl(fi).folder,'su_id2reg.csv'),'file')~=2
                continue
            end
            do=detectImportOptions(fullfile(fl(fi).folder,'su_id2reg.csv'));
            do=setvartype(do,do.VariableNames,{'double','double','char','double','char','char','char','char','char','char'});
            regtbl=readtable(fullfile(fl(fi).folder,'su_id2reg.csv'),do);
            regsel=ismember(regtbl.index,suids);
            su_meta.allcid=[su_meta.allcid;suids];
            su_meta.sess=[su_meta.sess;repmat(wtsessidx,numel(suids),1)];
            wtsessidx=wtsessidx+1;
            su_meta.reg_tree=[su_meta.reg_tree,table2cell(regtbl(regsel,5:10)).'];
            pathsuffix=fl(fi).folder;
            su_meta.allpath=[su_meta.allpath;repmat({pathsuffix},numel(suids),1)];
        end
        su_meta.allcid=uint16(su_meta.allcid);
        su_meta.sess=int32(su_meta.sess);

        if opt.adjust_white_matter
            % if strcmp(opt.criteria,'WT')
            %     su_meta.reg_tree=ephys.get_adjusted_reg_tree('adjust_white_matter',opt.adjust_white_matter);
            % else
                warning("Missing adjusted file, fall back to previous alignment")
            % end
        end

        if opt.save_file
            blame=vcs.blame();
            save(fullfile('binary',opt.filename),'su_meta','blame');
        end
    end
    
end
opt_=opt;
su_meta_=su_meta;
end