function reg_conn_bz(opt)
arguments
%     opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO','MYWT'})}='neupix'
    opt.data (:,1) struct =[]
    opt.prefix (1,:) char = 'BZWT'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
    opt.inhibit (1,1) logical = false;
end

%TODO merge with load_sig_pair script
if isempty(opt.data)
    if ~strcmp(opt.criteria,'Learning')
        if opt.inhibit
            load('sums_conn_inhibit.mat','sums_conn_str')
        else
            load('sums_conn.mat','sums_conn_str')
        end
    else
        load('sums_conn_learning.mat','sums_conn_str')
    end

else
    sums_conn_str=opt.data;
end
fprintf('Total sessions %d\n',length(sums_conn_str));
for fidx=1:length(sums_conn_str)
    tic
    disp(fidx);
    fpath=sums_conn_str(fidx).folder; %session data folder

    if contains(fpath,'SPKINFO')
        pc_stem=replace(regexp(fpath,'(?<=SPKINFO/).*$','match','once'),'/',filesep());
    else
        pc_stem=replace(fpath,'/',filesep());
    end
    if ispc
        inputf=fullfile('K:','neupix','SPKINFO',pc_stem,'FR_All_1000.hdf5');
    elseif isunix
        inputf=fullfile('/home/zx/neupix/SPKINFO',pc_stem,'FR_All_1000.hdf5');
    end
    all_su=int32(h5read(inputf,'/SU_id'));


    sig_con=int32(sums_conn_str(fidx).sig_con); % significant functional coupling
    if numel(sig_con)==2, sig_con=reshape(sig_con,1,2);end
    pair_comb_one_dir=nchoosek(all_su,2); % all pairs combination
    [sig_meta,pair_meta]=bz.util.get_meta(sig_con,pair_comb_one_dir,pc_stem,'criteria',opt.criteria); % assign meta info
    
    %mirror unidirection pair data
    fields={'suid','reg','wf_good','waveid'};
    for fi=fields
        %TODO online genenrate session tag
%         sig.(fi{1})=cat(1,sig.(fi{1}),sig_meta.(fi{1}));
        pair_meta.(fi{1})=cat(1,pair_meta.(fi{1}),flip(pair_meta.(fi{1}),ndims(pair_meta.(fi{1}))));%uni-dir to bi-dir
%         pair.(fi{1})=cat(1,pair.(fi{1}),pair_meta.(fi{1}));
    end
    if opt.inhibit
        save(fullfile('bzdata',sprintf('%s_conn_w_reg_%d_inhibitory.mat',opt.prefix,fidx)),'sig_meta','pair_meta','pc_stem','-v7','-nocompression')
    else
        save(fullfile('bzdata',sprintf('%s_conn_w_reg_%d.mat',opt.prefix,fidx)),'sig_meta','pair_meta','pc_stem','-v7','-nocompression')
%         structWriteHDF5('fpath',fullfile('bzdata',sprintf('%s_conn_w_reg_%d.hdf5',opt.prefix,fidx)),...
%             'datasets',struct('sig_meta',sig_meta,'pair_meta',pair_meta),...
%             'pc_stem',pc_stem)
    end
    toc
end
end


function structWriteHDF5(opt)
arguments
    opt.fpath (1,:) char {mustBeNonempty}
    opt.datasets (1,1) struct {mustBeNonempty}
    opt.pc_stem (1,:) char {mustBeNonempty}
end
    if exist(opt.fpath,'file')
        delete(opt.fpath)
    end
    h5create(opt.fpath,'/pc_stem',1,'Datatype','string')
    h5write(opt.fpath,'/pc_stem',string(opt.pc_stem))
    
    dataset_name=fieldnames(opt.datasets);
    for ii=1:numel(dataset_name)
        field_name=fieldnames(opt.datasets.(dataset_name{ii}));
        for jj=1:numel(field_name)
            h5writeOne(opt.fpath,opt.datasets,dataset_name{ii},field_name{jj});
        end
    end
end

function h5writeOne(fpath,datasets,dataset_name,field_name)
h5create(fpath,...
    sprintf('/%s/%s',dataset_name,field_name),...
    size(datasets.(dataset_name).(field_name)),...
    'Datatype',class(datasets.(dataset_name).(field_name)))
h5write(fpath,...
    sprintf('/%s/%s',dataset_name,field_name),...
    datasets.(dataset_name).(field_name))
end

