function reg_conn_my_combined(prefix, opt)
arguments
    prefix (1,:) char
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix'})}='neupix'
    opt.data (:,1) struct =[]
    opt.bin_range (1,2) double = [1 2]
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
end

%TODO merge with load_sig_pair script
if isempty(opt.data)
    if strcmp(opt.criteria,'WT')
        load(sprintf('sums_conn_my_%d_%d_combined.mat',opt.bin_range(1),opt.bin_range(2)),'sums_conn_str')
    elseif strcmp(opt.criteria,'Learning')
        load(sprintf('sums_conn_my_Learning_%d_%d_combined.mat',opt.bin_range(1),opt.bin_range(2)),'sums_conn_str')
    end
else
    sums_conn_str=opt.data;
end

for fidx=1:length(sums_conn_str)
    tic
    disp(fidx);
    if isempty(sums_conn_str(fidx).sig_con)
        disp('No significant connection');
        continue
    end
    if contains(sums_conn_str(fidx).folder,'SPKINFO')
        pc_stem=replace(regexp(sums_conn_str(fidx).folder,'(?<=SPKINFO/).*$','match','once'),'/','\');
    else
        pc_stem=replace(sums_conn_str(fidx).folder,'/','\');
    end
    sig_con=int32(sums_conn_str(fidx).sig_con); % significant functional coupling
    [sig_meta,~]=bz.util.get_meta(sig_con,[],pc_stem,'type',opt.type,'criteria',opt.criteria); % assign meta info
    pair_meta=[];
    blame=vcs.blame();
    if strcmp(opt.criteria,'WT')
        save(fullfile('mydata',sprintf('%s_conn_w_reg_my_%d.mat',prefix,fidx)),'sig_meta','pair_meta','pc_stem','blame','-v7','-nocompression')
    elseif strcmp(opt.criteria,'Learning')
        save(fullfile('mydata',sprintf('%s_conn_w_reg_my_learning_%d.mat',prefix,fidx)),'sig_meta','pair_meta','pc_stem','blame','-v7','-nocompression')
    end
    toc
end
end