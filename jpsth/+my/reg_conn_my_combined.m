function reg_conn_my_combined(opt)
arguments
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix'})}='neupix'
    opt.data (:,1) struct =[]
    opt.prefix (1,:) char = '0517'
    opt.bin_range (1,2) double = [1 2]
end

%TODO merge with load_sig_pair script
if isempty(opt.data)
    if strcmp(opt.type,'neupix')
        load('sums_conn_my_combined.mat','sums_conn_str')
    end
else
    sums_conn_str=opt.data;
end

for fidx=11:length(sums_conn_str)
    tic
    disp(fidx);
    if isempty(sums_conn_str(fidx).sig_con)
        disp('No significant connection');
        continue
    end
    pc_stem=replace(regexp(sums_conn_str(fidx).folder,'(?<=SPKINFO/).*$','match','once'),'/','\');
    sig_con=int32(sums_conn_str(fidx).sig_con); % significant functional coupling
    [sig_meta,~]=bz.util.get_meta(sig_con,[],pc_stem,'type',opt.type); % assign meta info
    pair_meta=[];
    blame=vcs.blame();
    save(fullfile('bzdata',sprintf('%s_conn_w_reg_%d.mat',opt.prefix,fidx)),'sig_meta','pair_meta','pc_stem','blame','-v7','-nocompression')
    toc
end
end