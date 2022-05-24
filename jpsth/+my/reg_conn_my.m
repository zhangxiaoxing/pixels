function reg_conn_my(opt)
arguments
    opt.type (1,:) char {mustBeMember(opt.type,{'MYWT'})}='MYWT'
    opt.data (:,1) struct =[]
    opt.prefix (1,:) char = '0516'
    opt.bin_range (1,2) double = [-2 7]
end

%TODO merge with load_sig_pair script
if isempty(opt.data)
    if strcmp(opt.type,'MYWT')
        load(sprintf('sums_conn_my_%d_%d.mat',opt.bin_range(1),opt.bin_range(2)),'sums_conn_str')
    end
else
    sums_conn_str=opt.data;
end

for fidx=1:length(sums_conn_str)
    tic
    disp(fidx);
    fpath=sums_conn_str(fidx).folder; %session data folder
    if strcmp(opt.type,'MYWT')
        pc_stem=replace(regexp(fpath,'(?<=SPKINFO/).*$','match','once'),'/','\');
%         inputf=fullfile('K:','neupix','SPKINFO',pc_stem,'FR_All_1000.hdf5');
%         all_su=int32(h5read(inputf,'/SU_id'));
    end
    pair_comb_one_dir=[];
    sig_meta=struct();
    if ~isempty(sums_conn_str(fidx).sums_conn_s1)
        sig_con_s1=int32(sums_conn_str(fidx).sums_conn_s1); % significant functional coupling
        [sig_meta.s1,~]=bz.util.get_meta(sig_con_s1,pair_comb_one_dir,pc_stem,'type',opt.type); % assign meta info
    end
    if ~isempty(sums_conn_str(fidx).sums_conn_s2)   
        sig_con_s2=int32(sums_conn_str(fidx).sums_conn_s2); % significant functional coupling
        [sig_meta.s2,~]=bz.util.get_meta(sig_con_s2,pair_comb_one_dir,pc_stem,'type',opt.type); % assign meta info
    end
    blame=vcs.blame();
    save(fullfile('mydata',sprintf('%s_conn_w_reg_my_%d_%d_%d.mat',opt.prefix,fidx,opt.bin_range(1),opt.bin_range(2))),'sig_meta','pc_stem','blame','-v7','-nocompression')
    toc
end
end