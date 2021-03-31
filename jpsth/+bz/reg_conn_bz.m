%TODO merge with load_sig_pair script

load('sums_conn.mat','sums_conn_str')
% 
% sig=struct(); % for significant connection
% sig.suid=[]; % cluster id assigned by kilosort, 2nd+ probe prefixed by probe# 
% sig.reg=[]; % brain region tree
% sig.mem_type=[]; % Wilcoxon rank summation p value
% sig.sess=[];
% 
% pair=sig; % for all pairs

for fidx=1:length(sums_conn_str)
    disp(fidx);
    fpath=sums_conn_str(fidx).folder; %session data folder
    pc_stem=replace(regexp(fpath,'(?<=SPKINFO/).*$','match','once'),'/','\');
    inputf=fullfile('K:','neupix','SPKINFO',pc_stem,'FR_All_1000.hdf5');
    all_su=int32(h5read(inputf,'/SU_id'));
    sig_con=int32(sums_conn_str(fidx).sig_con); % significant functional coupling
    pair_comb_one_dir=nchoosek(all_su,2); % all pairs combination
    [sig_meta,pair_meta]=bz.util.get_meta(sig_con,pair_comb_one_dir,pc_stem); % assign meta info
    
    fields={'suid','reg','wrsp','selec','mem_type'};
    for fi=fields
        %TODO online genenrate session tag
%         sig.(fi{1})=cat(1,sig.(fi{1}),sig_meta.(fi{1}));
        pair_meta.(fi{1})=cat(1,pair_meta.(fi{1}),flip(pair_meta.(fi{1}),ndims(pair_meta.(fi{1}))));%uni-dir to bi-dir
%         pair.(fi{1})=cat(1,pair.(fi{1}),pair_meta.(fi{1}));
    end
    tic
    save(fullfile('bzdata',sprintf('0315_conn_w_reg_%d.mat',fidx)),'sig_meta','pair_meta','pc_stem','-v7','-nocompression')
    toc
end
