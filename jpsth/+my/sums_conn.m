% assume 'sums' is loaded in workspace. Otherwise load corresponding
% XCORR_delay_bin.mat file first
% A peak at a negative lag (i.e., AI>0,. i.e. L-R>0) for stat.xcorr(chan1,chan2,:) means that chan1 is leading
% chan2.
function [out,avail]=sums_conn(opt)
arguments
    opt.prefix (1,:) char = '0319'
    opt.bin_range (1,2) double = [-2 7]
    opt.ntrial (1,1) double  {mustBePositive,mustBeInteger} = 20
    opt.nspk_thres (1,1) double {mustBePositive,mustBeInteger} = 1000
%     opt.overwrite (1,1) logical = false
%     opt.to_plot=false;
%     opt.to_save=false;
%     opt.to_process_reverb=true;
end
fs=dir(fullfile('mydata',sprintf('%s_XCORR_duo_f*_6_%d_%d_2msbin.mat',opt.prefix,opt.bin_range(1),opt.bin_range(2))));
sums_conn_str=[];
for fidx=1:size(fs,1)
    disp(fidx);
    load(fullfile(fs(fidx).folder,fs(fidx).name),'sums');
    thresh=norminv(0.9975); %0.05 bonferroni corrected
    xc_s1=sums{2};
    xshuf_s1=sums{3};
    xc_s2=sums{4};
    xshuf_s2=sums{5};
    if isempty(xc_s1) || isempty(xc_s2)
        out=[];
        avail=false;
        continue
    end
    f=struct();
    [f.sums_conn_s1,f.meta_s1]=my.get_sig(xc_s1,xshuf_s1,'ntrial',opt.ntrial,'nspk_thres',opt.nspk_thres,'thresh',thresh);
    [f.sums_conn_s2,f.meta_s2]=my.get_sig(xc_s2,xshuf_s2,'ntrial',opt.ntrial,'nspk_thres',opt.nspk_thres,'thresh',thresh);
    f.folder=sums{1};
    sums_conn_str=[sums_conn_str;f];
end
blame=vcs.blame();
save(sprintf('sums_conn_my_%d_%d.mat',opt.bin_range(1),opt.bin_range(2)),'sums_conn_str','blame');
end
