% This code was seperated from original xcorr code to streamline pipeline
% Appends per-neuron infomation in the sums file
% input is list of files, output is cell of sums.

bin_range=[-2 -1];
rpt_workaround={'M23_20191109_g0';'191018-DPA-Learning5_28_g1';'191226_64_learning6_g0_imec1_cleaned'};

sus_trans=h5read('../transient_6.hdf5','/sus_trans');
reg_list=h5read('../transient_6.hdf5','/reg');
cid_list=h5read('../transient_6.hdf5','/cluster_id');
path_list=h5read('../transient_6.hdf5','/path');
if isunix
    path_list=replace(path_list,'\','/');
end

sust=find(sus_trans(:,1));

disp(length(cid_list));
% pause;

trans=find(sus_trans(:,2));
supool=[sust;trans]';
counter=[];
done=[];

singleTracks=dir('/home/zx/neupix/wyt/DataSum/singleProbe');
singleTrackList={singleTracks(3:end).name};
% error_list=cell(0);
sums=cell(0,8);
fs=dir(sprintf('/home/zx/pixels/jpsth/0623_selec_XCORR_duo*_delay_6_%d_%d_2msbin.mat',bin_range(1),bin_range(2)));
disp(length(fs))
keyboard
% pause
for i=1:length(fs)
    tracks=2;
    disp(i);
    fstr=load(fullfile(fs(i).folder,fs(i).name));
    if isunix
        dpath=replace(fstr.sums{2},'\','/');
    else
        dpath=fstr.sums{2};
    end
    
    dpath_stem=regexp(dpath,'([^\\/]*)','match','once');
    if any(startsWith(rpt_workaround,dpath_stem))
        disp('redundant work around folder:');
        disp(dpath);
        continue
    end
    if any(startsWith(singleTrackList,dpath_stem))
        disp('Single track:');
        disp(dpath);
        tracks=1;
    end
    if tracks==2
        dpath=replace(dpath,'imec1','imec0');
    end
    wffile=fullfile(['/home/zx/neupix/DataSum/',dpath],'wf_stats.hdf5');% posix
    if isfile(wffile)
        wfstats=h5read(wffile,'/wf');
        wfstats1=[];
        xc_s1=fstr.sums{5};
        for lblidx=1:size(xc_s1.label,1)
            currlbl=str2double(xc_s1.label{lblidx,1});
            if currlbl>=10000
                imecNo=1;
                currlbl=currlbl-10000;
                if isempty(wfstats1)
                    if isfile(replace(wffile,'imec0','imec1'))
                        wfstats1=h5read(replace(wffile,'imec0','imec1'),'/wf');
                    else
                        disp(replace(wffile,'imec0','imec1'))
                        keyboard
                    end
                end
            else
                imecNo=0;
            end
            if imecNo==0
                wfidx=find(wfstats(:,1)==currlbl);
                if ~isempty(wfidx)
                    xc_s1.label{lblidx,2}=wfstats(wfidx,:);
                    raw_wf_file=fullfile(['/home/zx/neupix/WF/neuropixel/',dpath],'waveform.mat'); %posix
                    raw_fstr=load(raw_wf_file);
                    raw_idx=find([raw_fstr.waveform{:,2}]==wfstats(wfidx,1));
                    xc_s1.label{lblidx,3}=raw_fstr.waveform{wfidx,4};
                    %% prefered sample, reg,
                    suid=find(startsWith(path_list,dpath) & cid_list==str2double(xc_s1.label{lblidx,1}));
                    % sust, transient, switched, unclassified, early_in_6s,
                    % late_in_6s, 7X prefer_s
                    prefered_sample=sus_trans(suid,7:end);
                    reg=regexp(reg_list{suid},'(\w|\d)+','match','once');
                    xc_s1.label{lblidx,4}=prefered_sample;
                    xc_s1.label{lblidx,5}=reg;
                    xc_s1.label{lblidx,6}=0;
                else
                    disp(wfidx);
                    keyboard
%                     pause
                end
            else %imec==1
                wfidx=find(wfstats1(:,1)==currlbl);
                if ~isempty(wfidx)
                    xc_s1.label{lblidx,2}=wfstats1(wfidx,:);
                    raw_wf_file=fullfile(['/home/zx/neupix/WF/neuropixel/',replace(dpath,'imec0','imec1')],'waveform.mat'); %posix
                    if isfile(raw_wf_file)
                        raw_fstr=load(raw_wf_file);
                        raw_idx=find([raw_fstr.waveform{:,2}]==wfstats1(wfidx,1));
                        xc_s1.label{lblidx,3}=raw_fstr.waveform{wfidx,4};
                        %% prefered sample, reg,
                        suid=find(startsWith(path_list,replace(dpath,'imec0','imec1')) & cid_list==currlbl);
                        % sust, transient, switched, unclassified, early_in_6s,
                        % late_in_6s, 7X prefer_s
                        prefered_sample=sus_trans(suid,7:end);
                        reg=regexp(reg_list{suid},'(\w|\d)+','match','once');
                        xc_s1.label{lblidx,4}=prefered_sample;
                        xc_s1.label{lblidx,5}=reg;
                        xc_s1.label{lblidx,6}=1;
                    else
                        disp(wfidx);
                        keyboard
                    end
                else
                    disp(wfidx)
                    keyboard
                end
            end
        end
        fstr.sums{5}=xc_s1;
    else
        disp(wffile);
        keyboard
        continue
    end
    sums(end+1,:)=fstr.sums;
end
disp('check the file name is correct!')
% keyboard
save(sprintf('0623_selec_XCORR_duo_sums_delay_6_%d_%d_2msbin.mat',bin_range(1),bin_range(2)),'sums','-v7.3')

