%TODO Error trial
function out=delay_dur_selective(opt)
arguments
    opt.window (1,:) double = 5:6 %
    opt.one_SU_showcase (1,1) logical = false % for the TCOM-FC joint showcase
end

meta_str=ephys.util.load_meta('type','neupix','delay',6);
homedir=ephys.util.getHomedir('type','raw');
fl=dir(fullfile(homedir,'**','FR_All_1000.hdf5'));
out=[];
for ii=1:size(fl,1)
    dpath=regexp(fl(ii).folder,'(?<=SPKINFO[\\/]).*$','match','once');
    fpath=fullfile(fl(ii).folder,fl(ii).name);
    disp(ii);
    fpath=replace(fpath,'\',filesep());
    pc_stem=replace(dpath,'/','\');
    sesssel=startsWith(meta_str.allpath,pc_stem);
    if ~any(sesssel), continue;end
    fr=h5read(fpath,'/FR_All');
    trial=h5read(fpath,'/Trials');
    suid=h5read(fpath,'/SU_id');
    if isempty(suid), continue; end
    sessid=ephys.path2sessid(pc_stem);
    s1_s3=(trial(:,5)==4 & trial(:,8)==3 & trial(:,9)>0 & trial(:,10)>0);
    s1_s6=(trial(:,5)==4 & trial(:,8)==6 & trial(:,9)>0 & trial(:,10)>0);
    s2_s3=(trial(:,5)==8 & trial(:,8)==3 & trial(:,9)>0 & trial(:,10)>0);
    s2_s6=(trial(:,5)==8 & trial(:,8)==6 & trial(:,9)>0 & trial(:,10)>0);

    s3e=(trial(:,8)==3 & trial(:,10)==0);
    s6e=(trial(:,8)==6 & trial(:,10)==0);

    flat=@(x) x(:);
    [s1p,s2p,sp,fr3,fr6,si,auc,sie]=deal(nan(size(suid)));
    for ss=1:numel(suid)
        s1p(ss)=ranksum(flat(fr(s1_s3,ss,opt.window)),flat(fr(s1_s6,ss,opt.window)));
        s2p(ss)=ranksum(flat(fr(s2_s3,ss,opt.window)),flat(fr(s2_s6,ss,opt.window)));
        sp(ss)=ranksum(flat(fr(s1_s3|s2_s3,ss,opt.window)),flat(fr(s1_s6|s2_s6,ss,opt.window)));
        fr3(ss)=mean(fr(s1_s3|s2_s3,ss,opt.window),'all');
        fr6(ss)=mean(fr(s1_s6|s2_s6,ss,opt.window),'all');

        fr3e(ss)=mean(fr(s3e,ss,opt.window),'all');
        fr6e(ss)=mean(fr(s6e,ss,opt.window),'all');
        
        [~,~,~,auc(ss)]=perfcurve([zeros(nnz(s1_s3|s2_s3),1);ones(nnz(s1_s6|s2_s6),1)],...
                    [sum(fr(s1_s3|s2_s3,ss,opt.window),3);sum(fr(s1_s6|s2_s6,ss,opt.window),3)], ...
                    fr6(ss)>fr3(ss));
        si(ss)=(fr3(ss)-fr6(ss))/(fr3(ss)+fr6(ss));
        sie(ss)=(fr3e(ss)-fr6e(ss))/(fr3e(ss)+fr6e(ss));
            
    end
    out=[out;ones(size(suid)).*double(ii),ones(size(suid)).*double(sessid),suid,s1p,s2p,sp,fr3,fr6,si,sie,auc];
end
end

