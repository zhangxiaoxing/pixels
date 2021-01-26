sus_trans=h5read('../transient_6.hdf5','/sus_trans');%4+2+7
allpath=h5read('../transient_6.hdf5','/path');
allcid=h5read('../transient_6.hdf5','/cluster_id');
reg_all=h5read('../transient_6.hdf5','/reg');

load reg_keep.mat
reg_sel=ismember(deblank(reg_all),reg_set([1:112,114,115]));
nnz(startsWith(allpath,'M48_20191205_g0') & any(sus_trans(:,[1 2 4]),2) & reg_sel)

auc=h5read('../transient_6.hdf5','/auc');
fr=h5read('../transient_6.hdf5','/fr');
sel=h5read('../transient_6.hdf5','/raw_selectivity');
wrsp=h5read('../transient_6.hdf5','/wrs_p');

sus_trans=h5read('../transient_6.hdf5','/sus_trans');
reg_list=h5read('../transient_6.hdf5','/reg');

nnz((sus_trans(:,3)) & reg_sel)


upath=unique(allpath);
length(upath);
spath=cell(0);
for i=1:length(upath)
    pathstub=regexp(upath{i},'^.*?(?=\\)','match','once');
    spath{end+1}=pathstub;
end
spath=sort(unique(spath))';
p=[]
for i=1:length(spath)
    p(end+1)=nnz(startsWith(allpath,spath{i}))
end

mean(p)
std(p)/sqrt(numel(p))




to0=[];
from0=[];
for i=3:7
    from0=[from0,nnz(su_var_t(:,i-1)==0 & su_var_t(:,i)>0)];
    to0=[to0,nnz(su_var_t(:,i-1)>0 & su_var_t(:,i)==0)];
end
mean(from0./size(su_var_t,1))
mean(to0./size(su_var_t,1))

cnt=0
df=diff(su_var_t(:,2:end)>0,1,2);
for i=1:size(df,1)
    if min(df(i,:))<0 && max(df(i,:))>0
        if find(df(i,:)<0,1,'first')<find(df(i,:)>0,1,'last')
            cnt=cnt+1;
        end
    end
end
cnt/size(su_var_t,1)



sus_trans=h5read('../transient_3.hdf5','/sus_trans');%4+2+7
reg_all=h5read('../transient_3.hdf5','/reg');

reg_set=unique(reg_all);
white_matter=find(startsWith(reg_set,'a'),1);
reg_set=reg_set(1:white_matter-1);
reg_set(strcmp(reg_set,'Unlabeled'))=[];

subtotal=nnz(ismember(deblank(reg_all),deblank(reg_set)))
trans=nnz(ismember(deblank(reg_all),deblank(reg_set))& any(sus_trans(:,[2 4]),2))
sus=nnz(ismember(deblank(reg_all),deblank(reg_set))& sus_trans(:,1))
trans/subtotal
sus/subtotal





