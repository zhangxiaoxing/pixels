%TODO brain region filter, olfaction filter.
function v=get_variance(opt)
arguments
    opt.keep_trials (1,1) logical = false
    opt.ctx (1,1) logical = false
    opt.seltype (1,:) char {mustBeMember(opt.seltype,{'olfactory','duration','odsel','any'})} = 'any'
end
persistent varmat opt_
if isempty(varmat) || ~isequaln(opt,opt_)
    [~,~,sessmap]=ephys.sessid2path(0);
    meta=ephys.util.load_meta();
    homedir=ephys.util.getHomedir();
    for ii=reshape(cell2mat(sessmap.keys()),1,[])
        disp(ii)
        fpath=fullfile(homedir,sessmap(ii),"FR_All_ 250.hdf5");
        fr=h5read(fpath,'/FR_All');
        trials=h5read(fpath,'/Trials');
        suid=h5read(fpath,'/SU_id');
        msel=true(size(suid));
        if ~strcmp(opt.seltype,'any')
            sesssel=meta.sess==ii;
            regsel=strcmp(meta.reg_tree(2,sesssel),'CTX');
            olfsel=meta.mem_type(sesssel)>0;
            switch opt.seltype
                case 'olfactory'
                    msel=olfsel;
            end
        end
        if opt.ctx
            msel=msel&regsel;
        end
        fr=fr(:,msel,:);
        suid=suid(msel);

        %     e1sel=find(trial(:,5)==4 & trial(:,8)==opt.delay & trial(:,10)==0);
        %     e2sel=find(trial(:,5)==8 & trial(:,8)==opt.delay & trial(:,10)==0);

        %     sess=['s',num2str(ii)];
        pca_str=struct();
        for ff=["d3s1","d3s2","d6s1","d6s2"]
            pca_str.(ff)=containers.Map('KeyType','int32','ValueType','any');
        end
        pca_str=per_su_process(suid,fr,trials,pca_str,opt);
        if opt.keep_trials
            varmat=[varmat;pca_str];
        else
            varmat=[varmat;cell2mat(pca_str.d3s1.values(num2cell(suid.'))).',...
                cell2mat(pca_str.d3s2.values(num2cell(suid.'))).',...
                cell2mat(pca_str.d6s1.values(num2cell(suid.'))).',...
                cell2mat(pca_str.d6s2.values(num2cell(suid.'))).'];
        end
    end
end
pcamat_=varmat;
opt_=opt;
end
function pca_str=per_su_process(suid,fr,trials,pca_str,opt)
% arguments
%     sess,suid,fr,trials,pca_str,opt
% end

d3s1sel=find(trials(:,5)==4 & trials(:,8)==3 & trials(:,9)>0 & trials(:,10)>0);
d3s2sel=find(trials(:,5)==8 & trials(:,8)==3 & trials(:,9)>0 & trials(:,10)>0);
d6s1sel=find(trials(:,5)==4 & trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0);
d6s2sel=find(trials(:,5)==8 & trials(:,8)==6 & trials(:,9)>0 & trials(:,10)>0);


for su=1:size(fr,2)
    if opt.keep_trials
        pca_str.d3s1(suid(su))=squeeze(fr(d3s1sel,su,:));
        pca_str.d3s2(suid(su))=squeeze(fr(d3s2sel,su,:));
        pca_str.d6s1(suid(su))=squeeze(fr(d6s1sel,su,:));
        pca_str.d6s2(suid(su))=squeeze(fr(d6s2sel,su,:));
    else
        pca_str.d3s1(suid(su))=squeeze(mean((fr(d3s1sel,su,:))));
        pca_str.d3s2(suid(su))=squeeze(mean((fr(d3s2sel,su,:))));
        pca_str.d6s1(suid(su))=squeeze(mean((fr(d6s1sel,su,:))));
        pca_str.d6s2(suid(su))=squeeze(mean((fr(d6s2sel,su,:))));
    end
end
end


