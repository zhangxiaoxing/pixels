% function out=duration_wave(opt)
% arguments
%     opt.plot (1,1) logical = false
%     opt.ctx (1,1) logical = false
%     opt.align_test (1,1) logical = false
%     opt.quick_merge (1,1) logical = false
% end
% persistent out_ opt_
% if (isempty(out_) || ~isequaln(opt,opt_))  && ~opt.quick_merge
    meta=ephys.util.load_meta();
    [~,~,sessmap]=ephys.sessid2path(0);
    homedir=ephys.util.getHomedir('type','raw');
    anovameta=struct();
    [anovameta.sess,anovameta.allcid,anovameta.anovap]=deal([]);
    for sesskey=reshape(cell2mat(sessmap.keys()),1,[])
        disp(sesskey)
        fpath=fullfile(homedir,sessmap(sesskey),"FR_All_1000.hdf5");
        fr=h5read(fpath,'/FR_All');
        trials=h5read(fpath,'/Trials');
        suid=h5read(fpath,'/SU_id');

        dur_resp=behav.tag_block(trials,'wt',false);
        block_meta=[trials,dur_resp(:,end)];
        
        anovanps=nan(size(fr,2),7); %1:Samp,2:Dur,3:Bin,4:Samp*Dur,5:Samp*Bin,6:Dur*Bin,7:Samp*Dur*Bin
        for su=1:size(fr,2)
            sumat=squeeze(fr(:,su,1:7));
            sufrvec=sumat(:);
            sampvec=repmat(trials(:,5),7,1);
            durvec=repmat(trials(:,8),7,1);
            binvec=reshape(repmat((1:7),size(fr,1),1),[],1);
            anovanps(su,:)=anovan(sufrvec,{sampvec,durvec,binvec},'model','full','display','off');
        end

        anovameta.sess=[anovameta.sess;repmat(sesskey,size(fr,2),1)];
        anovameta.allcid=[anovameta.allcid;suid];
        anovameta.anovap=[anovameta.anovap;anovanps];
    end


% end

% meta=ephys.util.load_meta();
% waveid=ephys.get_wave_id(meta.sess,meta.allcid);
% nnz(anovameta.anovap(:,7)<0.05)

meta=ephys.util.load_meta();
idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));

BSsel=strcmp(meta.reg_tree(1,:),'BS') & ~strcmp(meta.reg_tree(5,:),'');
CHsel=strcmp(meta.reg_tree(1,:),'CH') & ~strcmp(meta.reg_tree(5,:),'');
grey_regs=unique(meta.reg_tree(5,BSsel | CHsel));

cnt=cellfun(@(x) nnz(strcmp(meta.reg_tree(5,:),x)), grey_regs);
grey_regs=grey_regs(cnt>100);

dur_indep_sel=(anovameta.anovap(:,1)<0.05 | anovameta.anovap(:,5)<0.05);
dur_dep_sel=(anovameta.anovap(:,4)<0.05 | anovameta.anovap(:,7)<0.05) &~dur_indep_sel;
dur_only_sel=(anovameta.anovap(:,2)<0.05 | anovameta.anovap(:,6)<0.05) & ~dur_dep_sel;

dur_only_reg_map=containers.Map();
sens_dur_reg_map=containers.Map();
dur_any_reg_map=containers.Map();
for r=grey_regs
    cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}));
    dur_only_cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}).' & dur_only_sel);
    sens_dur_cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}).' & dur_dep_sel);
    dur_any_cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}).' & (dur_dep_sel|dur_only_sel));
    dur_only_reg_map(r{1})=[dur_only_cnt/cnt,dur_only_cnt,cnt];
    sens_dur_reg_map(r{1})=[sens_dur_cnt/cnt,sens_dur_cnt,cnt];
    dur_any_reg_map(r{1})=[dur_any_cnt/cnt,dur_any_cnt,cnt];
end



feat_reg_map=dur_any_reg_map;


% To Line 98, K:\code\jpsth\+hier\dur_sens_ordinal_hierachy.m
% for GLM go K:\code\jpsth\+wave\connectivity_proportion_GLM.m

