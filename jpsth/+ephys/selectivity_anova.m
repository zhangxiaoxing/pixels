%%CONST
seq_mdl=false;
%%
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
        if seq_mdl
            anovanps=nan(size(fr,2),15); 
            % 1:Samp,2:Dur,3:Bin,4:Seq,5:Samp*Dur,6:Samp*Bin,7:Samp*Seq,8:Dur*Bin,9:Dur*Seq,10:Bin*Seq,
            % 11:Samp*Dur*Bin, 12:Samp*Dur*Seq, 13:Samp*Bin*Seq,% 14:Dur*Bin*Seq,
            % 15:Samp*Dur*Bin*Seq

        else
            anovanps=nan(size(fr,2),7); %1:Samp,2:Dur,3:Bin,4:Samp*Dur,5:Samp*Bin,6:Dur*Bin,7:Samp*Dur*Bin
        end
        for su=1:size(fr,2)
            frmat=squeeze(fr(:,su,1:7));
            sufrvec=frmat(:);
            sampvec=repmat(trials(:,5),7,1);
            durvec=repmat(trials(:,8),7,1);
            binvec=reshape(repmat((1:7),size(fr,1),1),[],1);
            if seq_mdl
                seqvec=repmat(block_meta(:,11),7,1);
                anovanps(su,:)=anovan(sufrvec,{sampvec,durvec,binvec,seqvec},'model','full','display','off');
            else
                anovanps(su,:)=anovan(sufrvec,{sampvec,durvec,binvec},'model','full','display','off');
            end
        end

        anovameta.sess=[anovameta.sess;repmat(sesskey,size(fr,2),1)];
        anovameta.allcid=[anovameta.allcid;suid];
        anovameta.anovap=[anovameta.anovap;anovanps];
    end
keyboard()

% end


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

seq_any_reg_map=containers.Map();
seq_wo_sample_reg_map=containers.Map();
dur_any_reg_map=containers.Map();
for r=grey_regs
    cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}));
    seq_any_cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}).' & dur_only_sel);
    seq_wo_sample_cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}).' & dur_dep_sel);
    dur_any_cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}).' & (dur_dep_sel|dur_only_sel));
    seq_any_reg_map(r{1})=[seq_any_cnt/cnt,seq_any_cnt,cnt];
    seq_wo_sample_reg_map(r{1})=[seq_wo_sample_cnt/cnt,seq_wo_sample_cnt,cnt];
    dur_any_reg_map(r{1})=[dur_any_cnt/cnt,dur_any_cnt,cnt];
end



feat_reg_map=dur_any_reg_map;


% To Line 98, K:\code\jpsth\+hier\dur_sens_ordinal_hierachy.m
% for GLM go K:\code\jpsth\+wave\connectivity_proportion_GLM.m




%%
if false
    figure('Color','w')
bar(sum(anovameta.anovap<0.05)./size(anovameta.anovap,1))


% 1:Samp,2:Dur,3:Bin,4:Seq,5:Samp*Dur,6:Samp*Bin,7:Samp*Seq,8:Dur*Bin,9:Dur*Seq,10:Bin*Seq,
% 11:Samp*Dur*Bin, 12:Samp*Dur*Seq, 13:Samp*Bin*Seq,% 14:Dur*Bin*Seq,
% 15:Samp*Dur*Bin*Seq

set(gca(),'XTick',1:15,'XTickLabel',{'1-Samp','2-Dur','3-Bin','4-Seq','5-Samp*Dur','6-Samp*Bin','7-Samp*Seq','8-Dur*Bin','9-Dur*Seq','10-Bin*Seq',...
'11-Samp*Dur*Bin', '12-Samp*Dur*Seq','13-Samp*Bin*Seq','14-Dur*Bin*Seq','15-Samp*Dur*Bin*Seq'},...
'YTick',0:0.1:0.6,'YTickLabel',0:10:60)
ylabel('Fraction of all neurons (%)')


% per reg, only consider dur and seq
meta=ephys.util.load_meta();
idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));

BSsel=strcmp(meta.reg_tree(1,:),'BS') & ~strcmp(meta.reg_tree(5,:),'');
CHsel=strcmp(meta.reg_tree(1,:),'CH') & ~strcmp(meta.reg_tree(5,:),'');
grey_regs=unique(meta.reg_tree(5,BSsel | CHsel));

cnt=cellfun(@(x) nnz(strcmp(meta.reg_tree(5,:),x)), grey_regs);
grey_regs=grey_regs(cnt>100);


seq_any_sel=any(anovameta.anovap(:,[4 7 9 10 12:15])<0.05,2);
seq_wo_sample_sel=seq_any_sel & ~any(anovameta.anovap(:,[1,5:7,11:13,15])<0.05,2);

samp_any_sel=any(anovameta.anovap(:,[1,5:7,11:13,15])<0.05,2);
samp_wo_seq_sel=samp_any_sel & ~any(anovameta.anovap(:,[4 7 9 10 12:15])<0.05,2);

seq_any_reg_map=containers.Map();
seq_wo_sample_reg_map=containers.Map();

samp_any_reg_map=containers.Map();
samp_wo_seq_reg_map=containers.Map();

for r=grey_regs
    cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}));
    seq_any_cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}).' & seq_any_sel);
    seq_wo_sample_cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}).' & seq_wo_sample_sel);

    samp_any_cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}).' & samp_any_sel);
    samp_wo_seq_cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}).' & samp_wo_seq_sel);

    seq_any_reg_map(r{1})=[seq_any_cnt/cnt,seq_any_cnt,cnt];
    seq_wo_sample_reg_map(r{1})=[seq_wo_sample_cnt/cnt,seq_wo_sample_cnt,cnt];

    samp_any_reg_map(r{1})=[samp_any_cnt/cnt,samp_any_cnt,cnt];
    samp_wo_seq_reg_map(r{1})=[samp_wo_seq_cnt/cnt,samp_wo_seq_cnt,cnt];
end

%% any
figure('Color','w')
sampmat=cell2mat(samp_any_reg_map.values(grey_regs).')
seqmat=cell2mat(seq_any_reg_map.values(grey_regs).')
scatter(sampmat(:,1),seqmat(:,1))
[r,p]=corr(sampmat(:,1),seqmat(:,1))
%% exclusive

sampmat=cell2mat(samp_wo_seq_reg_map.values(grey_regs).')
seqmat=cell2mat(seq_wo_sample_reg_map.values(grey_regs).')
[r,p]=corr(sampmat(:,1),seqmat(:,1))
figure('Color','w')
hold on
for rr=reshape(grey_regs,1,[])
    c=ephys.getRegColor(rr{1},'large_area',true);
    ss=samp_wo_seq_reg_map(rr{1});
    qq=seq_wo_sample_reg_map(rr{1});
    scatter(ss(1),qq(1),9,c,'filled','o')
    text(ss(1),qq(1),rr{1},'HorizontalAlignment','center','VerticalAlignment','top','Color',c);
end

feat_reg_map=seq_wo_sample_reg_map;

% To Line 98, K:\code\jpsth\+hier\dur_sens_ordinal_hierachy.m
% for GLM go K:\code\jpsth\+wave\connectivity_proportion_GLM.m



end




%% dur v seq



% per reg, only consider dur and seq
meta=ephys.util.load_meta();
idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));

BSsel=strcmp(meta.reg_tree(1,:),'BS') & ~strcmp(meta.reg_tree(5,:),'');
CHsel=strcmp(meta.reg_tree(1,:),'CH') & ~strcmp(meta.reg_tree(5,:),'');
grey_regs=unique(meta.reg_tree(5,BSsel | CHsel));

cnt=cellfun(@(x) nnz(strcmp(meta.reg_tree(5,:),x)), grey_regs);
grey_regs=grey_regs(cnt>100);


seq_any_sel=any(anovameta.anovap(:,[4 7 9 10 12:15])<0.05,2);
dur_any_sel=any(anovameta.anovap(:,[2 5 8 9 11 12 14 15])<0.05,2);

seq_wo_dur_sel=seq_any_sel & ~dur_any_sel;
dur_wo_seq_sel=samp_any_sel & ~seq_any_sel;

seq_any_reg_map=containers.Map();
seq_wo_dur_reg_map=containers.Map();

dur_any_reg_map=containers.Map();
dur_wo_seq_reg_map=containers.Map();

for r=grey_regs
    cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}));
    seq_any_cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}).' & seq_any_sel);
    seq_wo_dur_cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}).' & seq_wo_dur_sel);

    dur_any_cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}).' & dur_any_sel);
    dur_wo_seq_cnt=nnz(strcmp(meta.reg_tree(5,:),r{1}).' & dur_wo_seq_sel);

    seq_any_reg_map(r{1})=[seq_any_cnt/cnt,seq_any_cnt,cnt];
    seq_wo_dur_reg_map(r{1})=[seq_wo_dur_cnt/cnt,seq_wo_dur_cnt,cnt];

    dur_any_reg_map(r{1})=[dur_any_cnt/cnt,dur_any_cnt,cnt];
    dur_wo_seq_reg_map(r{1})=[dur_wo_seq_cnt/cnt,dur_wo_seq_cnt,cnt];
end

%% any
figure('Color','w')
durmat=cell2mat(dur_any_reg_map.values(grey_regs).')
seqmat=cell2mat(seq_any_reg_map.values(grey_regs).')
scatter(durmat(:,1),seqmat(:,1))
[r,p]=corr(durmat(:,1),seqmat(:,1))

%% exclusive

durmat=cell2mat(dur_wo_seq_reg_map.values(grey_regs).')
seqmat=cell2mat(seq_wo_dur_reg_map.values(grey_regs).')
[r,p]=corr(durmat(:,1),seqmat(:,1))
figure('Color','w')
hold on
for rr=reshape(grey_regs,1,[])
    c=ephys.getRegColor(rr{1},'large_area',true);
    ss=dur_wo_seq_reg_map(rr{1});
    qq=seq_wo_dur_reg_map(rr{1});
    scatter(ss(1),qq(1),9,c,'filled','o')
    text(ss(1),qq(1),rr{1},'HorizontalAlignment','center','VerticalAlignment','top','Color',c);
end