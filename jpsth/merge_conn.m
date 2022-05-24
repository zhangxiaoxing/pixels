if ~exist('bin_range','var')
    disp('missing bin_range');
    return
end
load(sprintf('0831_conn_chain_duo_6s_%d_%d.mat',bin_range(1),bin_range(2)),'pair_reg','pair_chain')
bz_conn_chain_S1=[];
bz_conn_reg_S1=[];
bz_spk_count_S1=[];
bz_conn_chain_S2=[];
bz_conn_reg_S2=[];
bz_spk_count_S2=[];
for sess=1:114
    sesbias=sess*100000;
    if ~exist(sprintf('0831_full_BZ_XCORR_duo_f%d_delay_6_%d_%d_2msbin.mat', ...
        sess,bin_range(1),bin_range(2)),'file')
        continue
    end
    fbz=load(sprintf('0831_full_BZ_XCORR_duo_f%d_delay_6_%d_%d_2msbin.mat', ...
        sess,bin_range(1),bin_range(2)));
    if ~isempty(fbz.sumsbz{5})
        s1=fbz.sumsbz{5}.sig_con;
        bz_conn_chain_S1=[bz_conn_chain_S1;reshape(fbz.sumsbz{5}.LUT(s1)+sesbias,[],2)];
        s1cc=s1;
        s1c=fbz.sumsbz{6};
        for i=1:size(s1,1)
            prec=s1c(s1c(:,1)==fbz.sumsbz{5}.LUT(s1(i,1)),2);
            postc=s1c(s1c(:,1)==fbz.sumsbz{5}.LUT(s1(i,2)),2);
            s1cc(i,:)=[prec,postc];
        end
        bz_spk_count_S1=[bz_spk_count_S1;s1cc];
        reg_S1=s1;
        for i=1:size(s1,1)
            reg1=find(pair_chain==fbz.sumsbz{5}.LUT(s1(i,1))+sesbias,1);
            reg2=find(pair_chain==fbz.sumsbz{5}.LUT(s1(i,2))+sesbias,1);
            if isempty(reg1) | isempty(reg2)
               reg_S1(i,:)=[0,0];
            else
               reg_S1(i,:)=[pair_reg(reg1),pair_reg(reg2)];
            end
        end
        bz_conn_reg_S1=[bz_conn_reg_S1;reg_S1];

    end
    if ~isempty(fbz.sumsbz{7})
        s2=fbz.sumsbz{7}.sig_con;
        bz_conn_chain_S2=[bz_conn_chain_S2;reshape(fbz.sumsbz{7}.LUT(s2)+sesbias,[],2)];
        s2cc=s2;
        s2c=fbz.sumsbz{8};
        for i=1:size(s2,1)
            prec=s2c(s2c(:,1)==fbz.sumsbz{7}.LUT(s2(i,1)),2);
            postc=s2c(s2c(:,1)==fbz.sumsbz{7}.LUT(s2(i,2)),2);
            s2cc(i,:)=[prec,postc];
        end
        bz_spk_count_S2=[bz_spk_count_S2;s2cc];
        reg_S2=s2;
        for i=1:size(s2,1)
           reg1=find(pair_chain==fbz.sumsbz{7}.LUT(s2(i,1))+sesbias,1);
           reg2=find(pair_chain==fbz.sumsbz{7}.LUT(s2(i,2))+sesbias,1);
            if isempty(reg1) | isempty(reg2)
               reg_S2(i,:)=[0,0];
            else
               reg_S2(i,:)=[pair_reg(reg1),pair_reg(reg2)];
           end
        end
        bz_conn_reg_S2=[bz_conn_reg_S2;reg_S2];
    end
end
sel_reg_s1=bz_conn_reg_S1(:,1)>0;
sel_reg_s2=bz_conn_reg_S2(:,1)>0;
bz_conn_chain_S1=bz_conn_chain_S1(sel_reg_s1,:);
bz_conn_chain_S2=bz_conn_chain_S2(sel_reg_s2,:);
bz_spk_count_S1=bz_spk_count_S1(sel_reg_s1,:);
bz_spk_count_S2=bz_spk_count_S2(sel_reg_s2,:);
bz_conn_reg_S1=bz_conn_reg_S1(sel_reg_s1,:);
bz_conn_reg_S2=bz_conn_reg_S2(sel_reg_s2,:);

blame=datetime();
save(sprintf('0831_conn_chain_duo_6s_%d_%d.mat',bin_range(1),bin_range(2)),...
    'bz_conn_chain_S1','bz_conn_chain_S2','bz_spk_count_S1','bz_spk_count_S2',...
    'bz_conn_reg_S1','bz_conn_reg_S2','blame','-append');



function re_merge_conn()
for bin=1:6
    bin_range=[bin,bin+1]
    load(sprintf('0831_conn_chain_duo_6s_%d_%d.mat',bin_range(1),bin_range(2)),'bz_conn_chain_S1','bz_conn_chain_S2','pair_chain','pref_pair');
    pchain=[pair_chain(:,1);pair_chain(:,2)];
    ppref=[pref_pair(:,1:6);pref_pair(:,7:12)];
    bz_pref_S1=repmat(bz_conn_chain_S1(:,1),1,6);
    for i=1:size(bz_conn_chain_S1,1)
        bz_pref_S1(i,:)=ppref(find(pchain==bz_conn_chain_S1(i,1),1),:);
    end
    bz_pref_S2=repmat(bz_conn_chain_S2(:,1),1,6);
    for i=1:size(bz_conn_chain_S2,1)
        bz_pref_S2(i,:)=ppref(find(pchain==bz_conn_chain_S2(i,1),1),:);
    end

    blame=datetime();
    save(sprintf('0831_conn_chain_duo_6s_%d_%d.mat',bin_range(1),bin_range(2)),'bz_pref_S1','bz_pref_S2','blame','-append');
end
end
