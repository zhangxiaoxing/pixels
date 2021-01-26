load('reg_keep.mat');
reg_set=reg_set(1:115);

load('0906_nonsel_conn_chain_duo_6s_1_2.mat','pair_chain','pair_reg','pref_pair');
pref=[max(pref_pair(:,1:6),[],2),max(pref_pair(:,7:end),[],2)];
dim=length(reg_set);
motif_candi=zeros(dim);

pair_mat=zeros(length(reg_set),length(reg_set),6);

for sess=1:115
    disp('#####')
    disp(sess)
    lbound=sess*100000;
    ubound=(sess+1)*100000;
    sess_sel=pair_chain(:,1)>=lbound & pair_chain(:,1)<ubound;
    disp(nnz(sess_sel))
    [uids_all,UA]=unique(pair_chain(sess_sel,:));
    %S1 S2
    ses_reg=pair_reg(sess_sel,:);
    regs_all=ses_reg(UA);
    reg_sel=regs_all<116;
    ses_pref=pref(sess_sel,:);
    prefs=ses_pref(UA);
    for onepref=1:2
        uids=uids_all(reg_sel & prefs==onepref);
        regs=regs_all(reg_sel & prefs==onepref);
        sess_reg=unique(regs);
        if numel(sess_reg)<2
            continue
        end
        for i=1:(length(sess_reg)-1)
            for j=(i+1):length(sess_reg)
                n1=nnz(regs==sess_reg(i));
                n2=nnz(regs==sess_reg(j));
                nx=numel(uids)-n1-n2;
                if nx<0
                    keyboard()
                end
                sums=0;
                %1,1,c x 2
                if nx>=2
                    sums=sums+n1*n2*nchoosek(nx,2);
                end
                %1,2,x
                if n2>=2
                    sums=sums+n1*nchoosek(n2,2)*nx;
                end
                %1,3
                if n2>=3
                    sums=sums+n1*nchoosek(n2,3);
                end
                %2,1,x
                if n1>=2
                    sums=sums+n2*nchoosek(n1,2)*nx;
                end
                %3,1
                if n1>=3
                    sums=sums+n2*nchoosek(n1,3);
                end
                %2,2
                if n1>=2 && n2>=2
                    sums=sums+nchoosek(n1,2)*nchoosek(n2,2);
                end
                motif_candi(sess_reg(i),sess_reg(j))=motif_candi(sess_reg(i),sess_reg(j))+sums;
            end
        end
    end
end

csvcell=cell(0,3);
for i=1:(length(motif_candi)-1)
    for j=(i+1):length(motif_candi)
        if motif_candi(i,j)>0
            csvcell(end+1,:)={reg_set{i},reg_set{j},motif_candi(i,j)};
        end
    end
end
writecell(csvcell,'motif_candi.csv');
