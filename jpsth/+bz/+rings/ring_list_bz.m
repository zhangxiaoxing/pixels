function ring_list_bz(type)
arguments
    type (1,:) char {mustBeMember(type,{'delay','delay_shuf'})}
end
%bzthres=250;  %TODO filter by spike number % not really necessary when using full-length data
if strcmp(type,'delay')
    rings=cell(max(sig.sess),3);
    for sess=1:max(sig.sess)
        disp(sess);
        for ring_size=3:5
            sess_sel = sig.sess==sess;
            sess_rings=bz.rings.find_rings_bz(sig.suid(sess_sel,:),ring_size);
            rings{sess,ring_size-2}=unique(bz.rings.flexsort(sess_rings),'rows');
        end
    end
    save(fullfile('bzdata','rings_bz.mat'),'rings');
end

%TODO anything after is not looked into as of Mar 22 2021

if strcmp(type,'delay_shuf')
    shufrpt=10;
    rings_shuf=cell(shufrpt,max(sig.sess),3);
    for rpt=1:shufrpt
        disp(rpt)
        for sess=1:max(sig.sess)
            sess_sel = sig.sess==sess;
            bz.rings.shuffle_conn_chain(sig.suid(sess_sel,:));
            [shufchainS2,shufregS2,shufprefS2]=bz.rings.shuffle_conn_chain(fstr{bin}.bz_conn_chain_S2,fstr{bin}.pair_chain,fstr{bin}.pair_reg,fstr{bin}.pref_pair);
            for ring_size=3:5
                if nnz(sel21)>0
                    sess_rings=count_motif_congru(shufchainS1(sel21,:),shufregS1(sel21,:),shufprefS1(sel21,:),bin,msize,1);
                    sess_rings=unique(flexsort(sess_rings),'rows');
                    rings_shuf(rpt,in_motif_idx,sess,bin,:)={sess_rings,onering2};
                end
            end
        end
    end
    
    if exist('delay_shuf','var') && delay_shuf
        save('rings_bz.mat','rings_shuf','-append');
    end
end
end