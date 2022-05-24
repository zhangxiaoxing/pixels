type='individual';
for bin=1:6
load(sprintf('0114_selec_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
for pref=1:2
    sel=find(pref_chain_S1(:,bin)==pref & pref_chain_S1(:,bin+1)==0 & pref_chain_S1(:,bin+6)==pref & pref_chain_S1(:,bin+7)==pref);
    sel=reshape(sel,1,[]);
    fidAll=idivide(int32(conn_chain_S1(sel,1)),int32(100000));
    uid1All=rem(conn_chain_S1(sel,1),100000);
    uid2All=rem(conn_chain_S1(sel,2),100000);
    for i=1:length(fidAll)
        fid=fidAll(i);
        uid1=uid1All(i);
        uid2=uid2All(i);
        close all;
        try
            plot_showcase_ft;
            keyboard
%             print('-dpng',sprintf('enhance_showcase_%d_%d_%d_%d.png',bin,fid,uid1,uid2),'-r300');
        catch
        end
    end

end
end

