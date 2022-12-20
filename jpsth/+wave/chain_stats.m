% chains_uf=wave.COM_chain(sel_meta);
% chains_uf_rev=wave.COM_chain(sel_meta,'reverse',true);
% blame=vcs.blame();
% save('chains_mix.mat','chains_uf','chains_uf_rev','blame')
load('chains_mix.mat','chains_uf','chains_uf_rev');

% Forward and reverse count

len=cellfun(@(x) numel(x),chains_uf.cids);
len_rev=cellfun(@(x) numel(x),chains_uf_rev.cids);

len_hist=histcounts(len,2.5:1:12.5);
len_rev_hist=histcounts(len_rev,2.5:1:12.5);

xx=3:12;

figure()
hold on;
fwdh=plot(xx,len_hist,'r-');
revh=plot(xx,len_rev_hist,'k-');
legend([fwdh,revh],{'Wave direction','Anti-wave direction'});
xlabel('Number of linked neuron')
ylabel('Total occurance')


% vs shuffle
% jpsth/+bz/+rings/shuffle_conn_bz_alt.m
% wave.COM_chain_shuf(wrs_mux_meta);
load('chains_shuf.mat','shuf_chains')
% TODO match region in recording data
su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
chains_uf.reg=cell(size(chains_uf.cids));
for sess=reshape(unique(chains_uf.sess),1,[])
    sesscid=su_meta.allcid(su_meta.sess==sess);
    sessreg=su_meta.reg_tree(5,su_meta.sess==sess);
    for cnid=reshape(find(chains_uf.sess==sess),1,[])
        [~,supos]=ismember(chains_uf.cids{cnid},sesscid);
        chains_uf.reg{cnid}=sessreg(supos);
    end

end %sess