% chains_uf=wave.COM_chain(sel_meta);
% chains_uf_rev=wave.COM_chain(sel_meta,'reverse',true);
% blame=vcs.blame();
% save('chains_mix.mat','chains_uf','chains_uf_rev','blame')

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