[sig,~]=bz.load_sig_sums_conn_file('pair',false);
fcsel=(sig.reg(:,5,1)==idmap.reg2ccfid('TT') & sig.reg(:,5,2)==idmap.reg2ccfid('DP'));
sess=sig.sess(fcsel);
cids=sig.suid(fcsel,:);
load('sums_conn_10.mat','sums_conn_str');

ccgsess=cellfun(@(x) ephys.path2sessid(x),{sums_conn_str.folder});

for ii=1:numel(sess)
    [spkID,~,~,~,~,~]=ephys.getSPKID_TS(sess(ii),'keep_trial',false);
    lead_cnt=nnz(spkID==cids(ii));

    ccgsel=ccgsess==sess(ii);
    ccgsuid=sums_conn_str(ccgsel).sig_con;
    susel=ccgsuid(:,1)==cids(ii,1) & ccgsuid(:,2)==cids(ii,2);
    figure()
    plot(sums_conn_str(ccgsel).ccg_sc(susel,:)./lead_cnt./0.0004)
    xlim(251+[-31,60])
    xline(251,'--r')
    title("S"+sess(ii)+"LN"+cids(ii,1)+"FN"+cids(ii,2))
    ylabel('Normalized FR (Hz)')
    grh=groot;
    if numel(grh.Children)>19
        keyboard()
    end
end

