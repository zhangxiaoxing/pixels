if false
    load('sums_conn.mat','sums_conn_str');
    sig=bz.load_sig_pair();
    idmap=load(fullfile('K:','code','align','reg_ccfid_map.mat'));
    com_map=wave.get_com_map('peak',false,'per_sec_stats',false);
    %caution, parallel index methods
    sumcon_sess=cellfun(@(x) ephys.path2sessid(...
        replace(regexp(x,'(?<=SPKINFO/).*$','match','once'),'/',filesep())...
        ),{sums_conn_str.folder}.');
    sumcon_allsess=cell2mat(arrayfun(@(x) repmat(sumcon_sess(x),...
        1,size(sums_conn_str(x).sig_con,1)),1:numel(sums_conn_str),'UniformOutput',false)).';
    sumcon_ccg=cell2mat({sums_conn_str.ccg_sc}.');
    sumcon_suids=cell2mat({sums_conn_str.sig_con}.');
    sumcon_qc=cell2mat({sums_conn_str.qc}.');
end
strict_sel=sumcon_qc(:,1)>0 & sumcon_qc(:,4)<25 & sumcon_qc(:,3)<4 & qc(:,2)<276 & qc(:,2)>251;
for sumcon_idx=find(strict_sel).'
    %transient, congruent, COM shift, cross region, good_wf
    sumconn_currsess=sumcon_allsess(sumcon_idx);
    sig_idx=sig.sess==sumconn_currsess & sig.suid(:,1)==sumcon_suids(sumcon_idx,1) & sig.suid(:,2)==sumcon_suids(sumcon_idx,2);
    if nnz(sig_idx)~=1, continue; end
    if ~all(sig.wf_good(sig_idx,:)), continue; end
    if ~all(ismember(sig.reg(sig_idx,1,:),[567,343])), continue; end %grey matter
    if any(sig.reg(sig_idx,5,:)==0) || diff(sig.reg(sig_idx,5,:))==0,continue;end
    if all(sig.mem_type(sig_idx,:)==2)
        samp='s1';
        rsamp=4;
    elseif all(sig.mem_type(sig_idx,:)==4)
        samp='s2';
        rsamp=8;
    else
        continue
    end
    %     keyboard
    if com_map.(['s',num2str(sumconn_currsess)]).(samp)(sumcon_suids(sumcon_idx,2))...
            -com_map.(['s',num2str(sumconn_currsess)]).(samp)(sumcon_suids(sumcon_idx,1))<0
        continue
    end
    fh=figure('Color','w');
    subplot(2,4,5)
    hold on;
    plot(sumcon_ccg(sumcon_idx,:),'-r');
    xlim([1,501])
    text(min(xlim()),max(ylim()),idmap.ccfid2reg(sig.reg(sig_idx,5,1)),'HorizontalAlignment','left','VerticalAlignment','top')
    text(max(xlim()),max(ylim()),idmap.ccfid2reg(sig.reg(sig_idx,5,2)),'HorizontalAlignment','right','VerticalAlignment','top')
    arrayfun(@(x) xline(x,'--k'),[251,251-25,251+25]);
    title(sprintf('%d,%d,%d,%d,%s',sumcon_idx,sumcon_allsess(sumcon_idx),sumcon_suids(sumcon_idx,1),sumcon_suids(sumcon_idx,2),samp));
    subplot(2,1,1)
    hold on
    bz.raster(sumconn_currsess,sumcon_suids(sumcon_idx,:),rsamp);
    waitfor(fh);
end



