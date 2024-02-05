global_init;
wt_sig=bz.load_sig_sums_conn_file('pair',false,'inhibit',false,'criteria','WT');
ln_sig=bz.load_sig_sums_conn_file('pair',false,'inhibit',false,'criteria','Learning');
nv_sig=bz.load_sig_sums_conn_file('pair',false,'inhibit',false,'criteria','Naive');

critset={'WT','Learning','Naive'};
% 
% categorical({'WT'},critset)
% 
% tmptbl=table('Size',[0,3],'VariableTypes',{'categorical','string','categorical'});
% for ii=1:10000
%     str=critset(rem(ii,3)+1);
%     tmptbl=[tmptbl;cell2table({categorical(str,critset),str{1},categorical(str,{'WT','Learning','Naive'})})];
% end
% 
% aa=tmptbl.Var1;
% bb=tmptbl.Var2;
% cc=tmptbl.Var3;
% whos aa bb cc



outtbl=table('Size',[0,4],'VariableTypes',{'categorical','int32','int32','int32'});
usess=unique(wt_sig.sess);
for sidx=reshape(usess,1,[])
    if rem(sidx,10)==0
        disp("WT "+num2str(sidx))
    end
    ssel=wt_sig.sess==sidx;
    ssuid=wt_sig.suid(ssel,:);
    gh=digraph(table(string(ssuid),'VariableNames',{'EndNodes'}));
    [~,binsize]=gh.conncomp('Type','weak');
    for ii=1:numel(binsize)
        outtbl=[outtbl;cell2table({categorical({'WT'},critset),sidx,ii,binsize(ii)})];
    end
end

usess=unique(ln_sig.sess);
for sidx=reshape(usess,1,[])
    if rem(sidx,10)==0
        disp("Learning "+num2str(sidx))
    end
    ssel=ln_sig.sess==sidx;
    ssuid=ln_sig.suid(ssel,:);
    gh=digraph(table(string(ssuid),'VariableNames',{'EndNodes'}));
    [~,binsize]=gh.conncomp('Type','weak');
    for ii=1:numel(binsize)
        outtbl=[outtbl;cell2table({categorical({'Learning'},critset),sidx,ii,binsize(ii)})];
    end
end

usess=unique(nv_sig.sess);
for sidx=reshape(usess,1,[])
    if rem(sidx,10)==0
        disp("Naive "+num2str(sidx))
    end
    ssel=nv_sig.sess==sidx;
    ssuid=nv_sig.suid(ssel,:);
    gh=digraph(table(string(ssuid),'VariableNames',{'EndNodes'}));
    [~,binsize]=gh.conncomp('Type','weak');
    for ii=1:numel(binsize)
        outtbl=[outtbl;cell2table({categorical({'Naive'},critset),sidx,ii,binsize(ii)})];
    end
end

outtbl.Properties.VariableNames={'Criteria','Session','ConnCompNo','Size'};
[n,edges]=histcounts(outtbl.Size(outtbl.Criteria=='Naive'),0:20:200)
[n,edges]=histcounts(outtbl.Size(outtbl.Criteria=='WT'),0:20:200)
[n,edges]=histcounts(outtbl.Size(outtbl.Criteria=='Learning'),0:20:200)


%%%%%%%%%%%%%%%%%%%%%%%%%%%% sessions %%%%%%%%%%%
wt_su_meta=ephys.util.load_meta("save_file",false,"adjust_white_matter",true,"criteria","WT","load_file",false,"skip_stats",true);
ln_su_meta=ephys.util.load_meta("save_file",false,"adjust_white_matter",true,"criteria","Learning","load_file",false,"skip_stats",true);
nv_su_meta=ephys.util.load_meta("save_file",false,"adjust_white_matter",true,"criteria","Naive","load_file",false,"skip_stats",true);

wt_sel_meta=ephys.get_a2_meta('load_file',false,'save_file',false,'criteria','WT');
ln_sel_meta=ephys.get_a2_meta('load_file',false,'save_file',false,'criteria','Learning');
nv_sel_meta=ephys.get_a2_meta('load_file',false,'save_file',false,'criteria','Naive');

critset={'WT','Learning','Naive'};
outtbl=table('Size',[0,3],'VariableTypes',{'categorical','int32','int32'});
usess=unique(wt_su_meta.sess);
for sidx=reshape(usess,1,[])
    if rem(sidx,10)==0
        disp("WT "+num2str(sidx))
    end
    cnt=nnz(wt_su_meta.sess==sidx);
    outtbl=[outtbl;cell2table({categorical({'WT'},critset),sidx,cnt})];
end

usess=unique(ln_su_meta.sess);
for sidx=reshape(usess,1,[])
    if rem(sidx,10)==0
        disp("Learning "+num2str(sidx))
    end
    cnt=nnz(ln_su_meta.sess==sidx);
    outtbl=[outtbl;cell2table({categorical({'Learning'},critset),sidx,cnt})];
end

usess=unique(nv_su_meta.sess);
for sidx=reshape(usess,1,[])
    if rem(sidx,10)==0
        disp("Naive "+num2str(sidx))
    end
    cnt=nnz(nv_su_meta.sess==sidx);
    outtbl=[outtbl;cell2table({categorical({'Naive'},critset),sidx,cnt})];
end

% wt/ln/nv per session -> sampling su -> build essemble
% sample 100 for session > 200 su
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SAMP_SU=100;
MIN_SESS_CNT=200;
PER_SESS_RPT=5;

critset={'WT','Learning','Naive'};
outtbl=table('Size',[0,4],'VariableTypes',{'categorical','int32','int32','cell'});
usess=intersect(unique(wt_su_meta.sess), unique(wt_sig.sess));
% TODO: investigate sessions with su but no SC

for sidx=reshape(usess,1,[])
    if rem(sidx,10)==0
        disp("WT "+num2str(sidx))
    end
    sesssel=wt_su_meta.sess==sidx;
    cnt=nnz(sesssel);
    if cnt>MIN_SESS_CNT
        for rpt=1:PER_SESS_RPT
            samp_su=randsample(wt_su_meta.allcid(sesssel),SAMP_SU,false);
            sig_sel=wt_sig.sess==sidx;
            sesssig=wt_sig.suid(sig_sel,:);
            samp_sig_sel=all(ismember(sesssig,samp_su),2);
            gh=digraph(string(sesssig(samp_sig_sel,1)),string(sesssig(samp_sig_sel,2)));
            [~,binsize]=gh.conncomp('Type','weak');
            outtbl=[outtbl;cell2table({categorical({'WT'},critset),sidx,rpt,{binsize}})];
        end
    end
end

usess=intersect(unique(ln_su_meta.sess), unique(ln_sig.sess));
for sidx=reshape(usess,1,[])
    if rem(sidx,10)==0
        disp("Learning "+num2str(sidx))
    end
    sesssel=ln_su_meta.sess==sidx;
    cnt=nnz(sesssel);
    if cnt>MIN_SESS_CNT
        for rpt=1:PER_SESS_RPT
            samp_su=randsample(ln_su_meta.allcid(sesssel),SAMP_SU,false);
            sig_sel=ln_sig.sess==sidx;
            sesssig=ln_sig.suid(sig_sel,:);
            samp_sig_sel=all(ismember(sesssig,samp_su),2);
            gh=digraph(string(sesssig(samp_sig_sel,1)),string(sesssig(samp_sig_sel,2)));
            [~,binsize]=gh.conncomp('Type','weak');
            outtbl=[outtbl;cell2table({categorical({'Learning'},critset),sidx,rpt,{binsize}})];
        end
    end
end

usess=intersect(unique(nv_su_meta.sess), unique(nv_sig.sess));
for sidx=reshape(usess,1,[])
    if rem(sidx,10)==0
        disp("Naive "+num2str(sidx))
    end
    sesssel=nv_su_meta.sess==sidx;
    cnt=nnz(sesssel);
    if cnt>MIN_SESS_CNT
        for rpt=1:PER_SESS_RPT
            samp_su=randsample(nv_su_meta.allcid(sesssel),SAMP_SU,false);
            sig_sel=nv_sig.sess==sidx;
            sesssig=nv_sig.suid(sig_sel,:);
            samp_sig_sel=all(ismember(sesssig,samp_su),2);
            gh=digraph(string(sesssig(samp_sig_sel,1)),string(sesssig(samp_sig_sel,2)));
            [~,binsize]=gh.conncomp('Type','weak');
            outtbl=[outtbl;cell2table({categorical({'Naive'},critset),sidx,rpt,{binsize}})];
        end
    end
end
outtbl.Properties.VariableNames={'Criteria','Session','rpt','Size'};

ranksum(cell2mat(outtbl.Size(outtbl.Criteria=='WT').'),cell2mat(outtbl.Size(outtbl.Criteria=='Learning').'))
ranksum(cell2mat(outtbl.Size(outtbl.Criteria=='WT').'),cell2mat(outtbl.Size(outtbl.Criteria=='Naive').'))
ranksum(cell2mat(outtbl.Size(outtbl.Criteria=='Naive').'),cell2mat(outtbl.Size(outtbl.Criteria=='Learning').'))


[h,p]=ttest2(cell2mat(outtbl.Size(outtbl.Criteria=='WT').'),cell2mat(outtbl.Size(outtbl.Criteria=='Learning').'))
[h,p]=ttest2(cell2mat(outtbl.Size(outtbl.Criteria=='WT').'),cell2mat(outtbl.Size(outtbl.Criteria=='Naive').'))
[h,p]=ttest2(cell2mat(outtbl.Size(outtbl.Criteria=='Naive').'),cell2mat(outtbl.Size(outtbl.Criteria=='Learning').'))

ranksum(cellfun(@(x) sum(x), outtbl.Size(outtbl.Criteria=='WT')),cellfun(@(x) sum(x), outtbl.Size(outtbl.Criteria=='Learning')))
ranksum(cellfun(@(x) sum(x), outtbl.Size(outtbl.Criteria=='WT')),cellfun(@(x) sum(x), outtbl.Size(outtbl.Criteria=='Naive')))
ranksum(cellfun(@(x) sum(x), outtbl.Size(outtbl.Criteria=='Learning')),cellfun(@(x) sum(x), outtbl.Size(outtbl.Criteria=='Naive')))


ranksum(cellfun(@(x) numel(x), outtbl.Size(outtbl.Criteria=='WT')),cellfun(@(x) numel(x), outtbl.Size(outtbl.Criteria=='Learning')))
ranksum(cellfun(@(x) numel(x), outtbl.Size(outtbl.Criteria=='WT')),cellfun(@(x) numel(x), outtbl.Size(outtbl.Criteria=='Naive')))
ranksum(cellfun(@(x) numel(x), outtbl.Size(outtbl.Criteria=='Learning')),cellfun(@(x) numel(x), outtbl.Size(outtbl.Criteria=='Naive')))


mean(cellfun(@(x) numel(x), outtbl.Size(outtbl.Criteria=='Naive')))
mean(cellfun(@(x) numel(x), outtbl.Size(outtbl.Criteria=='Learning')))
mean(cellfun(@(x) numel(x), outtbl.Size(outtbl.Criteria=='WT')))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

usess=unique(ln_su_meta.sess);
for sidx=reshape(usess,1,[])
    if rem(sidx,10)==0
        disp("Learning "+num2str(sidx))
    end
    cnt=nnz(ln_su_meta.sess==sidx);
    outtbl=[outtbl;cell2table({categorical({'Learning'},critset),sidx,cnt})];
end

usess=unique(nv_su_meta.sess);
for sidx=reshape(usess,1,[])
    if rem(sidx,10)==0
        disp("Naive "+num2str(sidx))
    end
    cnt=nnz(nv_su_meta.sess==sidx);
    outtbl=[outtbl;cell2table({categorical({'Naive'},critset),sidx,cnt})];
end




