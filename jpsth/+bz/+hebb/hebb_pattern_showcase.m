if ~exist('sums_all','var')
    load(fullfile('binary','sums_ring_stats_all.mat'),'sums_all');
end
load(fullfile('bzdata','rings_bz_vs_shuf.mat'),'rings');
if ~exist('sums_conn_str','var')
    load(fullfile("binary","sums_conn_10.mat"),'sums_conn_str');
    [sig,~]=bz.load_sig_sums_conn_file('pair',false);
    idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
    %caution, parallel index methods
    sumcon_sess=cellfun(@(x) ephys.path2sessid(...
        replace(regexp(x,'(?<=SPKINFO/).*$','match','once'),'/',filesep())...
        ),{sums_conn_str.folder}.');
    sumcon_allsess=cell2mat(arrayfun(@(x) repmat(sumcon_sess(x),...
        1,size(sums_conn_str(x).sig_con,1)),1:numel(sums_conn_str),'UniformOutput',false)).';
    sumcon_ccg=cell2mat({sums_conn_str.ccg_sc}.');
    sumcon_suids=cell2mat({sums_conn_str.sig_con}.');
    sumcon_qc=cell2mat({sums_conn_str.qc}.');
    strict_sel=sumcon_qc(:,2)>=252 & sumcon_qc(:,2)<=265 & sumcon_qc(:,4)>=2 & sumcon_qc(:,4)<=40 & sumcon_qc(:,5)>248;
    sumcon_suids=sumcon_suids(strict_sel,:);
    sumcon_allsess=sumcon_allsess(strict_sel,:);
    %1:Polar 2:Time of peak 3:Noise peaks 4:FWHM 5:rising edge
    %6:falling edge
end
load('hebb_pattern_bz.mat','hebbPattern');
load(fullfile('binary','su_meta.mat'),'su_meta');
fstr=load(fullfile('binary','wrs_mux_meta.mat'));
sel_meta=fstr.wrs_mux_meta;
clear fstr

usess=unique(cell2mat(hebbPattern(:,1)));
id2reg=@(currid,allid,reg) arrayfun(@(x) reg(allid==x),currid);
ring2list=@(in) [cell2mat(arrayfun(@(x) in(x:x+1),(1:numel(in)-1).','UniformOutput',false));in(end),in(1)];

greys=ephys.getGreyRegs('range','grey','mincount',0);
congruset={[1 5],[2 5],[3 6],[4 6]};

for sess=18%reshape(usess,1,[])
    % sufficient loop activity
    activ3=cell2mat(sums_all{1}(cell2mat(sums_all{1}(:,1))==sess,3));
    tags3=sums_all{1}(cell2mat(sums_all{1}(:,1))==sess,5);
    activ4=cell2mat(sums_all{2}(cell2mat(sums_all{2}(:,1))==sess,3));
    tags4=sums_all{2}(cell2mat(sums_all{2}(:,1))==sess,5);


    for selidx=1:4
        seltype=congruset{selidx};
        %% CTX & Congru
        ctx_congru_sel=su_meta.sess==sess & ismember(sel_meta.wave_id,seltype) & ismember(su_meta.reg_tree(5,:).',greys);
        sess_cid=su_meta.allcid(ctx_congru_sel);
        sess_reg=su_meta.reg_tree(5,ctx_congru_sel);
        sessplist=hebbPattern(cell2mat(hebbPattern(:,1))==sess,2:4);

        sess_ring_sel=arrayfun(@(x) all(ismember(cell2mat(sessplist(x,:)),sess_cid),'all'),1:size(sessplist,1));
        sessplist=sessplist(sess_ring_sel,:);
        if isempty(sessplist)
            continue
        end
        %% multi-region

        ring_reg=arrayfun(@(x) unique(id2reg(cell2mat(sessplist(x,:)),sess_cid,sess_reg)),1:size(sessplist,1),'UniformOutput',false);
        multi_reg_sel=cellfun(@(x) numel(x)>2, ring_reg);
        sessplist=sessplist(multi_reg_sel,:);
        ring_reg=ring_reg(multi_reg_sel);
        %% sufficient loop activity
        for ii=1:size(sessplist,1)
            disp([selidx,ii,size(sessplist,1)])
            disp(sessplist(ii,:))
            disp(subsref(arrayfun(@(x) id2reg(cell2mat(sessplist(x,:)),sess_cid,sess_reg),ii,'UniformOutput',false),substruct('{}',{':'})))
            r1idx=all(ismember(activ3,sessplist{ii,1}),2);
            r2idx=all(ismember(activ4,sessplist{ii,2}),2);
            r3idx=all(ismember(activ4,sessplist{ii,3}),2);
            if any(r1idx) && any(r2idx) && any(r3idx)
                patt_su=unique(cell2mat(sessplist(ii,1:3)));
% for generalized patterns                
%                 all3_sel=all(ismember(activ3,patt_su),2);
%                 all4_sel=all(ismember(activ4,patt_su),2);
% %% specific hebb
                all3_sel=r1idx;
                all4_sel=r2idx | r3idx;

                bz.hebb.plotone(sess,activ3,tags3,all3_sel,activ4,tags4,all4_sel,selidx,ii,ring_reg{ii});
            end
        end
    end
end



function ccg_selec(sessplist)
%% typical ccg waveform % diffcult to meet criteria
out=cell(0);
for ii=1:size(sessplist,1)
    conn=[ring2list(sessplist{ii,1});ring2list(sessplist{ii,2});ring2list(sessplist{ii,3})];
    typical_ccg=arrayfun(@(x) any(sumcon_allsess==sess & sumcon_suids(:,1)==conn(x,1) & sumcon_suids(:,2)==conn(x,2)),1:size(conn,1));
    disp([nnz(typical_ccg),numel(typical_ccg)]);
    if nnz(typical_ccg)<numel(typical_ccg)
        continue;
    end
    out=[out;sessplist(ii,:)];
end
end
