function fc_com_consistency_stepwise(opt)
arguments
    opt.level (1,1) double {mustBeInteger} = 5 %allen ccf structure detail level
    opt.decision (1,1) logical = false % return statistics of decision period, default is delay period
    opt.hiermap (1,:) char {mustBeMember(opt.hiermap,{'pvsst','OBM1','AON'})} = 'AON'
    opt.mem_type (1,:) char {mustBeMember(opt.mem_type,{'congru','incong','any'})} = 'congru' %TCOM for NONMEM undefined
end


% [~,com_meta]=wave.per_region_COM('decision',opt.decision,'wave',opt.wave);
% com_map=list2map(com_meta);
if strcmp(opt.hiermap,'AON')
    com_map_6=wave.get_com_map('wave','any6','delay',6);
    com_map_3=wave.get_com_map('wave','any3','delay',3);
else
    %TODO: LAT
end

sig=bz.load_sig_pair('type','neupix','prefix','BZWT','criteria','WT','load_waveid',true);
[~,is_same,LAT2AON,AON2LAT,hierv]=bz.util.diff_at_level(sig.reg,'hierarchy',true,'range','grey','hiermap',opt.hiermap,'mincount',0);
% hiersu=sort(hierv(hierv>0));
uv=unique(hierv(all(hierv>0,2)));
minv=min(hierv(all(hierv>0,2)),[],2);
maxv=max(hierv(all(hierv>0,2)),[],2);
opti=[];
for t1=numel(uv)-1:-1:3
    for t2=t1-1:-1:2
        block1=nnz(minv>uv(t1));
        block2=nnz(minv>uv(t2) &maxv<=uv(t1));
        block3=nnz(maxv<=uv(t2));
        opti=[opti;t1,t2,block1,block2,block3,block1.*block2.*block3];
    end
end
[~,loc]=max(opti(:,6));

l_bdry=uv(opti(loc,1));
h_bdry=uv(opti(loc,2));
disp(nnz(is_same | AON2LAT | LAT2AON))
fc_com_pvsst_stats=[];
usess=unique(sig.sess);
for sii=reshape(usess,1,[]) %iterate through sessions
    sesssel=sig.sess==sii;
    suid=sig.suid(sesssel,:);
    waveid=sig.waveid(sesssel,:);
    comsess_6=nan(size(suid));
    if isfield(com_map_6,['s',num2str(sii)])
        s1keys=int32(cell2mat(com_map_6.(['s',num2str(sii)]).s1.keys()));
        s1sel=ismember(suid,s1keys);
        comsess_6(s1sel)=cell2mat(com_map_6.(['s',num2str(sii)]).s1.values(num2cell(suid(s1sel))));
        s2keys=int32(cell2mat(com_map_6.(['s',num2str(sii)]).s2.keys()));
        s2sel=ismember(suid,s2keys);
        comsess_6(s2sel)=cell2mat(com_map_6.(['s',num2str(sii)]).s2.values(num2cell(suid(s2sel))));
    end
    comsess_3=nan(size(suid));
    if isfield(com_map_3,['s',num2str(sii)])
        s1keys=int32(cell2mat(com_map_3.(['s',num2str(sii)]).s1.keys()));
        s1sel=ismember(suid,s1keys);
        comsess_3(s1sel)=cell2mat(com_map_3.(['s',num2str(sii)]).s1.values(num2cell(suid(s1sel))));
        s2keys=int32(cell2mat(com_map_3.(['s',num2str(sii)]).s2.keys()));
        s2sel=ismember(suid,s2keys);
        comsess_3(s2sel)=cell2mat(com_map_3.(['s',num2str(sii)]).s2.values(num2cell(suid(s2sel))));
    end

    regsess=squeeze(sig.reg(sesssel,5,:));
    fc_com_pvsst_stats=[fc_com_pvsst_stats;double(sii).*ones(size(suid(:,1))),double(suid),comsess_6,comsess_3,double(regsess),double(waveid)];
    %========================================================1===================2:3==========4:5=======6:7=========8:9==============10:11====
end
save('fc_com_pvsst_stats.mat','fc_com_pvsst_stats');
switch opt.mem_type
    case 'congru'
        congrusel=all(ismember(fc_com_pvsst_stats(:,10:11),[1 3 5]),2) |all(ismember(fc_com_pvsst_stats(:,10:11),[2 4 6]),2);
    case 'incong'
        congrusel=any(ismember(fc_com_pvsst_stats(:,10:11),[1 3 5]),2) & any(ismember(fc_com_pvsst_stats(:,10:11),[2 4 6]),2);
end

% COM vs. hier
[mmsame,mmAON2LAT,mmLAT2AON,psame,pAON2LAT,pLAT2AON]=deal(nan(3,3));
[samecnt,AON2LATcnt,LAT2AONcnt,sameci,AON2LATci,LAT2AONci]=deal(nan(3,3,2));

boundry=[max(hierv,[],'all'),l_bdry,h_bdry,0];

% same
for idx=1:3
    stepsel=max(hierv,[],2)<=boundry(idx) & min(hierv,[],2)>boundry(idx+1);
    comsame=(fc_com_pvsst_stats(stepsel & is_same(:,opt.level) & congrusel,4:7)); % COM
    early2late_cnt=nnz([comsame(:,2);comsame(:,4)]>[comsame(:,1);comsame(:,3)]); % time early to late
    late2early_cnt=nnz([comsame(:,2);comsame(:,4)]<[comsame(:,1);comsame(:,3)]); % time late to early
    samecnt(idx,idx,:)=[early2late_cnt,late2early_cnt];
    [mmsame(idx,idx),sameci(idx,idx,:)]=binofit(early2late_cnt,early2late_cnt+late2early_cnt);
    psame(idx,idx)=2*binocdf(min(early2late_cnt,late2early_cnt),early2late_cnt+late2early_cnt,0.5);
end
%AON2LAT
for fromidx=1:3
    for toidx=1:3
%         disp([fromidx,toidx])
        from_sel=hierv(:,1)<=boundry(fromidx) & hierv(:,1)>boundry(fromidx+1);
        to_sel=hierv(:,2)<=boundry(toidx) & hierv(:,2)>boundry(toidx+1);
        
        comAON2LAT=(fc_com_pvsst_stats(from_sel & to_sel & AON2LAT(:,opt.level)  & congrusel,4:7));
        early2late_cnt=nnz([comAON2LAT(:,2);comAON2LAT(:,4)]>[comAON2LAT(:,1);comAON2LAT(:,3)]);
        late2early_cnt=nnz([comAON2LAT(:,2);comAON2LAT(:,4)]<[comAON2LAT(:,1);comAON2LAT(:,3)]);
        disp([early2late_cnt,late2early_cnt])
        AON2LATcnt(fromidx,toidx,:)=[early2late_cnt,late2early_cnt];
        [mmAON2LAT(fromidx,toidx),AON2LATci(fromidx,toidx,:)]=binofit(early2late_cnt,early2late_cnt+late2early_cnt);
        pAON2LAT(fromidx,toidx)=2*binocdf(min(early2late_cnt,late2early_cnt),early2late_cnt+late2early_cnt,0.5);   
    end
end

%LAT2AON
for fromidx=1:3
    for toidx=1:3
%         disp([fromidx,toidx])
        from_sel=hierv(:,1)<=boundry(fromidx) & hierv(:,1)>boundry(fromidx+1);
        to_sel=hierv(:,2)<=boundry(toidx) & hierv(:,2)>boundry(toidx+1);
        
        comLAT2AON=(fc_com_pvsst_stats(from_sel & to_sel & LAT2AON(:,opt.level)  & congrusel,4:7));
        early2late_cnt=nnz([comLAT2AON(:,2);comLAT2AON(:,4)]>[comLAT2AON(:,1);comLAT2AON(:,3)]);
        late2early_cnt=nnz([comLAT2AON(:,2);comLAT2AON(:,4)]<[comLAT2AON(:,1);comLAT2AON(:,3)]);
        disp([early2late_cnt,late2early_cnt])
        LAT2AONcnt(fromidx,toidx,:)=[early2late_cnt,late2early_cnt];
        [mmLAT2AON(fromidx,toidx),LAT2AONci(fromidx,toidx,:)]=binofit(early2late_cnt,early2late_cnt+late2early_cnt);
        pLAT2AON(fromidx,toidx)=2*binocdf(min(early2late_cnt,late2early_cnt),early2late_cnt+late2early_cnt,0.5);   
    end
end

keyboard()
