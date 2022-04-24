function fc_com_consistency_stepwise(opt)
arguments
    opt.level (1,1) double {mustBeInteger} = 5 %allen ccf structure detail level
    opt.decision (1,1) logical = false % return statistics of decision period, default is delay period
    opt.hiermap (1,:) char {mustBeMember(opt.hiermap,{'pvsst','OBM1','AON'})} = 'AON'
    opt.mem_type (1,:) char {mustBeMember(opt.mem_type,{'congru','incong','any'})} = 'congru' %TCOM for NONMEM undefined
    opt.new_data (1,1) logical = false
end



com_map_6=wave.get_com_map('wave','any6','delay',6);
com_map_3=wave.get_com_map('wave','any3','delay',3);
com_map_d=wave.get_dur_com_map();


sig=bz.load_sig_pair('type','neupix','prefix','BZWT','criteria','WT','load_waveid',true);
% hiersu=sort(hierv(hierv>0));
if opt.new_data
    fc_com_pvsst_stats=[];
    usess=unique(sig.sess);
    for sii=reshape(usess,1,[]) %iterate through sessions
        sesssel=sig.sess==sii;
        suid=sig.suid(sesssel,:);
        waveid=sig.waveid(sesssel,:);

        [comsess_sens_6,comsess_sens_3,comsess_dur_3,comsess_dur_6]=deal(nan(size(suid)));
        if isfield(com_map_6,['s',num2str(sii)])
            s1keys=int32(cell2mat(com_map_6.(['s',num2str(sii)]).s1.keys()));
            s1sel=ismember(suid,s1keys);
            comsess_sens_6(s1sel)=cell2mat(com_map_6.(['s',num2str(sii)]).s1.values(num2cell(suid(s1sel))));
            s2keys=int32(cell2mat(com_map_6.(['s',num2str(sii)]).s2.keys()));
            s2sel=ismember(suid,s2keys);
            comsess_sens_6(s2sel)=cell2mat(com_map_6.(['s',num2str(sii)]).s2.values(num2cell(suid(s2sel))));
        end

        if isfield(com_map_3,['s',num2str(sii)])
            s1keys=int32(cell2mat(com_map_3.(['s',num2str(sii)]).s1.keys()));
            s1sel=ismember(suid,s1keys);
            comsess_sens_3(s1sel)=cell2mat(com_map_3.(['s',num2str(sii)]).s1.values(num2cell(suid(s1sel))));
            s2keys=int32(cell2mat(com_map_3.(['s',num2str(sii)]).s2.keys()));
            s2sel=ismember(suid,s2keys);
            comsess_sens_3(s2sel)=cell2mat(com_map_3.(['s',num2str(sii)]).s2.values(num2cell(suid(s2sel))));
        end


        if isfield(com_map_d,['s',num2str(sii)])
            d3keys=int32(cell2mat(com_map_d.(['s',num2str(sii)]).d3.keys()));
            d3sel=ismember(suid,d3keys);
            comsess_dur_3(d3sel)=cell2mat(com_map_d.(['s',num2str(sii)]).d3.values(num2cell(suid(d3sel))));
            d6keys=int32(cell2mat(com_map_d.(['s',num2str(sii)]).d6.keys()));
            d6sel=ismember(suid,d6keys);
            comsess_dur_6(d6sel)=cell2mat(com_map_d.(['s',num2str(sii)]).d6.values(num2cell(suid(d6sel))));
        end

        regsess=squeeze(sig.reg(sesssel,5,:));
        fc_com_pvsst_stats=[fc_com_pvsst_stats;double(sii).*ones(size(suid(:,1))),double(suid),comsess_sens_6,comsess_sens_3,double(regsess),double(waveid),comsess_dur_6,comsess_dur_3];
        %========================================================1[sess]=============2:3============4:5==============6:7===========8:9============10:11========12:13=========14:15========
    end
    save('fc_com_pvsst_stats.mat','fc_com_pvsst_stats');
else
    load('fc_com_pvsst_stats.mat','fc_com_pvsst_stats');
end
switch opt.mem_type
    case 'congru'
        congrusel=all(ismember(fc_com_pvsst_stats(:,10:11),[1 3 5]),2) |all(ismember(fc_com_pvsst_stats(:,10:11),[2 4 6]),2);
    case 'incong'
        congrusel=any(ismember(fc_com_pvsst_stats(:,10:11),[1 3 5]),2) & any(ismember(fc_com_pvsst_stats(:,10:11),[2 4 6]),2);
end

idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
CTX_sel=squeeze(sig.reg(:,2,:)==idmap.reg2ccfid('CTX'));
CNU_sel=squeeze(sig.reg(:,2,:)==idmap.reg2ccfid('CNU'));
BS_sel=squeeze(sig.reg(:,1,:)==idmap.reg2ccfid('BS'));


% COM vs. hier
[mm_sens_3reg,p_sens_3reg,mm_dur_3reg,p_dur_3reg]=deal(nan(3,3));
[cnt_sens_3reg,ci_sens_3reg,cnt_dur_3reg,ci_dur_3reg]=deal(nan(3,3,2));

reg_cell={CTX_sel,CNU_sel,BS_sel};

% sensory
for fromidx=1:3
    for toidx=1:3
        %         disp([fromidx,toidx])
        from_sel=reg_cell{fromidx}(:,1);
        to_sel=reg_cell{toidx}(:,2);

        com_3reg=(fc_com_pvsst_stats(from_sel & to_sel  & congrusel,4:7));
        early2late_cnt=nnz([com_3reg(:,2);com_3reg(:,4)]>[com_3reg(:,1);com_3reg(:,3)]);
        late2early_cnt=nnz([com_3reg(:,2);com_3reg(:,4)]<[com_3reg(:,1);com_3reg(:,3)]);
        disp([early2late_cnt,late2early_cnt])
        cnt_sens_3reg(fromidx,toidx,:)=[early2late_cnt,late2early_cnt];
        [mm_sens_3reg(fromidx,toidx),ci_sens_3reg(fromidx,toidx,:)]=binofit(early2late_cnt,early2late_cnt+late2early_cnt);
        p_sens_3reg(fromidx,toidx)=2*binocdf(min(early2late_cnt,late2early_cnt),early2late_cnt+late2early_cnt,0.5);
    end
end

% duration
for fromidx=1:3
    for toidx=1:3
        %         disp([fromidx,toidx])
        from_sel=reg_cell{fromidx}(:,1);
        to_sel=reg_cell{toidx}(:,2);

        com_3reg=(fc_com_pvsst_stats(from_sel & to_sel,12:15));
        early2late_cnt=nnz([com_3reg(:,2);com_3reg(:,4)]>[com_3reg(:,1);com_3reg(:,3)]);
        late2early_cnt=nnz([com_3reg(:,2);com_3reg(:,4)]<[com_3reg(:,1);com_3reg(:,3)]);
        disp([early2late_cnt,late2early_cnt])
        cnt_dur_3reg(fromidx,toidx,:)=[early2late_cnt,late2early_cnt];
        [mm_dur_3reg(fromidx,toidx),ci_dur_3reg(fromidx,toidx,:)]=binofit(early2late_cnt,early2late_cnt+late2early_cnt);
        p_dur_3reg(fromidx,toidx)=2*binocdf(min(early2late_cnt,late2early_cnt),early2late_cnt+late2early_cnt,0.5);
    end
end

keyboard()
