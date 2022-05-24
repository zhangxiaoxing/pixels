function dec=get_fc_decoding(opt)
arguments
    opt.rpts (1,1) double = 2
    opt.trials (1,1) double {mustBeMember(opt.trials,10:10:50)} = 20
    opt.denovo (1,1) logical = true
    opt.type (1,:) char {mustBeMember(opt.type,{'all','none','both','pre','post'})} = 'all'
    opt.delay (1,1) double {mustBeMember(opt.delay,[3 6])} = 6
    opt.debug (1,1) logical = false
end

if opt.delay==6
    warning('Delay set to default 6')
end


if opt.denovo
    dec=struct();
    dec.s1=[];
    dec.s2=[];
    ubound=[opt.debug,~opt.debug]*[20;163];% reduce debug data size
    
    for fidx=1:ubound
        if ~rem(fidx,20), fprintf('Data %d of 163\n',fidx);end
        fpath=fullfile('K:','code','jpsth','fcdata',sprintf('fc_decoding_f%d.mat',fidx));
        if ~isfile(fpath)
            continue
        end
        load(fpath,'sums','trials','folder');
        trials=behav.procPerf(trials,'all'); % col 9 = WT, col 10 = correct
        
        wtsel=all(trials(:,9:10),2);
        
        su_sums=fc.dec.su_sel(sums,folder,opt.delay,'type',opt.type);
        if isempty(su_sums)
            continue
        end
        
        if opt.delay==6
            sel_S1=wtsel & trials(:,5)==4 & trials(:,8)== 6;
            sel_S2=wtsel & trials(:,5)==8 & trials(:,8)== 6;
        elseif opt.delay==3
            sel_S1=wtsel & trials(:,5)==4 & trials(:,8)== 3;
            sel_S2=wtsel & trials(:,5)==8 & trials(:,8)== 3;
        end
        if min([nnz(sel_S1),nnz(sel_S2)])<opt.trials
            if opt.debug, fprintf('Skip file %d\n',fidx); end
            continue
        end
        
        out=fc.dec.rnd_sample(su_sums,sel_S1,sel_S2,opt.rpts);
        dec.s1=cat(1,dec.s1,out.s1);
        dec.s2=cat(1,dec.s2,out.s2);
    end
    assignin('base','dec',dec);
else
    dec=evalin('base','dec');
end
end