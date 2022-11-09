% obsolete code
% abandoned as of 22.11.09

function dec=get_fc_decoding(sel_meta,opt)

arguments
    sel_meta
    opt.rpts (1,1) double = 2
    opt.trials (1,1) double {mustBeMember(opt.trials,10:10:50)} = 20
    opt.denovo (1,1) logical = true
    opt.type (1,:) char {mustBeMember(opt.type,{'all','none','both','pre','post'})} = 'all'
    opt.debug (1,1) logical = false
end

if opt.delay==6
    warning('Delay set to default 6')
end


if opt.denovo
    dec=struct();
    dec.s1=[];
    dec.s2=[];

    fcfl=dir(fullfile('fccoding','fc_coding_*.mat'));
    for fidx=1:numel(fcfl)
        disp(fidx);
        fpath=fullfile(fcfl(fidx).folder,fcfl(fidx).name);
        fstr=load(fpath,'onesess');
        trials=fstr.onesess.trials;
%         behav.procPerf(trials,'mode','all'); % col 9 = WT, col 10 = correct
%         wtsel=all(trials(:,9:10),2);
        
        sel_S1=wtsel & trials(:,5)==4 & trials(:,8)== 3;
        sel_S2=wtsel & trials(:,5)==8 & trials(:,8)== 3;

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