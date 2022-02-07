function out=duration_wave(opt)
arguments
    opt.plot (1,1) logical = true
    opt.ctx (1,1) logical = true
end
persistent out_ opt_
if isempty(out_) || ~isequaln(opt,opt_)
    meta=ephys.util.load_meta();
    [~,~,sessmap]=ephys.sessid2path(0);
    homedir=ephys.util.getHomedir('type','raw');
    for sesskey=reshape(cell2mat(sessmap.keys()),1,[])
        disp(sesskey)
        fpath=fullfile(homedir,sessmap(sesskey),"FR_All_ 250.hdf5");
        fr=h5read(fpath,'/FR_All');
        trials=h5read(fpath,'/Trials');
        suid=h5read(fpath,'/SU_id');

        dur_resp=behav.tag_block(trials,'wt',false);
        block_meta=[trials,dur_resp(:,end)];
        sesssel=meta.sess==sesskey;
        if opt.ctx
            regsel=strcmp(meta.reg_tree(2,sesssel),'CTX');
            fr=fr(:,regsel,:);
            suid=suid(regsel);
        end
        data.(sprintf('s%d',sesskey)).fr=fr;
        data.(sprintf('s%d',sesskey)).block_meta=block_meta;
        data.(sprintf('s%d',sesskey)).su_meta=[suid,suid];
        data.(sprintf('s%d',sesskey)).su_meta(:,1)=sesskey;
    end

    [d3fr,d6fr,wrsp,wave_meta]=deal([]);

    for sesskey=reshape(fieldnames(data),1,[])
        disp(sesskey)
        d3sel=data.(sesskey{1}).block_meta(:,8)==3 & all(data.(sesskey{1}).block_meta(:,9:10)==1,2);
        d6sel=data.(sesskey{1}).block_meta(:,8)==6 & all(data.(sesskey{1}).block_meta(:,9:10)==1,2);
        fr=data.(sesskey{1}).fr;
        for su=1:size(fr,2)
            [onewrsp]=deal(nan(1,size(fr,3)./4));
            for bin=1:4:(size(fr,3)-3)
                onebin3=mean(fr(d3sel,su,bin:bin+3),3);
                onebin6=mean(fr(d6sel,su,bin:bin+3),3);
                onewrsp((bin+3)./4)=ranksum(onebin3,onebin6);
                %             onesu3((bin+3)./4)=mean(onebin3);
                %             onesu6((bin+3)./4)=mean(onebin6);
            end
            d3fr=[d3fr;reshape(mean(fr(d3sel,su,:),1),1,[])];
            d6fr=[d6fr;reshape(mean(fr(d6sel,su,:),1),1,[])];
            wrsp=[wrsp;onewrsp];
        end
        wave_meta=[wave_meta;data.(sesskey{1}).su_meta];
    end
    %% selidx
    selidx=(d3fr-d6fr)./(d3fr+d6fr);
    opt_=opt;
    out_=struct();
    out_.wrsp=wrsp;
    out_.meta=wave_meta;
    out_.selidx=selidx;
end
out=out_;

% COM=(selidx*(1:14).')./sum(selidx,2);
% [~,cidx]=sort(COM);
% figure();
% imagesc(selidx(cidx,:),[-1,1])

if opt.plot
    gk = fspecial('gaussian', [3 3], 1);
    sigsel=any(wrsp(:,1:7)<0.05,2);
    figure()
    histogram(sum(wrsp(:,1:7)<0.05,2))

    %% preferred non preferred

    perf3=sigsel & mean(selidx(:,1:28),2)>0;
%     d3allfrn=normalize([d3fr(perf3,1:28),d6fr(perf3,1:28)],2,'range');
%     d3d3frn=d3allfrn(:,1:28);
%     d3d6frn=d3allfrn(:,29:end);
    C=mean([d3fr(perf3,1:28),d6fr(perf3,1:28)],2);
    S=max(abs([d3fr(perf3,1:28),d6fr(perf3,1:28)]-C),[],2);
    d3d3frn=(d3fr(perf3,1:28)-C)./S;
    d3d6frn=(d6fr(perf3,1:28)-C)./S;
    d3pos=d3d3frn;
    d3pos(d3pos<0)=0;
    COM3=(d3pos*(1:28).')./sum(d3pos,2);
    [~,cidx3]=sort(COM3);
    fh3=figure('Color','w','Position',[32,32,750,250]);
    subplot(1,2,1)
    imagesc(conv2(d3d3frn(cidx3,:),gk,'same'),[-1,1])
    set(gca(),'YDir','normal')
    colormap('turbo')
    ylabel('Neuron #')
    set(gca(),'XTick',(4:12:40)+0.5,'XTickLabel',-3:3:6)
    arrayfun(@(x) xline(x,'--k'),[12 16 28 32]+0.5);
    xlabel('Time (s)')
    title('In 3s trials')
    subplot(1,2,2)
    imagesc(conv2(d3d6frn(cidx3,:),gk,'same'),[-1,1])
    set(gca(),'YDir','normal')
    colormap('turbo')
    ylabel('Neuron #')
    set(gca(),'XTick',(4:12:40)+0.5,'XTickLabel',-3:3:6)
    arrayfun(@(x) xline(x,'--k'),[12 16 40 44]+0.5);
    xlabel('Time (s)')
    title('In 6s trials')
    sgtitle('Prefer 3s delay')
    xlim([1.5,28.5])
    exportgraphics(fh3,'duration_3s.pdf','ContentType','vector')

    perf6=sigsel & mean(selidx(:,1:28),2)<0;
%     d6allfrn=normalize([d3fr(perf6,1:28),d6fr(perf6,1:28)],2,'range');
%     d6d3frn=d6allfrn(:,1:28);
%     d6d6frn=d6allfrn(:,29:end);
    C=mean([d3fr(perf6,1:28),d6fr(perf6,1:28)],2);
    S=max(abs([d3fr(perf6,1:28),d6fr(perf6,1:28)]-C),[],2);
    d6d3frn=(d3fr(perf6,1:28)-C)./S;
    d6d6frn=(d6fr(perf6,1:28)-C)./S;
    d6pos=d6d6frn;
    d6pos(d6pos<0)=0;
    COM6=(d6pos*(1:28).')./sum(d6pos,2);

    [~,cidx6]=sort(COM6);
    fh6=figure('Color','w','Position',[32,32,750,250]);
    subplot(1,2,1)
    imagesc(conv2(d6d6frn(cidx6,:),gk,'same'),[-1,1])
    set(gca(),'YDir','normal')
    colormap('turbo')
    ylabel('Neuron #')
    set(gca(),'XTick',(4:12:40)+0.5,'XTickLabel',-3:3:6)
    arrayfun(@(x) xline(x,'--k'),[12 16 40 44]+0.5);
    xlabel('Time (s)')
    title('In 6s trials')
    subplot(1,2,2)
    imagesc(conv2(d6d3frn(cidx6,:),gk,'same'),[-1,1])
    set(gca(),'YDir','normal')
    colormap('turbo')
    set(gca(),'XTick',(4:12:40)+0.5,'XTickLabel',-3:3:6)
    arrayfun(@(x) xline(x,'--k'),[12 16 28 32]+0.5);
    xlabel('Time (s)')
    ylabel('Neuron #')
    title('In 3s trials')
    sgtitle('Prefer 6s delay')
    exportgraphics(fh6,'duration_6s.pdf','ContentType','vector')
end
end