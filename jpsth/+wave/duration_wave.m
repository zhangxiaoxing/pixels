
[~,~,sessmap]=ephys.sessid2path(0);
meta=ephys.util.load_meta();
homedir=ephys.util.getHomedir('type','raw');
for ii=reshape(cell2mat(sessmap.keys()),1,[])
    disp(ii)
    fpath=fullfile(homedir,sessmap(ii),"FR_All_ 250.hdf5");
    fr=h5read(fpath,'/FR_All');
    trials=h5read(fpath,'/Trials');
    suid=h5read(fpath,'/SU_id');

    dur_resp=behav.tag_block(trials,'wt',false);
    block_meta=[trials,dur_resp(:,end)];
    sesssel=meta.sess==ii;
    regsel=strcmp(meta.reg_tree(2,sesssel),'CTX');

    fr=fr(:,regsel,:);

    out.(sprintf('s%d',ii)).fr=fr;
    out.(sprintf('s%d',ii)).block_meta=block_meta;
end

[d3fr,d6fr,wrsp,meta]=deal([]);

for sesskey=reshape(fieldnames(out),1,[])
    disp(sesskey)
    d3sel=out.(sesskey{1}).block_meta(:,8)==3 & all(out.(sesskey{1}).block_meta(:,9:10)==1,2);
    d6sel=out.(sesskey{1}).block_meta(:,8)==6 & all(out.(sesskey{1}).block_meta(:,9:10)==1,2);
    fr=out.(sesskey{1}).fr;
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
        meta=[meta;str2double(regexp(sesskey{1},'\d*','match','once')),su];
    end
end
keyboard()
%% selidx
selidx=(d3fr-d6fr)./(d3fr+d6fr);
% COM=(selidx*(1:14).')./sum(selidx,2);
% [~,cidx]=sort(COM);
% figure();
% imagesc(selidx(cidx,:),[-1,1])

gk = fspecial('gaussian', [3 3], 1);
sigsel=any(wrsp(:,1:7)<0.05,2);

%% preferred non preferred
perf3=sigsel & mean(selidx(:,1:28),2)>0;
d3allfrn=normalize([d3fr(perf3,:),d6fr(perf3,:)],2,'range');
d3d3frn=d3allfrn(:,1:56);
d3d6frn=d3allfrn(:,57:end);
d3pos=d3d3frn;
d3pos(d3pos<0)=0;
COM3=(d3pos*(1:56).')./sum(d3pos,2);
[~,cidx3]=sort(COM3);
figure('Color','w');
subplot(1,2,1)
imagesc(conv2(d3d3frn(cidx3,:),gk,'same'),[0,0.8])
set(gca(),'YDir','normal')
colormap('turbo')
ylabel('Neuron #')
set(gca(),'XTick',(4:12:40)+0.5,'XTickLabel',-3:3:6)
arrayfun(@(x) xline(x,'--k'),[12 16 28 32]+0.5);
xlabel('Time (s)')
title('In 3s trials')
subplot(1,2,2)
imagesc(conv2(d3d6frn(cidx3,:),gk,'same'),[0,0.8])
set(gca(),'YDir','normal')
colormap('turbo')
ylabel('Neuron #')
set(gca(),'XTick',(4:12:40)+0.5,'XTickLabel',-3:3:6)
arrayfun(@(x) xline(x,'--k'),[12 16 40 44]+0.5);
xlabel('Time (s)')
title('In 6s trials')
sgtitle('Prefer 3s delay')

perf6=sigsel & mean(selidx(:,1:28),2)<0;
d6allfrn=normalize([d3fr(perf6,:),d6fr(perf6,:)],2,'range');
d6d3frn=d6allfrn(:,1:56);
d6d6frn=d6allfrn(:,57:end);
d6pos=d6d6frn;
d6pos(d6pos<0)=0;
COM6=(d6pos*(1:56).')./sum(d6pos,2);
[~,cidx6]=sort(COM6);
figure('Color','w');
subplot(1,2,1)
imagesc(conv2(d6d6frn(cidx6,:),gk,'same'),[0,0.8])
set(gca(),'YDir','normal')
colormap('turbo')
ylabel('Neuron #')
set(gca(),'XTick',(4:12:40)+0.5,'XTickLabel',-3:3:6)
arrayfun(@(x) xline(x,'--k'),[12 16 40 44]+0.5);
xlabel('Time (s)')
title('In 6s trials')
subplot(1,2,2)
imagesc(conv2(d6d3frn(cidx6,:),gk,'same'),[0,0.8])
set(gca(),'YDir','normal')
colormap('turbo')
set(gca(),'XTick',(4:12:40)+0.5,'XTickLabel',-3:3:6)
arrayfun(@(x) xline(x,'--k'),[12 16 28 32]+0.5);
xlabel('Time (s)')
ylabel('Neuron #')
title('In 3s trials')
sgtitle('Prefer 6s delay')