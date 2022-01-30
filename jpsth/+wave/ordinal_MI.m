%TODO brain region filter, olfaction filter.
% function out=get_ordinal_mat(opt)
%% gen data
[~,~,sessmap]=ephys.sessid2path(0);
meta=ephys.util.load_meta();
homedir=ephys.util.getHomedir('type','raw');
mi_all=[];
anovanp=[];
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

    for suidx=1:size(fr,2)
        [mi,ap]=deal(nan(1,14));
        for bin=1:14
            win=(bin*4)-3:bin*4;
            [xx,yy]=deal([]);
            for ord=int32([31:34,61:64])
                dur=idivide(ord,10);
                dur_ord=rem(ord,10);
                trl_sel=trials(:,8)==dur & dur_resp(:,3)==dur_ord & all(trials(:,9:10)>0,2);
                dx=mean(fr(trl_sel,suidx,win),3);
                xx=[xx;dx];
                yy=[yy;double(ord).*ones(size(dx))];
            end
            mi(bin)=discrete_continuous_info_fast(yy,xx,5,2);
            ap(bin)=anovan(xx,{yy},'display','off');
        end
        mi_all=[mi_all;mi];
        anovanp=[anovanp;ap];
    end
end
keyboard();
sig_sel=any(anovanp(:,1:7)<0.05,2);
mi_sel=mi_all(sig_sel,:);
COM=mi_sel(:,1:7)*((1:7).')./sum(mi_sel(:,1:7),2);
[~,cidx]=sort(COM);
gk = fspecial('gaussian', [3 3], 1);
fh=figure('Color','w');
imagesc(conv2(mi_sel(cidx,1:7),gk,'same'),[0,3])
colormap('turbo');
set(gca(),'YDir','normal')
ColorbarWithAxis([0,3],'Mutual information (bit)')
set(gca(),'XTick',1.5:3:7.5,'XTickLabel',-3:3:3)
xlabel('Time (s)')
ylabel('Neuron #')
arrayfun(@(x) xline(x,'--w'),[3.5,4.5]);
exportgraphics(fh,'Ordinal_MI.pdf','ContentType','vector')


