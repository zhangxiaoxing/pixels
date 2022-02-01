%TODO brain region filter, olfaction filter.
% function out=get_ordinal_mat(opt)
%% gen data
if false %gendata
    addpath('K:\Lib\MutualInfo\');
[~,~,sessmap]=ephys.sessid2path(0);
meta=ephys.util.load_meta();
homedir=ephys.util.getHomedir('type','raw');
[mi_all,anovanp,su_meta]=deal([]);

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
                trl_sel=trials(:,8)==dur & dur_resp(:,3)==dur_ord;
                dx=mean(fr(trl_sel,suidx,win),3);
                xx=[xx;dx];
                yy=[yy;double(ord).*ones(size(dx))];
            end
            mi(bin)=discrete_continuous_info_fast(yy,xx,5,2);
            ap(bin)=anovan(xx,{yy},'display','off');
        end
        mi_all=[mi_all;mi];
        anovanp=[anovanp;ap];
        su_meta=[su_meta;ii,suidx,suid(suidx)];
    end
end
end
if false
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
end
%% SU showcase
su_sel=any(anovanp(:,1:7)<0.01 & mi_all(:,1:7)>2,2);
su_plot=su_meta(su_sel,:);
mi_plot=mi_all(su_sel,1:7);
anp_plot=anovanp(su_sel,1:7);

for ii=1:size(su_plot,1)
    reg_tree=meta.reg_tree(:,meta.sess==su_plot(ii,1) & meta.allcid==su_plot(ii,3));
    disp([anp_plot(ii,:);mi_plot(ii,:)])
%     disp(reg_tree)
%     if ~strcmp(reg_tree{2},'CTX')
%         continue
%     end
%     [~,~,~,~,~,FT_SPIKE]=ephys.getSPKID_TS(ii,"keep_trial",true,"suids",[su_meta(ii,3)]);
%     dur_resp=behav.tag_block(FT_SPIKE.trialinfo,'wt',false);
    fpath=fullfile(homedir,sessmap(su_plot(ii,1)),"FR_All_ 250.hdf5");
    fr=h5read(fpath,'/FR_All');
    trials=h5read(fpath,'/Trials');
    suid=h5read(fpath,'/SU_id');

    dur_resp=behav.tag_block(trials,'wt',false);
    cmap=[winter(4);autumn(4)];

    fh=figure('Color','w');
    subplot(2,1,1)
    hold on;
    ph=nan(8,1);
    for ord=int32([31:34,61:64])
        dur=idivide(ord,10);
        dur_ord=rem(ord,10);
        trl_sel=trials(:,8)==dur & dur_resp(:,3)==dur_ord;
        cidx=dur_ord;
        if dur==6
            cidx=cidx+4;
        end
        ph(cidx)=plot(smooth(squeeze(mean(fr(trl_sel,su_plot(ii,2),1:28),1)),3),'Color',cmap(cidx,:));
    end
    arrayfun(@(x) xline(x,'--k'),[12,16]+0.5)
    title(strjoin(reg_tree([2,5]),'-'))
    legend(ph,{'3s No.1','3s No.2','3s No.3','3s No.4','6s No.1','6s No.2','6s No.3','6s No.4'},'location','eastoutside')
    
    bins=find(anp_plot(ii,:)<0.05 & mi_plot(ii,:)>2);
    for bb=1:numel(bins)
        subplot(numel(bins)*2,1,bb+numel(bins))
        hold on
        for ord=int32([31:34,61:64])
            dur=idivide(ord,10);
            dur_ord=rem(ord,10);
            trl_sel=trials(:,8)==dur & dur_resp(:,3)==dur_ord & all(trials(:,9:10)>0,2);
            cidx=dur_ord;
            if dur==6
                cidx=cidx+4;
            end
            mm=mean(fr(trl_sel,su_plot(ii,2),bins(bb)*4-3:bins(bb)*4),'all');
            sem=std(reshape(fr(trl_sel,su_plot(ii,2),bins(bb)*4-3:bins(bb)*4),1,[]))./sqrt((nnz(trl_sel).*4));
            bar(cidx,mm,'grouped','FaceColor',cmap(cidx,:));
            errorbar(cidx,mm,sem,'k.')
        end
        title(sprintf('sess#%d|suid#%d|bin #%d',su_plot(ii,1),su_plot(ii,3),bins(bb)))
    end
    keyboard()
    close(fh)
end



