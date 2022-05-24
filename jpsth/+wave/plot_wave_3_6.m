function out=plot_wave_3_6(opt)
arguments
    opt.sess (1,1) double {mustBeMember(opt.sess,1:116)} = 102
    opt.sortby (1,:) char {mustBeMember(opt.sortby,{'3s','6s'})} = '6s'
    opt.exportgraphics (1,1) logical = false
    opt.plot_global (1,1) logical = true
    opt.plot_session (1,1) logical = false
    opt.comb_set (1,:) double {mustBeInteger,mustBePositive} = 1:3
    opt.bootrpt (1,1) double {mustBeInteger,mustBePositive} = 3;
    opt.plot_2d_corr (1,1) logical = true
    opt.plot_corr_dist (1,1) logical = true
end
%% global
out=struct();
out.corr2=struct();
out.corr_dist=struct();
persistent com_map_3_sel com_map_3_alt com_map_6_sel com_map_6_alt com_map_3_sel_si com_map_3_alt_si com_map_6_sel_si com_map_6_alt_si

if isempty(com_map_3_sel)
end
for alt_comb=opt.comb_set
    switch alt_comb
        case 1
            com_map_3=wave.get_com_map('curve',true,'rnd_half',false,'delay',3,'wave','both');
            com_map_6=wave.get_com_map('curve',true,'rnd_half',false,'delay',6,'wave','both');
            tag='Selective in both 3s and 6s trials';
            fn='corr_3_6_both.pdf';
        case 2
            com_map_3=wave.get_com_map('curve',true,'rnd_half',false,'delay',3,'wave','only6');
            com_map_6=wave.get_com_map('curve',true,'rnd_half',false,'delay',6,'wave','only6');%             com_map_6_si=com_map_6_sel_si;
            tag='Selective only in 6s trials';
            fn='corr_3_6_only6.pdf';
        case 3
            com_map_3=wave.get_com_map('curve',true,'rnd_half',false,'delay',3,'wave','only3');
            com_map_6=wave.get_com_map('curve',true,'rnd_half',false,'delay',6,'wave','only3');
            tag='Selective only in 3s trials';
            fn='corr_3_6_only3.pdf';
    end
    fss=intersect(fieldnames(com_map_6),fieldnames(com_map_3)).';

    [immata,immatb,immat_anti_a,immat_anti_b,com_a,com_b]=deal([]);
    for fs1=fss
        fs=char(fs1);
        for pref=["s1","s2"]
            curr_key=intersect(cell2mat(com_map_3.(fs).(pref).keys),cell2mat(com_map_6.(fs).(pref).keys));
            heat3=cell2mat(com_map_3.(fs).([char(pref),'curve']).values(num2cell(curr_key.')));
            heat6=cell2mat(com_map_6.(fs).([char(pref),'curve']).values(num2cell(curr_key.')));
            anti3=cell2mat(com_map_3.(fs).([char(pref),'anticurve']).values(num2cell(curr_key.')));
            anti6=cell2mat(com_map_6.(fs).([char(pref),'anticurve']).values(num2cell(curr_key.')));
            COM3=cell2mat(values(com_map_3.(fs).(pref),num2cell(curr_key)));
            COM6=cell2mat(values(com_map_6.(fs).(pref),num2cell(curr_key)));
            %         keyboard()
            immata=[immata;heat3];
            immatb=[immatb;heat6];
            immat_anti_a=[immat_anti_a;anti3];
            immat_anti_b=[immat_anti_b;anti6];
            com_a=[com_a;COM3.'];
            com_b=[com_b;COM6.'];
        end
    end

    r=corr(com_a,com_b);

    currscale=max(abs([immata,immatb,immat_anti_a,immat_anti_b]),[],2);
    immata=immata./currscale;
    immatb=immatb./currscale;

    flat_mat=@(x) x(:);
    scale_mat=@(x) (x(:,1:2:size(x,2))+x(:,2:2:size(x,2)))./2;

    r2s=corr(flat_mat(immata),flat_mat(scale_mat(immatb)));
    r2sci=bootci(opt.bootrpt,@(x,y) corr(x,y),flat_mat(immata),flat_mat(scale_mat(immatb)));
    r2e=corr(flat_mat(immata),flat_mat(immatb(:,1:12)));
    r2eci=bootci(opt.bootrpt,@(x,y) corr(x,y),flat_mat(immata),flat_mat(immatb(:,1:12)));
    r2l=corr(flat_mat(immata),flat_mat(immatb(:,13:end)));
    r2lci=bootci(opt.bootrpt,@(x,y) corr(x,y),flat_mat(immata),flat_mat(immatb(:,13:end)));
    typestr=(replace(fn,'.pdf',''));
    out.corr2.(typestr)=[r2s,r2e,r2l];
    out.corr2ci.(typestr)=[r2sci,r2eci,r2lci];

    out.corr_dist.(typestr)=[arrayfun(@(x) corr(immata(x,:).',scale_mat(immatb(x,:)).'),1:size(immata,1));...
        arrayfun(@(x) corr(immata(x,:).',immatb(x,1:12).'),1:size(immata,1))];

    immat_anti_a=immat_anti_a./currscale;
    immat_anti_b=immat_anti_b./currscale;
    if opt.plot_global
        %% similarity confusion matrix
        ct=[immata,immatb,immat_anti_a,immat_anti_b];
        confumat=nan(size(ct,2));
        for ii=1:size(ct,2)
            for jj=1:size(ct,2)
                confumat(ii,jj)=corr(ct(:,ii),ct(:,jj));
            end
        end



        fch=figure('Color','w','Position',[32,32,850,750]);
        imagesc(confumat,[-1,1]);
        set(gca,'YDir','normal')
        colormap('turbo')
        cbh=colorbar();
        set(gca,'XTick',[0.5,12.5,24.5,36.5,48.5,60.5,72.5],'XTickLabel',{'0','3/0','3','6/0','3/0','3','6'})
        set(gca,'YTick',[0.5,12.5,24.5,36.5,48.5,60.5,72.5],'YTickLabel',{'0','3/0','3','6/0','3/0','3','6'})
        xlabel('Time(s) L->R: 3s-pref; 6s-pref; 3s nonpref; 6s nonpref')
        ylabel('Time(s) L->R: 3s-pref; 6s-pref; 3s nonpref; 6s nonpref')
        title(tag)
        cbh.Label.String='Pearson correlation r';


        %% wave
        fh=figure('Color','w','Position',[32,32,1400,800]);
        if alt_comb==3
            [~,sortidx]=sort(com_a);
        else
            [~,sortidx]=sort(com_b);
        end
        plotOne(2,immata(sortidx,:),com_a(sortidx),'sub_dim',[2,2],'title','FR, preferred, 3s','cmap','turbo');
        plotOne(1,immatb(sortidx,:),com_b(sortidx),'sub_dim',[2,2],'title','FR, preferred, 6s','cmap','turbo');
        plotOne(4,immat_anti_a(sortidx,:),com_a(sortidx),'sub_dim',[2,2],'title','FR, non-preferred, 3s','cmap','turbo');
        plotOne(3,immat_anti_b(sortidx,:),com_b(sortidx),'sub_dim',[2,2],'title','FR, non-preferred, 6s','cmap','turbo');
        sgtitle(sprintf('%s, FRTC r=%0.3f',tag,r));
        exportgraphics(fh,fn,'ContentType','vector')
    end
end
%     exportgraphics(fh,'wave_3s_6s_global.pdf','ContentType','vector');
save('wave_3_6_corr_data.mat','out');

%% 2d corr
if opt.plot_2d_corr
    fh=figure('Color','w','Position',[32,32,330,150]);
    hold on
    bh=bar([out.corr2.corr_3_6_both;out.corr2.corr_3_6_only6;out.corr2.corr_3_6_only3],'FaceColor','w','EdgeColor','k');
    bh(2).FaceColor=ones(1,3)./2;
    bh(3).FaceColor='k';

    xmat=cell2mat({bh.XEndPoints}.');
    fn=fieldnames(out.corr2ci);
    for fidx=1:numel(fn)
        errorbar(xmat(:,fidx),out.corr2.(fn{fidx}),...
            out.corr2ci.(fn{fidx})(1,:)-out.corr2.(fn{fidx}),...
            out.corr2ci.(fn{fidx})(2,:)-out.corr2.(fn{fidx}),'k.')
    end

    set(gca,'XTick',[])
    ylabel('FR correlation (r)')
    legend(bh,{'Scaled','Early half','Late half'},'Location','northoutside','Orientation','horizontal')
    exportgraphics(fh,'wave_3_6_FR_corr.pdf','ContentType','vector')
end
%% per su corr
if opt.plot_corr_dist
    fhb=figure('Color','w');
    subplot(1,3,1)
    scatter(out.corr_dist.corr_3_6_both(1,:),out.corr_dist.corr_3_6_both(2,:),4,'ko','filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.5)
    xlabel('3s correlation to scaled 6s')
    ylabel('3s correlation to early 6s')
    title('Both 3 and 6s')
    xlim([0,1])
    ylim([0,1])
    subplot(1,3,2)
    scatter(out.corr_dist.corr_3_6_only3(1,:),out.corr_dist.corr_3_6_only3(2,:),4,'ko','filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.5)
    xlabel('3s correlation to scaled 6s')
    ylabel('3s correlation to early 6s')
    title('Only 3s')
    xlim([0,1])
    ylim([0,1])
    subplot(1,3,3)
    scatter(out.corr_dist.corr_3_6_only6(1,:),out.corr_dist.corr_3_6_only6(2,:),4,'ko','filled','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.5)
    xlabel('3s correlation to scaled 6s')
    ylabel('3s correlation to early 6s')
    title('Only 6s')
    xlim([0,1])
    ylim([0,1])
    exportgraphics(fhb,'corr_6s_scale_vs_6s_early.pdf','ContentType','vector')
end
%% session
if opt.plot_session
    for sess=opt.sess
        %     while true
        disp(sess)
        fs=sprintf('s%d',sess);
        com_map_3=wave.get_com_map('onepath',['SPKINFO/',ephys.sessid2path(sess)],'curve',true,'rnd_half',false,'delay',3);
        if strcmp(opt.sortby,'3s')
            com_map_6=wave.get_com_map('onepath',['SPKINFO/',ephys.sessid2path(sess)],'curve',true,'rnd_half',false,'delay',6,'partial','early3in6');
        else
            com_map_6=wave.get_com_map('onepath',['SPKINFO/',ephys.sessid2path(sess)],'curve',true,'rnd_half',false,'delay',6);
        end
        if ~isfield(com_map_6,fs) || ~isfield(com_map_3,fs)
            continue
        end

        samp_key_S1=intersect(cell2mat(com_map_3.(fs).s1.keys),cell2mat(com_map_6.(fs).s1.keys));
        samp_key_S2=intersect(cell2mat(com_map_3.(fs).s2.keys),cell2mat(com_map_6.(fs).s2.keys));
        switch opt.sortby
            case '3s'
                COMS1=cell2mat(values(com_map_3.(fs).s1,num2cell(samp_key_S1)));
                COMS2=cell2mat(values(com_map_3.(fs).s2,num2cell(samp_key_S2)));
                if numel(COMS1)+numel(COMS2)<50
                    continue;
                end
                r=corr([COMS1,COMS2].',...
                    [cell2mat(values(com_map_6.(fs).s1,num2cell(samp_key_S1))),...
                    cell2mat(values(com_map_6.(fs).s2,num2cell(samp_key_S2)))].');
                plotidces=1:2;
            case '6s'
                COMS1=cell2mat(values(com_map_6.(fs).s1,num2cell(samp_key_S1)));
                COMS2=cell2mat(values(com_map_6.(fs).s2,num2cell(samp_key_S2)));
                if numel(COMS1)+numel(COMS2)<50
                    continue;
                end
                r=corr([COMS1,COMS2].',...
                    [cell2mat(values(com_map_3.(fs).s1,num2cell(samp_key_S1))),...
                    cell2mat(values(com_map_3.(fs).s2,num2cell(samp_key_S2)))].');
                plotidces=2:-1:1;

        end

        if r<0.78
            continue
        end


        sortmat=[ones(size(COMS1)),2*ones(size(COMS2));...
            double(samp_key_S1),double(samp_key_S2);...
            COMS1,COMS2].';

        sortmat=sortrows(sortmat,3);

        immata=[];
        com_a=[];
        for ri=1:size(sortmat,1)
            one_heat=com_map_3.(fs).(sprintf('s%dcurve',sortmat(ri,1)))(sortmat(ri,2));
            immata=[immata;one_heat];
            com_a=[com_a;com_map_3.(fs).(sprintf('s%d',sortmat(ri,1)))(sortmat(ri,2))];
        end
        immatb=[];
        com_b=[];
        for ri=1:size(sortmat,1)
            one_heat=com_map_6.(fs).(sprintf('s%dcurve',sortmat(ri,1)))(sortmat(ri,2));
            immatb=[immatb;one_heat];
            com_b=[com_b;com_map_6.(fs).(sprintf('s%d',sortmat(ri,1)))(sortmat(ri,2))];
        end
        immata=immata./max(abs([immata,immatb]),[],2);
        immatb=immatb./max(abs([immata,immatb]),[],2);

        fh=figure('Color','w','Position',[32,32,1080,140]);
        plotOne(plotidces(1),immata,com_a);
        plotOne(plotidces(2),immatb,com_b);
        shufidx=randsample(size(immatb,1),size(immatb,1));
        plotOne(3,immatb(shufidx,:),com_b(shufidx));
        sgtitle(num2str([sess,r]));
        %         exportgraphics(fh,sprintf('wave_half_half_%d.png',sess));
        keyboard()
        exportgraphics(fh,sprintf('wave_3s_6s_%d.pdf',sess));
        %         close all
        %         waitfor(fh)

        %     end
    end
end
end

function plotOne(subidx,imdata,comdata,opt)
arguments
    subidx (1,1) double {mustBeInteger,mustBePositive}
    imdata (:,:) double
    comdata = []
    opt.sub_dim (1,2) = [1,3]
    opt.plot_com (1,1) logical = false
    opt.scale (1,2) double = [-1,1]
    opt.title (1,:) char = []
    opt.cmap (1,:) char = 'turbo'

end
subplot(opt.sub_dim(1),opt.sub_dim(2),subidx);
hold on
colormap(opt.cmap);
gk = fspecial('gaussian', [3 3], 1);
% if size(imdata,2)<20
%     cpos=get(axes,'Position');
%     axes(sh,'Position',[cpos(1),cpos(2),0.5*(diff(cpos([1 3]))),cpos(4)])
% end
imagesc(conv2(imdata,gk,'same'),opt.scale)
if opt.plot_com && exist('comdata','var') && ~isempty(comdata)
    scatter(comdata,1:numel(comdata),2,'o','MarkerFaceColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
end
if size(imdata,2)>20
    set(gca(),'XTick',[0.5,20.5],'XTickLabel',[0,5]);
    xline(12.5,'--w','LineWidth',1);
    xlim([0.5,size(imdata,2)+0.5])
else
    xlim([0.5,size(imdata,2)*2+0.5])
    set(gca(),'XTick',[0.5,12.5],'XTickLabel',[0,3]);
end
colorbar();
ylim([0.5,size(imdata,1)+0.5])
if numel(opt.title)>0
    title(opt.title)
end

end