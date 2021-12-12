function plot_wave_3_6(opt)
arguments
    opt.sess (1,1) double {mustBeMember(opt.sess,1:116)} = 102
    opt.sortby (1,:) char {mustBeMember(opt.sortby,{'3s','6s'})} = '6s'
    opt.exportgraphics (1,1) logical = false
    opt.plot_global (1,1) logical = true
    opt.plot_session (1,1) logical = false
end
%% global
persistent com_map_3_sel com_map_3_alt com_map_6_sel com_map_6_alt com_map_3_sel_si com_map_3_alt_si com_map_6_sel_si com_map_6_alt_si
if opt.plot_global

    if isempty(com_map_3_sel)
        com_map_3_sel=wave.get_com_map('curve',true,'rnd_half',false,'delay',3);
        com_map_3_alt=wave.get_com_map('curve',true,'rnd_half',false,'delay',3,'alt_3_6',true);
        com_map_6_sel=wave.get_com_map('curve',true,'rnd_half',false,'delay',6);
        com_map_6_alt=wave.get_com_map('curve',true,'rnd_half',false,'delay',6,'alt_3_6',true);
    
        com_map_3_sel_si=wave.get_com_map('curve',true,'rnd_half',false,'delay',3,'selidx',true);
        com_map_3_alt_si=wave.get_com_map('curve',true,'rnd_half',false,'delay',3,'selidx',true,'alt_3_6',true);
        com_map_6_sel_si=wave.get_com_map('curve',true,'rnd_half',false,'delay',6,'selidx',true);
        com_map_6_alt_si=wave.get_com_map('curve',true,'rnd_half',false,'delay',6,'selidx',true,'alt_3_6',true);
    end
    for alt_comb=1:3
        switch alt_comb
            case 1
                com_map_6=com_map_6_sel;
                com_map_3=com_map_3_sel;
                com_map_6_si=com_map_6_sel_si;
                com_map_3_si=com_map_3_sel_si;
                tag='Selective in both 3s and 6s trials';
                fn='corr_3_6_both.pdf';
            case 2
                com_map_6=com_map_6_sel;
                com_map_3=com_map_3_alt;
                com_map_6_si=com_map_6_sel_si;
                com_map_3_si=com_map_3_alt_si;
                tag='Selective only in 6s trials';
                fn='corr_3_6_only6.pdf';
            case 3
                com_map_6=com_map_6_alt;
                com_map_3=com_map_3_sel;
                com_map_6_si=com_map_6_alt_si;
                com_map_3_si=com_map_3_sel_si;
                tag='Selective only in 3s trials';
                fn='corr_3_6_only3.pdf';
        end
        fss=intersect(fieldnames(com_map_6),fieldnames(com_map_3)).';

        [immata,immatb,immat_anti_a,immat_anti_b,immat_si_a,immat_si_b,com_a,com_b]=deal([]);
        for fs1=fss
            fs=char(fs1);
            for pref=["s1","s2"]
                curr_key=intersect(cell2mat(com_map_3.(fs).(pref).keys),cell2mat(com_map_6.(fs).(pref).keys));
                heat3=cell2mat(com_map_3.(fs).([char(pref),'curve']).values(num2cell(curr_key.')));
                heat6=cell2mat(com_map_6.(fs).([char(pref),'curve']).values(num2cell(curr_key.')));
                anti3=cell2mat(com_map_3.(fs).([char(pref),'anticurve']).values(num2cell(curr_key.')));
                anti6=cell2mat(com_map_6.(fs).([char(pref),'anticurve']).values(num2cell(curr_key.')));
                si3=cell2mat(com_map_3_si.(fs).([char(pref),'curve']).values(num2cell(curr_key.')));
                si6=cell2mat(com_map_6_si.(fs).([char(pref),'curve']).values(num2cell(curr_key.')));
                COM3=cell2mat(values(com_map_3.(fs).(pref),num2cell(curr_key)));
                COM6=cell2mat(values(com_map_6.(fs).(pref),num2cell(curr_key)));
                %         keyboard()
                immata=[immata;heat3];
                immatb=[immatb;heat6];
                immat_anti_a=[immat_anti_a;anti3];
                immat_anti_b=[immat_anti_b;anti6];
                immat_si_a=[immat_si_a;si3];
                immat_si_b=[immat_si_b;si6];
                com_a=[com_a;COM3.'];
                com_b=[com_b;COM6.'];
            end
        end

        r=corr(com_a,com_b);
        keyboard()

        currscale=max(abs([immata,immatb,immat_anti_a,immat_anti_b]),[],2);
        immata=immata./currscale;
        immatb=immatb./currscale;
        immat_anti_a=immat_anti_a./currscale;
        immat_anti_b=immat_anti_b./currscale;        
        fh=figure('Color','w','Position',[32,32,1400,800]);
        if alt_comb==3
            [~,sortidx]=sort(com_a);
        else
            [~,sortidx]=sort(com_b);
        end
        plotOne(2,immata(sortidx,:),com_a(sortidx),'sub_dim',[3,2],'title','FR, preferred, 3s','cmap','turbo');
        plotOne(1,immatb(sortidx,:),com_b(sortidx),'sub_dim',[3,2],'title','FR, preferred, 6s','cmap','turbo');
        plotOne(4,immat_anti_a(sortidx,:),com_a(sortidx),'sub_dim',[3,2],'title','FR, non-preferred, 3s','cmap','turbo');
        plotOne(3,immat_anti_b(sortidx,:),com_b(sortidx),'sub_dim',[3,2],'title','FR, non-preferred, 6s','cmap','turbo');
        plotOne(6,immat_si_a(sortidx,:),com_a(sortidx),'sub_dim',[3,2],'scale',[-0.6,0.6],'title','Selectivit index, 3s','cmap','turbo');
        plotOne(5,immat_si_b(sortidx,:),com_b(sortidx),'sub_dim',[3,2],'scale',[-0.6,0.6],'title','Selectivit index, 6s','cmap','turbo');
%         shufidx=randsample(size(immatb,1),size(immatb,1));
%         plotOne(3,immatb(shufidx,:),com_b(shufidx));
        sgtitle(sprintf('%s, FRTC r=%0.3f',tag,r));
        exportgraphics(fh,fn,'ContentType','vector')
    end
end
%     exportgraphics(fh,'wave_3s_6s_global.pdf','ContentType','vector');


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