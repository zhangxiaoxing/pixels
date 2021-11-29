function plot_wave_3_6(opt)
arguments
    opt.sess (1,1) double {mustBeMember(opt.sess,1:116)} = 102
    opt.sortby (1,:) char {mustBeMember(opt.sortby,{'3s','6s'})} = '6s'
    opt.exportgraphics (1,1) logical = false
    opt.plot_global (1,1) logical = true
    opt.plot_session (1,1) logical = false
end
%% global
if opt.plot_global
    com_map_3=wave.get_com_map('curve',true,'rnd_half',false,'delay',3);
    com_map_6=wave.get_com_map('curve',true,'rnd_half',false,'delay',6);

    %iterate
    fss=intersect(fieldnames(com_map_6),fieldnames(com_map_3)).';
    [immata,immatb,com_a,com_b]=deal([]);
    for fs1=fss
        fs=char(fs1);
        for pref=["s1","s2"]
            samp_key=intersect(cell2mat(com_map_3.(fs).(pref).keys),cell2mat(com_map_6.(fs).(pref).keys));
            heat3=cell2mat(com_map_3.(fs).([char(pref),'curve']).values(num2cell(samp_key.')));
            heat6=cell2mat(com_map_6.(fs).([char(pref),'curve']).values(num2cell(samp_key.')));
            COM3=cell2mat(values(com_map_3.(fs).(pref),num2cell(samp_key)));
            COM6=cell2mat(values(com_map_6.(fs).(pref),num2cell(samp_key)));
            %         keyboard()
            immata=[immata;heat3];
            immatb=[immatb;heat6];
            com_a=[com_a;COM3.'];
            com_b=[com_b;COM6.'];
        end
    end

    r=corr(com_a,com_b);
    [com_b,sortidx]=sort(com_b);

    immata=immata./max(abs([immata,immatb]),[],2);
    immatb=immatb./max(abs([immata,immatb]),[],2);

    fh=figure('Color','w','Position',[32,32,1080,140]);
    plotOne(2,immata(sortidx,:),com_a(sortidx));
    plotOne(1,immatb(sortidx,:),com_b);
    shufidx=randsample(size(immatb,1),size(immatb,1));
    plotOne(3,immatb(shufidx,:),com_b(shufidx));
    sgtitle(sprintf('r=%0.3f',r));
    exportgraphics(fh,'wave_3s_6s_global.pdf','ContentType','vector');
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

function plotOne(subidx,imdata,comdata)
subplot(1,3,subidx);
hold on
colormap('turbo');
gk = fspecial('gaussian', [3 3], 1);
imagesc(conv2(imdata,gk,'same'),[-1 1])
scatter(comdata,1:numel(comdata),2,'o','MarkerFaceColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
if size(imdata,2)>20
    set(gca(),'XTick',[0.5,20.5],'XTickLabel',[0,5]);
    xline(12.5,'--w','LineWidth',1);
else
    set(gca(),'XTick',[0.5,12.5],'XTickLabel',[0,3]);
end
colorbar();
ylim([0.5,size(imdata,1)+0.5])
xlim([0.5,size(imdata,2)+0.5])
end