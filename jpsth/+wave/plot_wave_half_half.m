% WIP, fix before use

function fhp=plot_wave_half_half(sel_meta,opt)
arguments
    sel_meta
    opt.sess_id (1,1) double {mustBePositive,mustBeInteger} = 102
    opt.multi_mod (1,1) logical = true
    opt.minr (1,1) double = 0.85
end
for sess=opt.sess_id
    fhp=[];
    %     while true
    fs=sprintf('s%d',sess);
    if opt.multi_mod
        [com_map_1h,com_map_2h]=wave.get_pct_com_map(sel_meta,'curve',true,'early_smooth',false,'rnd_half',true,'one_sess',102);
        samp_key_S1=intersect(cell2mat(com_map_1h.(fs).olf_s1.com6.keys),cell2mat(com_map_2h.(fs).olf_s1.com6.keys));
        samp_key_S2=intersect(cell2mat(com_map_1h.(fs).olf_s2.com6.keys),cell2mat(com_map_2h.(fs).olf_s2.com6.keys));
        COMS1=cell2mat(com_map_1h.(fs).olf_s1.com6.values(num2cell(samp_key_S1)));
        COMS2=cell2mat(com_map_1h.(fs).olf_s2.com6.values(num2cell(samp_key_S2)));
    else
        com_map=wave.get_com_map(sel_meta,'onepath',['SPKINFO/',ephys.sessid2path(sess)],'curve',true,'rnd_half',true,'delay',6,'wave','anyContext2');
        if ~isfield(com_map,fs)
            disp([sess,0,0])
            continue
        end
        samp_key_S1=intersect(cell2mat(com_map.(fs).c1a.keys),cell2mat(com_map.(fs).c1b.keys));
        samp_key_S2=intersect(cell2mat(com_map.(fs).c2a.keys),cell2mat(com_map.(fs).c2b.keys));
        COMS1=cell2mat(values(com_map.(fs).c1a,num2cell(samp_key_S1)));
        COMS2=cell2mat(values(com_map.(fs).c2a,num2cell(samp_key_S2)));
    end

    if numel([samp_key_S1,samp_key_S2])<100
        disp([sess,numel([samp_key_S1,samp_key_S2]),0])
        continue
    end
    if opt.multi_mod
        r=corr([COMS1,COMS2].',...
            [cell2mat(values(com_map_2h.(fs).olf_s1.com6,num2cell(samp_key_S1))),...
            cell2mat(values(com_map_2h.(fs).olf_s2.com6,num2cell(samp_key_S2)))].');
    else
        r=corr([COMS1,COMS2].',...
            [cell2mat(values(com_map.(fs).c1b,num2cell(samp_key_S1))),...
            cell2mat(values(com_map.(fs).c2b,num2cell(samp_key_S2)))].');
    end
    if r<opt.minr
        disp([sess,numel([samp_key_S1,samp_key_S2]),r])
        continue
    end

    if opt.multi_mod
        %% s1
        immata1=cell2mat(com_map_1h.(fs).olf_s1.s1d6.values(num2cell(samp_key_S1.')));
        immatb1=cell2mat(com_map_2h.(fs).olf_s1.s1d6.values(num2cell(samp_key_S1.')));
        immat_anti1=cell2mat(com_map_1h.(fs).olf_s1.s2d6.values(num2cell(samp_key_S1.')));
        TCOMb1=cell2mat(com_map_2h.(fs).olf_s1.com6.values(num2cell(samp_key_S1.')));
        %         fh1=plot_multi(immata1,immatb1,immat_anti1,COMS1,TCOMb1);
        %         sgtitle(num2str([sess,r,1]));
        %% s2
        immata2=cell2mat(com_map_1h.(fs).olf_s2.s2d6.values(num2cell(samp_key_S2.')));
        immatb2=cell2mat(com_map_2h.(fs).olf_s2.s2d6.values(num2cell(samp_key_S2.')));
        immat_anti2=cell2mat(com_map_1h.(fs).olf_s2.s1d6.values(num2cell(samp_key_S2.')));
        TCOMb2=cell2mat(com_map_2h.(fs).olf_s2.com6.values(num2cell(samp_key_S2.')));
        %         fh2=plot_multi(immata2,immatb2,immat_anti2,COMS2,TCOMb2);
        %         sgtitle(num2str([sess,r,2]));
    else
        %% s1
        immata1=cell2mat(com_map.(fs).c1acurve.values(num2cell(samp_key_S1.')));
        immatb1=cell2mat(com_map.(fs).c1bcurve.values(num2cell(samp_key_S1.')));
        immat_anti1=cell2mat(com_map.(fs).c1aanticurve.values(num2cell(samp_key_S1.')));
        TCOMb1=cell2mat(com_map.(fs).c1b.values(num2cell(samp_key_S1.')));
        %         fh1=plot_multi(immata1,immatb1,immat_anti1,COMS1,TCOMb1);
        %         sgtitle(num2str([sess,r,1]));
        %% s2
        immata2=cell2mat(com_map.(fs).c2acurve.values(num2cell(samp_key_S2.')));
        immatb2=cell2mat(com_map.(fs).c2bcurve.values(num2cell(samp_key_S2.')));
        immat_anti2=cell2mat(com_map.(fs).c2aanticurve.values(num2cell(samp_key_S2.')));
        TCOMb2=cell2mat(com_map.(fs).c2b.values(num2cell(samp_key_S2.')));
        %         fh2=plot_multi(immata2,immatb2,immat_anti2,COMS2,TCOMb2);
        %         sgtitle(num2str([sess,r,2]));        
    end
    fhp=plot_multi([immata1;immata2],[immatb1;immatb2],[immat_anti1;immat_anti2],[COMS1,COMS2],[TCOMb1;TCOMb2]);
    sgtitle(num2str([sess,r,2]));

%     keyboard()
    %         exportgraphics(fh1,sprintf('wave_half_half_%d_S1.pdf',sess));
    %         exportgraphics(fh2,sprintf('wave_half_half_%d_S2.pdf',sess));
%     exportgraphics(fhp,sprintf('wave_half_half_%d_prefered.pdf',sess));

    %         close all
    %         waitfor(fh)

    %     end
end
end
function fh=plot_multi(immata,immatb,immat_anti,TCOMa,TCOMb)
[sortedTCOM,sortidx]=sort(TCOMa);
normscale=max([immata,immatb,immat_anti],[],2);
fh=figure('Color','w','Position',[32,32,1000,150]);
plotOne(1,immata(sortidx,:)./normscale(sortidx),sortedTCOM);
plotOne(2,immatb(sortidx,:)./normscale(sortidx),TCOMb(sortidx));
shufidx=randsample(size(immatb,1),size(immatb,1));
plotOne(3,immatb(shufidx,:)./normscale(sortidx),TCOMb(shufidx));
plotOne(4,immat_anti(sortidx,:)./normscale(sortidx))
end
function plotOne(subidx,imdata,comdata)
subplot(1,4,subidx);
hold on
colormap('turbo');
% colormap('gray');
gk = fspecial('gaussian', [3 3], 1);
imagesc(conv2(imdata,gk,'same'),[0 0.7])
if exist('comdata','var') && ~isempty(comdata)
    scatter(comdata,1:numel(comdata),4,'o','MarkerFaceColor','k','MarkerFaceAlpha',1,'MarkerEdgeColor','w');
end
set(gca(),'XTick',[0.5,20.5],'XTickLabel',[0,5]);
colorbar();
ylim([0.5,size(imdata,1)+0.5])
xlim([0.5,24.5])
end