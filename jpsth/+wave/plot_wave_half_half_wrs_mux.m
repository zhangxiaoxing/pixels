function fhp=plot_wave_half_half_wrs_mux(sel_meta,opt)
arguments
    sel_meta
    opt.sess_id (1,1) double {mustBePositive,mustBeInteger} = 102
end
for sess=opt.sess_id
    %     while true
    fs=sprintf('s%d',sess);
    com_map=wave.get_com_map(sel_meta,'onepath',['SPKINFO/',ephys.sessid2path(sess)],'curve',true,'rnd_half',true,'delay',6,'wave','anyContext2');
    if ~isfield(com_map,fs)
        disp([sess,0,0])
        continue
    end
    samp_key_S1=intersect(cell2mat(com_map.(fs).c1a.keys),cell2mat(com_map.(fs).c1b.keys));
    samp_key_S2=intersect(cell2mat(com_map.(fs).c2a.keys),cell2mat(com_map.(fs).c2b.keys));
    COMS1=cell2mat(values(com_map.(fs).c1a,num2cell(samp_key_S1)));
    COMS2=cell2mat(values(com_map.(fs).c2a,num2cell(samp_key_S2)));
    if numel([samp_key_S1,samp_key_S2])<100
        disp([sess,numel([samp_key_S1,samp_key_S2]),0])
        continue
    end
    r=corr([COMS1,COMS2].',...
        [cell2mat(values(com_map.(fs).c1b,num2cell(samp_key_S1))),...
        cell2mat(values(com_map.(fs).c2b,num2cell(samp_key_S2)))].');
    if r<0.85
        disp([sess,numel([samp_key_S1,samp_key_S2]),r])
        continue
    end
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
normscale=max(abs([immata,immatb,immat_anti]),2);
fh=figure('Color','w','Position',[32,32,1080,140]);
plotOne(1,immata(sortidx,:)./normscale(sortidx).',sortedTCOM);
plotOne(2,immatb(sortidx,:)./normscale(sortidx).',TCOMb(sortidx));
shufidx=randsample(size(immatb,1),size(immatb,1));
plotOne(3,immatb(shufidx,:)./normscale(sortidx).',TCOMb(shufidx));
plotOne(4,immat_anti(sortidx,:)./normscale(sortidx).')
end
function plotOne(subidx,imdata,comdata)
subplot(1,4,subidx);
hold on
colormap('turbo');
gk = fspecial('gaussian', [3 3], 1);
imagesc(conv2(imdata,gk,'same'),[0 1])
if exist('comdata','var') && ~isempty(comdata)
    scatter(comdata,1:numel(comdata),2,'o','MarkerFaceColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
end
set(gca(),'XTick',[0.5,20.5],'XTickLabel',[0,5]);
colorbar();
ylim([0.5,size(imdata,1)+0.5])
xlim([0.5,24.5])
end