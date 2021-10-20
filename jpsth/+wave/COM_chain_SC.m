function COM_chain_SC(opt)
arguments
    opt.strict (1,1) logical = true %strict FC criteria
%     opt.min_diff (1,1)
    opt.screen (1,1) logical = true
end
load('sums_conn.mat','sums_conn_str');
meta_str=ephys.util.load_meta('type','neupix');
% warning('partial iteration for illustration')
figidx=1;
gk = fspecial('gaussian', [3 3], 1);
for fidx=1:numel(sums_conn_str)
    if opt.strict
        ccgqc=sums_conn_str(fidx).qc; %reference quality control parameter
        strict_sel=ccgqc(:,2)>=252 & ccgqc(:,4)>=2 & ccgqc(:,4)<=40 & ccgqc(:,5)>248;
        %1:Polar 2:Time of peak 3:Noise peaks 4:FWHM 5:rising edge
        %6:falling edge
        oneccg=sums_conn_str(fidx).ccg_sc(strict_sel,:); %ccg
        onecon=sums_conn_str(fidx).sig_con(strict_sel,:); %jitter controlled significant functional coupling
    disp([fidx,nnz(strict_sel),numel(strict_sel)]);    
    else % full input from English, Buzsaki code
        oneccg=sums_conn_str(fidx).ccg_sc;
        onecon=sums_conn_str(fidx).sig_con;
    end
    %TODO nonmem,incongruent
    onecom=wave.get_com_map('onepath',sums_conn_str(fidx).folder,'curve',true); %per su center of mass, normalized FR
    skey=fieldnames(onecom);
    for samp=["s1","s2"]
        mapkeys=cell2mat(onecom.(skey{1}).(samp).keys); %pre-selected transient selective su
        typesel=all(ismember(int32(onecon(:,:)),mapkeys),2);
        typesigcon=onecon(typesel,:);
        con_com_diff=arrayfun(@(x) onecom.(skey{1}).(samp)(x),int32(onecon(typesel,:)));
        dirsel=con_com_diff(:,2)-con_com_diff(:,1)>(10/250); % Assuming 250ms bin
        dirsigcon=typesigcon(dirsel,:);
        upre=unique(dirsigcon(:,1)).';
        
        chains=cell(0);
        for i=upre
            onechain=cell(0);
            cpre=i;
            while true %many to many wave iteration
                newpair=dirsigcon(ismember(dirsigcon(:,1),cpre),:);
                if isempty(newpair)
                    if numel(onechain)>2
                        chains{end+1}=onechain;
                    end
                    break
                else
                    onechain{end+1}={...
                        newpair,...
                        arrayfun(@(x) onecom.(skey{1}).(samp)(x),newpair)};
                    cpre=newpair(:,2);
                end
            end
        end
        if numel(chains)>0
            for ii=1:numel(chains)
                split_chains=wave.recursive_chain(chains{ii},[]); %one su per order
                for pp=1:2:size(split_chains,1)
                    onechain=split_chains(pp:pp+1,:);
                    if onechain(2,end)-onechain(2,1)<4
                        continue
                    end
                    fh=figure('Color','w','Position',[100,100,500,500]); %illustration
                    subplot(2,2,3);
                    hold on;
                    % envelop
                    yyaxis left
                    suids=onechain(1,:);
                    hm=[];
                    for mm=1:numel(suids)
                        oneheat=conv2(onecom.(skey{1}).(sprintf('%sheat',samp))(suids(mm)),gk,'same');
                        hm=[hm;oneheat;zeros(5,size(oneheat,2))];
                    end
                    colormap(flip(colormap('gray')));
                    imagesc(hm);
                    set(gca(),'XTick',0:8:24,'XTickLabel',0:2:6,'YTick',[])
                    xlim([0,24])
                    ylim([0,size(hm,1)-5]);
                    xlabel('Delay time (s)');
                    ylabel('Example trials')
%                     colorbar()
                    %FC
                    yyaxis right
                    for jj=1:size(onechain,2)
                        text(onechain(2,jj),jj,num2str(onechain(1,jj)));
                        plot(onechain(2,jj),jj,'ko','MarkerFaceColor','k');
                    end
                    for jj=1:size(onechain,2)-1
                        plot(onechain(2,jj:jj+1),[jj;jj+1],'r-');
                    end
                    ylim([0.5,size(onechain,2)+0.5]);
                    xlim([0,24]);
                    set(gca(),'XTick',[0,8,16,24],'XTickLabel',[0,2,4,6],'YTick',0:numel(onechain)+1)
                    xlabel('COM time in delay (s)')
                    ylabel('Functional coupling order in wave ->');
                    title(sprintf('Session %s, sample %s',skey{1},samp));
                    
                    
                    for kk=1:size(onechain,2)-1
                        conn=onechain(1,kk:kk+1);
                        subplot(2,size(onechain,2)-1,kk);
                        ccgsel=onecon(:,1)==conn(1) & onecon(:,2)==conn(2);
                        plot(oneccg(ccgsel,:),'-r');
                        arrayfun(@(x) xline(x,'--k'),[226,251,276]);
                        xlim([191,311]);
                        set(gca(),'XTick',201:25:301,'XTickLabel',-20:10:20)
                        fpath=replace(regexp(sums_conn_str(fidx).folder,'(?<=SPKINFO/).*','match','once'),'/','\');
                        reg1=meta_str.reg_tree{5,startsWith(meta_str.allpath,fpath)...
                            & meta_str.allcid==uint16(conn(1))};
                        reg2=meta_str.reg_tree{5,startsWith(meta_str.allpath,fpath)...
                            & meta_str.allcid==uint16(conn(2))};
                        
                        text(191,max(ylim()),{num2str(conn(1)),reg1},'HorizontalAlignment','left','VerticalAlignment','top');
                        text(301,max(ylim()),{num2str(conn(2)),reg2},'HorizontalAlignment','right','VerticalAlignment','top');
                        
                    end
                    subplot(2,2,4)
                    hold on
                    suids=onechain(1,:);
                    if false
                        cmap=colormap('cool');
                        deltac=floor(size(cmap,1)./numel(suids));
                        for mm=1:numel(suids)
                            plot((onecom.(skey{1}).(sprintf('%scurve',samp))(suids(mm))),...
                                '-','Color',cmap((mm-1)*deltac+1,:));
                            xline(onecom.(skey{1}).(samp)(suids(mm)),'--','Color',cmap((mm-1)*deltac+1,:));
                        end
                        set(gca(),'XTick',0:8:24,'XTickLabel',0:2:6)
                        xlabel('Delay time (s)');
                        ylabel('Normalized firing rate [0,1]');
                    end
                    sgtitle(sprintf('fidx %d chain %d-%d',fidx,ii,pp));
                    if opt.screen
                        exportgraphics(fh,sprintf('SC\\SC%05d.png',figidx),'Resolution',300);
                        keyboard()
                    else
                        keyboard()
                        exportgraphics(fh,'SC\SC_select.pdf');
                    end
                    figidx=figidx+1;
                    close(fh);
                end
            end
        else
            disp('Empty Chain');
        end
    end
end
end



