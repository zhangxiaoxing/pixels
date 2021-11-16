
function session_graph(opt)
arguments
    opt.gen_data (1,1) logical = false
    opt.plot_scatter (1,1) logical = false % in x-y correlation style
    opt.plot_bar (1,1) logical = false % per-session pair-wise connected dots
    opt.plot_shuf_bar (1,1) logical = true
    opt.exportgraph (1,1) logical = true
    opt.save_preview (1,1) logical = false
end
if opt.gen_data
    [sig,~]=bz.load_sig_pair('pair',false);
    meta=ephys.util.load_meta();
    ctxsel=all(squeeze(sig.reg(:,2,:))==688,2);
    for sess=reshape(unique(sig.sess),1,[])
        for rtype=["within","cross","all"]
            switch rtype
                case "within"
                    rsel=all(squeeze(sig.reg(:,5,:))>0,2) & squeeze(sig.reg(:,5,1)==sig.reg(:,5,2));
                case "cross"
                    rsel=all(squeeze(sig.reg(:,5,:))>0,2) & squeeze(sig.reg(:,5,1)~=sig.reg(:,5,2));
                case "all"
                    rsel=all(squeeze(sig.reg(:,5,:))>0,2);
            end

            %meta select data region size

            [s1c,s2c,incongc,nonc,anyc]=deal({'Source','Target'});
            [s1n,s2n,nonn,anyn]=deal({'Id','Label'});

            sess_sel=find(sig.sess==sess & ctxsel & rsel);

            %same s1, same s2, nonmem, any
            % {source,target}
            for si=reshape(sess_sel,1,[])
                anyc=[anyc;{sig.suid(si,1),sig.suid(si,2)}];
                if any(ismember(sig.mem_type(si,:),1:2),'all') && any(ismember(sig.mem_type(si,:),3:4),'all')
                    incongc=[incongc;{sig.suid(si,1),sig.suid(si,2)}];
                end
                if all(sig.mem_type(si,:)==0,'all') % non
                    nonc=[nonc;{sig.suid(si,1),sig.suid(si,2)}];
                elseif all(ismember(sig.mem_type(si,:),1:2),'all')
                    s1c=[s1c;{sig.suid(si,1),sig.suid(si,2)}];
                elseif all(ismember(sig.mem_type(si,:),3:4),'all')
                    s2c=[s2c;{sig.suid(si,1),sig.suid(si,2)}];
                end
            end
            writecell(s1c,fullfile('gephidata',sprintf('s1_edge_%03d_%s.csv',sess,rtype)));
            writecell(s2c,fullfile('gephidata',sprintf('s2_edge_%03d_%s.csv',sess,rtype)));
            writecell(incongc,fullfile('gephidata',sprintf('incong_edge_%03d_%s.csv',sess,rtype)));
            writecell(nonc,fullfile('gephidata',sprintf('non_edge_%03d_%s.csv',sess,rtype)));
            writecell(anyc,fullfile('gephidata',sprintf('any_edge_%03d_%s.csv',sess,rtype)));

            s1n=[s1n;num2cell(repmat(meta.allcid(meta.sess==sess & ismember(meta.mem_type,1:2).' & strcmp(meta.reg_tree(2,:),'CTX').'),1,2))];
            s2n=[s2n;num2cell(repmat(meta.allcid(meta.sess==sess & ismember(meta.mem_type,3:4).' & strcmp(meta.reg_tree(2,:),'CTX').'),1,2))];
            incongn=[s1n;s2n(2:end,:)];
            nonn=[nonn;num2cell(repmat(meta.allcid(meta.sess==sess & (meta.mem_type==0).' & strcmp(meta.reg_tree(2,:),'CTX').'),1,2))];
            anyn=[anyn;num2cell(repmat(meta.allcid(meta.sess==sess & (meta.mem_type>=0).' & strcmp(meta.reg_tree(2,:),'CTX').'),1,2))];

            writecell(s1n,fullfile('gephidata',sprintf('s1_node_%03d_%s.csv',sess,rtype)));
            writecell(s2n,fullfile('gephidata',sprintf('s2_node_%03d_%s.csv',sess,rtype)));
            writecell(incongn,fullfile('gephidata',sprintf('incong_node_%03d_%s.csv',sess,rtype)));
            writecell(nonn,fullfile('gephidata',sprintf('non_node_%03d_%s.csv',sess,rtype)));
            writecell(anyn,fullfile('gephidata',sprintf('any_node_%03d_%s.csv',sess,rtype)));

            %% Node-matching surrogate nonmem and incong
            %             keyboard()
            %             match s1 then s2
            for rpt=1:100
                if size(nonn,1)>size(s1n,1)
                    shufn1=nonn([1;randsample(size(nonn,1)-1,size(s1n,1)-1)+1],:);
                    edgesel1=find(all(ismember(cell2mat(nonc(2:end,:)),unique(int32(cell2mat(shufn1(2:end,:))))),2))+1;
                    shufc1=nonc([1;edgesel1],:);
                    writecell(shufc1,fullfile('gephidata',sprintf('shuf1_non_edge_%03d_%s_%d.csv',sess,rtype,rpt)));
                    writecell(shufn1,fullfile('gephidata',sprintf('shuf1_non_node_%03d_%s_%d.csv',sess,rtype,rpt)));
                end
                if size(nonn,1)>size(s2n,1)
                    shufn2=nonn([1;randsample(size(nonn,1)-1,size(s2n,1)-1)+1],:);
                    edgesel2=find(all(ismember(cell2mat(nonc(2:end,:)),unique(int32(cell2mat(shufn2(2:end,:))))),2))+1;
                    shufc2=nonc([1;edgesel2],:);
                    writecell(shufc2,fullfile('gephidata',sprintf('shuf2_non_edge_%03d_%s_%d.csv',sess,rtype,rpt)));
                    writecell(shufn2,fullfile('gephidata',sprintf('shuf2_non_node_%03d_%s_%d.csv',sess,rtype,rpt)));
                end
            end
        end
    end
end
if opt.plot_scatter || opt.plot_bar
    for rtype=["within","cross","all"]
        congru_S1=readmatrix(fullfile('gephidata',sprintf('s1_%s_gephi_graph_sums.csv',rtype)));
        congru_S2=readmatrix(fullfile('gephidata',sprintf('s2_%s_gephi_graph_sums.csv',rtype)));
        incong=readmatrix(fullfile('gephidata',sprintf('incong_%s_gephi_graph_sums.csv',rtype)));
        nonmem=readmatrix(fullfile('gephidata',sprintf('non_%s_gephi_graph_sums.csv',rtype)));

        %     any=readmatrix(fullfile('gephidata',sprintf('any_%s_gephi_graph_sums.csv',rtype)));
%         if opt.plot_scatter
%             plotOne(congru_S1(node_sel,[1 3]),nonmem(node_sel,[1 3]),'Connected comp.');
%             plotOne(congru_S1(node_sel,[1 4]),nonmem(node_sel,[1 4]),'Avg. cluster coef.');
%             plotOne(congru_S1(node_sel,[1 5]),nonmem(node_sel,[1 5]),'Density');
%             plotOne(congru_S1(node_sel,[1 6]),nonmem(node_sel,[1 6]),'Avg. degree');
%             plotOne(congru_S1(node_sel,[1 7]),nonmem(node_sel,[1 7]),'Avg. path length');
%         end
        if opt.plot_bar
            %% nonmem
            node_sel1=congru_S1(:,1)>20 & nonmem(:,1)>20;
            node_sel2=congru_S2(:,1)>20 & nonmem(:,1)>20;
            fh=figure('Color','w','Position',[100,100,950,200]);
            plotOneBar([congru_S1(node_sel1,1);congru_S2(node_sel2,1)],[nonmem(node_sel1,1);nonmem(node_sel2,1)],[1,7],'Node','con_non',rtype);
            plotOneBar([congru_S1(node_sel1,3);congru_S2(node_sel2,3)],[nonmem(node_sel1,3);nonmem(node_sel2,3)],[2,7],'Connected comp.','con_non',rtype);
            plotOneBar([congru_S1(node_sel1,4);congru_S2(node_sel2,4)],[nonmem(node_sel1,4);nonmem(node_sel2,4)],[3,7],'Avg. cluster coef.','con_non',rtype);
            plotOneBar([congru_S1(node_sel1,5);congru_S2(node_sel2,5)],[nonmem(node_sel1,5);nonmem(node_sel2,5)],[4,7],'Density','con_non',rtype);
            plotOneBar([congru_S1(node_sel1,6);congru_S2(node_sel2,6)],[nonmem(node_sel1,6);nonmem(node_sel2,6)],[5,7],'Avg. degree','con_non',rtype);
            plotOneBar([congru_S1(node_sel1,7);congru_S2(node_sel2,7)],[nonmem(node_sel1,7);nonmem(node_sel2,7)],[6,7],'Avg. path length','con_non',rtype);
            plotOneBar([congru_S1(node_sel1,8);congru_S2(node_sel2,8)],[nonmem(node_sel1,8);nonmem(node_sel2,8)],[7,7],'Diameter','con_non',rtype);
            if opt.exportgraph
                exportgraphics(fh,fullfile('gephidata',filenamerule(sprintf('func_conn_stats_con_non_%s.pdf',rtype))),'ContentType','vector');
            end
            if opt.save_preview
                exportgraphics(fh,fullfile('gephidata',filenamerule(sprintf('func_conn_stats_con_non_%s.png',rtype))),'ContentType','image');
            end
            
            %% incong    
            node_sel1=congru_S1(:,1)>20 & incong(:,1)>20;
            node_sel2=congru_S2(:,1)>20 & incong(:,1)>20;
            fh=figure('Color','w','Position',[100,100,950,200]);
            plotOneBar([congru_S1(node_sel1,1);congru_S2(node_sel2,1)],[incong(node_sel1,1);incong(node_sel2,1)],[1,7],'Node','con_incong',rtype);
            plotOneBar([congru_S1(node_sel1,3);congru_S2(node_sel2,3)],[incong(node_sel1,3);incong(node_sel2,3)],[2,7],'Connected comp.','con_incong',rtype);
            plotOneBar([congru_S1(node_sel1,4);congru_S2(node_sel2,4)],[incong(node_sel1,4);incong(node_sel2,4)],[3,7],'Avg. cluster coef.','con_incong',rtype);
            plotOneBar([congru_S1(node_sel1,5);congru_S2(node_sel2,5)],[incong(node_sel1,5);incong(node_sel2,5)],[4,7],'Density','con_incong',rtype);
            plotOneBar([congru_S1(node_sel1,6);congru_S2(node_sel2,6)],[incong(node_sel1,6);incong(node_sel2,6)],[5,7],'Avg. degree','con_incong',rtype);
            plotOneBar([congru_S1(node_sel1,7);congru_S2(node_sel2,7)],[incong(node_sel1,7);incong(node_sel2,7)],[6,7],'Avg. path length','con_incong',rtype);
            plotOneBar([congru_S1(node_sel1,8);congru_S2(node_sel2,8)],[incong(node_sel1,8);incong(node_sel2,8)],[7,7],'Diameter','con_incong',rtype);
            if opt.exportgraph
                exportgraphics(fh,fullfile('gephidata',filenamerule(sprintf('func_conn_stats_con_incong_%s.pdf',rtype))),'ContentType','vector');
            end
            if opt.save_preview
                exportgraphics(fh,fullfile('gephidata',filenamerule(sprintf('func_conn_stats_con_incong_%s.png',rtype))),'ContentType','image');
            end
        end
    end
end
%% node matching shuffles
if opt.plot_shuf_bar
    for rtype=["within","cross","all"]
        congru_S1=readmatrix(fullfile('gephidata',sprintf('s1_%s_gephi_graph_sums.csv',rtype)));
        congru_S2=readmatrix(fullfile('gephidata',sprintf('s2_%s_gephi_graph_sums.csv',rtype)));
        [shuf1,shuf2]=deal([]);
        for rpt=1:100
            shuf1=cat(3,shuf1,readmatrix(fullfile('gephidata',sprintf('shuf1_non_%s_%d_gephi_graph_sums.csv',rtype,rpt))));
            shuf2=cat(3,shuf2,readmatrix(fullfile('gephidata',sprintf('shuf2_non_%s_%d_gephi_graph_sums.csv',rtype,rpt))));
        end
        shufsum1=nanmean(shuf1,3);
        shufsum2=nanmean(shuf2,3);
        node_sel1=congru_S1(:,1)>20 & shufsum1(:,1)>20;
        node_sel2=congru_S2(:,1)>20 & shufsum2(:,1)>20;
        fh=figure('Color','w','Position',[100,100,950,200]);
        plotOneBar([congru_S1(node_sel1,1);congru_S2(node_sel2,1)],[shufsum1(node_sel1,1);shufsum2(node_sel2,1)],[1,7],'Node','con_shuf',rtype);
        plotOneBar([congru_S1(node_sel1,3);congru_S2(node_sel2,3)],[shufsum1(node_sel1,3);shufsum2(node_sel2,3)],[2,7],'Connected comp.','con_shuf',rtype);
        plotOneBar([congru_S1(node_sel1,4);congru_S2(node_sel2,4)],[shufsum1(node_sel1,4);shufsum2(node_sel2,4)],[3,7],'Avg. cluster coef.','con_shuf',rtype);
        plotOneBar([congru_S1(node_sel1,5);congru_S2(node_sel2,5)],[shufsum1(node_sel1,5);shufsum2(node_sel2,5)],[4,7],'Density','con_shuf',rtype);
        plotOneBar([congru_S1(node_sel1,6);congru_S2(node_sel2,6)],[shufsum1(node_sel1,6);shufsum2(node_sel2,6)],[5,7],'Avg. degree','con_shuf',rtype);
        plotOneBar([congru_S1(node_sel1,7);congru_S2(node_sel2,7)],[shufsum1(node_sel1,7);shufsum2(node_sel2,7)],[6,7],'Avg. path length','con_shuf',rtype);
        plotOneBar([congru_S1(node_sel1,8);congru_S2(node_sel2,8)],[shufsum1(node_sel1,8);shufsum2(node_sel2,8)],[7,7],'Diameter','con_shuf',rtype);
        if opt.exportgraph
            exportgraphics(fh,fullfile('gephidata',filenamerule(sprintf('func_conn_stats_con_shuf_%s.pdf',rtype))),'ContentType','vector');
        end
        if opt.save_preview
            exportgraphics(fh,fullfile('gephidata',filenamerule(sprintf('func_conn_stats_con_shuf_%s.png',rtype))),'ContentType','image');
        end
    end
end
for rtype=["within","cross","all"]
    %% raw path
    congru_path=[readmatrix(fullfile('gephidata',sprintf('s1_%s_gephi_component_path.csv',rtype)));...
        readmatrix(fullfile('gephidata',sprintf('s2_%s_gephi_component_path.csv',rtype)))];
    %         incong_path=readmatrix(fullfile('gephidata',sprintf('incong_%s_gephi_component_path.csv',rtype)));
    nonmem_path=readmatrix(fullfile('gephidata',sprintf('non_%s_gephi_component_path.csv',rtype)));
    fh=figure();
    hold on;
    scatter(congru_path(:,1),congru_path(:,2),4,'ro')
    scatter(nonmem_path(:,1),nonmem_path(:,2),4,'ko')
    xlim([0,50])

    E=[0:10:120];
    [congruY,congruE]=discretize(congru_path(:,1),E);
    [nonY,nonE]=discretize(nonmem_path(:,1),E);
    figure()
    hold on;
    for i=1:max(nonY)
        plot(i,mean(congru_path(congruY==i,2)),'ro')
        plot(i,mean(nonmem_path(nonY==i,2)),'ko')
    end

    %% shuf path
    congru_path=[readmatrix(fullfile('gephidata',sprintf('s1_%s_gephi_component_path.csv',rtype)));...
        readmatrix(fullfile('gephidata',sprintf('s2_%s_gephi_component_path.csv',rtype)))];
    shuf_path=[];
    for rpt=1:100
        shuf_path=cat(1,shuf_path,readmatrix(fullfile('gephidata',sprintf('shuf1_non_%s_%d_gephi_component_path.csv',rtype,rpt))));
        shuf_path=cat(1,shuf_path,readmatrix(fullfile('gephidata',sprintf('shuf2_non_%s_%d_gephi_component_path.csv',rtype,rpt))));
    end
%     fh=figure();
%     hold on;
%     scatter(congru_path(:,1),congru_path(:,2),4,'ro')
%     scatter(shuf_path(:,1),shuf_path(:,2),4,'ko')
%     xlim([5,100])

    E=[0:5:20,30:10:50];
    congruY=discretize(congru_path(:,1),E);
    shufY=discretize(shuf_path(:,1),E);
    discreMat=[];
    for i=1:max(shufY)
        cpath=congru_path(congruY==i,2);
        spath=shuf_path(shufY==i,2);
        cdiam=congru_path(congruY==i,3);
        sdiam=shuf_path(shufY==i,3);
        cradi=congru_path(congruY==i,4);
        sradi=shuf_path(shufY==i,4);
%% path
        discreMat=[discreMat;mean(cpath),... % path length
            bootci(500,@(x) mean(x),cpath).',...
            mean(spath),...
            bootci(500,@(x) mean(x),spath).',...
            mean(cdiam),... % diameter
            bootci(500,@(x) mean(x),cdiam).',...
            mean(sdiam),...
            bootci(500,@(x) mean(x),sdiam).',...
            mean(cradi),... % radius
            bootci(500,@(x) mean(x),cradi).',...
            mean(sradi),...
            bootci(500,@(x) mean(x),sradi).'...            
            ];
    end
    XX=diff(E,1,2)./2+E(1:end-1);
    pos=1:6:13;
    lbls={'Average shortest path length','Component diamemter','Component radius'};
    fns={'path','diameter','radius'};
    for ii=1:3
        pltIdx=pos(ii);
        altIdx=pltIdx+3;
        fh=figure('Color','w','Position',[100,100,150,200]);
        hold on
        plot(XX,discreMat(:,pltIdx),'-r.');
        plot(XX,discreMat(:,altIdx),'-k.');
        errorbar(XX,discreMat(:,pltIdx),discreMat(:,pltIdx+1)-discreMat(:,pltIdx),discreMat(:,pltIdx+2)-discreMat(:,pltIdx),'r.');
        errorbar(XX,discreMat(:,altIdx),discreMat(:,altIdx+1)-discreMat(:,altIdx),discreMat(:,altIdx+2)-discreMat(:,altIdx),'k.');
        xlabel('Nodes')
        ylabel(lbls{ii})
        ylim([0,max(ylim())]);
        set(gca(),'XTick',0:50:100)
        xlim([0,100])
        anovan([congru_path(:,ii+1);shuf_path(:,ii+1)],{[congruY;shufY],[zeros(size(congruY));ones(size(shufY))]})
        keyboard();
        if opt.exportgraph
            exportgraphics(fh,fullfile('gephidata',sprintf('func_conn_%s_con_shuf_%s.pdf',fns{ii},rtype)),'ContentType','vector');
        end
        if opt.save_preview
            exportgraphics(fh,fullfile('gephidata',sprintf('func_conn_%s_con_shuf_%s.png',fns{ii},rtype)),'ContentType','image');
        end
    end
end
%%
end

function plotOne(exp,ref,lbl)
figure()
hold on;
scatter(exp(:,1),exp(:,2),4,'r');
scatter(ref(:,1),exp(:,2),4,'b');
title(lbl);
end


function fh=plotOneBar(congru,incongru,subidx,lbl1,lbl2,lbl3)
arguments
    congru
    incongru
    subidx (1,2) double {mustBeInteger,mustBePositive}
    lbl1 (1,:) char
    lbl2 (1,:) char
    lbl3 (1,:) char
end
subplot(1,subidx(2),subidx(1));
hold on;
if size(incongru,1)*2==size(congru,1)
    incongru=repmat(incongru,2,1);
end
ydata=[incongru,congru];
xdata=ones(size(ydata));
if contains(lbl1,'comp.')
    xdata(:,1)=xdata(:,1)+rand(size(ydata,1),1)*0.3;
    xdata(:,2)=xdata(:,2)+1-rand(size(ydata,1),1)*0.3;
else
    xdata(:,1)=xdata(:,1)+rand(size(ydata,1),1)*0.2;
    xdata(:,2)=xdata(:,2)+1-rand(size(ydata,1),1)*0.2;
end
plot(xdata',ydata','-k')
plot(xdata(:,1),ydata(:,1),'ko','MarkerSize',4,'MarkerFaceColor','none','MarkerEdgeColor','k')
plot(xdata(:,2),ydata(:,2),'ro','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','k')
errorbar(0.75,mean(ydata(:,1)),std(ydata(:,1))/sqrt(size(ydata,1)),'ko','MarkerSize',6);
errorbar(2.25,mean(ydata(:,2)),std(ydata(:,2))/sqrt(size(ydata,1)),'ro','MarkerSize',6);
xlim([0.5,2.5])
ylim([0,max(ylim())])
ylabel(strjoin({lbl2,lbl3,lbl1}),'Interpreter','none');
xlabel(sprintf('%0.4f',signrank(congru,incongru)));
set(gca,'XTick',[])
end

function out=filenamerule(s)
out=replace(replace(s,' ','_'),'.','_');
out(end-3)='.';
end
