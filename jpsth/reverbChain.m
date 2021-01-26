fstr=cell(1,6);
for bin=1:6
    fstr{bin}=load(sprintf('0116_memory_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
end
onesample(fstr,true,1);
onesample(fstr,true,2);
% onesample(2,fstr);


function onesample(fstr,to_plot,prefidx)
if ~exist('to_plot','var')
    to_plot=false;
end

if ~exist('prefidx','var')
    prefidx=1;
end

% for bin=1:6
%     disp(length(fstr{bin}.conn_chain_S1))
%     disp(length(fstr{bin}.conn_chain_S1)/length(fstr{bin}.pair_chain))
% end


% prefsamp=1;
    prefsamp=1;
    bin1sel1=(fstr{1}.(sprintf('pref_chain_S%d',prefidx))(:,1)==prefsamp & all(fstr{1}.(sprintf('pref_chain_S%d',prefidx))(:,7:8)==prefsamp,2));
    bin1postu=unique(fstr{1}.(sprintf('conn_chain_S%d',prefidx))(bin1sel1,2));

    bin2sel1=ismember(fstr{2}.(sprintf('conn_chain_S%d',prefidx))(:,1),bin1postu) & all(fstr{2}.(sprintf('pref_chain_S%d',prefidx))(:,8:9)==prefsamp,2);
    bin2postu=unique(fstr{2}.(sprintf('conn_chain_S%d',prefidx))(bin2sel1,2));

    bin3sel1=ismember(fstr{3}.(sprintf('conn_chain_S%d',prefidx))(:,1),bin2postu) & all(fstr{3}.(sprintf('pref_chain_S%d',prefidx))(:,9:10)==prefsamp,2);
    bin3postu=unique(fstr{3}.(sprintf('conn_chain_S%d',prefidx))(bin3sel1,2));

    bin4sel1=ismember(fstr{4}.(sprintf('conn_chain_S%d',prefidx))(:,1),bin3postu) & all(fstr{4}.(sprintf('pref_chain_S%d',prefidx))(:,10:11)==prefsamp,2);
    bin4postu=unique(fstr{4}.(sprintf('conn_chain_S%d',prefidx))(bin4sel1,2));

    bin5sel1=ismember(fstr{5}.(sprintf('conn_chain_S%d',prefidx))(:,1),bin4postu) & all(fstr{5}.(sprintf('pref_chain_S%d',prefidx))(:,11:12)==prefsamp,2);
    
    prefsamp=2;
    bin1sel2=(fstr{1}.(sprintf('pref_chain_S%d',prefidx))(:,1)==prefsamp & all(fstr{1}.(sprintf('pref_chain_S%d',prefidx))(:,7:8)==prefsamp,2));
    bin1postu=unique(fstr{1}.(sprintf('conn_chain_S%d',prefidx))(bin1sel2,2));

    bin2sel2=ismember(fstr{2}.(sprintf('conn_chain_S%d',prefidx))(:,1),bin1postu) & all(fstr{2}.(sprintf('pref_chain_S%d',prefidx))(:,8:9)==prefsamp,2);
    bin2postu=unique(fstr{2}.(sprintf('conn_chain_S%d',prefidx))(bin2sel2,2));

    bin3sel2=ismember(fstr{3}.(sprintf('conn_chain_S%d',prefidx))(:,1),bin2postu) & all(fstr{3}.(sprintf('pref_chain_S%d',prefidx))(:,9:10)==prefsamp,2);
    bin3postu=unique(fstr{3}.(sprintf('conn_chain_S%d',prefidx))(bin3sel2,2));

    bin4sel2=ismember(fstr{4}.(sprintf('conn_chain_S%d',prefidx))(:,1),bin3postu) & all(fstr{4}.(sprintf('pref_chain_S%d',prefidx))(:,10:11)==prefsamp,2);
    bin4postu=unique(fstr{4}.(sprintf('conn_chain_S%d',prefidx))(bin4sel2,2));

    bin5sel2=ismember(fstr{5}.(sprintf('conn_chain_S%d',prefidx))(:,1),bin4postu) & all(fstr{5}.(sprintf('pref_chain_S%d',prefidx))(:,11:12)==prefsamp,2);
    
%%
if to_plot
    cmap=[1,1,1;1,0,0];
    fh=figure('Color','w','Position',[100,100,200,1000]);
    subplot(6,1,1)
    imagesc(fstr{1}.(sprintf('pref_chain_S%d',prefidx))(bin1sel1|bin1sel2,1:6))
    colormap(cmap);
    ax1=gca();
    subplot(6,1,2)
    imagesc(fstr{1}.(sprintf('pref_chain_S%d',prefidx))(bin1sel1|bin1sel2,7:12))
    colormap(cmap);
    ax2=gca();
    subplot(6,1,3)
    imagesc(fstr{2}.(sprintf('pref_chain_S%d',prefidx))(bin2sel1|bin2sel2,7:12))
    colormap(cmap);
    ax3=gca();
    subplot(6,1,4)
    imagesc(fstr{3}.(sprintf('pref_chain_S%d',prefidx))(bin3sel1|bin3sel2,7:12))
    colormap(cmap);
    ax4=gca();
    subplot(6,1,5)
    imagesc(fstr{4}.(sprintf('pref_chain_S%d',prefidx))(bin4sel1|bin4sel2,7:12))
    colormap(cmap);
    ax5=gca();
    subplot(6,1,6)
    imagesc(fstr{5}.(sprintf('pref_chain_S%d',prefidx))(bin5sel1|bin5sel2,7:12))
    colormap(cmap);
    ax6=gca();
    % linkaxes([ax1 ax2 ax3 ax4 ax5 ax6])
    for ax=[ax1 ax2 ax3 ax4 ax5 ax6]
        ax.XTick=[];
        ax.YTick=[1,max(ax.YTick)];
    end
    ax6.XTick=1:6;
    ax6.XLabel.String='Time bin (sec)';
%     keyboard
    exportgraphics(fh,sprintf('prolonged_showcase_%d.pdf',prefidx),'ContentType','vector');
end
end

% conn_pref_mat=[fstr{1}.pref_chain_S1(bin1sel,1:6);...
%     nan(25,6);...
%     fstr{1}.pref_chain_S1(bin1sel,7:12);...
%     nan(25,6);...
%     fstr{2}.pref_chain_S1(bin2sel,7:12);...
%     nan(25,6);...
%     fstr{3}.pref_chain_S1(bin3sel,7:12);...
%     nan(25,6);...
%     fstr{4}.pref_chain_S1(bin4sel,7:12);...
%     nan(25,6);...
%     fstr{5}.pref_chain_S1(bin5sel,7:12)];
% % close all
% fh=figure();
% imagesc(conn_pref_mat);
% % exportgraphics(fh,'enhance_showcase_v01.pdf','ContentType','vector')
