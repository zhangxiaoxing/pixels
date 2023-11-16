% generate data to feed gephi for the chained loops illustration
% 2023.2.6

% chain + loop
global_init;
sig=bz.load_sig_sums_conn_file('pair',false);
su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
greys=ephys.getGreyRegs('range','grey');
wrs_mux_meta=ephys.get_wrs_mux_meta();
com_map=wave.get_pct_com_map(wrs_mux_meta,'early_smooth',false);


load(fullfile('bzdata','sums_ring_stats_all.mat'));% 1X3
rstats=bz.rings.rings_reg_pie(sums_all,'plot',false);% 1X3
clear sums_all
rstats(:,5)=num2cell(cellfun(@(x) all(ismember(x,greys),"all"),rstats(:,8)));

load(fullfile('bzdata','chains_mix.mat'),'chains_uf');
chains=chains_uf;
clear chains_uf;
chains.len=cellfun(@(x) numel(x),chains.cids);
chains.reg=cell(size(chains.sess));
for sess=reshape(unique(chains.sess),1,[])
    sesscid=su_meta.allcid(su_meta.sess==sess);
    sessreg=su_meta.reg_tree(5,su_meta.sess==sess);
    sess_reg_map=containers.Map(num2cell(sesscid),sessreg);
    for cc=reshape(find(chains.sess==sess),1,[])
        chains.reg{cc}=sess_reg_map.values(num2cell(chains.cids{cc}));
    end
end
chains.regsel=cellfun(@(x) all(ismember(x,greys),"all"),chains.reg);


if false % for y axis brain region order, used selective percentage instead
    tcom3_maps=cell(1,3);
    for typeidx=1:3
        type=subsref(["mixed","olf","dur"],struct(type='()',subs={{typeidx}}));
        [fcom3.(type).collection,fcom3.(type).com_meta]=wave.per_region_COM(...
            com_map,'sel_type',type,'com_field','com3');
        ureg=intersect(ephys.getGreyRegs('range','grey'),...
            fcom3.(type).collection(:,2));
        [~,tcidx]=ismember(ureg,fcom3.(type).collection(:,2));
        tcom3_maps{typeidx}=containers.Map(...
            ureg,num2cell(cellfun(@(x) x/4, fcom3.(type).collection(tcidx,1))));
    end


    tcom6_maps=cell(1,3);
    for typeidx=1:3
        type=subsref(["mixed","olf","dur"],struct(type='()',subs={{typeidx}}));
        [fcom6.(type).collection,fcom6.(type).com_meta]=wave.per_region_COM(...
            com_map,'sel_type',type,'com_field','com6');
        ureg=intersect(ephys.getGreyRegs('range','grey'),...
            fcom6.(type).collection(:,2));
        [~,tcidx]=ismember(ureg,fcom6.(type).collection(:,2));
        tcom6_maps{typeidx}=containers.Map(...
            ureg,num2cell(cellfun(@(x) x/4, fcom6.(type).collection(tcidx,1))));
    end
end

[map_cells,~]=ephys.pct_reg_bars(wrs_mux_meta,'skip_plot',true); % only need map_cells for tcom-frac corr
reg_prop=subsref(cell2mat(map_cells{2}.values(ureg)),struct(type={'()'},subs={{':',1}}));
[reg_prop,tcomidx]=sort(reg_prop);

% showcase
% with both wave

wtags=["s1d3","s1d6","s2d3","s2d6"];
for wid=1:4
    cwsel=chains.wave==wtags(wid);
    rwsel=cellfun(@(x) ismember(wid,x),rstats(:,7));
    chain_sess_all=unique(chains.sess(cwsel));
    ring_sess_all=unique([rstats{rwsel,1}]);
    sess_all=intersect(chain_sess_all,ring_sess_all);
    for sess=reshape(sess_all,1,[])

        sesscid=su_meta.allcid(su_meta.sess==sess);
        sessreg=su_meta.reg_tree(5,su_meta.sess==sess);
        sess_reg_map=containers.Map(num2cell(sesscid),sessreg);

        if nnz(cell2mat(rstats(:,1))==sess & rwsel)>3
            % TODO join congurent olf & dur wave
            ring_id=unique([rstats{cell2mat(rstats(:,1))==sess & rwsel,3}]);
            chain_id=unique([chains.cids{chains.sess==sess & cwsel & chains.len>4 & chains.regsel}]);
            
            sessid=unique([ring_id,chain_id]);
            if numel(sessid)<30
                continue
            end

            disp(wid)
            disp(sess)

            sigsel=sig.sess==sess & all(ismember(sig.suid,sessid),2);

            edges=sig.suid(sigsel,:);
            edges(:,3:5)=0;
            %export to gephi
%             csvcell={'Source','Target','A','B','C'};
            
            for ii=reshape(find(cell2mat(rstats(:,1))==sess & rwsel),1,[]) %ring
                onering=rstats{ii,3};
                for nn=1:numel(onering)-1
                    edges(edges(:,1)==onering(nn) & edges(:,2)==onering(nn+1),3)=1;
                end
                edges(edges(:,1)==onering(end) & edges(:,2)==onering(1),3)=1;
            end

            for ii=reshape(find(chains.sess==sess & cwsel & chains.len>4),1,[]) %chain
                onechain=chains.cids{ii};
                for nn=1:numel(onechain)-1
                    edges(edges(:,1)==onechain(nn) & edges(:,2)==onechain(nn+1),4)=1;
                end
            end
            edges(:,5)=edges(:,3)+2*edges(:,4);
            csvcell=[{'Source','Target','Loop','Chain','Joint'};num2cell(edges(edges(:,5)>0,:))];
            writecell(csvcell,fullfile('bzdata',sprintf('conn4gephi_s%dw_%d.csv',sess,wid)));
            if ismember(wid,[1 3])
                sesscommap=[com_map.("s"+num2str(sess)).(wtags(wid)).com3;...
                com_map.("s"+num2str(sess)).olf_s1.com3;...
                com_map.("s"+num2str(sess)).olf_s2.com3;...
                com_map.("s"+num2str(sess)).dur_d3.com3];
            else
                sesscommap=[com_map.("s"+num2str(sess)).(wtags(wid)).com6;...
                com_map.("s"+num2str(sess)).olf_s1.com6;...
                com_map.("s"+num2str(sess)).olf_s2.com6;...
                com_map.("s"+num2str(sess)).dur_d3.com6];
            end

            regtcom=cell2mat(tcom3_maps{2}.values(sess_reg_map.values(num2cell(sessid))));
            
%             keyboard()
            nodecell=[{'Id','Label','XX','YY'};...
            num2cell([sessid;sessid;cell2mat(sesscommap.values(num2cell(sessid)));regtcom.*8].')];
            writecell(nodecell,fullfile('bzdata',sprintf('node4gephi_s%dw_%d.csv',sess,wid)));
        end
    end
end

