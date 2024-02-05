MAX_PATH_NUM=1000;

global_init;
wt_sig=bz.load_sig_sums_conn_file('pair',false,'inhibit',false,'criteria','WT');
usess=unique(wt_sig.sess);
if read_chain
    error("N/A yet")
else
    all_chain=cell(0);
    for sidx=reshape(usess,1,[])
        if rem(sidx,10)==0
            disp(sidx)
        end
        ssel=wt_sig.sess==sidx;
        ssuid=wt_sig.suid(ssel,:);
        gh=digraph(table(string(ssuid),'VariableNames',{'EndNodes'}));
        % plot(gh,'Layout','circle')
        gnodes=gh.Nodes.Name;
        pairs_unidir=nchoosek(gnodes,2);
        pairs=[pairs_unidir;fliplr(pairs_unidir)];
        for ii=1:height(pairs)
            path=gh.allpaths(pairs(ii,1),pairs(ii,2),'MinPathLength',3,'MaxNumPaths',MAX_PATH_NUM);
            pnum=numel(path);
            if ~isempty(path)
                % all_chain.add({sidx,path,numel(path)==MAX_PATH_NUM});
                for jj=1:height(path)
                    all_chain=[all_chain;{sidx,ii,path{jj},pnum==MAX_PATH_NUM}];
                end
            end
        end
    end
end
blame=vcs.blame();
save(fullfile('binary','digraph_chains.mat'),'all_chain','blame')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for sidx=reshape(usess,1,[])
    if rem(sidx,10)==0
        disp(sidx)
    end
    ssel=wt_sig.sess==sidx;
    ssuid=wt_sig.suid(ssel,:);
    gh=digraph(table(string(ssuid),'VariableNames',{'EndNodes'}));
    keyboard()
    if gh.hascycles
        keyboard()
    end
end


plist=cell(all_chain.toArray());

