if false
    if isempty(su_meta)
        load(fullfile('binary','su_meta.mat'))
    end

    if isempty(sel_meta)
        fstr=load(fullfile('binary','wrs_mux_meta.mat'));
        sel_meta=fstr.wrs_mux_meta;
        clear fstr
    end

    reg_com_maps=wave.get_reg_com_maps(sel_meta);


    load(fullfile('binary','sums_ring_stats_all.mat'),'sums_all');
    usess=unique(cell2mat([sums_all{1}(:,1);sums_all{2}(:,1);sums_all{3}(:,1)]));
end
loopssums=[];

for sess=reshape(usess,1,[])
    sesssel=su_meta.sess==sess;
    sessreg=su_meta.reg_tree(5,sesssel);
    sesswave=sel_meta.wave_id(sesssel);
    sesscid=uint16(su_meta.allcid(sesssel));
    wavedict=dictionary(sesscid,sesswave);
    regdict=dictionary(sesscid,sessreg.');
    for delay=[3 6]
        disp(sess)
        for rsize=3:5
            col=sums_all{rsize-2};
            rsel=cell2mat(col(:,1))==sess;
            for ridx=reshape(find(rsel),1,[])
                currreg=regdict(col{ridx,3});
                currwave=wavedict(col{ridx,3});

                %% MEM OR NONMEM
                if all(currwave==0)
                    currmem="nonmem";
                elseif delay==3 && (all(ismember(currwave,[1 5]),'all')...
                        ||all(ismember(currwave,[3 6]),'all'))
                    currmem="mem3s";
                elseif delay==6 && (all(ismember(currwave,[2 5]),'all')...
                        ||all(ismember(currwave,[4 6]),'all'))
                    currmem="mem6s";
                end

                %% WITHIN, CONSIST, INCON, MISSING
                if any(ismissing(currreg),'all')
                    currrelation="missing";
                elseif numel(unique(currreg))==1
                    currrelation="within";
                elseif delay==3 && all(reg_com_maps.tcom3_maps.odor_only.isKey(currreg),'all')
                    tcoms=cell2mat(reg_com_maps.tcom3_maps.odor_only.values(currreg));
                    congrucount=[arrayfun(@(x) tcoms(x)>tcoms(x-1),2:numel(tcoms)),tcoms(1)>tcoms(end)];
                    incongcount=[arrayfun(@(x) tcoms(x)<tcoms(x-1),2:numel(tcoms)),tcoms(1)<tcoms(end)];
                    if nnz(congrucount)>nnz(incongcount)
                        currrelation="consist";
                    elseif nnz(congrucount)<nnz(incongcount)
                        currrelation="inconsist";
                    else
                        currrelation="balanced";
                    end
                elseif delay==6 && all(reg_com_maps.tcom6_maps.odor_only.isKey(currreg),'all')
                    tcoms=cell2mat(reg_com_maps.tcom6_maps.odor_only.values(currreg));
                    congrucount=[arrayfun(@(x) tcoms(x)>tcoms(x-1),2:numel(tcoms)),tcoms(1)>tcoms(end)];
                    incongcount=[arrayfun(@(x) tcoms(x)<tcoms(x-1),2:numel(tcoms)),tcoms(1)<tcoms(end)];
                    if nnz(congrucount)>nnz(incongcount)
                        currrelation="consist";
                    elseif nnz(congrucount)<nnz(incongcount)
                        currrelation="inconsist";
                    else
                        currrelation="balanced";
                    end
                else
                    currrelation="undefined";
                end
                loopssums=[loopssums;cell2table(...
                    {currmem,delay,sess,rsize,col(ridx,3),{currwave},{currreg},currrelation},...
                    "VariableNames",{'memory','delay','session','size','cid','wave','reg','relation'})];...
            end
        end
    end
end

nonmemcnt=[nnz(loopssums.memory=="nonmem" & loopssums.relation=="within"),...
nnz(loopssums.memory=="nonmem" & loopssums.relation=="consist"),...
nnz(loopssums.memory=="nonmem" & loopssums.relation=="inconsist"),...
nnz(loopssums.memory=="nonmem" & loopssums.relation=="balanced")];

memcnt=[nnz(loopssums.memory~="nonmem" & loopssums.relation=="within"),...
nnz(loopssums.memory~="nonmem" & loopssums.relation=="consist"),...
nnz(loopssums.memory~="nonmem" & loopssums.relation=="inconsist"),...
nnz(loopssums.memory~="nonmem" & loopssums.relation=="balanced")];

figure()
tiledlayout(1,2)
nexttile
bar(nonmemcnt)
title('Non-memory')
ylabel("Number of occurrence")
set(gca,'XTick',1:4,'XTickLabel',{'Within','Consistent','Inconsistent','Balanced'})
nexttile
bar(memcnt)
title('Memory')
ylabel("Number of occurrence")
set(gca,'XTick',1:4,'XTickLabel',{'Within','Consistent','Inconsistent','Balanced'})
sgtitle('Loops SC direction vs. wave')



