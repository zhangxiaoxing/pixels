function fh=COM_half_half_wrs_mux(opt)
arguments
    opt.filename (1,:) char = 'com_halfs_100.mat'
    opt.plot (1,1) logical = true
end
fstr=load(opt.filename);
r_olf_dur_mux=stats_file(fstr.com_halfs);

gg=[repmat(1,1,100),repmat(2,1,100),repmat(3,1,100),repmat(4,1,300)];

fh=figure('Color','w','Position',[32,32,275,235]);
hold on
boxplot(reshape(r_olf_dur_mux,1,[]),gg);
ylabel('TCOM Pearson''s r')
set(gca(),'XTick',1:4,'XTickLabel',{'OLF','DUR','Mux','Shuffle'},'YTick',0:0.5:1);
yline(0,'-k')
xlim([0.5,4.5]);
ylim([-0.1,1])

end

function r_olf_dur_mux=stats_file(com_halfs)
r_olf_dur_mux=nan(size(com_halfs,1),6);
for rpti=1:size(com_halfs)
    com_1h=com_halfs{rpti,1};
    com_2h=com_halfs{rpti,2};
    sesses=reshape(intersect(fieldnames(com_1h),fieldnames(com_2h)),1,[]);
    olf_com_mat=[];
    dur_com_mat=[];
    mux_com_mat=[];
    for sess=sesses
        fs=sess{1};
        %% olf, dur, mix
        for stype=["olf_s1","olf_s2"]
            if isfield(com_1h.(fs),stype) && isfield(com_2h.(fs),stype)
                sukeys=num2cell(intersect(cell2mat(com_1h.(fs).(stype).com.keys),...
                    cell2mat(com_2h.(fs).(stype).com.keys))); % TODO possible intersect error trials 
                olf_com_mat=[olf_com_mat;cell2mat(com_1h.(fs).(stype).com.values(sukeys)).',...
                    cell2mat(com_2h.(fs).(stype).com.values(sukeys)).'];
            end
        end
        for stype=["dur_d3","dur_d6"]
            if isfield(com_1h.(fs),stype) && isfield(com_2h.(fs),stype)
                sukeys=num2cell(intersect(cell2mat(com_1h.(fs).(stype).com.keys),...
                    cell2mat(com_2h.(fs).(stype).com.keys))); % TODO possible intersect error trials 
                dur_com_mat=[dur_com_mat;cell2mat(com_1h.(fs).(stype).com.values(sukeys)).',...
                    cell2mat(com_2h.(fs).(stype).com.values(sukeys)).'];
            end
        end
        for stype=["s1d3","s1d6","s2d3","s2d6"]
            if isfield(com_1h.(fs),stype) && isfield(com_2h.(fs),stype)
                sukeys=num2cell(intersect(cell2mat(com_1h.(fs).(stype).com.keys),...
                    cell2mat(com_2h.(fs).(stype).com.keys))); % TODO possible intersect error trials 
                mux_com_mat=[mux_com_mat;cell2mat(com_1h.(fs).(stype).com.values(sukeys)).',...
                    cell2mat(com_2h.(fs).(stype).com.values(sukeys)).'];
            end
        end
    end

    r_olf_dur_mux(rpti,1)=corr(olf_com_mat(:,1),olf_com_mat(:,2));
    r_olf_dur_mux(rpti,2)=corr(dur_com_mat(:,1),dur_com_mat(:,2));
    r_olf_dur_mux(rpti,3)=corr(mux_com_mat(:,1),mux_com_mat(:,2));
    r_olf_dur_mux(rpti,4)=corr(olf_com_mat(:,1),randsample(olf_com_mat(:,2),size(olf_com_mat,1)));
    r_olf_dur_mux(rpti,5)=corr(dur_com_mat(:,1),randsample(dur_com_mat(:,2),size(dur_com_mat,1)));
    r_olf_dur_mux(rpti,6)=corr(mux_com_mat(:,1),randsample(mux_com_mat(:,2),size(mux_com_mat,1)));
end
end

