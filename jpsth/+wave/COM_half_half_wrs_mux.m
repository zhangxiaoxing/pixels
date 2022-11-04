function [fh,r_olf_dur_mux]=COM_half_half_wrs_mux(opt)
arguments
    opt.filename (1,:) char = 'com_halfs_100.mat'
    opt.err_filename (1,:) char = 'com_error.mat'
    opt.plot (1,1) logical = true
end

fstr=load(opt.filename);
estr=load(opt.err_filename);
r_olf_dur_mux=stats_file(fstr.com_halfs,estr.com_map_err); %correct shuffle error

gg=[repmat(1,1,100),repmat(3,1,100),repmat(5,1,100),...
    repmat(7,1,300),...
    repmat(2,1,100),repmat(4,1,100),repmat(6,1,100)]; % group vector for box plot

fh=figure('Color','w','Position',[32,32,275,235]);
hold on
boxplot(reshape(r_olf_dur_mux,1,[]),gg);
ylabel('TCOM Pearson''s r')
set(gca(),'XTick',1:7,'XTickLabel',{'OlfC','OlfE','DurC','DurE','MuxC','MuxE','Shuffle'},'YTick',0:0.5:1);
yline(0,'-k')
xlim([0.5,7.5]);
ylim([-0.1,1.25])

pwrs=arrayfun(@(x) ranksum(r_olf_dur_mux(:,x),r_olf_dur_mux(:,x+6)), 1:3);%correct error
panova=anovan(reshape(r_olf_dur_mux,[],1),gg.','display','off'); % cross all grp
title([sprintf('wrs%.3f,',pwrs),sprintf('anova%.3f',panova)])
% fh=figure('Color','w','Position',[32,32,275,235]);
% hold on
% swarmchart(gg,reshape(r_olf_dur_mux,1,[]),1);
end

function r_olf_dur_mux=stats_file(com_halfs,com_err)
r_olf_dur_mux=nan(size(com_halfs,1),9);
for rpti=1:size(com_halfs)
    com_1h=com_halfs{rpti,1};
    com_2h=com_halfs{rpti,2};
    sesses=reshape(intersect(fieldnames(com_1h),fieldnames(com_2h)),1,[]);
    olf_com_mat=[];
    dur_com_mat=[];
    mux_com_mat=[];

    err_olf_com_mat=[];
    err_dur_com_mat=[];
    err_mux_com_mat=[];

    for sess=sesses
        fs=sess{1};
        %% olf, dur, mix
        for stype=["olf_s1","olf_s2"]
            if isfield(com_1h.(fs),stype) && isfield(com_2h.(fs),stype)
                sukeys=num2cell(intersect(cell2mat(com_1h.(fs).(stype).com.keys),...
                    cell2mat(com_2h.(fs).(stype).com.keys))); 
                olf_com_mat=[olf_com_mat;cell2mat(com_1h.(fs).(stype).com.values(sukeys)).',...
                    cell2mat(com_2h.(fs).(stype).com.values(sukeys)).'];
            end
            %% error trial
            % 1h<>err, 2h<>err, (skip or placeholder?)
            if isfield(com_err,fs) && isfield(com_1h.(fs),stype) && isfield(com_err.(fs),stype)
                sukeys=num2cell(intersect(cell2mat(com_1h.(fs).(stype).com.keys),...
                    cell2mat(com_err.(fs).(stype).com.keys)));
                err_olf_com_mat=[err_olf_com_mat;cell2mat(com_1h.(fs).(stype).com.values(sukeys)).',...
                    cell2mat(com_err.(fs).(stype).com.values(sukeys)).'];
            end

        end
        for stype=["dur_d3","dur_d6"]
            if isfield(com_1h.(fs),stype) && isfield(com_2h.(fs),stype)
                sukeys=num2cell(intersect(cell2mat(com_1h.(fs).(stype).com.keys),...
                    cell2mat(com_2h.(fs).(stype).com.keys))); 
                dur_com_mat=[dur_com_mat;cell2mat(com_1h.(fs).(stype).com.values(sukeys)).',...
                    cell2mat(com_2h.(fs).(stype).com.values(sukeys)).'];
            end

            if isfield(com_err,fs) && isfield(com_1h.(fs),stype) && isfield(com_err.(fs),stype)
                sukeys=num2cell(intersect(cell2mat(com_1h.(fs).(stype).com.keys),...
                    cell2mat(com_err.(fs).(stype).com.keys)));
                err_dur_com_mat=[err_dur_com_mat;cell2mat(com_1h.(fs).(stype).com.values(sukeys)).',...
                    cell2mat(com_err.(fs).(stype).com.values(sukeys)).'];
            end

        end
        for stype=["s1d3","s1d6","s2d3","s2d6"]
            if isfield(com_1h.(fs),stype) && isfield(com_2h.(fs),stype)
                sukeys=num2cell(intersect(cell2mat(com_1h.(fs).(stype).com.keys),...
                    cell2mat(com_2h.(fs).(stype).com.keys))); % TODO possible intersect error trials 
                mux_com_mat=[mux_com_mat;cell2mat(com_1h.(fs).(stype).com.values(sukeys)).',...
                    cell2mat(com_2h.(fs).(stype).com.values(sukeys)).'];
            end
            if isfield(com_err,fs) && isfield(com_1h.(fs),stype) && isfield(com_err.(fs),stype)
                sukeys=num2cell(intersect(cell2mat(com_1h.(fs).(stype).com.keys),...
                    cell2mat(com_err.(fs).(stype).com.keys))); % TODO possible intersect error trials 
                err_mux_com_mat=[err_mux_com_mat;cell2mat(com_1h.(fs).(stype).com.values(sukeys)).',...
                    cell2mat(com_err.(fs).(stype).com.values(sukeys)).'];
            end

        end
    end

    r_olf_dur_mux(rpti,1)=corr(olf_com_mat(:,1),olf_com_mat(:,2));
    r_olf_dur_mux(rpti,2)=corr(dur_com_mat(:,1),dur_com_mat(:,2));
    r_olf_dur_mux(rpti,3)=corr(mux_com_mat(:,1),mux_com_mat(:,2));
    r_olf_dur_mux(rpti,4)=corr(olf_com_mat(:,1),randsample(olf_com_mat(:,2),size(olf_com_mat,1)));
    r_olf_dur_mux(rpti,5)=corr(dur_com_mat(:,1),randsample(dur_com_mat(:,2),size(dur_com_mat,1)));
    r_olf_dur_mux(rpti,6)=corr(mux_com_mat(:,1),randsample(mux_com_mat(:,2),size(mux_com_mat,1)));
    r_olf_dur_mux(rpti,7)=corr(err_olf_com_mat(:,1),err_olf_com_mat(:,2));
    r_olf_dur_mux(rpti,8)=corr(err_dur_com_mat(:,1),err_dur_com_mat(:,2));
    r_olf_dur_mux(rpti,9)=corr(err_mux_com_mat(:,1),err_mux_com_mat(:,2));


end
end

