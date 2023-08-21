function [per_sess_struct,per_sess_mat]=stats_replay_sess(stats_all,opt)
arguments
    stats_all cell
    opt.feat_sel = []
    % opt.two_or_more (1,1) logical = false
end

per_sess_struct=struct();
for ii=1:numel(stats_all)
    motif_set=stats_all{ii};
    if isempty(opt.feat_sel)
        featsel=true(size(motif_set.count(:,1)));
    else
        featsel=opt.feat_sel;
    end
    %
    % if opt.two_or_more
    %     [gc,gr]=groupcounts(motif_set.condition.');
    %     gcmap=containers.Map(gr,num2cell(gc));
    % end

    for jj=1:size(motif_set.count,2)
        onekey=motif_set.condition{jj};
        % if opt.two_or_more && gcmap(onekey)<=2
        %     continue
        % end
        if ~isfield(per_sess_struct,onekey)
            per_sess_struct.(onekey)=cell2struct({motif_set.count(featsel,jj).';motif_set.time(featsel,jj).';motif_set.tag(jj)},{'count';'time';'tag'});
        else
            if isequaln(per_sess_struct.(onekey).time,motif_set.time(featsel,jj).')
                per_sess_struct.(onekey).count=[per_sess_struct.(onekey).count;motif_set.count(featsel,jj).'];
                per_sess_struct.(onekey).tag=[per_sess_struct.(onekey).tag;motif_set.tag(jj)];
            else
                disp("time does not match")
                keyboard()
            end
        end
    end
end

outcell=struct2cell(per_sess_struct);
per_sess_mat=cell2mat(cellfun(@(x) sum(x.count,1)./x.time,outcell,'UniformOutput',false));
end
