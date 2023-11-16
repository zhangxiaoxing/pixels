% only evaluate SU-waveid, spike-time independent

% nonmem,overlap,independent,within,cross,ringsize,session,cids,[nan]
function sums=rings_wave(sel_meta,opt)
arguments
    sel_meta
    opt.shufid double {mustBeScalarOrEmpty} = 0
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
    opt.odor_only (1,1) logical
end

% persistent meta rings_shuf

switch opt.criteria
    case 'WT'
        load(fullfile('binary','su_meta.mat'),'su_meta');
    case 'Learning'
        su_meta=ephys.util.load_meta("save_file",false,"adjust_white_matter",true,"criteria","Learning","load_file",false,"skip_stats",true);
    case 'any'
        error("Unfinished")
end


if isempty(sel_meta)
    switch opt.criteria
        case 'WT'
            fstr=load(fullfile('binary','wrs_mux_meta.mat'));
            sel_meta=fstr.wrs_mux_meta;
            clear fstr
        case 'Learning'
            sel_meta=ephys.get_wrs_mux_meta('load_file',false,'save_file',false,'criteria','Learning','extend6s',true);
        case 'any'
            error("Unfinished")
    end
end
switch opt.criteria
    case 'WT'
        fstr=load(fullfile('binary','rings_bz_vs_shuf.mat'));
    case 'Learning'
        fstr=load(fullfile('binary','LN_rings_bz_vs_shuf.mat'));
    otherwise
        error("Unfinished");
end

if opt.shufid==0
    rings=fstr.rings;
else
    rings=fstr.rings_shuf{opt.shufid};
end
% TODO actual data
meta_sess=-1;
sums=cell(0,7);
for sess=1:size(rings,1)
    for rsidx=1:3
        if isempty(rings{sess,rsidx})
            continue
        end
        % TODO: join_wave, join_reg
        if meta_sess~=sess
            meta_sess=sess;
            sesscid=su_meta.allcid(su_meta.sess==meta_sess);
            sessreg=su_meta.reg_tree(5,su_meta.sess==meta_sess);
            sess_reg_map=containers.Map(num2cell(sesscid),sessreg);
            sesswaveid=sel_meta.wave_id(su_meta.sess==meta_sess);
            sess_wave_map=containers.Map(num2cell(sesscid),num2cell(sesswaveid));

        end
        onesess=rings{sess,rsidx};
        for ridx=1:size(onesess,1)
%             disp([sess,rsidx,ridx,size(sums,1)]);
            curr_waveid=cell2mat(sess_wave_map.values(num2cell(onesess(ridx,:))));
            curr_reg=sess_reg_map.values(num2cell(onesess(ridx,:)));
            if all(~ismissing(curr_reg),"all")
                if numel(unique(curr_reg))==1
                    reg_type='same';
                else
                    reg_type='cross';
                end
            else
                reg_type='missing';
            end
            [rwid,sel_type]=bz.rings.ring_wave_type(curr_waveid,'odor_only',opt.odor_only);
            sums=[sums;{sess,rsidx+2,onesess(ridx,:),rwid,sel_type,curr_waveid,curr_reg,reg_type}];
        end
%         within_sel=reg_sel & all(reg_tagged,2) & all(reg_tagged(:,2:end)==reg_tagged(:,1),2);
%         cross_sel=reg_sel & all(reg_tagged,2) & any(reg_tagged(:,2:end)~=reg_tagged(:,1),2);
%         supp_stats=[supp_stats;ctx_sel,...%ctx
%             reg_sel & any(wave_ids>0,2) & any(wave_ids==0,2),...%mem-nonmem
%             reg_sel & any(ismember(wave_ids,[1 3 5]),2) & any(ismember(wave_ids,[2 4 6]),2) & all(wave_ids~=0,2),... % different sample
%             reg_sel & (all(ismember(wave_ids,[1 3 5]),2) | all(ismember(wave_ids,[2 4 6]),2)),...% same sample, if early = false;
%             reg_sel & (all(wave_ids==5,2) | all(wave_ids==6,2))];% both
%         % congru-odor congru-dur incongru nonmem
%         % x within, cross
%         % TODO: per region 
% 
%         % nonmem,overlap,independent,within,cross,ringsize,session,cids,[nan]
%         out=[out;...
%             all(wave_ids==0,2) & reg_sel,...
%             reg_sel & all(ismember(wave_ids,[1,5]),2) | all(ismember(wave_ids,[2,6]),2) | all(ismember(wave_ids,[3,5]),2) | all(ismember(wave_ids,[4,6]),2),...
%             reg_sel & (any(wave_ids==1,2) & any(wave_ids==3,2) & all(ismember(wave_ids,[1 3]),2)) | (any(wave_ids==2,2) & any(wave_ids==4,2) & all(ismember(wave_ids,[2 4]),2)),...
%             within_sel,...
%             cross_sel,...
%             (rsidx+2)*ones(size(within_sel)),...
%             sess*ones(size(within_sel)),...
%             rings{sess,rsidx},...
%             nan(numel(within_sel),3-rsidx)];
    end
end

