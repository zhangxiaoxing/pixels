% Generate a list that contains cyclic FCs. Regardless of
% how frequently the associated structure appears in actual recordings. 

% in use as of 2023.04.13

% TODO: streamline multiple dataset output.
% TODO: Option to load rings from file.

% Pending delete as of 2023.10.29

function rings=ring_list_bz(opt)
arguments
    opt.shufid {mustBeScalarOrEmpty} =[]
    opt.wave (1,1) logical = false
end
%bzthres=250;  %TODO filter by spike number % not really necessary when using full-length data
if isempty(opt.shufid)
    sig=bz.load_sig_sums_conn_file('pair',false);
    fname='rings_bz.mat';
else
    load('bz_ring_shufs.mat','shufs')
    sig=shufs{opt.shufid};
    fname=sprintf('rings_bz_shuf_%d.mat',opt.shufid);
end
% if
% load(fullfile('bzdata','rings_bz.mat'))
% else
rings=cell(max(sig.sess),3);
for sess=1:max(sig.sess)
    disp(sess);
    for ring_size=3:5
        sess_sel = sig.sess==sess;
        if nnz(sess_sel)<3, continue;end
        sess_rings=bz.rings.find_rings_bz(sig.suid(sess_sel,:),ring_size);
        rings{sess,ring_size-2}=unique(bz.rings.flexsort(sess_rings),'rows');
    end
end
% end

if opt.wave
    % tag wave for rings
    su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
    wrs_mux_meta=ephys.get_wrs_mux_meta();
    rings_wave=cell2struct({[],[],[],cell(0)},{'sess','dur','wave','cids'},2);

    % Separated for historical-compatibility reasons. Consider merge when
    % ready
    rings_incon=cell2struct({[],[],[],cell(0)},{'sess','dur','wave','cids'},2);
    rings_nonmem=cell2struct({[],[],[],cell(0)},{'sess','dur','wave','cids'},2);

    for sidx=1:size(rings,1)
        sesscid=su_meta.allcid(su_meta.sess==sidx);
        sesswid=wrs_mux_meta.wave_id(su_meta.sess==sidx);
        for rSizeIdx=1:3
            for ridx=1:size(rings{sidx,rSizeIdx},1)
                rcids=rings{sidx,rSizeIdx}(ridx,:);
                [~,supos]=ismember(rcids,sesscid);
                rwids=sesswid(supos);
%                 keyboard();end;end;end
                [rwid,seltype]=bz.rings.ring_wave_type(rwids);
                if strcmp(rwid,'others')
                    continue
                end

                if strcmp(rwid,'nonmem')
                    for dd=[3 6]
                        rings_nonmem.sess=[rings_nonmem.sess;sidx];
                        rings_nonmem.cids=[rings_nonmem.cids;{rcids}];
                        rings_nonmem.wave=[rings_nonmem.wave;'NA'];
                        rings_nonmem.dur=[rings_nonmem.dur;dd];
                    end
                    continue
                end

                if strcmp(rwid,'incongru')
                    for dd=[3 6]
                        rings_incon.sess=[rings_incon.sess;sidx];
                        rings_incon.cids=[rings_incon.cids;{rcids}];
                        rings_incon.wave=[rings_incon.wave;'NA'];
                        rings_incon.dur=[rings_incon.dur;dd];
                    end
                    continue
                end


                if ~all(ismember(rwids,1:8))
                    % should not happen
                    keyboard()
                end
                if all(rwids==5,'all')
                    wid="olf_s1";
                    dur=[3 6];
                elseif all(rwids==6,'all')
                    wid="olf_s2";
                    dur=[3 6];
                elseif all(rwids==7,'all')
                    wid="dur_d3";
                    dur=3;
                elseif all(rwids==8,'all')
                    wid="dur_d6";
                    dur=6;
                elseif all(rwids==5 | rwids==1,'all') ...
                        || all(rwids==7 | rwids==1,'all')
                    wid="s1d3";
                    dur=3;
                elseif all(rwids==5 | rwids==2,'all')...
                        || all(rwids==8 | rwids==2,'all')
                    wid="s1d6";
                    dur=6;
                elseif all(rwids==6 | rwids==3,'all')...
                        || all(rwids==7 | rwids==3,'all')
                    wid="s2d3";
                    dur=3;
                elseif all(rwids==6 | rwids==4,'all')...
                        || all(rwids==8 | rwids==4,'all')
                    wid="s2d6";
                    dur=6;
                else
                    % unclassified group. should not happen
                    keyboard()
                end

                for dd=dur
                    rings_wave.sess=[rings_wave.sess;sidx];
                    rings_wave.cids=[rings_wave.cids;{rcids}];
                    rings_wave.wave=[rings_wave.wave;wid];
                    rings_wave.dur=[rings_wave.dur;dd];
                end
            end
        end
    end
    blame=vcs.blame();
    save(fullfile('bzdata','rings_bz_wave.mat'),'rings_wave','rings_incon','rings_nonmem','blame');
end

blame=vcs.blame();
save(fullfile('bzdata',fname),'rings','blame');
end