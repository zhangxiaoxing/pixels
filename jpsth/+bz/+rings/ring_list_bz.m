% Generate a list that contains cyclic FCs. Regardless of
% how frequently the associated structure appears in actual recordings. 
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

if opt.wave
    % TODO: tag wave
    su_meta=ephys.util.load_meta('skip_stats',true,'adjust_white_matter',true);
    wrs_mux_meta=ephys.get_wrs_mux_meta();
    rings_wave=struct();
    [rings_wave.sess,rings_wave.dur,rings_wave.wave]=deal([]);
    rings_wave.cids=cell(0);
    for sidx=1:size(rings,1)
        sesscid=su_meta.allcid(su_meta.sess==sidx);
        sesswid=wrs_mux_meta.wave_id(su_meta.sess==sidx);
        for rSizeIdx=1:3
            for ridx=1:size(rings{sidx,rSizeIdx},1)
                rcids=rings{sidx,rSizeIdx}(ridx,:);
                [~,supos]=ismember(rcids,sesscid);
                rwids=sesswid(supos);
                if ~all(ismember(rwids,1:8))
                    continue;
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
                    continue;% incongruent
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
    save(fullfile('bzdata','rings_bz_wave.mat'),'rings_wave','blame');

end

blame=vcs.blame();
save(fullfile('bzdata',fname),'rings','blame');
end