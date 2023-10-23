% The preferred shuffle method as of 2022.12.05

% The distinction between the original shuffle approach and this one is
% whether the FC rate difference of memory and non-memory neurons
% should be taken into account (original) or not (This). 

function shufs=shuffle_conn_bz_alt(opt)
arguments
    opt.poolsize (1,1) double {mustBeInteger,mustBePositive} = 2
    opt.rpt (1,1) double {mustBeInteger,mustBePositive} = 2
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
    opt.fn (1,:) char = "bz_ring_shufs.mat"
end
warning("Using "+opt.poolsize+" workers");

poolo=parpool(opt.poolsize);
shufs=cell(opt.rpt,1);
parfor rpt=1:opt.rpt
    shufs{rpt}=shuffle_conn_onerpt(opt);
end
blame=vcs.blame();
blame.author_tag=['shuffled （structural) connection data for further motif analysis,' ...
    ' Not distinguishing F.C. rate between memory and non-memory neurons' ...
    ' and brain areas '];
save(fullfile("binary",opt.fn),'shufs','blame','opt')
delete(poolo);
end

function shuf=shuffle_conn_onerpt(opt)

shuf=struct();
shuf_list=[];
% cnt=[];
global_init;
[sig,pair]=bz.load_sig_sums_conn_file('pair',true,'criteria',opt.criteria);

sig_reg_sel=all(ismember(sig.reg(:,1,:),[343,567]),3);
pair_reg_sel=all(ismember(pair.reg(:,1,:),[343,567]),3);

usess=unique(sig.sess);
%per session
for si=1:numel(usess)
    curr_sess=usess(si);
    sig_sess_sel=sig.sess==curr_sess;
    pair_sess_sel=pair.sess==curr_sess;
    curr_cnt=nnz(sig_sess_sel & sig_reg_sel);
    if curr_cnt==0, continue;end
    pair_curr_sel=pair_reg_sel & pair_sess_sel;
    shuf_list=[shuf_list;randsample(find(pair_curr_sel),curr_cnt)];
end

for fi=["suid","sess"]
    shuf.(fi)=pair.(fi)(shuf_list,:);
end
shuf.reg=pair.reg(shuf_list,:,:);
% random switch direction
toswitch=randsample(numel(shuf_list),floor(numel(shuf_list)./2));
shuf.suid(toswitch,:)=shuf.suid(toswitch,[2,1]);
shuf.reg(toswitch,:,:)=shuf.reg(toswitch,:,[2,1]);
end