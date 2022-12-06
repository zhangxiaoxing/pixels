% The preferred shuffle method as of 2022.12.05

% The distinction between the original shuffle approach and this one is
% whether or not the FC rate difference of memory and non-memory neurons
% should be taken into account. 

function shufs=shuffle_conn_bz_alt(opt)
arguments
    opt.poolsize (1,1) double {mustBeInteger,mustBePositive} = 2
    opt.rpt (1,1) double {mustBeInteger,mustBePositive} = 2
end
poolo=parpool(opt.poolsize);
shufs=cell(opt.rpt,1);
parfor rpt=1:opt.rpt
    shufs{rpt}=shuffle_conn_onerpt(rpt);
end
blame=vcs.blame();
blame.author_tag=['shuffled （structural) connection data for further loop analysis,' ...
    ' Not distinguishing F.C. rate between memory and non-memory neurons' ...
    ' and brain areas '];
save('bz_ring_shufs.mat','shufs','blame')
delete(poolo);
end

function shuf=shuffle_conn_onerpt(rptid,opt)
arguments
    rptid (1,1) double {mustBeInteger,mustBePositive}
    opt.reg_branch (1,1) double {mustBeMember(opt.reg_branch,1:5)}=5;
%     opt.mem_type (1,:) char
%     {mustBeMember(opt.mem_type,{'congru','nonmem'})}='congru'  Will not
%     respect memory type in alternative method
end
shuf=struct();
shuf_list=[];
% cnt=[];
global_init;
[sig,pair]=bz.load_sig_sums_conn_file('pair',true);

sig_reg_sel=all(ismember(sig.reg(:,1,:),[343,567]),3);
pair_reg_sel=all(ismember(pair.reg(:,1,:),[343,567]),3);

usess=unique(sig.sess);
%per session
for si=1:numel(usess)
%     fprintf('%d,%d\n',rptid,si);
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