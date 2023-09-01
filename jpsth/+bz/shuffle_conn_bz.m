function shuffle_conn_bz(opt)
arguments
    opt.poolsize (1,1) double {mustBeInteger,mustBePositive} = 2
    opt.rpt (1,1) double {mustBeInteger,mustBePositive} = 5
end
parpool(opt.poolsize)
shufs=cell(opt.rpt,1);
parfor rpt=1:opt.rpt
    shufs{rpt}=shuffle_conn_onerpt(rpt);
end
save('bz_ring_shufs.mat','shufs')
end

function shuf=shuffle_conn_onerpt(rptid,opt)
arguments
    rptid (1,1) double {mustBeInteger,mustBePositive}
    opt.reg_branch (1,1) double {mustBeMember(opt.reg_branch,1:5)}=5;
    opt.mem_type (1,:) char {mustBeMember(opt.mem_type,{'congru','nonmem'})}='congru'
end
shuf=struct();
shuf_list=[];
% cnt=[];
[sig,pair]=bz.load_sig_pair('pair',true);
sig_mem_sel1=all(ismember(sig.mem_type,1:2),2);
sig_mem_sel2=all(ismember(sig.mem_type,3:4),2);
pair_mem_sel1=all(ismember(pair.mem_type,1:2),2);
pair_mem_sel2=all(ismember(pair.mem_type,3:4),2);
% pair_dir_sel=pair.suid(:,1)<pair.suid(:,2);


usess=unique(sig.sess);
%congruent,
%per session
for si=1:numel(usess)
    fprintf('%d,%d\n',rptid,si);
    curr_sess=usess(si);
    sig_sess_sel=sig.sess==curr_sess;
    pair_sess_sel=pair.sess==curr_sess;
    ureg=unique(pair.reg(pair.sess==curr_sess,opt.reg_branch,:));
    ureg=ureg(ureg~=0);
    %shuffle between
    %per region pair?
    if numel(ureg)>=2
        reg_comb=nchoosek(ureg,2);
    else
        reg_comb=[];
    end
    %     sesscnt=0;
    seq=[1,2;2,1];
    for rci=1:size(reg_comb,1)
        for dir=1:2
            %fwd,rev
            sig_reg_sel=sig.reg(:,opt.reg_branch,1)==reg_comb(rci,seq(dir,1)) ...
                & sig.reg(:,opt.reg_branch,2)==reg_comb(rci,seq(dir,2));
            if strcmp(opt.mem_type,'congru')
                curr_cnt=nnz(sig_sess_sel & sig_reg_sel & (sig_mem_sel1 | sig_mem_sel2));
                if curr_cnt==0, continue;end
                pair_reg_sel=pair.reg(:,opt.reg_branch,1)==reg_comb(rci,seq(dir,1)) ...
                    & pair.reg(:,opt.reg_branch,2)==reg_comb(rci,seq(dir,2));
                pair_subs_sel=(pair_sess_sel & pair_reg_sel...
                    &(pair_mem_sel1|pair_mem_sel2));
                shuf_list=[shuf_list;randsample(find(pair_subs_sel),curr_cnt)];
                %                 sesscnt=sesscnt+curr_cnt;
            end
        end
    end
    %     cnt=[cnt;curr_sess,sesscnt];
    %shuffle within
    %per region
    
    for ri=1:numel(ureg)
        curr_reg=ureg(ri);
        sig_reg_sel=sig.reg(:,opt.reg_branch,1)==curr_reg ...
            & sig.reg(:,opt.reg_branch,2)==curr_reg;
        if strcmp(opt.mem_type,'congru')
            curr_cnt=nnz(sig_sess_sel & sig_reg_sel & (sig_mem_sel1 | sig_mem_sel2));
            if curr_cnt==0, continue;end
            %get pair substrate
            pair_reg_sel=pair.reg(:,opt.reg_branch,1)==curr_reg ...
                & pair.reg(:,opt.reg_branch,2)==curr_reg;
            pair_subs_sel=(pair_sess_sel & pair_reg_sel...
                &(pair_mem_sel1|pair_mem_sel2));
            shuf_list=[shuf_list;randsample(find(pair_subs_sel),curr_cnt)];
        end
    end
end
for fi=["suid","sess","mem_type"]
    shuf.(fi)=pair.(fi)(shuf_list,:);
end
shuf.reg=pair.reg(shuf_list,:,:);
end


% assignin('base','shuf_list',shuf_list)
% nnz(sig_mem_sel1 | sig_mem_sel2)
% nnz(sig.reg(:,5,1)>0 & sig.reg(:,5,2)>0 & (sig_mem_sel1 | sig_mem_sel2))
%
% shuf_mem_sel1=all(ismember(shuf.mem_type,1:2),2);
% shuf_mem_sel2=all(ismember(shuf.mem_type,3:4),2);
% nnz(shuf_mem_sel1 | shuf_mem_sel2)
% nnz(shuf.reg(:,5,1)==shuf.reg(:,5,2) & (shuf_mem_sel1 | shuf_mem_sel2))

%%non-mem