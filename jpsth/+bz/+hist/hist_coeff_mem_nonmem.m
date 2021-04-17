function hist_coeff_mem_nonmem(sess,type,opt)
% fit history effect/short term plasticity/time constant between pairs of
% neuron with a general linear model
%
% use \jpsth\+bz\+hist\plot_file_data.m for visualization
arguments
    sess (1,1) double {mustBeInteger,mustBePositive,mustBeNonempty}
    type (1,:) char {mustBeMember(type,{'congru','incongru','non-mem'})}
    opt.prefix (1,:) char = '0331'
    opt.tsbin_size (1,1) double = 600
    opt.postspike (1,1) logical = true
    opt.fc_effi (1,1) logical = false
    opt.fc_prob (1,1) logical = false
end
sig=bz.load_sig_pair();
switch type
    case 'congru'
        typesel=sig.sess==sess & (all(ismember(sig.mem_type,1:2),2) | all(ismember(sig.mem_type,3:4),2));
    case 'incongru'
        typesel=sig.sess==sess ...
            & ((ismember(sig.mem_type(:,1),1:2) & ismember(sig.mem_type(:,2),3:4)) ...
            |(ismember(sig.mem_type(:,1),3:4) & ismember(sig.mem_type(:,2),1:2)));
    case 'non-mem'
        typesel=sig.sess==sess & all(sig.mem_type==0,2);
end
dim=nnz(typesel);
if dim==0
    return
end
idces=find(typesel);
sess_suids=nan(dim,2);
postspk=nan(dim,11); %sess,pre-cid,post-cid,incept+coefficients
fc_eff=nan(dim,11);
fc_prob=nan(dim,11);
maxiter=false(dim,3);
sidx=1;
for i=reshape(idces,1,[])
    fprintf('%d of %d\n',sidx,dim);
    sess_suids(sidx,:)=sig.suid(i,:);
    [postspk(sidx,:),fc_eff(sidx,:),fc_prob(sidx,:),maxiter(sidx,:)]=...
        bz.hist.history_coeff(sig.sess(i),sig.suid(i,:),...
        'tsbin_size',opt.tsbin_size,...
        'postspike',opt.postspike,...
        'fc_effi',opt.fc_effi,...
        'fc_prob',opt.fc_prob);
    %maxiter->[SPK,FC_EFF,FC_PROB] 
    sidx=sidx+1;
end
%TODO return if whatever-empty
save(sprintf('%s_stp_%s_%d_%d.mat',opt.prefix,type,sess,opt.tsbin_size),'fc_eff','fc_prob','postspk','sess_suids','maxiter','sess','type');
end