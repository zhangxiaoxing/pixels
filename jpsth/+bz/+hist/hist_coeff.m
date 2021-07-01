function hist_coeff(sess,opt)
% fit history effect/short term plasticity/time constant between pairs of
% neuron with a general linear model
%
% use \jpsth\+bz\+hist\plot_file_data.m for visualization

arguments
    sess (1,1) double {mustBeInteger,mustBePositive,mustBeNonempty}
    opt.prefix (1,:) char = 'BZWT'
    opt.tsbin_size (1,1) double = 600
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO','MY'})}='neupix'
    opt.laser (1,:) char {mustBeMember(opt.laser,{'on','off','any'})} = 'any'
    opt.epoch (1,:) char {mustBeMember(opt.epoch,{'delay','ITI','any'})} = 'any'
    opt.criteria (1,:) char {mustBeMember(opt.criteria,{'Learning','WT','any'})} = 'WT'
    opt.correct_error (1,:) char {mustBeMember(opt.correct_error,{'correct','error','any'})} = 'any'
end

sig=bz.load_sig_pair('type',opt.type,'prefix',opt.prefix,'criteria',opt.criteria);
typesel=sig.sess==sess;
dim=nnz(typesel);
if dim==0
    return
end
idces=find(typesel);
sess_suids=nan(dim,2);
postspk=nan(dim,11); %incept+coefficients
skip=false(dim,1);
[pp,rsq]=deal(nan(dim,1));
sidx=1;
for i=reshape(idces,1,[])
    fprintf('%d of %d\n',sidx,dim);
    sess_suids(sidx,:)=sig.suid(i,:);
    [postspk(sidx,:),skip(sidx),pp(sidx),rsq(sidx)]=...
        bz.hist.history_coeff_one(sig.sess(i),sig.suid(i,:),...
        'tsbin_size',opt.tsbin_size,...
        'type',opt.type,...
        'laser',opt.laser,...
        'epoch',opt.epoch,...
        'criteria',opt.criteria,...
        'correct_error',opt.correct_error);
    %maxiter->[SPK,FC_EFF,FC_PROB] 
    sidx=sidx+1;
end
%TODO return if whatever-empty
switch opt.laser
    case 'on',laser_suffix='_laserOn';
    case 'off',laser_suffix='_laserOff';
    case 'any',laser_suffix='';
end

switch opt.epoch
    case 'ITI',epoch_suffix='_ITI';
    case 'delay',epoch_suffix='_delay';
    case 'any',epoch_suffix='';
end

switch opt.criteria
    case 'Learning',epoch_suffix='_Learning';
end

blame=vcs.blame();
save(sprintf('%s_stp_%s_%d_%d%s%s.mat',...
opt.prefix,opt.correct_error,sess,opt.tsbin_size,laser_suffix,epoch_suffix),...
'postspk','pp','rsq','sess_suids','skip','sess','blame');
end
