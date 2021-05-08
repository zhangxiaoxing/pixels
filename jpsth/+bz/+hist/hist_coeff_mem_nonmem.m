function hist_coeff_mem_nonmem(sess,mtype,opt)
% fit history effect/short term plasticity/time constant between pairs of
% neuron with a general linear model
%
% use \jpsth\+bz\+hist\plot_file_data.m for visualization

arguments
    sess (1,1) double {mustBeInteger,mustBePositive,mustBeNonempty}
    mtype (1,:) char {mustBeMember(mtype,{'congru','incongru','non-mem'})}
    opt.prefix (1,:) char = '0331'
    opt.tsbin_size (1,1) double = 600
    opt.type (1,:) char {mustBeMember(opt.type,{'neupix','AIOPTO'})}='neupix'
    opt.laser (1,:) char {mustBeMember(opt.laser,{'on','off','any'})} = 'any'
    opt.epoch (1,:) char {mustBeMember(opt.epoch,{'delay','ITI','any'})} = 'any'
end

sig=bz.load_sig_pair('type',opt.type,'prefix',opt.prefix);
switch mtype
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
postspk=nan(dim,11); %incept+coefficients
skip=false(dim);
sidx=1;
for i=reshape(idces,1,[])
    fprintf('%d of %d\n',sidx,dim);
    sess_suids(sidx,:)=sig.suid(i,:);
    [postspk(sidx,:),skip(sidx)]=...
        bz.hist.history_coeff(sig.sess(i),sig.suid(i,:),...
        'tsbin_size',opt.tsbin_size,...
        'type',opt.type,...
        'laser',opt.laser,...
        'epoch',opt.epoch);
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

save(sprintf('%s_stp_%s_%d_%d%s%s.mat',...
opt.prefix,mtype,sess,opt.tsbin_size,laser_suffix,epoch_suffix),...
'postspk','sess_suids','skip','sess','mtype');
end