function [mem_type,per_bin]=get_mem_type(stat,selec,opt)
% 0=NM,1=S1 sust, 2=S1 trans, 3=S2 sust, 4=S2 trans,-1=switched
arguments
    stat (:,:) double % assuming 1 sec bin
    selec (14,:) double % assuming 1 sec bin
    opt.alpha (1,1) double = 0.05
    opt.delay (1,1) double {mustBeMember(opt.delay,[3,6])}
    opt.fdr (1,1) logical = false
end
assert(ismember(size(stat,1),[14,5,8]))

mem_type=nan(1,size(selec,2));
per_bin=nan(opt.delay,size(selec,2));
for i=1:size(selec,2)
    if ~opt.fdr
        sel_bin=find(stat(5:(opt.delay+4),i)<opt.alpha);
    else
        sel_bin=find(stat(2:(opt.delay+1),i)<opt.alpha);
    end
    if isempty(sel_bin)
        mem_type(i)=0;
        per_bin(:,i)=0;
    else
        ssign=sign(selec(sel_bin+4,i));
        if all(ssign>=0) || all(ssign<=0) % non-switched 
            if nnz(ssign)==opt.delay  % sustained
                if ssign(1)==1
                    mem_type(i)=1;
                    per_bin(:,i)=1;
                else
                    mem_type(i)=3;
                    per_bin(:,i)=2;
                end
            else % transient
                per_bin(:,i)=0;
                if any(ssign==1)
                    mem_type(i)=2;
                    per_bin(sel_bin,i)=1;
                else
                    mem_type(i)=4;
                    per_bin(sel_bin,i)=2;
                end
            end
        else %switched
            mem_type(i)=-1;
            per_bin(:,i)=-1;
        end
    end
end
mem_type=int32(mem_type);
per_bin=int32(per_bin);

% TODO: parameter controlled file op
% homedir=ephys.util.getHomedir();
% h5create(fullfile(homedir,'transient_6.hdf5'),'/mem_type',size(mem_type),'Datatype','int32')
% h5write(fullfile(homedir,'transient_6.hdf5'), '/mem_type', mem_type);


end