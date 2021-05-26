function [mem_type,per_bin]=get_mem_type(wrs_p,selec,opt)
% 0=NM,1=S1 sust, 2=S1 trans, 3=S2 sust, 4=S2 trans,-1=switched
arguments
    wrs_p (14,:) double % assuming 1 sec bin
    selec (14,:) double % assuming 1 sec bin
    opt.alpha (1,1) double = 0.05
end
mem_type=nan(1,size(selec,2));
per_bin=nan(6,size(selec,2));
for i=1:size(selec,2)
    sel_bin=find(wrs_p(5:10,i)<opt.alpha);
    if isempty(sel_bin)
        mem_type(i)=0;
        per_bin(:,i)=0;
    else
        ssign=sign(selec(sel_bin+4,i));
        if numel(unique(ssign))==1 % non-switched 
            if numel(ssign)==6  % sustained
                if ssign(1)==1
                    mem_type(i)=1;
                    per_bin(:,i)=1;
                else
%                     if ssign(1)~=-1,keyboard;end % one time check
                    mem_type(i)=3;
                    per_bin(:,i)=2;
                end
            else % transient
                per_bin(:,i)=0;
                if ssign(1)==1
                    mem_type(i)=2;
                    per_bin(sel_bin,i)=1;
                else
%                     if ssign(1)~=-1,keyboard;end % one time check
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
% h5create(fullfile(opt.homedir,'transient_6.hdf5'),'/mem_type',size(mem_type),'Datatype','int32')
% h5write(fullfile(opt.homedir,'transient_6.hdf5'), '/mem_type', mem_type);


end