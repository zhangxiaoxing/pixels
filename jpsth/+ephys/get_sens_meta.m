function out_=get_sens_meta(opt)
arguments
    opt.permutation (1,1) logical = true
%     opt.merge_bin (1,1) logical = true
    opt.load_file (1,1) logical = true
    opt.save_file (1,1) logical = false
    opt.perm_repeat (1,1) double {mustBePositive,mustBeInteger} = 1000
    opt.uneven_duration (1,1) logical = true % use all 6-sec delay
end

persistent out opt_
if isempty(out) || ~isequaln(opt,opt_)
    if opt.uneven_duration
        bin3s=5:7;
        bin6s=5:10;
    else
        bin3s=5:7;
        bin6s=5:7;
    end


    if opt.load_file
        load('perm_sens.mat','sens_meta')
        out=sens_meta;
    else

        [~,~,sessmap]=ephys.sessid2path(0);
        homedir=ephys.util.getHomedir('type','raw');
        sesskeys=cell2mat(sessmap.keys());
        [out.wrs_p_d3,out.wrs_p_d6,out.fdr_d3,out.fdr_d6,out.selec_d3,out.selec_d6,...
            out.mem_type_d3,out.per_bin_d3,out.mem_type_d6,out.per_bin_d6,out.wave_id]=deal([]);
        for sessid=sesskeys
            disp(sessid);
            fpath=fullfile(homedir,sessmap(sessid),"FR_All_1000.hdf5");
            fr=h5read(fpath,'/FR_All');
            trials=h5read(fpath,'/Trials');
%             suid=h5read(fpath,'/SU_id');
            s2sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,5)==8;
            s1sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,5)==4;

            d3sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,8)==3;
            d6sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,8)==6;
%             if opt.merge_bin, bins=1; else, bins=4; end
            [wrs_p_d3,fdr_d3,selec_d3]=deal(nan(size(fr,2),numel(bin3s)));
            [wrs_p_d6,fdr_d6,selec_d6]=deal(nan(size(fr,2),numel(bin6s)));
            for suidx=1:size(fr,2)
%                 if opt.merge_bin
%                     if opt.permutation
%                         wrs_p_d3(suidx,:)=permutation_test_1d(mean(fr(d3sel & s1sel,suidx,bin3s),3),mean(fr(d3sel & s2sel,suidx,bin3s),3),opt.perm_repeat);
%                          %TODO: delay-dependent
%                         wrs_p_d6(suidx,:)=permutation_test_1d(mean(fr(d6sel & s1sel,suidx,bin6s),3),mean(fr(d6sel & s2sel,suidx,bin6s),3),opt.perm_repeat);
%                     else
%                         wrs_p_d3(suidx,:)=ranksum(mean(fr(d3sel & s1sel,suidx,5:7),3),mean(fr(d3sel & s2sel,suidx,5:7),3));
%                          %TODO: delay-dependent
%                         wrs_p_d6(suidx,:)=ranksum(mean(fr(d6sel & s1sel,suidx,5:7),3),mean(fr(d6sel & s2sel,suidx,5:7),3));
%                     end
%                     selec_d3(suidx,:)=sel_idx(fr(d3sel & s1sel,suidx,5:7),fr(d3sel & s2sel,suidx,5:7));
%                      %TODO: delay-dependent
%                     selec_d6(suidx,:)=sel_idx(fr(d6sel & s1sel,suidx,5:7),fr(d6sel & s2sel,suidx,5:7));
%                 else
                    if opt.permutation
                         %TODO: delay-dependent
                        wrs_p_d3(suidx,:)=arrayfun(@(x) permutation_test_1d(fr(d3sel & s1sel,suidx,x),fr(d3sel & s2sel,suidx,x),opt.perm_repeat),bin3s);
                        wrs_p_d6(suidx,:)=arrayfun(@(x) permutation_test_1d(fr(d6sel & s1sel,suidx,x),fr(d6sel & s2sel,suidx,x),opt.perm_repeat),bin6s);
                    else
                         %TODO: delay-dependent
                        wrs_p_d3(suidx,:)=arrayfun(@(x) ranksum(fr(d3sel & s1sel,suidx,x),fr(d3sel & s2sel,suidx,x)),bin3s);
                        wrs_p_d6(suidx,:)=arrayfun(@(x) ranksum(fr(d6sel & s1sel,suidx,x),fr(d6sel & s2sel,suidx,x)),bin6s);
                    end
                    fdr_d3(suidx,:)=mafdr(wrs_p_d3(suidx,:),'BHFDR',true);
                    fdr_d6(suidx,:)=mafdr(wrs_p_d6(suidx,:),'BHFDR',true);

                    selec_d3(suidx,:)=sel_idx(fr(d3sel & s1sel,suidx,5:7),fr(d3sel & s2sel,suidx,bin3s));
                    selec_d6(suidx,:)=sel_idx(fr(d6sel & s1sel,suidx,5:7),fr(d6sel & s2sel,suidx,bin6s));

                     %TODO: delay-dependent
%                     selec_s1(suidx,:)=arrayfun(@(x) sel_idx(fr(d & c3sel,suidx,x),fr(cS1sel & c6sel,suidx,x)),4:7);
%                     selec_s2(suidx,:)=arrayfun(@(x) sel_idx(fr(cS2sel & c3sel,suidx,x),fr(cS2sel & c6sel,suidx,x)),4:7);
%                 end
            end
%             if opt.merge_bin
%                 [mem_type_d3,per_bin_d3]=get_mem_type(wrs_p_d3,selec_d3);
%                 [mem_type_d6,per_bin_d6]=get_mem_type(wrs_p_d6,selec_d6);
%             else
                [mem_type_d3,per_bin_d3]=get_mem_type(fdr_d3,selec_d3);
                [mem_type_d6,per_bin_d6]=get_mem_type(fdr_d6,selec_d6);
%             end
            wave_id=get_wave_id(mem_type_d3,mem_type_d6);
            out.wrs_p_d3=[out.wrs_p_d3;wrs_p_d3];
            out.wrs_p_d6=[out.wrs_p_d6;wrs_p_d6];
            out.fdr_d3  =[out.fdr_d3  ;fdr_d3  ];
            out.fdr_d6  =[out.fdr_d6  ;fdr_d6  ];
            out.selec_d3=[out.selec_d3;selec_d3];
            out.selec_d6=[out.selec_d6;selec_d6];
            out.mem_type_d3=[out.mem_type_d3;mem_type_d3];
            out.mem_type_d6=[out.mem_type_d6;mem_type_d6];
            out.per_bin_d3 =[out.per_bin_d3 ;per_bin_d3 ];
            out.per_bin_d6 =[out.per_bin_d6 ;per_bin_d6 ];
            out.wave_id=[out.wave_id;wave_id];
            if opt.save_file
                sens_meta=out;
                save('perm_sens.mat','sens_meta');
            end
        end
    end
end
opt_=opt;
out_=out;
end

%%
function si=sel_idx(A,B)
if any([A(:);B(:)])
    si=(mean(A(:))-mean(B(:)))./(mean(A(:))+mean(B(:)));
else
    si=0;
end
end

function [mem_type,per_bin]=get_mem_type(stat,selec,opt)
arguments
    stat
    selec
    opt.alpha (1,1) double = 0.05
end
mem_type=nan(size(selec,1),1);
if size(stat,2)>1
    per_bin=nan(size(selec,1),size(stat,2));
    for i=1:size(selec,1)
         %TODO: delay-dependent
        sel_bin=find(stat(i,:)<opt.alpha);
        if isempty(sel_bin)
            mem_type(i)=0;
            per_bin(i,:)=0;
        else
            ssign=sign(selec(i,sel_bin));
            if all(ssign>=0) || all(ssign<=0) % non-switched
                per_bin(i,:)=0;
                if numel(sel_bin)==size(stat,2)  %sust, TODO replace magic number
                    m_types=[1,3];
                else % transient
                    m_types=[2,4];
                end
                if any(ssign>0) %favor d6, sust=1,transient=2
                    mem_type(i)=m_types(1);
                    per_bin(i,sel_bin)=1;
                else    %favor d3, sust=3, transient=4
                    mem_type(i)=m_types(2);
                    per_bin(i,sel_bin)=2;
                end
            else %switched
                mem_type(i)=-1;
                per_bin(i,:)=-1;
            end
        end
    end
else
    mem_type(stat>=opt.alpha)=0;
    mem_type(stat<opt.alpha & selec>0)=2;
    mem_type(stat<opt.alpha & selec<0)=4;
    per_bin=[];
end
end

function wave_id=get_wave_id(mem_type_d3,mem_type_d6)
wave_id=zeros(size(mem_type_d3))-1;
wave_id(mem_type_d3==0 & mem_type_d6==0)=0;
wave_id(ismember(mem_type_d3,1:2) & mem_type_d6==0)=1; %d3d6 d6/
wave_id(ismember(mem_type_d3,3:4) & mem_type_d6==0)=2; %d3d3 d6/
wave_id(ismember(mem_type_d6,1:2) & mem_type_d3==0)=3; %d6d6 d3/
wave_id(ismember(mem_type_d6,3:4) & mem_type_d3==0)=4; %d6d3 d3/
wave_id(ismember(mem_type_d3,1:2) & ismember(mem_type_d6,1:2))=5; %d3d6 d6
wave_id(ismember(mem_type_d3,3:4) & ismember(mem_type_d6,3:4))=6; %d3d6 d3
end