function out_=get_dur_meta(opt)
arguments
    opt.permutation (1,1) logical = true
    opt.merge_bin (1,1) logical = false
    opt.load_file (1,1) logical = true
    opt.save_file (1,1) logical = false
    opt.perm_repeat (1,1) double {mustBePositive,mustBeInteger} = 1000
end

persistent out opt_
if isempty(out) || ~isequaln(opt,opt_)
    if opt.load_file
        load('perm_dur.mat','dur_meta')
        out=dur_meta;
    else
        [~,~,sessmap]=ephys.sessid2path(0);
        homedir=ephys.util.getHomedir();
        sesskeys=cell2mat(sessmap.keys());
        [out.wrs_p_s1,out.wrs_p_s2,out.fdr_s1,out.fdr_s2,out.selec_s1,out.selec_s2,...
            out.mem_type_s1,out.per_bin_s1,out.mem_type_s2,out.per_bin_s2,out.wave_id]=deal([]);
        for sessid=sesskeys
            disp(sessid);
            %         if false
            %             fpath=fullfile(homedir,sessmap(sessid),"FR_All_ 250.hdf5");
            %             fr25=h5read(fpath,'/FR_All');
            %         end
            fpath=fullfile(homedir,sessmap(sessid),"FR_All_1000.hdf5");
            fr=h5read(fpath,'/FR_All');
            trials=h5read(fpath,'/Trials');
            suid=h5read(fpath,'/SU_id');
            c6sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,8)==6;
            c3sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,8)==3;

            cS1sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,5)==4;
            cS2sel=trials(:,9)~=0 & trials(:,10)~=0 & trials(:,5)==8;
            if opt.merge_bin, bins=1; else, bins=4; end
            [wrs_p_s1,wrs_p_s2,fdr_s1,fdr_s2,selec_s1,selec_s2]=deal(nan(size(fr,2),bins));
            for suidx=1:size(fr,2)
                if opt.merge_bin
                    if opt.permutation
                        wrs_p_s1(suidx,:)=permutation_test_1d(mean(fr(cS1sel & c3sel,suidx,5:7),3),mean(fr(cS1sel & c6sel,suidx,5:7),3),opt.perm_repeat);
                        wrs_p_s2(suidx,:)=permutation_test_1d(mean(fr(cS2sel & c3sel,suidx,5:7),3),mean(fr(cS2sel & c6sel,suidx,5:7),3),opt.perm_repeat);
                    else
                        wrs_p_s1(suidx,:)=ranksum(mean(fr(cS1sel & c3sel,suidx,5:7),3),mean(fr(cS1sel & c6sel,suidx,5:7),3));
                        wrs_p_s2(suidx,:)=ranksum(mean(fr(cS2sel & c3sel,suidx,5:7),3),mean(fr(cS2sel & c6sel,suidx,5:7),3));
                    end
                    selec_s1(suidx,:)=sel_idx(fr(cS1sel & c3sel,suidx,5:7),fr(cS1sel & c6sel,suidx,5:7));
                    selec_s2(suidx,:)=sel_idx(fr(cS2sel & c3sel,suidx,5:7),fr(cS2sel & c6sel,suidx,5:7));
                else
                    if opt.permutation
                        wrs_p_s1(suidx,:)=arrayfun(@(x) permutation_test_1d(fr(cS1sel & c3sel,suidx,x),fr(cS1sel & c6sel,suidx,x),opt.perm_repeat),4:7);
                        wrs_p_s2(suidx,:)=arrayfun(@(x) permutation_test_1d(fr(cS2sel & c3sel,suidx,x),fr(cS2sel & c6sel,suidx,x),opt.perm_repeat),4:7);
                    else
                        wrs_p_s1(suidx,:)=arrayfun(@(x) ranksum(fr(cS1sel & c3sel,suidx,x),fr(cS1sel & c6sel,suidx,x)),4:7);
                        wrs_p_s2(suidx,:)=arrayfun(@(x) ranksum(fr(cS2sel & c3sel,suidx,x),fr(cS2sel & c6sel,suidx,x)),4:7);
                    end
                    fdr_s1(suidx,:)=mafdr(wrs_p_s1(suidx,:),'BHFDR',true);
                    fdr_s2(suidx,:)=mafdr(wrs_p_s2(suidx,:),'BHFDR',true);
                    selec_s1(suidx,:)=arrayfun(@(x) sel_idx(fr(cS1sel & c3sel,suidx,x),fr(cS1sel & c6sel,suidx,x)),4:7);
                    selec_s2(suidx,:)=arrayfun(@(x) sel_idx(fr(cS2sel & c3sel,suidx,x),fr(cS2sel & c6sel,suidx,x)),4:7);
                end
            end
            if opt.merge_bin
                [mem_type_s1,per_bin_s1]=get_mem_type(wrs_p_s1,selec_s1);
                [mem_type_s2,per_bin_s2]=get_mem_type(wrs_p_s2,selec_s2);
            else
                [mem_type_s1,per_bin_s1]=get_mem_type(fdr_s1,selec_s1);
                [mem_type_s2,per_bin_s2]=get_mem_type(fdr_s2,selec_s2);
            end
            wave_id=get_wave_id(mem_type_s1,mem_type_s2);
            out.wrs_p_s1=[out.wrs_p_s1;wrs_p_s1];
            out.wrs_p_s2=[out.wrs_p_s2;wrs_p_s2];
            out.fdr_s1  =[out.fdr_s1  ;fdr_s1  ];
            out.fdr_s2  =[out.fdr_s2  ;fdr_s2  ];
            out.selec_s1=[out.selec_s1;selec_s1];
            out.selec_s2=[out.selec_s2;selec_s2];
            out.mem_type_s1=[out.mem_type_s1;mem_type_s1];
            out.mem_type_s2=[out.mem_type_s2;mem_type_s2];
            out.per_bin_s1 =[out.per_bin_s1 ;per_bin_s1 ];
            out.per_bin_s2 =[out.per_bin_s2 ;per_bin_s2 ];
            out.wave_id=[out.wave_id;wave_id];
            if opt.save_file
                dur_meta=out;
                save('perm_dur.mat','dur_meta');
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
    per_bin=nan(size(selec,1),3);
    for i=1:size(selec,1)
        sel_bin=find(stat(i,2:4)<opt.alpha);
        if isempty(sel_bin)
            mem_type(i)=0;
            per_bin(i,:)=0;
        else
            ssign=sign(selec(i,sel_bin+1));
            if all(ssign>=0) || all(ssign<=0) % non-switched
                per_bin(i,:)=0;
                if numel(sel_bin)==3  %sust, TODO replace magic number
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


function wave_id=get_wave_id(mem_type_s1,mem_type_s2)
wave_id=zeros(size(mem_type_s1))-1;
wave_id(mem_type_s1==0 & mem_type_s2==0)=0;
wave_id(ismember(mem_type_s1,1:2) & mem_type_s2==0)=1; %s1d6 s2/
wave_id(ismember(mem_type_s1,3:4) & mem_type_s2==0)=2; %s1d3 s2/
wave_id(ismember(mem_type_s2,1:2) & mem_type_s1==0)=3; %s2d6 s1/
wave_id(ismember(mem_type_s2,3:4) & mem_type_s1==0)=4; %s2d3 s1/
wave_id(ismember(mem_type_s1,1:2) & ismember(mem_type_s2,1:2))=5; %s1s2 d6
wave_id(ismember(mem_type_s1,3:4) & ismember(mem_type_s2,3:4))=6; %s1s2 d3
end