%TODO merge with reg_conn_bz script
function [sig_]=load_sig_pair(opt)
arguments
    opt.type (1,:) char {mustBeMember(opt.type,{'MYWT'})}='MYWT'
    opt.prefix (1,:) char = '0516'
    opt.bin_range (1,2) double = [-2 7]
end
persistent sig  type_
if isempty(sig) || ~strcmp(opt.type,type_)
    fl=dir(fullfile('mydata',sprintf('%s_conn_w_reg_my_*_%d_%d.mat',opt.prefix,opt.bin_range(1),opt.bin_range(2))));
    if size(fl,1)<9, warning('Files not found');return; end
    sig=struct(); % for significant connect
    for s=["s1","s2"]
        sig.(s).suid=cell(0); % cluster id assigned by kilosort, 2nd+ probe prefixed by probe#
        sig.(s).reg=cell(0); % brain region tree
        %     sig.wrsp=cell(0); % Wilcoxon rank summation p value
        %     sig.selec=cell(0); % selectivity index
        sig.(s).sess=cell(0);
        sig.(s).mem_type=cell(0); % 0=NM,1=S1 sust, 2=S1 trans, 3=S2 sust, 4=S2 trans,-1=switched, see epys.get_mem_type
        
        fields={'suid','reg','mem_type','per_bin'};
        for fidx=1:size(fl,1)
            disp(fidx);
            fstr=load(fullfile(fl(fidx).folder,fl(fidx).name));
            if isfield(fstr.sig_meta,'s1') && isfield(fstr.sig_meta,'s2')
                for fi=fields
                    sig.(s).(fi{1}){fidx}=fstr.sig_meta.(s).(fi{1});
                end
                sig.(s).sess{fidx}=repmat(ephys.path2sessid(fstr.pc_stem,'type','neupix'),size(fstr.sig_meta.(s).suid,1),1);
            else
                for fi=fields
                    sig.(s).(fi{1}){fidx}=int32([]);
                end
                sig.(s).sess{fidx}=int32([]);
            end
        end
        for fi=[fields,{'sess'}]
            sig.(s).(fi{1})=cell2mat(reshape(sig.(s).(fi{1}),[],1));
        end
    end
end
sig_=sig;
type_=opt.type;
end