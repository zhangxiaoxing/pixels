function [sig_same, sig_h2l, sig_l2h,pair_same, pair_h2l, pair_l2h]=get_reg_hier_relation(opt)
arguments
    opt.load_int_data (1,1) logical = true
    opt.save_int_data (1,1) logical = false
    opt.hiermap (1,:) char {mustBeMember(opt.hiermap,{'pvsst','OBM1','AON','LAT'})} = 'AON'
    opt.sig_reg
    opt.pair_reg
end
if opt.load_int_data
    hier_str=load('conn_prob_bars_hier.mat');
    if ~isempty(opt.sig_reg)
        assert(isequaln(opt.sig_reg,hier_str.sig_reg))
    end
    if ~isempty(opt.pair_reg)
        assert(isequaln(opt.pair_reg,hier_str.pair_reg))
    end

    sig_same=hier_str.sig_same;
    sig_h2l=hier_str.sig_h2l ;
    sig_l2h=hier_str.sig_l2h  ;
    pair_same=hier_str.pair_same;
    pair_h2l=hier_str.pair_h2l;
    pair_l2h=hier_str.pair_l2h ;

else
    [~,sig_same,sig_h2l,sig_l2h]=bz.util.diff_at_level(opt.sig_reg,'hierarchy',true,'hiermap',opt.hiermap);
    [~,pair_same,pair_h2l,pair_l2h]=bz.util.diff_at_level(opt.pair_reg,'hierarchy',true,'hiermap',opt.hiermap);
    sig_reg=opt.sig_reg;
    pair_reg=opt.pair_reg;
    %save intermediate data
    if opt.save_int_data
        save('conn_prob_bars_hier.mat',...
            'pair_h2l',...
            'pair_l2h',...
            'pair_same',...
            'sig_h2l',...
            'sig_l2h',...
            'sig_same','pair_reg','sig_reg')
    end
end


