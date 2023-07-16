global gather_config
gather_config=struct();
gather_config.fc_win=10;
gather_config.adjust_white_matter=true;
gather_config.corr_type='Pearson';
gather_config.fnsuffix='_10ms_adj_pearson';
gather_config.FC_thresh=100;
gather_config.sel_su_per_region=20;
gather_config.odpath=fullfile(getenv('userprofile'),'OneDrive','Neupix');

if strcmp(gather_config.corr_type,'Pearson')
    gather_config.corr_ln_log='PearsonLinearLog';
    gather_config.corr_log_log='PearsonLogLog';
else
    gather_config.corr_ln_log='Spearman';
    gather_config.corr_log_log='Spearman';
end