function wave_stay_disappear(wrs_mux_meta)

global_init();

% wave

% selectivity? activity?
% odor-only? merge both?
osel=ismember(wrs_mux_meta.wave_id,5:6);

% wave-stay
stay_sel=osel & any(wrs_mux_meta.p_olf6<0.05,2);

% wave-disappear
disap_sel=osel & all(wrs_mux_meta.p_olf6>=0.05,2);

% F.C.
[sig,pair]=bz.load_sig_sums_conn_file('pair',true);

% stay s1:15 stays2:16 disapps1:25 disapps2:26 others:unchanged

SDid=wrs_mux_meta.wave_id;
SDid(stay_sel & SDid==5)=15;
SDid(stay_sel & SDid==6)=16;
SDid(disap_sel & SDid==5)=25;
SDid(disap_sel & SDid==6)=26;

sig=bz.join_fc_waveid(sig,SDid);
pair=bz.join_fc_waveid(pair,SDid);

sig_generic=bz.join_fc_waveid(sig,wrs_mux_meta.wave_id);
pair_generic=bz.join_fc_waveid(pair,wrs_mux_meta.wave_id);

% same region, cross region
sig_same_reg_sel=all(ismember(sig.reg(:,1,:),[343,567]),3) ...
    & all(sig.reg(:,5,:)>0,3) ...
    & sig.reg(:,5,1)==sig.reg(:,5,2);

sig_diff_reg_sel=all(ismember(sig.reg(:,1,:),[343,567]),3) ...
    & all(sig.reg(:,5,:)>0,3) ...
    & sig.reg(:,5,1)~=sig.reg(:,5,2);

pair_same_reg_sel=all(ismember(pair.reg(:,1,:),[343,567]),3) ...
    & all(pair.reg(:,5,:)>0,3) ...
    & pair.reg(:,5,1)==pair.reg(:,5,2);

pair_diff_reg_sel=all(ismember(pair.reg(:,1,:),[343,567]),3) ...
    & all(pair.reg(:,5,:)>0,3) ...
    & pair.reg(:,5,1)~=pair.reg(:,5,2);


%% congru to, congru from, incongru to, incongru from.
stay_samereg=[nnz(sig_same_reg_sel & sig.waveid(:,1)==15 & sig.waveid(:,2)==15),nnz(pair_same_reg_sel &pair.waveid(:,1)==15 & pair.waveid(:,2)==15);...
nnz(sig_same_reg_sel & sig.waveid(:,1)==16 & sig.waveid(:,2)==16),nnz(pair_same_reg_sel &pair.waveid(:,1)==16 & pair.waveid(:,2)==16)];

disapp_samereg=[nnz(sig_same_reg_sel & sig.waveid(:,1)==25 & sig.waveid(:,2)==25),nnz(pair_same_reg_sel &pair.waveid(:,1)==25 & pair.waveid(:,2)==25);...
nnz(sig_same_reg_sel & sig.waveid(:,1)==26 & sig.waveid(:,2)==26),nnz(pair_same_reg_sel &pair.waveid(:,1)==26 & pair.waveid(:,2)==26)];

[cong_same_hat,cong_same_ci]=binofit(sum(stay_samereg(:,1)),sum(stay_samereg(:,2)))
[disapp_same_hat,disapp_same_ci]=binofit(sum(disapp_samereg(:,1)),sum(disapp_samereg(:,2)))

[~,~,psame]=crosstab([zeros(sum(stay_samereg(:,2)),1);ones(sum(disapp_samereg(:,2)),1)],...
    [1:sum(stay_samereg(:,2))>sum(stay_samereg(:,1)),1:sum(disapp_samereg(:,2))>sum(disapp_samereg(:,1))])
%========================
stay_diffreg=[nnz(sig_diff_reg_sel & sig.waveid(:,1)==15 & sig.waveid(:,2)==15),nnz(pair_diff_reg_sel &pair.waveid(:,1)==15 & pair.waveid(:,2)==15);...
nnz(sig_diff_reg_sel & sig.waveid(:,1)==16 & sig.waveid(:,2)==16),nnz(pair_diff_reg_sel &pair.waveid(:,1)==16 & pair.waveid(:,2)==16)];

disapp_diffreg=[nnz(sig_diff_reg_sel & sig.waveid(:,1)==25 & sig.waveid(:,2)==25),nnz(pair_diff_reg_sel &pair.waveid(:,1)==25 & pair.waveid(:,2)==25);...
nnz(sig_diff_reg_sel & sig.waveid(:,1)==26 & sig.waveid(:,2)==26),nnz(pair_diff_reg_sel &pair.waveid(:,1)==26 & pair.waveid(:,2)==26)];

[cong_diff_hat,cong_diff_ci]=binofit(sum(stay_diffreg(:,1)),sum(stay_diffreg(:,2)))
[disapp_diff_hat,disapp_diff_ci]=binofit(sum(disapp_diffreg(:,1)),sum(disapp_diffreg(:,2)))

[~,~,pdiff]=crosstab([zeros(sum(stay_diffreg(:,2)),1);ones(sum(disapp_diffreg(:,2)),1)],...
    [1:sum(stay_diffreg(:,2))>sum(stay_diffreg(:,1)),1:sum(disapp_diffreg(:,2))>sum(disapp_diffreg(:,1))])


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% merge all odor neuron
% from; to;

stay_samereg=[nnz(sig_same_reg_sel & sig_generic.waveid(:,1)==5 & sig.waveid(:,2)==15),nnz(pair_same_reg_sel &pair_generic.waveid(:,1)==5 & pair.waveid(:,2)==15);...
nnz(sig_same_reg_sel & sig_generic.waveid(:,1)==6 & sig.waveid(:,2)==16),nnz(pair_same_reg_sel &pair_generic.waveid(:,1)==6 & pair.waveid(:,2)==16)];

disapp_samereg=[nnz(sig_same_reg_sel & sig_generic.waveid(:,1)==5 & sig.waveid(:,2)==25),nnz(pair_same_reg_sel &pair_generic.waveid(:,1)==5 & pair.waveid(:,2)==25);...
nnz(sig_same_reg_sel & sig_generic.waveid(:,1)==6 & sig.waveid(:,2)==26),nnz(pair_same_reg_sel &pair_generic.waveid(:,1)==6 & pair.waveid(:,2)==26)];

[cong_same_hat,cong_same_ci]=binofit(sum(stay_samereg(:,1)),sum(stay_samereg(:,2)))
[disapp_same_hat,disapp_same_ci]=binofit(sum(disapp_samereg(:,1)),sum(disapp_samereg(:,2)))

[~,~,psame]=crosstab([zeros(sum(stay_samereg(:,2)),1);ones(sum(disapp_samereg(:,2)),1)],...
    [1:sum(stay_samereg(:,2))>sum(stay_samereg(:,1)),1:sum(disapp_samereg(:,2))>sum(disapp_samereg(:,1))])
%========================
stay_diffreg=[nnz(sig_diff_reg_sel & sig_generic.waveid(:,1)==5 & sig.waveid(:,2)==15),nnz(pair_diff_reg_sel &pair_generic.waveid(:,1)==5 & pair.waveid(:,2)==15);...
nnz(sig_diff_reg_sel & sig_generic.waveid(:,1)==6 & sig.waveid(:,2)==16),nnz(pair_diff_reg_sel &pair_generic.waveid(:,1)==6 & pair.waveid(:,2)==16)];

disapp_diffreg=[nnz(sig_diff_reg_sel & sig_generic.waveid(:,1)==5 & sig.waveid(:,2)==25),nnz(pair_diff_reg_sel &pair_generic.waveid(:,1)==5 & pair.waveid(:,2)==25);...
nnz(sig_diff_reg_sel & sig_generic.waveid(:,1)==6 & sig.waveid(:,2)==26),nnz(pair_diff_reg_sel &pair_generic.waveid(:,1)==6 & pair.waveid(:,2)==26)];

[cong_diff_hat,cong_diff_ci]=binofit(sum(stay_diffreg(:,1)),sum(stay_diffreg(:,2)))
[disapp_diff_hat,disapp_diff_ci]=binofit(sum(disapp_diffreg(:,1)),sum(disapp_diffreg(:,2)))

[~,~,pdiff]=crosstab([zeros(sum(stay_diffreg(:,2)),1);ones(sum(disapp_diffreg(:,2)),1)],...
    [1:sum(stay_diffreg(:,2))>sum(stay_diffreg(:,1)),1:sum(disapp_diffreg(:,2))>sum(disapp_diffreg(:,1))])




