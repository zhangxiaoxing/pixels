load('opti.mat','both_map','exclu_map','sink_src_mat',"sink_ccfid","src_ccfid")
idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
grey_regs=ephys.getGreyRegs();
% prob=optimproblem('ObjectiveSense','min')
% b=optimvar('b',2,["bt","bl"])


sensv=subsref(cell2mat(both_map.values(grey_regs).'), struct(type='()',subs={{':',1}}));
durv=subsref(cell2mat(exclu_map.values(grey_regs).'), struct(type='()',subs={{':',1}}));
[~,kloc]=ismember(cell2mat(idmap.reg2ccfid.values(grey_regs)),sink_ccfid);
TTv=sink_src_mat(kloc,src_ccfid==idmap.reg2ccfid('TT'));
LATv=sink_src_mat(kloc,src_ccfid==idmap.reg2ccfid('LAT'));
[r,p]=corr(TTv,sensv,'type','Spearman')
