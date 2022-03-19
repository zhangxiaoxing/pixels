%featii,epochii,regii,double(glmxmeta(regii,:)),r,p
load('one_reg_corr_list.mat','one_reg_corr_list')
idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
uregs=unique(one_reg_corr_list(:,6));
for feat=1:2:5
    disp(feat)
    for uii=reshape(uregs,1,[])
%         disp(uii)
        for io=1:2 %afferent, efferent
%             disp(io)
            sel=one_reg_corr_list(:,6)==uii & one_reg_corr_list(:,4)==io & one_reg_corr_list(:,1)==feat;
            if nnz(sel)==0
                continue
            end
            r=(one_reg_corr_list(sel,7));
            if max((r))<0.5
                continue
            end
            fh=figure('Color','w');
            plot(one_reg_corr_list(sel,2),one_reg_corr_list(sel,7))
            title(idmap.ccfid2full(one_reg_corr_list(find(sel,1,'first'),6)));
            ylim([-0.75,0.75]);
            set(gca(),'XTick',1:6,'XTickLabel',{'Sample','EarlyD','LateD','Test','PostRw','PreSamp'})
%             disp(feat)
%             disp(io)
            waitfor(fh)
        end
    end
end