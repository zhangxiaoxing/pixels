load('selec_duo_conn_chain_duo_6s_6_7.mat','reg_chain_list','pref_chain_list','conn_chain_list');
load('reg_keep.mat','ckeep');
reg_set=find(ckeep);
curr_pref=2;

done=[];
ccb1=conn_chain_list{1}(all(ismember(reg_chain_list{1},reg_set),2) & ...
    pref_chain_list{1}(:,1)==curr_pref &...
    all(pref_chain_list{1}(:,7:8)==curr_pref,2),:);%connection_chain_bin1
chain_sums=nan(0,7);
for s1=1:length(ccb1)
    fprintf('%d of %d\n',s1,length(ccb1));
    source=ccb1(s1,1);
    if ~ismember(source,done)
        step2=ccb1(ccb1(:,1)==source,2);
        ccb2=conn_chain_list{2}(ismember(conn_chain_list{2}(:,1),step2) &...
            ismember(reg_chain_list{2}(:,2),reg_set) &...
            all(pref_chain_list{2}(:,8:9)==curr_pref,2),:);%connection_chain_bin1
        if ~isempty(ccb2)
            s2set=unique(ccb2(:,1));
            for s2=1:length(s2set)
                step3=ccb2(ccb2(:,1)==s2set(s2),2);
                ccb3=conn_chain_list{3}(ismember(conn_chain_list{3}(:,1),step3) &...
                    ismember(reg_chain_list{3}(:,2),reg_set) &...
                    all(pref_chain_list{3}(:,9:10)==curr_pref,2),:);%connection_chain_bin1
                if ~isempty(ccb3)
                    s3set=unique(ccb3(:,1));
                    for s3=1:length(s3set)
                        step4=ccb3(ccb3(:,1)==s3set(s3),2);
                        ccb4=conn_chain_list{4}(ismember(conn_chain_list{4}(:,1),step4) &...
                            ismember(reg_chain_list{4}(:,2),reg_set) &...
                            all(pref_chain_list{4}(:,10:11)==curr_pref,2),:);%connection_chain_bin1
                        if ~isempty(ccb4)
                            s4set=unique(ccb4(:,1));
                            for s4=1:length(s4set)
                                step5=ccb4(ccb4(:,1)==s4set(s4),2);
                                ccb5=conn_chain_list{5}(ismember(conn_chain_list{5}(:,1),step5) &...
                                    ismember(reg_chain_list{5}(:,2),reg_set) &...
                                    all(pref_chain_list{5}(:,11:12)==curr_pref,2),:);%connection_chain_bin1
                                if ~isempty(ccb5)
                                    s5set=unique(ccb5(:,1));
                                    for s5=1:length(s5set)
                                        step6=ccb5(ccb5(:,1)==s5set(s5),2);
                                        ccb6=conn_chain_list{6}(ismember(conn_chain_list{6}(:,1),step6) &...
                                            ismember(reg_chain_list{6}(:,2),reg_set) &...
                                            all(pref_chain_list{6}(:,12)==curr_pref,2),:);%connection_chain_bin1
                                        if ~isempty(ccb6)
                                            s6set=unique(ccb6(:,1));
                                            for s6=1:length(s6set)
                                                step7=ccb6(ccb6(:,1)==s6set(s6),2);
                                                for s7=1:length(step7)
                                                    chain_sums(end+1,:)=[source,s2set(s2),s3set(s3),s4set(s4),s5set(s5),s6set(s6),step7(s7)];
                                                    
                                                    if rem(chain_sums,1000000)==0
                                                        save('cross_bin_chain.mat','sums');
                                                    end
                                                end
                                            end
                                        end
                                        
                                    end
                                end
                                
                            end
                        end
                        
                    end
                end
            end
        end
        
        done=[done,source];
    end
end
chainSumsUniq=nan(0,7);
for i=1:length(chain_sums)
    if length(unique(chain_sums(i,:)))==7
        chainSumsUniq(end+1,:)=chain_sums(i,:);
    end
end



    %% export for gephi.

    % csvcell=cell(1,4);
    % csvcell(1,:)={'source','target','type','timeset'};

    % for bin=1:2
    %     conn_chain=conn_chain_list{bin};
    %     reg_chain=reg_chain_list{bin};
    %     
    %     for j=1:length(conn_chain)
    %         if all(ismember(reg_chain(j,:),reg_set))
    %             csvcell(end+1,:)={conn_chain(j,1),conn_chain(j,2),reg_chain(j,1),sprintf('<[%d,%d]>',bin,bin+1)};
    %         end
    %     end
    % end    
    % writecell(csvcell,sprintf('conn_chain.csv',bin));

    
    
    
return

%% 2nd order stats
disp('2nd order stats');
keyboard
stats2order=false;
if stats2order
    load('reverb2nd.mat')
    load('reg_keep.mat')
    arrow2set=unique(reg_2nd,'rows');
    ckeepsel=ismember(arrow2set(:,1),find(ckeep)) & ismember(arrow2set(:,2),find(ckeep));
    arrow2set=arrow2set(ckeepsel,: );
    arrow2count=zeros(length(arrow2set),1);
    for i=1:length(arrow2set)
        arrow2count(i)=nnz(all(reg_2nd==arrow2set(i,:),2));
    end
    load('pair_mat_duo_6s_1_2.mat')
    arrow2pair=zeros(length(arrow2set),1);
    for i=1:length(arrow2set)
        arrow2pair(i)=pair_mat(arrow2set(i,2),arrow2set(i,1));
    end
    pair_sel=arrow2pair>40;
    arrow2ratio=arrow2count(pair_sel)./arrow2pair(pair_sel);
    [~,I]=sort(arrow2ratio,'descend');
    arrow2set=arrow2set(pair_sel,:);
    result=cell(length(arrow2set),2);
    for i=1:length(arrow2set)
        result{i,1}=sprintf('%s->%s',reg_set{arrow2set(I(i),1)},reg_set{arrow2set(I(i),2)});
    end
    result(:,2)=mat2cell(arrow2ratio(I),ones(length(arrow2set),1));
end

%load all 6 bins
