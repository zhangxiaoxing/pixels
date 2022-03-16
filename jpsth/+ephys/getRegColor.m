function c=getRegColor(reg,opt)
arguments
    reg
    opt.large_area (1,1) logical = false
end
persistent idmap
if isempty(idmap)
    idmap=load(fullfile('..','align','reg_ccfid_map.mat'));
end
if opt.large_area

    if any(strcmp(idmap.reg2tree(reg),'PAL'))
       c=[255,0,255]./255;
    elseif any(strcmp(idmap.reg2tree(reg),'MB'))
       c=[0,255,255]./255;
    elseif any(strcmp(idmap.reg2tree(reg),'STR'))
        c=[133,102,40]./255;
    elseif any(strcmp(idmap.reg2tree(reg),'HPF'))
        c=[0,0,255]./255;        
    elseif any(strcmp(idmap.reg2tree(reg),'OLF'))    
        c=[255,0,0]./255;
    elseif any(strcmp(idmap.reg2tree(reg),'TH'))        
        c=[179,34,117]./255;
    elseif any(strcmp(idmap.reg2tree(reg),'Isocortex'))
        c=[0,0,0];
    elseif any(strcmp(idmap.reg2tree(reg),'CTXsp'))
        c=[99,89,160]./255;
    end
else
    if ismember(reg,{'SS','MO'})
        c=[64,98,159]./255;
    elseif ismember(reg,{'ORB','RSP','ACA','PTLp','FRP'})
        c=[220,169,88]./255;
    elseif ismember(reg,{'VIS','AUD'})
        c=[133,102,40]./255;
    elseif ismember(reg,{'AI','VISC','GU','PERI','ECT','TEa'})
        c=[179,34,117]./255;
    elseif ismember(reg,{'PL','ILA'})
        c=[99,89,160]./255;
    elseif ismember(reg,{'AON','DP','TT','PIR'})
        c=[255,0,0]./255;
    elseif ismember(reg,{'BST','SI','GPe','GPi'})
        c=[0,0,255]./255;
    elseif ismember(reg,{'CEA','ACB','LS','CP'})
        c=[255,0,255]./255;
    elseif ismember(reg,{'HIP','RHP'})
        c=[0,0,0];
    elseif ismember(reg,{'EPd'})
        c=[0,255,255]./255;
    else
        c=[127,127,127]./255;
        disp('Missing grouping data')
        disp(reg)
    end
end

end