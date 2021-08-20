function c=getRegColor(reg)
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