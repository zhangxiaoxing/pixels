function figtable(fh,th,tbl,opt)
arguments
    fh % figure handle
    th % tile handle
    tbl % content
    opt.title (1,:) char =[]
end
th.Visible='off';
if isa(tbl,"double")
    txtstr=[string(opt.title);string(num2str(tbl,'%.3f || '))];
elseif isa(tbl,"cell")
    if ~isempty(opt.title)
        tbl=[opt.title;tbl];
    end
    txtstr=tbl;
end
txth=annotation(fh,'textbox','String',txtstr,'Interpreter','none');
txth.Units=th.Units;
txth.Position=th.Position;
txth.Position(3)=1-th.Position(1)-0.02;
end