function figtable(fh,th,tbl)
uith=uitable(fh,'Data',tbl);
th.Visible='off';
uith.Units=th.Units;
uith.Position=th.Position;
uith.Position(3)=1-th.Position(1)-0.02;
uith.ColumnWidth={uith.Position(3)*fh.Position(3)-2};
uith.RowName=[];
uith.ColumnName=[];
end