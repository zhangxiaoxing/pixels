function appendfig(opt)
arguments
    opt.datetime (1,1) logical = true
    opt.tag (1,:) char = ''
    opt.close (1,1) logical = false
    opt.fn (1,:) char = 'C:\Users\Libra\OneDrive\Neupix\pct_decoding.pdf'
end

if exist(opt.fn,"file")
    [fid,errmsg] = fopen(opt.fn, 'a');
    if ~isempty(errmsg)
        warning("Error open file "+string(opt.fn));
        return
    else
        fclose(fid);
    end
end


if opt.datetime
    datestr=subsref(char(datetime()),struct('type',{'()'},'subs',{{[1:6,12:17]}}));
    annostr=strjoin({opt.tag,datestr});
else
    annostr=opt.tag;
end
fhandles=get(groot(),'Children');
for hc=reshape(fhandles,1,[])
    ah=annotation(hc,'textbox',[0 0.9 1 0.1], ...
        'String',annostr,'EdgeColor','none','Interpreter','none');

    exportgraphics(hc,opt.fn,'ContentType','vector','Append',true);
    delete(ah);
    if opt.close
        close(hc)
    end
end
fprintf(['[\bUpdated file ',replace(opt.fn,'\','/'),']\b\n'])
