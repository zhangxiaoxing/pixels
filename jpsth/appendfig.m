function appendfig(opt)
arguments
    opt.datetime (1,1) logical = true
    opt.tag (1,:) char = ''
    opt.close (1,1) logical = false
    opt.path (1,:) = ''
    opt.fn (1,:) char = 'pct_decoding.pdf'
    opt.multi (1,1) logical = false
end

if isempty(opt.path)
    opt.path=fullfile(getenv('userprofile'),'OneDrive','Neupix');
end

fpath=fullfile(opt.path,opt.fn);

if exist(fpath,"file")
    [fid,errmsg] = fopen(fpath, 'a');
    if ~isempty(errmsg)
        warning("Error open file "+string(fpath));
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
if numel(fhandles)>1 && ~opt.multi
    warning("more than 1 figures, continue?")
    keyboard()
end
for hc=reshape(fhandles,1,[])
    ah=annotation(hc,'textbox',[0 0.9 1 0.1], ...
        'String',annostr,'EdgeColor','none','Interpreter','none');

    exportgraphics(hc,fpath,'ContentType','vector','Append',true);
    delete(ah);
    if opt.close
        close(hc)
    end
end
fprintf(['[\bUpdated file ',replace(fpath,'\','/'),']\b\n'])
