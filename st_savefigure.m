function cfg = st_savefigure(cfg, fh)
ttic = tic;
mtic = memtic;
functionname = getfunctionname();
fprintf([functionname ' function started\n']);

if ~isfield(cfg,'functionname')
    cfg.functionname = functionname;
end

if ~isfield(cfg,'subfolder')
    cfg.subfolder = 'unknown';
end

dt = now;


cfg.subfolder               = ft_getopt(cfg, 'subfolder', 'unknown');
cfg.figureoutputformat      = ft_getopt(cfg, 'figureoutputformat', 'png');
cfg.figureoutputresolution  = ft_getopt(cfg, 'figureoutputresolution', 300);
cfg.timestamp               = ft_getopt(cfg, 'timestamp', 'yes');
cfg.folderstructure         = ft_getopt(cfg, 'folderstructure', 'yes');

timestampfix = '';
    if istrue(cfg.timestamp)
        timestampfix = ['_' datestr(dt,'yyyy-mm-dd-HH-MM-SS-FFF')];
    end
    
    subfolderpath = '';
    if istrue(cfg.folderstructure)
        subfolderpath = ['res' filesep];
        if ~isdir([subfolderpath cfg.functionname])
            mkdir([subfolderpath cfg.functionname]);
        end
        if ~isdir([subfolderpath cfg.functionname filesep cfg.subfolder])
            mkdir([subfolderpath cfg.functionname filesep cfg.subfolder]);
        end
        subfolderpath = [subfolderpath cfg.functionname filesep cfg.subfolder];
        [path filename ext] = fileparts(cfg.figureoutputfile);
        postpath = subfolderpath;
        if ~isempty(postpath)
            postpath = [postpath filesep];
        else
            postpath = [postpath];
        end
        cfg.figureoutputfile = [postpath filename timestampfix ext];
    else
        [path filename ext] = fileparts(cfg.figureoutputfile);
        postpath = path;
        if ~isempty(postpath)
            postpath = [postpath filesep];
        else
            postpath = [postpath];
        end
        cfg.figureoutputfile = [postpath filename timestampfix ext];
    end
            
    switch cfg.figureoutputformat
        case 'fig'
            [path filename ext] = fileparts(cfg.figureoutputfile);
            if ~strcomp(ext,['.' cfg.figureoutputformat])
                cfg.figureoutputfile = [cfg.figureoutputfile  '.fig'];
            end
            saveas(fh, [cfg.figureoutputfile  '.fig']);
        case 'eps'
            print(fh,['-d' 'epsc'],['-r' num2str(cfg.figureoutputresolution)],[cfg.figureoutputfile]);
        otherwise
            print(fh,['-d' cfg.figureoutputformat],['-r' num2str(cfg.figureoutputresolution)],[cfg.figureoutputfile]);
    end

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end