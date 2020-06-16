function montage = st_read_montage(cfg, filename)
% ST_READ_MONTAGE read a montage file
%
% Use as
%   montage = st_read_montage(cfg, filename)
%
% Configutation parameter can be empty, e.g. cfg = []
%
% Optional configuration parameters are
%   cfg.datatype         = string, either 'columns' (e.g. *.tsv, *.csv, *.txt)
%                          or 'xml' (e.g. *.xml) (default =
%                          'columns')
%   cfg.columndelimimter = string, of the column delimiter, must be either
%                          ',', ' ', '|' or '\t' (a tab) (default = ',')
%   cfg.skiplines        = scalar, number of lines to skip in file (default = 0)%
%   cfg.fileencoding     = string, with encoding e.g. 'UTF-8', see matlab help of
%                          READTABLE for FileEncoding, (default = '', try system specific)
%
% See also ST_READ_SCORING, ST_SLEEPONSET

timerVal = tic;
memtic
st = dbstack;
functionname = st.name;

fprintf([functionname ' function started\n']);

cfg.datatype           = ft_getopt(cfg, 'datatype', 'columns');
cfg.columndelimimter   = ft_getopt(cfg, 'columndelimimter', ',');
cfg.skiplines          = ft_getopt(cfg, 'skiplines', 0);
cfg.fileencoding       = ft_getopt(cfg, 'fileencoding', '');


if ~strcmp(cfg.datatype, 'columns')
    % TODO implement the use of tree structure in XML files, e.g. with
    % XPath, or xslt
    ft_error('The datatype parameter only supports the option ''columns'' for now.');
end


parampairs = {};
parampairs = [parampairs, {'ReadVariableNames',false}];
parampairs = [parampairs, {'HeaderLines',cfg.skiplines}];

if ~isempty(cfg.columndelimimter)
    parampairs = [parampairs, {'Delimiter',cfg.columndelimimter}];
end

if ~isempty(cfg.fileencoding)
    parampairs = [parampairs, {'FileEncoding',cfg.fileencoding}];
end

montageTable = readtable(filename,parampairs{:});

montage = [];
montage.labelold = table2cell(montageTable(1,2:end));
montage.labelnew  = table2cell(montageTable(2:end,1))';
montage.tra = double(cellfun(@str2num,table2cell(montageTable(2:end,2:end))));


fprintf([functionname ' function finished\n']);
toc(timerVal)
memtoc
end



