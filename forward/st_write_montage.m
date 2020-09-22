function filename = st_write_montage(cfg, montage)
% ST_WRITE_MONTAGE write a montage file
%
% Use as
%   filepath = st_write_montage(cfg, montage)
%
% Configutation parameters
%   cfg.filename         = string, path to the new montage (default =
%   'montage.txt')
%
% Optional configuration parameters are
%   cfg.datatype         = string, either 'columns' (e.g. *.tsv, *.csv, *.txt)
%                          (default = 'columns')
%   cfg.columndelimimter = string, of the column delimiter, must be either
%                          ',', ' ', '|' or '\t' (a tab) (default = ',')
%   cfg.fileencoding     = string, with encoding e.g. 'UTF-8', see matlab help of
%                          READTABLE for FileEncoding, (default = '', try system specific)
%
% See also ST_READ_SCORING, ST_SLEEPONSET

% Copyright (C) 2019-, Frederik D. Weber
%
% This file is part of SleepTrip, see http://www.sleeptrip.org
% for the documentation and details.
%
%    SleepTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    SleepTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    SleepTrip is a branch of FieldTrip, see http://www.fieldtriptoolbox.org
%    and adds funtionality to analyse sleep and polysomnographic data.
%    SleepTrip is under the same license conditions as FieldTrip.
%
%    You should have received a copy of the GNU General Public License
%    along with SleepTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

timerVal = tic;
mtic = memtic;
st = dbstack;
functionname = st.name;

fprintf([functionname ' function started\n']);

hasFilename = false;
if isfield(cfg, 'filename')
    hasFilename = true;
end
cfg.filename           = ft_getopt(cfg, 'filename', 'montage');
cfg.datatype           = ft_getopt(cfg, 'datatype', 'columns');
cfg.columndelimimter   = ft_getopt(cfg, 'columndelimimter', ',');
cfg.skiplines          = ft_getopt(cfg, 'skiplines', 0);
cfg.fileencoding       = ft_getopt(cfg, 'fileencoding', '');


if ~strcmp(cfg.datatype, 'columns')
    % TODO implement the use of tree structure in XML files, e.g. with
    % XPath, or xslt
    ft_error('The datatype parameter only supports the option ''columns'' for now.');
end


if ~hasFilename
    switch cfg.columndelimimter
        case {',', '|', ' '}
            filename = [cfg.filename '.txt'];
        case {'\t'}
            filename = [cfg.filename '.tsv'];
        otherwise
            ft_error('cfg.columndelimimter = %s is not supported',cfg.columndelimimter)
    end
else
    filename = cfg.filename;
end

fid = fopen(filename,'Wt'); % the capital W means opening in nonflush mode for faster writing!!!

%fid = fopen([filename],'Wt'); % the capital W means opening in nonflush
%mode for faster writing!!! t is for text mode.
%fid = fopen([filename],'wt'); % the capital W means opening in nonflush
%mode for faster writing!!! t is for textmode
%write header of ouptufile

% use "old/new" instead of "org/new"
montage = fixoldorg(montage);

headerwritestr = ['new_channels'];
for iCol = 1:numel(montage.labelold)
    headerwritestr = [headerwritestr cfg.columndelimimter '%s'];
end
fprintf(fid,[headerwritestr],montage.labelold{:});
fprintf(fid,['\n']);


for iRow = 1:numel(montage.labelnew)
    rowwritestr = montage.labelnew{iRow};
    for iCol = 1:numel(montage.labelold)
        rowwritestr = [rowwritestr cfg.columndelimimter '%f'];
    end
    row = num2cell(montage.tra(iRow,:));
    fprintf(fid,[rowwritestr],row{:});
    fprintf(fid,['\n']);

    %fprintf(fid,rowwritestr,table{1,:});
end

filepath = fid;
fclose(fid);



fprintf([functionname ' function finished\n']);
toc(timerVal)
memtoc(mtic)
end



