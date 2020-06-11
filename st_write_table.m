function [filename] = st_write_table(table, filename, delimiter)

% ST_WRITE_TABLE write out table fast
%
% Use as
%   [filename] = st_write_table(table, filename)
%   [filename] = st_write_table(table, filename, delimiter)
%
% See also ST_WRITE_RES, ST_APPEND_DATA

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

%tic
if nargin < 3
    delimiter = ',';
end

ft_progress('init', 'text', ['Initializing writing wait...']);

colnames = table.Properties.VariableNames;
rowwritestr = '';
headerwritestr = '';
for iCol = 1:numel(colnames)
    if iCol == 1
        deli = '';
    else
        deli = delimiter;
    end
    headerwritestr = [headerwritestr deli '%s'];
    
    if isnumeric(table.(colnames{iCol}))
        if isfloat(table.(colnames{iCol}))
            if all(floor(table.(colnames{iCol})) == table.(colnames{iCol}))
                rowwritestr = [rowwritestr deli '%d'];
            else
                rowwritestr = [rowwritestr deli '%f'];
            end
        else
            rowwritestr = [rowwritestr deli '%d'];
        end
    else
        rowwritestr = [rowwritestr deli '%s'];
    end
end

%headerwritestr
%rowwritestr
%filename = 'test.csv';

fid = fopen([filename],'Wt'); % the capital W means opening in nonflush mode for faster writing!!!
%fid = fopen([filename],'Wt'); % the capital W means opening in nonflush
%mode for faster writing!!! t is for text mode.
%fid = fopen([filename],'wt'); % the capital W means opening in nonflush
%mode for faster writing!!! t is for textmode
%write header of ouptufile
fprintf(fid,[headerwritestr '\n'],colnames{:});

%checkCell = true;
table = table2cell(table);
rowwritestr = [rowwritestr '\n'];
nRows = size(table,1);
for iRow = 1:nRows
    if mod(iRow,100) == 1
        ft_progress(iRow/nRows, 'writing row %d of %d (%d percent)', iRow, nRows, fix(100*iRow/nRows));
    end
   % if checkCell
   %     thecells  = cellfun(@iscell,table(iRow,:))
   % end
    fprintf(fid,rowwritestr,table{iRow,:});
    %fprintf(fid,rowwritestr,table{1,:});
end
ft_progress(nRows/nRows, 'writing row %d of %d (%d percent)', nRows, nRows, fix(100*nRows/nRows));


fclose(fid);
ft_progress('close');
%toc

%tic
%writetable(table,filename)
%toc

