function [filename] = st_write_table(table, filename)

% ST_WRITE_TABLE write out table fast
%
% Use as
%   [filename] = st_write_table(table, filename)
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

%table = res_slowwaves_events{1}.table;
%table = cat(1,table,table,table,table,table,table,table,table,table,table,table);
%table = cat(1,table,table,table,table,table,table,table,table,table,table,table);

%tic
delimiter = ',';

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

fid = fopen([filename],'wt');
%write header of ouptufile
fprintf(fid,[headerwritestr '\n'],colnames{:});

table = table2cell(table);
for iRow = 1:size(table,1)
    fprintf(fid,[rowwritestr '\n'],table{iRow,:});
end

fclose(fid)
%toc

%tic
%writetable(table,filename)
%toc

