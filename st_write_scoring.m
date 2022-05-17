function filename = st_write_scoring(cfg, scoring)
% ST_WRITE_SCORING write a scoring out in a file
%
% Use as
%   filepath = st_write_scoring(cfg, scoring)
%
% Configutation parameters
%   cfg.filename         = string, path to the new montage (default =
%   'scoring')
%
% Optional configuration parameters are
%   cfg.datatype         = string, either 'sleeptrip' (i.e. *.mat) or
%                          'numbersincolumns' (e.g. *.tsv, *.csv, *.txt) in the format readable by
%                          SpiSOP or 'spisop'
%                          (default = 'sleeptrip')
%   cfg.columndelimimter = string, of the column delimiter, must be either
%                          ',', ' ', '|' or '\t' (a tab) (default = '\t')
%   cfg.fileencoding     = string, with encoding e.g. 'UTF-8', see Matlab helpf of
%                          fopen for InEncoding, (default = '', try system specific)
%   cfg.to               = string, if it is set it will convert to a known standard
%                          see ST_SCORINGCONVERT for details
%   cfg.scoremap          = a scoremap in case scoring standard is 'custom'
%                          see ST_SCORINGCONVERT for details
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

ttic = tic;
mtic = memtic;
functionname = getfunctionname();
fprintf([functionname ' function started\n']);

hasFilename = false;
if isfield(cfg, 'filename')
    hasFilename = true;
end
cfg.filename           = ft_getopt(cfg, 'filename', 'scoring');
cfg.datatype           = ft_getopt(cfg, 'datatype', 'numbersincolumns');
cfg.columndelimimter   = ft_getopt(cfg, 'columndelimimter', '\t');
cfg.fileencoding       = ft_getopt(cfg, 'fileencoding', '');


if ~strcmp(cfg.datatype, 'numbersincolumns') && ~strcmp(cfg.datatype, 'sleeptrip') && ~strcmp(cfg.datatype, 'spisop')
    % TODO implement the use of tree structure in XML files, e.g. with
    % XPath, or xslt
    ft_error('The datatype parameter only supports the option ''sleeptrip'' or ''numbersincolumns'' or ''spisop'' for now.');
end

if isfield(cfg,'to')
    cfg_sc = [];
    cfg_sc.to = cfg.to;
    if strcmp(scoring.standard,'custom')
        if ~isfield(cfg,'scoremap')
            ft_error('conversion with custom scoring standard requires cfg.scoremap to be defined. see ST_SCORINGCONVERT for details')
        end
        cfg_sc.scoremap = cfg.scoremap;
    end
    scoring = st_scoringconvert(cfg_sc, scoring);
end


switch cfg.datatype
    case 'sleeptrip'
        
    case {'numbersincolumns' 'spisop'}
        if strcmp(scoring.standard,'custom')
            ft_error('writing custom scoring standard to cfg.datatype = ''%s'' is not supported. Chose another data type or use the cfg.to option to convert to a non-custom scoring standard.',cfg.datatype)
        end
end

if ~hasFilename
    switch cfg.datatype
        case 'sleeptrip'
            filename = [cfg.filename '.mat'];
        case 'numbersincolumns'
            switch cfg.columndelimimter
                case {','}
                    filename = [cfg.filename '.csv'];
                case {'|', ' '}
                    filename = [cfg.filename '.txt'];
                case {'\t'}
                    filename = [cfg.filename '.tsv'];
                otherwise
                    ft_error('cfg.columndelimimter = %s is not supported',cfg.columndelimimter)
            end
        case 'spisop'
            % treat all the delimiters the same and write in a csv
            switch cfg.columndelimimter
                case {',', '|', ' ', '\t'}
                    %filename = [cfg.filename '.txt'];
                otherwise
                    ft_error('cfg.columndelimimter = %s is not supported',cfg.columndelimimter)
            end
        otherwise
            ft_error('cfg.datatype = %s is not supported',cfg.datatype)
    end
else
    switch cfg.datatype
        case 'sleeptrip'
            [pathstr,name,ext] = fileparts(cfg.filename);
            if ~strcmp(ext,'.mat')
                filename = [cfg.filename '.mat'];
            else
                filename = cfg.filename;
            end
        case 'numbersincolumns'
            [pathstr,name,ext] = fileparts(cfg.filename);
            switch cfg.columndelimimter
                case {','}
                    if ~strcmp(ext,'.csv')
                        filename = [cfg.filename '.csv'];
                    else
                        filename = cfg.filename;
                    end
                case {'|', ' '}
                    if ~strcmp(ext,'.txt')
                        filename = [cfg.filename '.txt'];
                    else
                        filename = cfg.filename;
                    end
                case {'\t'}
                    if ~strcmp(ext,'.tsv')
                        filename = [cfg.filename '.tsv'];
                    else
                        filename = cfg.filename;
                    end
                otherwise
                    ft_error('cfg.columndelimimter = %s is not supported',cfg.columndelimimter)
            end
        case 'spisop'
            [pathstr,name,ext] = fileparts(cfg.filename);
            switch cfg.columndelimimter
                case {',', '|', ' ', '\t'}
                    if ~strcmp(ext,'.txt')
                        filename = [cfg.filename '.txt'];
                    else
                        filename = cfg.filename;
                    end
                otherwise
                    ft_error('cfg.columndelimimter = %s is not supported',cfg.columndelimimter)
            end
        otherwise
            ft_error('cfg.datatype = %s is not supported',cfg.datatype)
    end
end



switch cfg.datatype
    case 'sleeptrip'
        save(filename,'scoring')
    case {'numbersincolumns' 'spisop'}
        cfg_tmp = [];
        cfg_tmp.to = 'number';
        temp_scoring = st_scoringconvert(cfg_tmp,scoring);
        if isfield(temp_scoring,'confidence')
            temp_hypn = [cellfun(@str2num,temp_scoring.epochs,'UniformOutput',true)' temp_scoring.excluded' temp_scoring.confidence'];
        else
            temp_hypn = [cellfun(@str2num,temp_scoring.epochs,'UniformOutput',true)' temp_scoring.excluded'];
        end
        writetable(array2table(temp_hypn),filename,'FileType','text','WriteVariableNames',false,'Delimiter',cfg.columndelimimter)
        
        if isfield(scoring,'arousals')
        	writetable(scoring.arousals,[filename '.arousals.tsv'],'FileType','text','WriteVariableNames',true,'Delimiter','\t')
        end
        
        if isfield(scoring,'artifacts')
        	writetable(scoring.artifacts,[filename '.artifacts.tsv'],'FileType','text','WriteVariableNames',true,'Delimiter','\t')
        end
        
        if isfield(scoring,'events')
        	writetable(scoring.events,[filename '.events.tsv'],'FileType','text','WriteVariableNames',true,'Delimiter','\t')
        end
    otherwise
        ft_error('cfg.datatype = %s is not supported',cfg.datatype)
end

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end
