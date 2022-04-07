function [filepaths] = st_write_res(cfg, varargin)

% ST_WRITE_RES write out result structures with result origins, type,
% appendition status and date stamp. muliple results can be written out at
% once
%
% Use as
%   [filepaths] = st_write_res(cfg, res1, ...)
%
% Configuration requires at least on of these parameters
% if none of these is defined and no time stamp is chosen, files might be
% overwritten.
%   cfg.prefix  = string added before filename (default = '')
%   cfg.infix   = string added inside filename (default = '')
%   cfg.postfix = string added after filename but before time stamp (default = '')
%
% Optional configuration parameters are
%   cfg.timestamp        = either 'yes' or 'no' if a time stamp should be
%                          added to filename (default = 'yes')
%   cfg.folderstructure  = either 'yes' or 'no' if a folder structure should
%                          be created with the result origin and type
%                          all results will be stored in "/res/..." (default = 'yes')
%   cfg.writetablemethod = either 'st_write_table' or 'writetable' if the
%                          awefully slow matlab 'writetable' should be used
%                          or the quick 'st_write_table' function that is
%                          an order of magnitude faster. (default = 'st_write_table')
%   cfg.delimiter        = string as the delimiter to separate row values ,
%                          save options are: ',' or '\t' (a tab) or '|' (a bar)
%                          maybe but avoid ';' or ' '
%                          if cfg.writetablemethod = 'st_write_table' there
%                          are more options and any string to be used.
%                          (default = ',')
%   cfg.fileending       = string of the file ending, e.g. '.csv' or '.tsv'
%                          default = '.csv' but '.tsv' (if cfg.delimiter =
%                          '\t')
%   
%
% See also ST_APPEND_DATA

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
dt = now;

ttic = tic;
mtic = memtic;
functionname = getfunctionname();
fprintf([functionname ' function started\n']);

% set the defaults
cfg.prefix                 = ft_getopt(cfg, 'prefix',  '');
cfg.infix                  = ft_getopt(cfg, 'infix',   '');
cfg.postfix                = ft_getopt(cfg, 'postfix', '');
cfg.timestamp              = ft_getopt(cfg, 'timestamp', 'yes');
cfg.folderstructure        = ft_getopt(cfg, 'folderstructure', 'yes');
cfg.writetablemethod       = ft_getopt(cfg, 'writetablemethod', 'st_write_table');
cfg.delimiter              = ft_getopt(cfg, 'delimiter', ',');
cfg.fileending             = ft_getopt(cfg, 'fileending', '.csv');

if strcmp(cfg.delimiter, '\t')
    cfg.fileending = '.tsv';
end

if all([isempty(cfg.prefix) isempty(cfg.infix) isempty(cfg.postfix)])
    if ~istrue(cfg.timestamp)
        ft_warning('FILES MIGHT BE OVERWRITTEN!\nNo timestamp and no prefix, no infix nor a postfix string for filenames were defined.')
    else
        ft_warning('No prefix, infix or postfix string for filenames were defined.')
    end
end

if nargin < 2
    ft_error('at least two arguemts need to be provided');
end
nRes = nargin-1;


filepaths = cell(nRes,1);
for iRes = 1:nRes
    res = varargin{iRes};
    appfix = '';
    if isfield(res,'appended')
        if res.appended
            appfix = ['_' 'appended'];
        end
    end
    timestampfix = '';
    if istrue(cfg.timestamp)
        timestampfix = ['_' datestr(dt,'yyyy-mm-dd-HH-MM-SS-FFF')];
    end
    prefix = '';
    if ~isempty(cfg.prefix)
        prefix = [cfg.prefix '_'];
    end
    infix = '';
    if ~isempty(cfg.infix)
        infix = ['_' cfg.infix];
    end
    postfix = '';
    if ~isempty(cfg.postfix)
        postfix = ['_' cfg.postfix];
    end
    
    subfolderpath = '';
    if istrue(cfg.folderstructure)
        subfolderpath = ['res' filesep];
        if ~isfolder([subfolderpath res.ori])
            mkdir([subfolderpath res.ori]);
        end
        if ~isfolder([subfolderpath res.ori filesep res.type])
            mkdir([subfolderpath res.ori filesep res.type]);
        end
        
        subfolderpath = [subfolderpath res.ori filesep res.type filesep];
        
    end
    
    %     postpostfix = '';
    %     if nRes > 1
    %             postpostfix = ['_resnum_' num2str(iRes)];
    %     end
    filepath = [subfolderpath prefix res.ori '_' res.type infix appfix postfix timestampfix cfg.fileending];
    fprintf('writing %s ... \n', filepath);
    switch cfg.writetablemethod
        case 'st_write_table'
            st_write_table(res.table,filepath,cfg.delimiter);%,...
            %'QuoteStrings',false...
            %);
        case 'writetable'
            writetable(res.table,filepath,...
                'FileType','text',...
                'WriteVariableNames',true,...
                'WriteRowNames',false,...
                'Delimiter',cfg.delimiter);%,...
            %'QuoteStrings',false...
            %);
        otherwise
            ft_error(['table write method ' cfg.writetablemethod ' unknown, see the help for all the options.']);
            
    end
    fprintf('... finished %s\n', filepath);
    filepaths{iRes} = filepath;
end


fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end
