function [res reschannelcolumnname] = st_res_filter(cfg, res)

% ST_RES FILTER fitler the properties and aggregate by channel of a 2-dimensional property (column of a table in the) result structure as
% retruned from SleepTrip functions (e.g. ST_SPINDLES, ST_POWER) 
%
% Use as
%   [res] = st_res_filter(cfg, res)
%   [res channelcolumnname] = st_res_filter(cfg, res)
%
% Configuration parameters can include:
%   cfg.channel       = the channels you want to select (default = 'all');
%   cfg.reschannelcolumnname = name of column in the res.table (default is
%                        chosen by the firs column that contains the string
%                        'channel';
%   cfg.filtercolumns = either a string or a cellstr of the columns you want to filter for, e.g.
%                        {'resnum', 'freq', 'channel'}
%   cfg.filtervalues  = either a string/value or cell of cells/cellstr, with the corresponding values to
%                       the columns defined cfg.filtercolumns, e.g. 
%                        {1:3, [1:20], {'C3', 'C4'}}
%
% See also ST_TOPOPLOTRES, ST_CHANNEL_EVENT_ERP, ST_CHANNEL_EVENT_TFR

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

% set defaults
cfg.channel  = ft_getopt(cfg, 'channel', 'all', 1);

if ~isfield(res,'table')
    ft_error('result structure needs to have a table with the properties in the columns')
end


idx_chan = find(~cellfun(@isempty,regexp(res.table.Properties.VariableNames,'channel','match','once')),1,'first');
if isempty(idx_chan)
    ft_error('result structure table needs to have a column containing called ''channel'' ')
end


cfg.reschannelcolumnname  = ft_getopt(cfg, 'reschannelcolumnname', res.table.Properties.VariableNames{idx_chan});

fprintf([functionname ' function initialized\n']);

%select the channels of interest
res.table = res.table(ismember(res.table.(cfg.reschannelcolumnname),ft_channelselection(cfg.channel,cellstr(unique(res.table.(cfg.reschannelcolumnname))))),:);

if isfield(cfg,'filtercolumns')
    if ~isfield(cfg,'filtervalues')
        ft_error(['if cfg.filtercolumns is defined then cfg.filtervalues also needs to be defined'])
    end
    if ~iscell(cfg.filtercolumns)
        cfg.filtercolumns = {cfg.filtercolumns};
    end
    if ~iscell(cfg.filtervalues)
        cfg.filtervalues = {cfg.filtervalues};
    end
    if numel(cfg.filtercolumns) ~= numel(cfg.filtervalues)
        ft_error(['cfg.filtercolumns needs to have the same amount of elements as cfg.filtervalues'])
    end
    
    for iCol = 1:numel(cfg.filtercolumns)
        fv = res.table{:,{cfg.filtercolumns{iCol}}};
        if ~isempty(fv)
            res.table = res.table(ismember(fv,cfg.filtervalues{iCol}),:);
        end
    end
end

reschannelcolumnname = cfg.reschannelcolumnname;
fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end

