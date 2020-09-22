function [scoring] = st_read_scoring(cfg,tableScoring)

% ST_READ_SCORING reads sleep scoring files and returns
% them in a well defined structure. It is a wrapper around a reader for 
% different importers.
%
% Use as
%   scoring = st_read_scoring(cfg)
%   scoring = st_read_scoring(cfg, tableScoring)
%
% 
% The configuration structure needs to specify
%   cfg.scoringfile      = string, the scoring file (and path)
%   cfg.scoringformat    = string, the scoring file format (and path) 
%
% optional paramters are
%   cfg.standard         = string, scoring standard either 'aasm' or AASM
%                          or 'rk' for Rechtschaffen&Kales or 'custom' for 
%                          which case a scoremap needs to be given (default = 'aasm')
%   cfg.epochlength      = scalar, epoch length in seconds, (default = 30)
%   cfg.dataoffset       = scalar, offest from data in seconds, 
%                          positive number for scoring starting after data
%                          negative number for scoring starting before data
%                          (default = 0)
%  cfg.fileencoding      = string, with encoding e.g. 'UTF-8', see matlab help of
%                          READTABLE for FileEncoding, (default = '', try system specific)
%
% Alternatively one can specify a more general data format with datatype
% with a configuration of onlye the following necessary optoins
%   cfg.scoringfile      = string, the scoring file (and path)
%   cfg.scoremap         = structure, a mapping from , see below
% 
% ...and additional options
%   cfg.datatype         = string, either 'columns' (e.g. *.tsv, *.csv, *.txt) 
%                          or 'xml' (e.g. *.xml), or 'spisop' for (SpiSOP) like input, or 'fasst' (for FASST toolbox
%                          export), (default =
%                          'columns')
%   cfg.columndelimimter = string, of the column delimiter, must be either
%                          ',', ' ', '|' or '\t' (a tab) (default = '\t')
%   cfg.skiplines        = scalar, number of lines to skip in file (default = 0)
%   cfg.skiplinesbefore  = string, lines skipped before line found matching string 
%   cfg.ignorelines      = Nx1 cell-array with strings that mark filtering/ignoring lines (default = {}, nothing ignored)
%   cfg.selectlines      = Nx1 cell-array with strings that should be selected (default = {}, not specified, all selected) 
%   cfg.datatype         = string, time in seconds
%   cfg.columnnum        = scalar, the column in which the scoring is stored (default = 1)
%   cfg.exclepochs       = string, if you want to read in a column with excluded epochs indicated either 'yes' or 'no' (default = 'no')
%   cfg.exclcolumnnum    = scalar, if cfg.exclepochs is 'yes' then this is the column in which the exclusion of epochs is stored (default = 2)
%   cfg.exclcolumnstr    = Nx1 cell-array with strings that mark exclusing of epochs, only if cfg.exclepochs is 'yes' this is relevant, (default = {'1', '2', '3'})
%
% Alternatively, if using the funciont like 
%   [scoring] = st_read_scoring(cfg, tableScoring)
% 
%   cfg.ignorelines      = Nx1 cell-array with strings that mark filtering/ignoring lines (default = {}, nothing ignored)
%   cfg.selectlines      = Nx1 cell-array with strings that should be selected (default = {}, not specified, all selected) 
%   cfg.datatype         = string, time in seconds
%   cfg.columnnum        = scalar, the column in which the scoring is stored (default = 1)
%   cfg.exclepochs       = string, if you want to read in a column with excluded epochs indicated either 'yes' or 'no' (default = 'no')
%   cfg.exclcolumnnum    = scalar, if cfg.exclepochs is 'yes' then this is the column in which the exclusion of epochs is stored (default = 2)
%   cfg.exclcolumnstr    = Nx1 cell-array with strings that mark exclusing of epochs, only if cfg.exclepochs is 'yes' this is relevant, (default = {'1', '2', '3'})
%
%
% A scoremap is specified as a structure with the fields
%   scoremap.labelold      = Nx1 cell-array of old labels in file
%   scoremap.labelnew      = Nx1 cell-array of new labels to be named
%   scoremap.unknown       = string, in case the occuring string to label is not
%                            covered in scoremap.old
%
% As an example, for an A SpiSOP to AASM scoring, 
%   scoremap = [];
%   scoremap.labelold  = {'0', '1',  '2',  '3',  '4',  '5', '8'};
%   scoremap.labelnew  = {'W', 'N1', 'N2', 'N3', 'N3', 'R', 'W'};
%   scoremap.unknown   = '?';
%
% See also ST_PREPROCESSING

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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% set the defaults
cfg.scoringformat      = ft_getopt(cfg, 'scoringformat', 'custom');
cfg.standard           = ft_getopt(cfg, 'standard', 'aasm');
cfg.datatype           = ft_getopt(cfg, 'datatype', 'columns');
cfg.columndelimimter   = ft_getopt(cfg, 'columndelimimter', '\t');
cfg.skiplines          = ft_getopt(cfg, 'skiplines', 0);
cfg.skiplinesbefore    = ft_getopt(cfg, 'skiplinesbefore', '');
cfg.ignorelines        = ft_getopt(cfg, 'ignorelines', {});
cfg.selectlines        = ft_getopt(cfg, 'selectlines', {});
cfg.columnnum          = ft_getopt(cfg, 'columnnum', 1);
cfg.exclepochs         = ft_getopt(cfg, 'exclepochs', 'no');
cfg.exclcolumnnum      = ft_getopt(cfg, 'exclcolumnnum', 2);          
cfg.exclcolumnstr      = ft_getopt(cfg, 'exclcolumnstr', {'1', '2', '3'});
cfg.epochlength        = ft_getopt(cfg, 'epochlength', 30);       
cfg.dataoffset         = ft_getopt(cfg, 'dataoffset', 0);  
cfg.fileencoding       = ft_getopt(cfg, 'fileencoding', '');  

% flag to determine which reading option to take
readoption = 'readtable'; % either 'readtable' or 'load'
if nargin > 1
    readoption = 'table';
end

if strcmp(cfg.standard,'custom') && ~isfield(cfg,'scoremap')
	ft_error('if the cfg.standard is set to ''custom'' it requires also a cfg.scoremap as parameter in the configuration.');
end

switch  cfg.scoringformat
    case 'custom'
        % do nothing
    case 'zmax'
        % ZMax exported csv
        scoremap = [];
        scoremap.labelold  = {'W', 'N1', 'N2', 'N3', 'R', ' U'};
        switch cfg.standard
            case 'aasm'
                scoremap.labelnew  = {'W', 'N1', 'N2', 'N3', 'R', '?'};
            case 'rk'
                ft_warning('the zmax data format is typically in AASM scoring, converting it to Rechtschaffen&Kales might distort results.');
                scoremap.labelnew  = {'W', 'S1', 'S2', 'S3', 'R', '?'};
        end
        scoremap.unknown   = '?';
        
        cfg.scoremap         = scoremap;
        cfg.columndelimimter = ',';
        cfg.ignorelines      = {'LOUT','LON'};
        cfg.columnnum        = 4;

    case 'somnomedics'
        % Somnomedics english version exported profile txt
        scoremap = [];
        scoremap.labelold  = {'Wake', 'N1', 'N2', 'N3', 'REM', 'A'};
        switch cfg.standard
            case 'aasm'
                scoremap.labelnew  = {'W',    'N1', 'N2', 'N3', 'R',   'A'};
            case 'rk'
                ft_warning('the somnomedics data format is typically in AASM scoring, converting it to Rechtschaffen&Kales might distort results.')
                scoremap.labelnew  = {'W',    'S1', 'S2', 'S3', 'R',   'A'};
        end
        scoremap.unknown   = '?';
        
        cfg.scoremap         = scoremap;
        cfg.columndelimimter = ';';
        cfg.skiplines        = 7;
        cfg.columnnum        = 2;
        
    case {'spisop' 'schlafaus' 'sleepin'}
        scoremap = [];
        scoremap.labelold  = {'0', '1',  '2',  '3',  '4',  '5', '8', '-1'};
        switch cfg.standard
            case 'aasm'
                scoremap.labelnew  = {'W', 'N1', 'N2', 'N3', 'N3', 'R', 'W',  '?'};
            case 'rk'
                scoremap.labelnew  = {'W', 'S1', 'S2', 'S3', 'S4', 'R', 'MT', '?'};
        end
        scoremap.unknown   = '?';
        
        cfg.scoremap         = scoremap;
        cfg.columnnum        = 1;
        cfg.exclepochs       = 'yes';
        cfg.exclcolumnnum    = 2;
        cfg.columndelimimter = '';
        cfg.exclcolumnstr = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'};
        readoption = 'load';
    case {'fasst'}
        scoremap = [];
        scoremap.labelold  = {'0', '1',  '2',  '3',  '4',  '5', '7'};
        switch cfg.standard
            case 'aasm'
                scoremap.labelnew  = {'W', 'N1', 'N2', 'N3', 'N3', 'R', '?'};
            case 'rk'
                scoremap.labelnew  = {'W', 'S1', 'S2', 'S3', 'S4', 'R', '?'};
        end
        scoremap.unknown   = '?';
        
        cfg.scoremap         = scoremap;
        cfg.columnnum        = 1;
        
    otherwise
end

if ~strcmp(cfg.datatype, 'columns') 
  % TODO implement the use of tree structure in XML files, e.g. with
  % XPath, or xslt
  ft_error('The datatype parameter only supports the option ''columns'' for now.');
end

% optionally get the data from the URL and make a temporary local copy
if nargin<2
filename = fetch_url(cfg.scoringfile);
if ~exist(filename, 'file')
    ft_error('The scoring file "%s" file was not found, cannot read in scoring information. No scoring created.', filename);
end
else
    cfg.scoringfile = [];
end

switch readoption
    case 'table';
        tableScoring = tableScoring;
    case 'readtable'
        parampairs = {};
        parampairs = [parampairs, {'ReadVariableNames',false}];
        parampairs = [parampairs, {'HeaderLines',cfg.skiplines}];
        
        if ~isempty(cfg.columndelimimter)
            parampairs = [parampairs, {'Delimiter',cfg.columndelimimter}];
        end
        
        if ~isempty(cfg.fileencoding)
            parampairs = [parampairs, {'FileEncoding',cfg.fileencoding}];
        end
        
        tableScoring = readtable(filename,parampairs{:});
        
    case 'load'
        hyp = load(filename);
        tableScoring = table(hyp(:,1),hyp(:,2));
        
    otherwise
        ft_error('the type %s to read scoring files is not handled. please choose a valid option', readoption);
end

tableScoringNcols = size(tableScoring,2);

if tableScoringNcols < cfg.columnnum
    ft_error('The scoring did contain only %d columns.\n The requested column number %d was not present.\n No epochs read in.', tableScoringNcols, cfg.columnnum);
end

if ~isempty(cfg.ignorelines) || ~isempty(cfg.selectlines)
    startline = tableScoring{:,1};
    if isfloat(startline)
        startline = cellstr(num2str(startline));
    end
    
    ignore = logical(zeros(size(tableScoring,1),1));
    if ~isempty(cfg.ignorelines)
        ignore = cellfun(@(x)  any(ismember(cfg.ignorelines, x)), startline, 'UniformOutput', 1);
    end
    
    select = logical(ones(size(tableScoring,1),1));
    if ~isempty(cfg.selectlines)
        select = cellfun(@(x)  any(ismember(cfg.selectlines, x)), startline, 'UniformOutput', 1);
    end
    % get only the rows that matter and update the new table dimension
    tableScoring = tableScoring((select & ~ignore),:);
    tableScoringNcols = size(tableScoring,2);
end


if ~isfield(cfg,'scoremap')
    ft_error('No scoremap was defined in the configuration file. Cannot translate scoring.');
else
    
end

if numel(cfg.scoremap.labelold) ~= numel(cfg.scoremap.labelnew)
    ft_error('Size of cfg.scoremap.labelold and cfg.scoremap.labelold does not match. Cannot translate scoring.');
end

if strcmp(cfg.exclepochs, 'yes')
    if tableScoringNcols < cfg.exclcolumnnum
        ft_warning('The scoring did contain only %d columns.\n The requested column number %d was not present.\n No epochs read for exclusion.', tableScoringNcols, cfg.exclcolumnnum);
    end
end

scoring = [];
scoring.ori = [];
scoring.ori.epochs = tableScoring{:,cfg.columnnum};
if strcmp(cfg.exclepochs, 'yes')
    scoring.ori.excluded = tableScoring{:,cfg.exclcolumnnum};
end
if isfloat(scoring.ori.epochs)
    scoring.ori.epochs = cellstr(arrayfun(@(x) sprintf('%d', x), scoring.ori.epochs, 'UniformOutput', false));
end

if strcmp(cfg.exclepochs, 'yes')
    if isfloat(scoring.ori.excluded)
        scoring.ori.excluded = cellstr(arrayfun(@(x) sprintf('%d', x), scoring.ori.excluded, 'UniformOutput', false));
    end
end

scoring.epochs = cell(1,numel(scoring.ori.epochs));
scoring.epochs(:) = {cfg.scoremap.unknown};
for iLabel = 1:numel(cfg.scoremap.labelold)
    old = cfg.scoremap.labelold{iLabel};
    new = cfg.scoremap.labelnew{iLabel};
    match = cellfun(@(x) strcmp(x, old), scoring.ori.epochs, 'UniformOutput', 1);
    scoring.epochs(match) = {new};
end

scoring.excluded = logical(zeros(1,numel(scoring.ori.epochs)));
if strcmp(cfg.exclepochs, 'yes')
    for iLabel = 1:numel(cfg.scoremap.labelold)
        match = cellfun(@(x)  any(ismember(cfg.exclcolumnstr, x)), scoring.ori.excluded, 'UniformOutput', 1);
        scoring.excluded(match) = true;
    end
end

scoring.ori.scoremap   = cfg.scoremap;
scoring.ori.scoringfile   = cfg.scoringfile;
scoring.ori.scoringformat   = cfg.scoringformat;
scoring.ori.table = tableScoring;

scoring.label = unique(cfg.scoremap.labelnew)';
scoring.label = ({scoring.label{:}})';%assure the vertical orientation

scoring.cfg = cfg;
scoring.epochlength = cfg.epochlength;
scoring.dataoffset = cfg.dataoffset;
scoring.standard = cfg.standard;
