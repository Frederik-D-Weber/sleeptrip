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
%   cfg.scoringformat    = string, the scoring file format
%                          either:
%                          'custom' with a scoremap file
%                          'zmax'   for hypnodyne corp Zmax exported scoring
%                                   files
%                          'somnomedics_english' or 'somnomedics' for
%                                   somnomedics exproted scoring files
%                          'spisop', 'schlafaus' or 'sleepin' for files
%                                   that are from Schlafaus/SpiSOP or sleepin
%                                   software
%                          'fasst' for scoring files exported from FASST
%                                  scoing software
%                          'u-sleep-30s' for scoring files exported as *.txt
%                                   from U-sleep in 30-s epochs
%                          'nin' .mat files from the NIN
%                          'compumedics_profusion_xml' for Compumedix Profusion software 
%                                   with .xml exported scorings and scoring
%                                   events (note that cfg.scoremap can be
%                                   defined to match the scoring)
%                          'nautus_remlogic_xml' or 'sleep_cure_solutions_xml' or 'remlogic_xml' or 'embla_xml' 
%                                   for .xml files exported from Nautus Embla® RemLogic™ PSG Software  with scorings and scoring events
%                                   (note that cfg.scoremap can be
%                                   defined to match the scoring)
%                          'sleeptrip' for scoring files exported SleepTrip as a .mat
%                                  containing a scoring structure
%                                  (named 'scoring')
%                          'rythm_dreem_json' for scoring files from Rythm
%                                  Dreem head wearable files
%
% optional paramters are
%   cfg.standard         = string, scoring standard either 'aasm' or AASM
%                          or 'rk' for Rechtschaffen&Kales or 'custom' for
%                          which case a scoremap needs to be given (default = 'aasm')
%   cfg.to               = string, if it is set it will convert to a known standard
%                          see ST_SCORINGCONVERT for details
%   cfg.scoringartifactssfile = string, path to the artifact file, by devault [cfg.scoringfile '.artifacts.tsv'] or if not present [cfg.scoringfile '.artifacts.csv] ;
%   cfg.scoringarousalsfile   = string, path to the arousal file, by devault  [cfg.scoringfile '.arousals.tsv'] or if not present [cfg.scoringfile '.arousals.csv] ;
%   cfg.scoringeventsfile     = string, path to the events file, by devault  [cfg.scoringfile '.events.tsv'] or if not present [cfg.scoringfile '.events.csv] ;
%   cfg.forceundefinedto = string,... and force the unsupported scoring labels to this string, see ST_SCORINGCONVERT for details
%   cfg.epochlength      = scalar, epoch length in seconds, (default = 30)
%   cfg.dataoffset       = scalar, offest from data in seconds,
%                          positive number for scoring starting after data
%                          negative number for scoring starting before data
%                          (default = 0)
%  cfg.fileencoding      = string, with encoding e.g. 'UTF-8', see matlab help of
%                          READTABLE for FileEncoding, (default = '', try system specific)
%  cfg.eventsoffset      = string, describing if the
%                          event times in the original scoring file refer
%                          either to 'scoringonset' or 'dataonset' (default
%                          = 'scoringonset'), note that events are then
%                          adapted with reference for dataonset always, see cfg.dataoffset.
%                          This will also apply for read in artifacts,
%                          arousals and all other events
%
% Alternatively one can specify a more general data format with datatype
% with a configuration of only the following necessary options
%   cfg.scoringfile      = string, the scoring file (and path)
%   cfg.scoremap         = structure, a mapping from , see below
%
% ...and additional options
%   cfg.datatype         = string, either 'columns' (e.g. *.tsv, *.csv, *.txt)
%                          or 'xml' (e.g. *.xml), or 'spisop' for (SpiSOP) like input, or 'fasst' (for FASST toolbox
%                          export), (default =
%                          'columns' (or overwritten in case the format is known to be xml))
%   cfg.columndelimimter = string, of the column delimiter, must be either
%                          ',', ' ', '|' or '\t' (a tab) (default = '\t')
%   cfg.skiplines        = scalar, number of lines to skip in file (default = 0)
%   cfg.skiplinesbefore  = string, lines skipped before line found matching string
%   cfg.ignorelines      = Nx1 cell-array with strings that mark filtering/ignoring lines (default = {}, nothing ignored)
%   cfg.selectlines      = Nx1 cell-array with strings that should be selected (default = {}, not specified, all selected)
%   cfg.columnnum        = scalar, the column in which the scoring is stored (default = 1)
%   cfg.exclepochs       = string, if you want to read in a column with excluded epochs indicated either 'yes' or 'no' (default = 'no')
%   cfg.exclcolumnnum    = scalar, if cfg.exclepochs is 'yes' then this is the column in which the exclusion of epochs is stored (default = 2)
%   cfg.exclcolumnstr    = Nx1 cell-array with strings that mark exclusing of epochs, only if cfg.exclepochs is 'yes' this is relevant, (default = {'1', '2', '3'})
%
% Alternatively, if using the function like
%   scoring = st_read_scoring(cfg, tableScoring)
%
%   cfg.ignorelines      = Nx1 cell-array with strings that mark filtering/ignoring lines (default = {}, nothing ignored)
%   cfg.selectlines      = Nx1 cell-array with strings that should be selected (default = {}, not specified, all selected)
%   cfg.datatype         = string, time in seconds
%   cfg.columnnum        = scalar, the column in which the scoring is stored (default = 1)
%   cfg.exclepochs       = string, if you want to read in a column with excluded epochs indicated either 'yes' or 'no' (default = 'no')
%   cfg.exclcolumnnum    = scalar, if cfg.exclepochs is 'yes' then this is the column in which the exclusion of epochs is stored (default = 2)
%   cfg.exclcolumnstr    = Nx1 cell-array with strings that mark exclusing of epochs, only if cfg.exclepochs is 'yes' this is relevant, (default = {'1', '2', '3'})
%
% A scoremap is specified as a structure with the fields
%   scoremap.labelold      = Nx1 cell-array of old labels in file
%   scoremap.labelnew      = Nx1 cell-array of new labels to be named
%   scoremap.unknown       = string, in case the occuring string to label is not
%                            covered in scoremap.old
%
% As an example, for an a SpiSOP to AASM scoring,
%   scoremap = [];
%   scoremap.labelold  = {'0', '1',  '2',  '3',  '4',  '5', '8', '-1'};
%   scoremap.labelnew  = {'W', 'N1', 'N2', 'N3', 'N3', 'R', 'W', '?'};
%   scoremap.unknown   = '?';
%
% As an example, for an a SpiSOP to Rechtschaffen&Kales scoring,
%   scoremap = [];
%   scoremap.labelold  = {'0', '1',  '2',  '3',  '4',  '5', '8', '-1'};
%   scoremap.labelnew  = {'W', 'S1', 'S2', 'S3', 'S4', 'R', 'MT', '?'};
%   scoremap.unknown   = '?';
%
%
%  a scoring can also contain a field arousals, artifacts and more general also
%  events, i.e. scoring.arousals, scoring.artifacts and scoring.events,
%  respectively
%  both are a table of the form

%  ---------------------------------------------
%  event    | start | stop  | duration | channel
%  ---------------------------------------------
%  arousal | 104.1 | 108.1 | 4.0      | EEG
%  arousal | 502.6 | 509.5 | 6.9      | EMG
%  ...     | ...   | ...   | ...      | channel
%  ---------------------------------------------
% 
% note that channel can also be used with placeholders like E*G and that
% stop is the optional field as start and duration are already sufficient
% columns/values (so the stop column is redundant)
%
%
% See also ST_PREPROCESSING

% Copyright (C) 2019-, Frederik D. Weber
% Thanks Roy Cox, Zsófia Zavecz for valuable suggestions and code snippets
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
cfg.eventsoffset       = ft_getopt(cfg, 'eventsoffset', 'scoringonset');

cfg.eventsoffset       = ft_getopt(cfg, 'eventsoffset', 'scoringonset');



% flag to determine which reading option to take
readoption = 'readtable'; % either 'readtable' or 'load'
if nargin > 1
    readoption = 'table';
end

hasArtifacts = false;
hasArousals = false;
hasEvents = false;
hasLightsOff = false;
hasLightsOn = false;

if isfield(cfg,'scoremap') && isfield(cfg,'standard')
    if ~strcmp(cfg.standard,'custom')
        ft_error('Using cfg.scoremap you need to set cfg.standard = ''custom''. To convert to a non-custom standard use the cfg.to option, e.g. cfg.to = ''aasm''')
    end
end

if isfield(cfg,'scoremap')
    if isfield(cfg,'standard') && ~strcmp(cfg.standard,'custom')
        ft_warning('setting cfg.standard = ''custom'' because a cfg.scoremap is defined in the configuration.');
        cfg.standard = 'custom';
    elseif ~isfield(cfg,'standard')
        ft_warning('setting cfg.standard = ''custom'' because a cfg.scoremap is defined in the configuration.');
        cfg.standard = 'custom';
    end
end

if isfield(cfg,'scoremap')
   % ft_warning('setting cfg.standard = ''custom'' because a cfg.scoremap is defined in the configuration.');
    if ~isfield(cfg,'to')
        ft_warning('You might want to define cfg.to as well to be explicit to which standard you want to convert to.');
    end
end

if strcmp(cfg.standard,'custom') && ~isfield(cfg,'scoremap')
    ft_error('if the cfg.standard is set to ''custom'' it requires also a cfg.scoremap as parameter in the configuration.');
end

% optionally get the data from the URL and make a temporary local copy

if nargin<2
    filename = fetch_url(cfg.scoringfile);
    if ~exist(filename, 'file')
        ft_error('The scoring file "%s" file was not found, cannot read in scoring information. No scoring created.', filename);
    end
    
    fileendings = {'.tsv','.csv'};
    
    for iFileending = 1:numel(fileendings)
        cfg.scoringarousalsfile = ft_getopt(cfg, 'scoringarousalsfile', [cfg.scoringfile '.arousals' fileendings{iFileending}] );
        filename_arousals = fetch_url(cfg.scoringarousalsfile);
        filename_arousals_ending = fileendings{iFileending};
        if ~exist(filename_arousals, 'file')
            ft_warning(['No scoring arousal file "%s" file with ' fileendings{iFileending} ' ending was not present.'], filename_arousals);
            filename_arousals = [];
        else
            break
        end
        
    end
    
    for iFileending = 1:numel(fileendings)
        cfg.scoringartifactfile = ft_getopt(cfg, 'scoringartifactfile', [cfg.scoringfile '.artifact' fileendings{iFileending}] );
        filename_artifacts = fetch_url(cfg.scoringartifactfile);
        filename_artifacts_ending = fileendings{iFileending};
        if ~exist(filename_artifacts, 'file')
            ft_warning(['No scoring arousal file "%s" file with ' fileendings{iFileending} ' ending was not present.'], filename_artifacts);
            filename_artifacts = [];
        else
            break
        end
    end
    
    for iFileending = 1:numel(fileendings)
        cfg.scoringeventsfile = ft_getopt(cfg, 'scoringeventsfile', [cfg.scoringfile '.events' fileendings{iFileending}] );
        filename_events = fetch_url(cfg.scoringeventsfile);
        filename_events_ending = fileendings{iFileending};
        if ~exist(filename_events, 'file')
            ft_warning(['No scoring events file "%s" file with ' fileendings{iFileending} ' ending was not present.'], filename_events);
            filename_events = [];
        else
            break
        end
    end
else
    cfg.scoringfile = [];
    cfg.scoringarousalsfile = [];
    cfg.scoringeventsfile = [];
    cfg.scoringartifactfile = [];
end

processTableStucture = true;

scoremap = [];
tableArousal = [];
tableArtifacts = [];
tableEvents = [];

switch  cfg.scoringformat
    case 'custom'
        % do nothing
        scoremap = cfg.scoremap;
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
        
        %cfg.scoremap         = scoremap;
        cfg.columndelimimter = ',';
        cfg.ignorelines      = {'LOUT','LON'};
        cfg.columnnum        = 4;
    case 'zmax_autoscored'
        % ZMax exported one column text from autoscoring
        scoremap = [];
        scoremap.labelold  = {'5', '1', '2', '3', '4', ' U'};
        switch cfg.standard
            case 'aasm'
                scoremap.labelnew  = {'W', 'N1', 'N2', 'N3', 'R', '?'};
            case 'rk'
                ft_warning('the zmax data format is typically in AASM scoring, converting it to Rechtschaffen&Kales might distort results.');
                scoremap.labelnew  = {'W', 'S1', 'S2', 'S3', 'R', '?'};
        end
        scoremap.unknown   = '?';
        
        %cfg.scoremap         = scoremap;
        cfg.columndelimimter = ',';
        %cfg.ignorelines      = {'LOUT','LON'};
        cfg.columnnum        = 1;
        
    case {'somnomedics_english', 'somnomedics'}
        % Somnomedics english version exported profile txt
        scoremap = [];
        scoremap.labelold  = {'Wake', 'N1', 'N2', 'N3', 'REM', 'A'};
        switch cfg.standard
            case 'aasm'
                scoremap.labelnew  = {'W',    'N1', 'N2', 'N3', 'R',   '?'};
            case 'rk'
                ft_warning('the somnomedics data format is typically in AASM scoring, converting it to Rechtschaffen&Kales might distort results.')
                scoremap.labelnew  = {'W',    'S1', 'S2', 'S3', 'R',   '?'};
        end
        scoremap.unknown   = '?';
        
        %cfg.scoremap         = scoremap;
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
        
        %cfg.scoremap         = scoremap;
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
        
        %cfg.scoremap         = scoremap;
        cfg.columnnum        = 1;
    case {'u-sleep-30s'}
        % the .txt version from the website
        scoremap = [];
        scoremap.labelold  = {'Wake', 'N1',  'N2',  'N3',  'N4',  'REM', '?'};
        switch cfg.standard
            case 'aasm'
                scoremap.labelnew  = {'W', 'N1', 'N2', 'N3', 'N3', 'R', '?'};
            case 'rk'
                ft_warning('the u-sleep data format (.txt) is typically in AASM scoring, converting it to Rechtschaffen&Kales might distort results.')
                scoremap.labelnew  = {'W', 'S1', 'S2', 'S3', 'S4', 'R', '?'};
        end
        scoremap.unknown   = '?';
        %cfg.scoremap         = scoremap;
        cfg.columnnum        = 1;
        cfg.skiplines        = 2;
    case {'nin'}
        % the .txt version from the website
        scoremap = [];
        scoremap.labelold  = {'0', '1',  '2',  '3',  'bla',  '5', '?'};
        switch cfg.standard
            case 'aasm'
                scoremap.labelnew  = {'W', 'N1', 'N2', 'N3', 'N3', 'R', '?'};
            case 'rk'
                ft_warning('the u-sleep data format (.txt) is typically in AASM scoring, converting it to Rechtschaffen&Kales might distort results.')
                scoremap.labelnew  = {'W', 'S1', 'S2', 'S3', 'S4', 'R', '?'};
        end
        scoremap.unknown   = '?';
        
        %cfg.scoremap         = scoremap;
        load(filename,'sleepscore')
        rawScores_NB = sleepscore(:,1); %column vector of Neurobit scores (values from [0 1 2 3 5])
        tableScoring = table(rawScores_NB);
        %cfg.to = 'aasm';
        readoption = 'table';
    case {'compumedics_profusion_xml'} % of Compumedix from the Profusion xml export
        if ~isfield(cfg, 'scoremap')
            scoremap = [];
            scoremap.labelold  = {'0', '1', '2', '3', '4', '5', '-1'};
            switch cfg.standard
                case 'aasm'
                    scoremap.labelnew  = {'W', 'N1', 'N2', 'N3', 'N3', 'R', '?'};
                case 'rk'
                    ft_warning('the compumedics_profusion_xml data format (.xml) is typically in AASM scoring, converting it to Rechtschaffen&Kales might distort results.')
                    scoremap.labelnew  = {'W', 'S1', 'S2', 'S3', 'S4', 'R', '?'};
            end
            scoremap.unknown   = '?';
        else
            scoremap = cfg.scoremap;
        end
        
       
      
        import javax.xml.xpath.*
        
        doc = xmlread(filename);
        factory = XPathFactory.newInstance;
        
        xpath = factory.newXPath;
        
        % <SleepStages>
        % <SleepStage>0</SleepStage>
        % <SleepStage>0</SleepStage>
        % <SleepStage>0</SleepStage>
                
        %eventtypes
        expression = xpath.compile('//SleepStages/SleepStage/text()');
        %expression = xpath.compile('//ScoredEvents/ScoredEvent[./Name/text()='Arousal (ASDA)']');

        eventStageValue = getValuesByExpression(doc,expression);
        
        expression = xpath.compile('//EpochLength/text()');
        eventEpochLengthValue = getValuesByExpression(doc,expression);
        if ~isempty(eventEpochLengthValue)
            epochlength = str2num(eventEpochLengthValue{1});
            if epochlength ~= cfg.epochlength
                ft_error('The requested cfg.epochlength = ''%d'' does not match the epoch length defined in the scoring file with %d seconds. Change cfg.epochlength accordingly.',cfg.epochlength, epochlength)
            end
        end
            
        tableScoring = table(eventStageValue,'VariableNames',{'stage'});
        

        % <Name>SpO2 artifact</Name>
        % <Start>17.3</Start>
        % <Duration>39.1</Duration>
        
        expression = xpath.compile('//ScoredEvents/ScoredEvent[contains(./Name/text(), ''Arousal'') or contains(./Name/text(), ''AROUSAL'')]/Name/text()');
        eventNames = getValuesByExpression(doc,expression);
        expression = xpath.compile('//ScoredEvents/ScoredEvent[contains(./Name/text(), ''Arousal'') or contains(./Name/text(), ''AROUSAL'')]/Start/text()');
        eventStarts = getValuesByExpression(doc,expression);
        expression = xpath.compile('//ScoredEvents/ScoredEvent[contains(./Name/text(), ''Arousal'') or contains(./Name/text(), ''AROUSAL'')]/Duration/text()');
        eventDurations = getValuesByExpression(doc,expression);
        expression = xpath.compile('//ScoredEvents/ScoredEvent[contains(./Name/text(), ''Arousal'') or contains(./Name/text(), ''AROUSAL'')]/Input/text()');
        eventInputs = getValuesByExpression(doc,expression);
        
        if ~isempty(eventNames)
            tableArousal = table(eventNames','VariableNames',{'event'});
            tableArousal = cat(2,tableArousal,table(cellfun(@str2num,eventStarts,'UniformOutput',true)',cellfun(@str2num,eventStarts,'UniformOutput',true)'+cellfun(@str2num,eventDurations,'UniformOutput',true)',cellfun(@str2num,eventDurations,'UniformOutput',true)','VariableNames',{'start','stop','duration'}));
            tableArousal = cat(2,tableArousal,table(eventInputs','VariableNames',{'channel'}));
        end
        
        expression = xpath.compile('//ScoredEvents/ScoredEvent[not(contains(./Name/text(), ''Arousal'') or contains(./Name/text(), ''AROUSAL''))]/Name/text()');
        eventNames = getValuesByExpression(doc,expression);
        expression = xpath.compile('//ScoredEvents/ScoredEvent[not(contains(./Name/text(), ''Arousal'') or contains(./Name/text(), ''AROUSAL''))]/Start/text()');
        eventStarts = getValuesByExpression(doc,expression);
        expression = xpath.compile('//ScoredEvents/ScoredEvent[not(contains(./Name/text(), ''Arousal'') or contains(./Name/text(), ''AROUSAL''))]/Duration/text()');
        eventDurations = getValuesByExpression(doc,expression);
        expression = xpath.compile('//ScoredEvents/ScoredEvent[not(contains(./Name/text(), ''Arousal'') or contains(./Name/text(), ''AROUSAL''))]/Input/text()');
        eventInputs = getValuesByExpression(doc,expression);
        
        if ~isempty(eventNames)
            tableEvents = table(eventNames','VariableNames',{'event'});
            tableEvents = cat(2,tableEvents,table(cellfun(@str2num,eventStarts,'UniformOutput',true)',cellfun(@str2num,eventStarts,'UniformOutput',true)'+cellfun(@str2num,eventDurations,'UniformOutput',true)',cellfun(@str2num,eventDurations,'UniformOutput',true)','VariableNames',{'start','stop','duration'}));
            tableEvents = cat(2,tableEvents,table(eventInputs','VariableNames',{'channel'}));
        end
        
        tableScoring = tableScoring;
        
        if ~isempty(tableArousal)
            tableArousal = tableArousal;
            hasArousals = true;
        end
%         if ~isempty(tableLightsOff)
%         	lightsoff_from_scoring_offset = tableLightsOff.start(end);
%             hasLightsOff = true;
%         end
%         if ~isempty(tableLightsOff)
%             lightson_from_scoring_offset = tableLightsOn.start(1);
%             hasLightsOn = true;
%         end
        if ~isempty(tableEvents)
            tableEvents = tableEvents;
            hasEvents = true;
        end
        
        cfg.columnnum        = 1;
        %cfg.exclepochs       = 'no';
        %cfg.exclcolumnnum    = 2;
        processTableStucture = true;
        
        readoption = 'xml';
        
    case {'nautus_remlogic_xml', 'sleep_cure_solutions_xml', 'remlogic_xml', 'embla_xml'} % of Nautus, Embla® RemLogic™ PSG Software
        if ~isfield(cfg, 'scoremap')
            scoremap = [];
            scoremap.labelold  = {'SLEEP-S0', 'SLEEP-S1', 'SLEEP-S2', 'SLEEP-S3', 'SLEEP-S4', 'SLEEP-REM', 'SLEEP-UNSCORED'};
            switch cfg.standard
                case 'aasm'
                    scoremap.labelnew  = {'W', 'N1', 'N2', 'N3', 'N3', 'R', '?'};
                case 'rk'
                    ft_warning('the nautus_remlogic_xml data format (.xml) is typically in AASM scoring, converting it to Rechtschaffen&Kales might distort results.')
                    scoremap.labelnew  = {'W', 'S1', 'S2', 'S3', 'S4', 'R', '?'};
            end
            scoremap.unknown   = '?';
        else
            scoremap = cfg.scoremap;
        end
        
        
        dtformat = 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSS';
        %datetimeDiffToOffsetInSeconds(offsetdts,dts,dtformat)
        %offsetdts = '2021-03-10T00:04:00.773124';
        %dts = '2021-03-10T00:05:01.884235';
      
        import javax.xml.xpath.*
        
        doc = xmlread(filename);
        factory = XPathFactory.newInstance;
        
        xpath = factory.newXPath;
        
        % <Events>
        % <Event>
        % <Type dt:dt="string">SLEEP-S0</Type>
        % <Location dt:dt="string">EEG-C3-M2</Location>
        % <StartTime dt:dt="string">2021-03-09T21:10:21.000000</StartTime>
        % <StopTime dt:dt="string">2021-03-09T21:10:51.000000</StopTime>
        
        
        
        eventtypes_stage = scoremap.labelold;
        %eventtypes_stage = {'SLEEP-S0', 'SLEEP-S1', 'SLEEP-S2', 'SLEEP-S3', 'SLEEP-REM', 'SLEEP-UNSCORED'};
        eventtypes_lighsoff = {'LIGHTS-OFF'};
        eventtypes_lighson = {'LIGHTS-ON'};
        eventtypes_arousal = { 'AROUSAL'};
        
        %eventtypes
        expression = xpath.compile('//EventTypes/EventType[*]/text()');
        eventTypes = getValuesByExpression(doc,expression);
        
        expression = xpath.compile('//Events/Event[*]/Type/text()');
        evenTypeValues = getValuesByExpression(doc,expression);
        
        expression = xpath.compile('//Events/Event[*]/Location/text()');
        evenTypeLocations = getValuesByExpression(doc,expression);
        
        expression = xpath.compile('//Events/Event[*]/StartTime/text()');
        evenTypeStartTimes = getValuesByExpression(doc,expression);
        
        expression = xpath.compile('//Events/Event[*]/StopTime/text()');
        evenTypeStopTimes = getValuesByExpression(doc,expression);
        
        
        first_sleep_stage = find(logical(ismember(evenTypeValues,eventtypes_stage)'),1,'first');
        if ~isempty(first_sleep_stage)
            offsetdts = evenTypeStartTimes{first_sleep_stage};
        end
        
        
        if (numel(evenTypeValues) ~= numel(evenTypeStartTimes)) || (numel(evenTypeValues) ~= numel(evenTypeStopTimes))
            ft_error('the Events are inconsistent with regard to the presence in their XML children: Type, StartTime, StopTime')
        end
        
        %offsetdts = '2021-03-10T00:04:00.773124';
        %dts = '2021-03-10T00:05:01.884235';
        dtformat = 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSS';
        tableScoring = table({},'VariableNames',{'stage'});
        tableScoring = cat(2,tableScoring,table([],[],[],'VariableNames',{'start','stop','duration'}));
        tableScoring = cat(2,tableScoring,table({},'VariableNames',{'channel'}));
        
        tableArousal = table({},'VariableNames',{'event'});
        tableArousal = cat(2,tableArousal,table([],[],[],'VariableNames',{'start','stop','duration'}));
        tableArousal = cat(2,tableArousal,table({},'VariableNames',{'channel'}));
        
        tableLightsOff = table({},'VariableNames',{'event'});
        tableLightsOff = cat(2,tableLightsOff,table([],[],[],'VariableNames',{'start','stop','duration'}));
        tableLightsOff = cat(2,tableLightsOff,table({},'VariableNames',{'channel'}));
        
        tableLightsOn = table({},'VariableNames',{'event'});
        tableLightsOn = cat(2,tableLightsOn,table([],[],[],'VariableNames',{'start','stop','duration'}));
        tableLightsOn = cat(2,tableLightsOn,table({},'VariableNames',{'channel'}));
        
        tableEvents = table({},'VariableNames',{'event'});
        tableEvents = cat(2,tableEvents,table([],[],[],'VariableNames',{'start','stop','duration'}));
        tableEvents = cat(2,tableEvents,table({},'VariableNames',{'channel'}));

        
        for iEvents = 1:numel(evenTypeValues)
            ev_name = strtrim(evenTypeValues{iEvents});
            
            ev_time_start = datetimeDiffToOffsetInSeconds(offsetdts,evenTypeStartTimes{iEvents},dtformat);
            ev_time_stop = datetimeDiffToOffsetInSeconds(offsetdts,evenTypeStopTimes{iEvents},dtformat);
            ev_time_duration =  ev_time_stop-ev_time_start;
            ev_channel = strtrim(evenTypeLocations{iEvents});

            
            switch ev_name
                case eventtypes_stage
                    tableScoring = cat(1,tableScoring,table({ev_name},ev_time_start,ev_time_stop,ev_time_duration,{ev_channel},'VariableNames',{'stage','start','stop','duration','channel'}));
                case eventtypes_lighsoff
                    tableLightsOff = cat(1,tableLightsOff,table({ev_name},ev_time_start,ev_time_stop,ev_time_duration,{ev_channel},'VariableNames',{'event','start','stop','duration','channel'}));
                case eventtypes_lighson
                    tableLightsOn = cat(1,tableLightsOn,table({ev_name},ev_time_start,ev_time_stop,ev_time_duration,{ev_channel},'VariableNames',{'event','start','stop','duration','channel'}));
                case eventtypes_arousal
                    tableArousal = cat(1,tableArousal,table({ev_name},ev_time_start,ev_time_stop,ev_time_duration,{ev_channel},'VariableNames',{'event','start','stop','duration','channel'}));
                otherwise
                	tableEvents = cat(1,tableEvents,table({ev_name},ev_time_start,ev_time_stop,ev_time_duration,{ev_channel},'VariableNames',{'event','start','stop','duration','channel'}));
                    ft_warning('an event with name name %s was not covered by the reading function',ev_name)
            end
        end
        
        tableScoring = tableScoring;
        
        if any(round(tableScoring.duration/cfg.epochlength,3) ~= 1)
            ft_error('cfg.epochlength = ''%d'' does not match the duration of all the epochs detected.',cfg.epochlength)
        end
        if any(round(tableScoring.start + tableScoring.duration,3) ~= round(tableScoring.stop,3))
        	ft_error('There is at least one epoch missing that is missing after detected epochs number ',strjoin(num2str(find(round(tableScoring.start + tableScoring.duration,3) ~= round(tableScoring.stop,3)))'))
        end
        
        if ~isempty(tableArousal)
            tableArousal = tableArousal;
            hasArousals = true;
        end
        if ~isempty(tableArtifacts)
            tableArtifacts = tableArtifacts;
            hasArtifacts = true;
        end
        if ~isempty(tableLightsOff)
        	lightsoff_from_scoring_offset = tableLightsOff.start(end);
            hasLightsOff = true;
        end
        if ~isempty(tableLightsOn)
            lightson_from_scoring_offset = tableLightsOn.start(1);
            hasLightsOn = true;
        end
        if ~isempty(tableEvents)
            tableEvents = tableEvents;
            hasEvents = true;
        end

        cfg.columnnum        = 1;
        %cfg.exclepochs       = 'no';
        %cfg.exclcolumnnum    = 2;
        processTableStucture = true;
        
        readoption = 'xml';
    case {'rythm_dreem_json'}
        % The sleep stages are as follows:
        % 0: wake, 
        % 1: N1,
        % 2: N2, 
        % 3: N3, 
        % 4:REM

        json = read_json(filename);
        
        scoremap = [];
        scoremap.labelold  = {'0','1', '2', '3', '4', '?'};
        switch cfg.standard
            case 'aasm'
                scoremap.labelnew  = {'W', 'N1', 'N2', 'N3', 'R', '?'};
            case 'rk'
                ft_warning('the dreem data format is typically in AASM scoring, converting it to Rechtschaffen&Kales might distort results.');
                scoremap.labelnew  = {'W', 'S1', 'S2', 'S3', 'R', '?'};
        end
        scoremap.unknown   = '?';
        
        json = ft_struct2string(json)';
        
        tableScoring = table(json,'VariableNames',{'stage'});

        cfg.columnnum        = 1;
        %cfg.exclepochs       = 'no';
        %cfg.exclcolumnnum    = 2;
        processTableStucture = true;
        
        readoption = 'json';
        
    case {'sleeptrip'}
        readoption = 'loadmat';
    otherwise
end

switch cfg.datatype 
    case 'columns'
    case 'xml'
        readoption = 'xml';
    case 'json'
        readoption = 'json';
    otherwise
    ft_error('cfg.datatype = ''%s'' unknown',cfg.datatype);
end


switch readoption
    case 'xml'
        
    case 'json'
        
    case 'table';
        tableScoring = tableScoring;
    case 'readtable'
        parampairs = {};
        parampairs = [parampairs, {'ReadVariableNames',false}];
        %parampairs = [parampairs, {'HeaderLines',cfg.skiplines}];
        
        if ~isempty(cfg.columndelimimter)
            parampairs = [parampairs, {'Delimiter',cfg.columndelimimter}];
        end
        
        if ~isempty(cfg.fileencoding)
            parampairs = [parampairs, {'FileEncoding',cfg.fileencoding}];
        end
        try % due to conflicting Matlab conventions in readtable parameters in different matlab versions
            parampairs2 = [parampairs, {'HeaderLines',cfg.skiplines}];
            tableScoring = readtable(filename,parampairs2{:});
        catch
            parampairs2 = [parampairs, {'NumHeaderLines',cfg.skiplines}];
            tableScoring = readtable(filename,parampairs2{:});
        end
    case 'load'
        hyp = load(filename);
        tableScoring = table(hyp(:,1),hyp(:,2));
    case 'loadmat'
        processTableStucture = false;
        scoring = load(filename, 'scoring');
    otherwise
        ft_error('the type %s to read scoring files is not handled. please choose a valid option', readoption);
end

if processTableStucture
    
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
    
    
    % if ~isfield(cfg,'scoremap')
    %     ft_error('No scoremap was defined in the configuration file. Cannot translate scoring.');
    % else
    %
    % end
    
    if numel(scoremap.labelold) ~= numel(scoremap.labelnew)
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
    scoring.epochs(:) = {scoremap.unknown};
    match_cum = logical(zeros(numel(scoring.ori.epochs),1));
    for iLabel = 1:numel(scoremap.labelold)
        old = scoremap.labelold{iLabel};
        new = scoremap.labelnew{iLabel};
        match = cellfun(@(x) strcmp(x, old), scoring.ori.epochs, 'UniformOutput', 1);
        match_cum = match_cum | logical(match(:));
        scoring.epochs(match) = {new};
    end
    
    if any(~match_cum)
        ft_warning('The sleep stages ''%s'' in the original/raw scoring were not covered in the scoremap and have thus been set to ''%s''',strjoin(unique(scoring.ori.epochs(~match_cum))),scoremap.unknown)
    end
    
    scoring.excluded = logical(zeros(1,numel(scoring.ori.epochs)));
    if strcmp(cfg.exclepochs, 'yes')
        %for iLabel = 1
        match = cellfun(@(x)  any(ismember(cfg.exclcolumnstr, x)), scoring.ori.excluded, 'UniformOutput', 1);
        scoring.excluded(match) = true;
        %end
    end
    
    if ~isempty(scoremap)
        scoring.ori.scoremap   = scoremap;
    end
    scoring.ori.scoringfile   = cfg.scoringfile;
    scoring.ori.scoringformat   = cfg.scoringformat;
    scoring.ori.table = tableScoring;
    
    scoring.label = unique(scoremap.labelnew)';
    scoring.label = ({scoring.label{:}})';%assure the vertical orientation
    
    scoring.cfg = cfg;
    scoring.epochlength = cfg.epochlength;
    scoring.dataoffset = cfg.dataoffset;
    scoring.standard = cfg.standard;
    if hasArousals
        switch cfg.eventsoffset
            case 'scoringonset'
                tableArousal.start = tableArousal.start + scoring.dataoffset;
                tableArousal.stop = tableArousal.stop + scoring.dataoffset;
        end
        scoring.arousals = tableArousal;
    end
    if hasLightsOff
        switch cfg.eventsoffset
            case 'scoringonset'
                scoring.lightsoff = lightsoff_from_scoring_offset + scoring.dataoffset;
            otherwise
                scoring.lightsoff = lightsoff_from_scoring_offset;
        end
    end
    if hasLightsOn
        switch cfg.eventsoffset
            case 'scoringonset'
                scoring.lightson = lightson_from_scoring_offset + scoring.dataoffset;
            otherwise
                scoring.lightson = lightson_from_scoring_offset;
        end
    end
    if hasArtifacts
        switch cfg.eventsoffset
            case 'scoringonset'
        tableArtifacts.start = tableArtifacts.start + scoring.dataoffset;
        tableArtifacts.stop = tableArtifacts.stop + scoring.dataoffset;
        end
        scoring.artifacts = tableArtifacts;
    end
    if hasEvents
        switch cfg.eventsoffset
            case 'scoringonset'
        tableEvents.start = tableEvents.start + scoring.dataoffset;
        tableEvents.stop = tableEvents.stop + scoring.dataoffset;
        end
        scoring.events = tableEvents;
    end
    
end

if nargin<2
if ~isempty(filename_arousals)
    
    switch filename_arousals_ending
        case '.tsv'
            tableArousal = readtable(filename_arousals,'FileType','text','ReadVariableNames',true,'Delimiter','\t');
        otherwise
            tableArousal = readtable(filename_arousals,'FileType','text','ReadVariableNames',true);
    end
    
    if ~ismember({'stop'}, tableArousal.Properties.VariableNames)
        tableArousal.stop = tableArousal.start + tableArousal.duration;
    end
    
    if ~ismember({'duration'}, tableArousal.Properties.VariableNames)
        tableArousal.duration = tableArousal.stop - tableArousal.start;
    end
    
    if ~ismember({'channel'}, tableArousal.Properties.VariableNames)
        tableArousal.channel = repmat('all',size(tableArousal,1),1);
    end
    
    switch cfg.eventsoffset
        case 'scoringonset'
            tableArousal.start = tableArousal.start + scoring.dataoffset;
            tableArousal.stop = tableArousal.stop + scoring.dataoffset;
    end
    if isfield(scoring, 'arousals')
        scoring.arousals = unique(cat(1,scoring.arousals,tableArousal));
    else
        scoring.arousals = tableArousal;
    end
    
end
end

if nargin<2
    if ~isempty(filename_artifacts)
        switch filename_artifacts_ending
            case '.tsv'
                tableArtifacts = readtable(filename_artifacts,'FileType','text','ReadVariableNames',true,'Delimiter','\t');
            otherwise
                tableArtifacts = readtable(filename_artifacts,'FileType','text','ReadVariableNames',true);
        end
        
        if ~ismember({'stop'}, tableArtifacts.Properties.VariableNames)
            tableArtifacts.stop = tableArtifacts.start + tableArtifacts.duration;
        end
        
        if ~ismember({'duration'}, tableArtifacts.Properties.VariableNames)
            tableArtifacts.duration = tableArtifacts.stop - tableArtifacts.start;
        end
        
        if ~ismember({'channel'}, tableArtifacts.Properties.VariableNames)
            tableArtifacts.channel = repmat('all',size(tableArtifacts,1),1);
        end
        
        switch cfg.eventsoffset
            case 'scoringonset'
                tableArtifacts.start = tableArtifacts.start + scoring.dataoffset;
                tableArtifacts.stop = tableArtifacts.stop + scoring.dataoffset;
        end
        if isfield(scoring, 'artifacts')
            scoring.artifacts = unique(cat(1,scoring.artifacts,tableArtifacts));
        else
            scoring.artifacts = tableArtifacts;
        end
        
    end
end


if nargin<2
    if ~isempty(filename_events)
        
        switch filename_events_ending
            case '.tsv'
                tableEvents = readtable(filename_events,'FileType','text','ReadVariableNames',true,'Delimiter','\t');
            otherwise
                tableEvents = readtable(filename_events,'FileType','text','ReadVariableNames',true);
        end
        
        if ~ismember({'stop'}, tableEvents.Properties.VariableNames)
            tableEvents.stop = tableEvents.start + tableEvents.duration;
        end
        
        if ~ismember({'duration'}, tableEvents.Properties.VariableNames)
            tableEvents.duration = tableEvents.stop - tableEvents.start;
        end
        
        if ~ismember({'channel'}, tableEvents.Properties.VariableNames)
            tableEvents.channel = repmat('all',size(tableEvents,1),1);
        end
        
        switch cfg.eventsoffset
            case 'scoringonset'
                tableEvents.start = tableEvents.start + scoring.dataoffset;
            tableEvents.stop = tableEvents.stop + scoring.dataoffset;
    end
    if isfield(scoring, 'events')
        scoring.events = unique(cat(1,scoring.events,tableEvents));
    else
        scoring.events = tableEvents;
    end
    
end
end

       

if isfield(cfg,'to')
    cfg_sc = [];
    cfg_sc.to = cfg.to;
    if isfield(cfg,'forceundefinedto')
        cfg_sc.forceundefinedto = cfg.forceundefinedto;
    end
    if strcmp(cfg.standard,'custom')
        cfg_sc.scoremap = scoremap;
    end
    scoring = st_scoringconvert(cfg_sc, scoring);
end
end

function diffsec = datetimeDiffToOffsetInSeconds(offsetdts,dts,dtformat)
%offsetdts = '2021-03-10T00:04:00.773124';
%dts = '2021-03-10T00:05:01.884235';
%dtformat = 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSS';
offsetdt = datetime(offsetdts,'InputFormat',dtformat);
dt = datetime(dts,'InputFormat',dtformat);
diffsec = datenum(dt-offsetdt)*24*3600;
end

function nodevalues = getValuesByExpression(doc,expression)
        import javax.xml.xpath.*

nodeList = expression.evaluate(doc,XPathConstants.NODESET);
nodevalues = {};
for iNode = 1:nodeList.getLength
    node = nodeList.item(iNode-1);
    nodevalues{iNode} = char(node.getNodeValue);
end
end

function json = read_json(filename)
ft_hastoolbox('jsonlab', 1);
json = loadjson(filename);
json = ft_struct2char(json); % convert strings into char-arrays
end