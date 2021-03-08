function [cfg] = st_scorebrowser(cfg, data)

% ST_SCOREBROWSER for scoring sleep data and plotting events
%
% Use as
%   [cfg scoring] = st_scorebrowser(cfg, data)
%   [cfg scoring] = st_scorebrowser(cfg) 
%       where the configuration structure contains the reference to the dataset
%       on your hard disk (see below), or use as
%
%   [cfg scoring] = st_scorebrowser()
%       where you will be asked to manually load the files step by step.
%
% where the input data is a structure as obtained from ST_PREPROCESSING or
% from FT_COMPONENTANALYSIS.
%
% If you want to browse data that is on disk, you have to specify
%   cfg.dataset                 = string with the filename
% Instead of specifying the dataset, you can also explicitely specify the
% name of the file containing the header information and the name of the
% file containing the data, using
%   cfg.datafile                = string with the filename
%   cfg.headerfile              = string with the filename
% Instead of specifiying a dataset, data header, a montage or a scoring you
% can also dot this interactively by setting
%   cfg.datainteractive         = 'yes' (default = 'no') 
%
% Optional configuration parameters are:
%   cfg.scoring  = structure provided by ST_READ_SCORING
%   cfg.standard = either 'aasm' for AASM scoring or 'rk' for
%                  Rechtschaffen&Kales scoring standard
%                  (default = 'aasm')
%   cfg.cuttoscoring = if the data should be cut according to scoring and
%                      dataoffset in case scoring is shorter than data. 
%                      This requires cfg.scoring to be set.
%                      (default = 'no')
%
% The following configuration options are supported:
%   cfg.startepoch              = number of epoch to start view in.
%   cfg.ylim                    = vertical scaling, can be 'maxmin', 'maxabs' or [ymin ymax] (default = 'maxabs')
%   cfg.zlim                    = color scaling to apply to component topographies, 'minmax', 'maxabs' (default = 'maxmin')
%   cfg.epochlength             = duration in seconds for scoring and
%                                 cutting the data up (default = 30)
%   cfg.trl                     = structure that defines the data segments of interest, only applicable for trial-based data
%   cfg.continuous              = 'yes' or 'no' whether the data should be interpreted as continuous or trial-based
%   cfg.channel                 = cell-array with channel labels, see FT_CHANNELSELECTION
%   cfg.plotlabels              = 'yes' (default), 'no', 'some'; whether
%                                 to plot channel labels in vertical viewmode ('some' plots one in every ten
%                                 labels; useful when plotting a large number of channels at a time)
%   cfg.ploteventlabels         = 'type=value', 'colorvalue' (default = 'type=value');
%   cfg.viewmode                = string, 'vertical', currently 'butterfly' or 'component' for visualizing components e.g. from an ICA are not supported (default is 'vertical')
%   cfg.artfctdef.xxx.artifact  = Nx2 matrix with artifact segments see FT_ARTIFACT_xxx functions
%   cfg.selectfeature           = string, name of feature to be selected/added (default = 'visual')
%   cfg.selectmode              = 'markartifact', 'markpeakevent', 'marktroughevent' (default = 'markartifact')
%   cfg.colorgroups             = 'sequential' 'allblack' 'jet' 'hsv' 'labelcharx' (x = xth character in label), 'chantype' or
%                                  vector with length(data/hdr.label) defining groups (default = 'sequential')
%   cfg.channelcolormap         = COLORMAP (default = customized lines map with 15 colors)
%   cfg.selfun                  = string, name of function which is evaluated using the right-click context menu
%                                  The selected data and cfg.selcfg are passed on to this function.
%   cfg.selcfg                  = configuration options for function in cfg.selfun
%   cfg.renderer                = string, 'opengl', 'zbuffer', 'painters', see MATLAB Figure Properties. If this function crashes, you should try 'painters'.
%   cfg.bgcolor                 = background color, either 'white' or
%                                 'dark' (default = 'white')
%   cfg.eegscale                = number, scaling to apply to the EEG channels prior to display
%   cfg.eogscale                = number, scaling to apply to the EOG channels prior to display
%   cfg.ecgscale                = number, scaling to apply to the ECG channels prior to display
%   cfg.emgscale                = number, scaling to apply to the EMG channels prior to display
%   cfg.megscale                = number, scaling to apply to the MEG channels prior to display
%   cfg.gradscale               = number, scaling to apply to the MEG gradiometer channels prior to display (in addition to the cfg.megscale factor)
%   cfg.magscale                = number, scaling to apply to the MEG magnetometer channels prior to display (in addition to the cfg.megscale factor)
%   cfg.mychanscale             = number, scaling to apply to the channels specified in cfg.mychan
%   cfg.mychan                  = Nx1 cell-array with selection of channels
%   cfg.chanscale               = Nx1 vector with scaling factors, one per channel specified in cfg.channel
%   cfg.compscale               = string, 'global' or 'local', defines whether the colormap for the topographic scaling is
%                                  applied per topography or on all visualized components (default 'global')
%
% In case of component viewmode, a layout is required. If no layout is
% give, an attempt is made to construct one from the sensor definition.
% EEG or MEG sensor positions can be present in the data or can be specified as
%   cfg.elec          = structure with electrode positions, see FT_DATATYPE_SENS
%   cfg.grad          = structure with gradiometer definition, see FT_DATATYPE_SENS
%   cfg.elecfile      = name of file containing the electrode positions, see FT_READ_SENS
%   cfg.gradfile      = name of file containing the gradiometer definition, see FT_READ_SENS
%
%
% The scaling to the EEG, EOG, ECG, EMG and MEG channels is optional and
% can be used to bring the absolute numbers of the different channel types
% in the same range (e.g. fT and uV). The channel types are determined from
% the input data using FT_CHANNELSELECTION.
%
% The "artifact" field in the output cfg is a Nx2 matrix comparable to the
% "trl" matrix of FT_DEFINETRIAL. The first column of which specifying the
% beginsamples of an artifact period, the second column contains the
% endsamples of the artifactperiods.
%
% NOTE for debugging: in case the databrowser crashes, use delete(gcf) to
% kill the figure.
%
% See also FT_PREPROCESSING, FT_REJECTARTIFACT, FT_ARTIFACT_EOG,
% FT_ARTIFACT_MUSCLE, FT_ARTIFACT_JUMP, FT_ARTIFACT_MANUAL,
% FT_ARTIFACT_THRESHOLD, FT_ARTIFACT_CLIP, FT_ARTIFACT_ECG,
% FT_COMPONENTANALYSIS
%
% See also FT_DATABROWSER, ST_READ_SCORING, FT_PREPROCESSING, FT_APPLY_MONTAGE

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


% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;


% do the general setup of the function
st_defaults


ft_preamble init
ft_preamble provenance
ft_preamble init
ft_preamble provenance
ft_preamble trackconfig
ft_preamble debug

hasdata = (nargin>1);
hascomp = hasdata && ft_datatype(data, 'comp');

% for backward compatibility
cfg = ft_checkconfig(cfg, 'unused',     {'comps', 'inputfile', 'outputfile'});
cfg = ft_checkconfig(cfg, 'renamed',    {'zscale', 'ylim'});
cfg = ft_checkconfig(cfg, 'renamedval', {'ylim', 'auto', 'maxabs'});
cfg = ft_checkconfig(cfg, 'renamedval', {'selectmode', 'mark', 'markartifact'});

% ensure that the preproc specific options are located in the cfg.preproc substructure
cfg = ft_checkconfig(cfg, 'createsubcfg',  {'preproc'});

% set the defaults
if ~isfield(cfg, 'ylim'),            cfg.ylim = 'maxabs';                 end
if ~isfield(cfg, 'artfctdef'),       cfg.artfctdef = struct;              end
if ~isfield(cfg, 'selectfeature'),   cfg.selectfeature = 'visual';        end % string or cell-array
if ~isfield(cfg, 'selectmode'),      cfg.selectmode = 'markartifact';     end
if ~isfield(cfg, 'blocksize'),       cfg.blocksize = [];                  end % now used for both continuous and non-continuous data, defaulting done below
if ~isfield(cfg, 'preproc'),         cfg.preproc = [];                    end % see preproc for options
if ~isfield(cfg, 'selfun'),          cfg.selfun = [];                     end % default functions: 'simpleFFT','multiplotER','topoplotER','topoplotVAR','movieplotER'
if ~isfield(cfg, 'selcfg'),          cfg.selcfg = [];                     end % defaulting done below, requires layouts/etc to be processed
if ~isfield(cfg, 'colorgroups'),     cfg.colorgroups = 'sequential';      end
if ~isfield(cfg, 'channelcolormap'), cfg.channelcolormap = [0.75 0 0;0 0 1;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];   end
if ~isfield(cfg, 'eegscale'),        cfg.eegscale = [];                   end
if ~isfield(cfg, 'eogscale'),        cfg.eogscale = [];                   end
if ~isfield(cfg, 'ecgscale'),        cfg.ecgscale = [];                   end
if ~isfield(cfg, 'emgscale'),        cfg.emgscale = [];                   end
if ~isfield(cfg, 'megscale'),        cfg.megscale = [];                   end
if ~isfield(cfg, 'magscale'),        cfg.magscale = [];                   end
if ~isfield(cfg, 'gradscale'),       cfg.gradscale = [];                  end
if ~isfield(cfg, 'chanscale'),       cfg.chanscale = [];                  end
if ~isfield(cfg, 'mychanscale'),     cfg.mychanscale = [];                end
if ~isfield(cfg, 'layout'),          cfg.layout = [];                     end
if ~isfield(cfg, 'plotlabels'),      cfg.plotlabels = 'yes';              end
if ~isfield(cfg, 'event'),           cfg.event = [];                      end % this only exists for backward compatibility and should not be documented
if ~isfield(cfg, 'continuous'),      cfg.continuous = [];                 end % the default is set further down in the code, conditional on the input data
if ~isfield(cfg, 'ploteventlabels'), cfg.ploteventlabels = 'type=value';  end
cfg.zlim           = ft_getopt(cfg, 'zlim',          'maxmin');
cfg.compscale      = ft_getopt(cfg, 'compscale',     'global');

cfg.standard        = ft_getopt(cfg, 'standard', 'aasm');
cfg.renderer        = ft_getopt(cfg, 'renderer');
cfg.epochlength     = ft_getopt(cfg, 'epochlength', 30);
cfg.bgcolor         = ft_getopt(cfg, 'bgcolor', 'white');
cfg.cuttoscoring    = ft_getopt(cfg, 'cuttoscoring', 'no');
cfg.viewmode        = ft_getopt(cfg, 'viewmode', 'vertical');
cfg.startepoch      = ft_getopt(cfg, 'startepoch', 1);
cfg.datainteractive = ft_getopt(cfg, 'datainteractive', 'no');


if istrue(cfg.datainteractive)
    ask_again = true;
    asks = 0;
    while ask_again
        try
            asks = asks+1;
            cfg_cd = [];
            cfg_cd.datafile = 'ask';
            cfg_cd.headerfile = 'ask';
            cfg_cd.montagefile = 'ask';
            cfg_cd.datagrammerfile = 'ask';
            cfg_cd.scoringfile = 'ask';
            
            [cfg_dhms] = st_choosedata(cfg_cd);
            
            %             if ~isempty(cfg_dhms.datafile) && ~isempty(cfg_dhms.headerfile)
            %                                     cfg.datafile = cfg_dhms.datafile;
            %                     cfg.headerfile = cfg_dhms.headerfile;
            %             else
            if isempty(cfg_dhms.datafile)
                return
            end
            if ~isempty(cfg_dhms.datafile)
                cfg_pp = [];
                if ~isempty(cfg_dhms.montagefile)
                    cfg_rm = [];
                    cfg_pp.montage = st_read_montage(cfg_rm, cfg_dhms.montagefile);
                end
                if ~isempty(cfg_dhms.datafile)
                    cfg_pp.datafile = cfg_dhms.datafile;
                    cfg_pp.headerfile = cfg_dhms.headerfile;
                else
                    cfg_pp.dataset = cfg_dhms.datafile;
                end
                data = st_preprocessing(cfg_pp);
                hasdata = true;
                
                if ~isempty(cfg_dhms.datagrammerfile)
                    cfg_rdg = [];
                    cfg.datagrammer = st_read_datagrammer(cfg_rdg, cfg_dhms.datagrammerfile);
                    %cfg_adg = [];
                    data = st_apply_datagrammer(data, cfg.datagrammer);
                end
            end
            if ~isempty(cfg_dhms.scoringfile)
                cfg_rs = [];
                cfg_rs.scoringfile = cfg_dhms.scoringfile;
                cfg_rs.scoringformat   = cfg_dhms.scoringformat;
                cfg_rs.standard = cfg_dhms.scoringstandard;
                cfg.scoring = st_read_scoring(cfg_rs);
                cfg.standard = cfg.scoring.standard;
            end
            
            prompt = {'Resample at sampling rate (Hz)'};
            title = 'Update sampling rate?';
            dims = [1 35];
            definput = cellstr(num2str([data.fsample]))';
            updated_samplerate = inputdlg(prompt,title,dims,definput);
            
            if ~isempty(updated_samplerate)
                updated_samplerate = str2num(updated_samplerate{1});
                
                if updated_samplerate ~= data.fsample
                    
                    cfg_rs = [];
                    cfg_rs.resamplefs = updated_samplerate;%frequency at which the data will be resampled (default = 256 Hz)
                    cfg_rs.detrend = 'no';
                    data = ft_resampledata(cfg_rs,data);
                end
            end
            
            ask_again = false;
       catch err
            
            if asks >= 5
                ask_again = false;
                ft_error('could not load any data interactively after ')
                return
            end
            answer_read = questdlg('FAILED to load data or setup. TRY AGAIN?', ...
                'Read in scoring?', ...
                'Yes','No','No');
    
            if ~istrue(answer_read)
                return
            end
        end
    end
end


cfg.channel  = ft_getopt(cfg, 'channel', 'all', 1);

if isfield(cfg,'scoring')
    cfg.blocksize = cfg.scoring.epochlength;
    cfg.epochlength = cfg.scoring.epochlength;   
end

cfg.blocksize  = cfg.epochlength;

cfg.yaxdisteqi = 'no'; % for the hypnogram equidistant places?

cfg.plot_stage_signatures = 'yes';

        cfg.doSleepScoring = 'yes';
        cfg.highlight_scoring_channels = 'yes';
        
        cfg.drawgrid = 'yes';
        cfg.drawgrid_seconds = [0.5 1 3];
        cfg.drawgrid_colors = {[0.9 0.9 0.9] [0.9 0.9 0.9] [0.5 0 0]};
        cfg.drawgrid_LineStyle = {':' '-' '-'};



%cfg_datbrow.artfctdef.EEG.artifact = [300 400; 600 900];
%cfg_datbrow.artfctdef.EMG.artifact = [200 500; 600 950];
cfg.artfctdef.EEG.artifact = [];
cfg.artfctdef.EMG.artifact = [];
cfg.artfctdef.EOG.artifact = [];
cfg.selectfeature = 'EEG';
cfg.selectmode = 'markartifact';
cfg.channel = 1:length(data.label);
cfg.chanscale = ones(1,length(data.label));
    
%     if strcmp(ApplyScalingSettings,'yes')
%         fileScalingSettings = listOfScalingSettingsFiles{iData};
%         curr_channel_scaling_settings_table = readtable([pathInputFolder filesep fileScalingSettings],'Delimiter',',');
%         for iChanEntry = 1:size(curr_channel_scaling_settings_table,1)
%             curr_chan_number = find(strcmp(curr_channel_scaling_settings_table.channel_label(iChanEntry),data.label));
%             if ~isempty(curr_chan_number)
%                 cfg.chanscale(curr_chan_number) = cfg.chanscale(curr_chan_number)*curr_channel_scaling_settings_table.zoom_factor(iChanEntry);
%             end
%         end
%     end
    
    %cfg_datbrow.channel = 6:8;
    %cfg_datbrow.chanscale = cfg_datbrow.chanscale(cfg_datbrow.channel);
    
    %cfg.colorgroups = 'jet';
    cfg.colorgroups       = ft_getopt(cfg, 'colorgroups', 'allblack');

    cfg.event_begin_end_color = [0 1 0];
    cfg.event_begin_end_color2 = [0 0 1];
    
    
%     if strcmp(ApplyEventmappingSettings,'yes')
%         fileEventmappingSettings = listOfEventmappingSettingsFiles{iData};
%         curr_channel_eventmapping_settings_table = readtable([pathInputFolder filesep fileEventmappingSettings],'Delimiter',',');
%     end
    
%     if strcmp(ApplyEventsSelection,'yes')
%         
%         eventsTargetPath = listOfEventsTarget1Paths{iData};
%         dsEventsTarget = dataset('File',eventsTargetPath,'Delimiter',',');
%         
%         if ~isempty(EventsTarget1FilterForColumn)
%             matchIndicator = zeros(size(dsEventsTarget,1),1);
%             for iComb = 1:length(EventsTarget1FilterValues)
%                 %iComp = 1
%                 tempCompTarget = EventsTarget1FilterValues{iComb};
%                 
%                 if iscell(dsEventsTarget.(EventsTarget1FilterForColumn))
%                     matchIndicator = matchIndicator | ( strcmp(dsEventsTarget.(EventsTarget1FilterForColumn), tempCompTarget) );
%                 else
%                     matchIndicator = matchIndicator | ( dsEventsTarget.(EventsTarget1FilterForColumn) ==  tempCompTarget);
%                 end
%             end
%             
%             dsEventsTarget = dsEventsTarget(matchIndicator,:);
%         end
%         
%         nEventsTarget = size(dsEventsTarget,1);
%         
%         cfg.begin_end_events = {};
%         for iCompChan = 1:length(data.label)
%             cfg.begin_end_events{iCompChan} = [];
%             curr_dat_label = data.label(iCompChan);
%             
%             curr_comp_channel = curr_dat_label;
%             
%             if strcmp(ApplyEventmappingSettings,'yes')
%                 curr_chan_map_number = find(strcmp(curr_dat_label,curr_channel_eventmapping_settings_table.channel_label));
%                 if ~isempty(curr_chan_map_number)
%                     curr_mapped_label = curr_channel_eventmapping_settings_table.event_channel_label(curr_chan_map_number);
%                     curr_comp_channel = curr_mapped_label;
%                 end
%             end
%             
%             
%             matchIndicator = ones(nEventsTarget,1);
%             
%             
%             for iComb = 1:length(EventsTarget1CompareColumns)
%                 %iComp = 1 % datasetnum
%                 %iComp = 2 % channel
%                 tempCompTarget = EventsTarget1CompareColumns{iComb};
%                 if iComb == 1
%                     curr_comp = iData;
%                 elseif iComb == 2
%                     curr_comp = curr_comp_channel;
%                 end
%                 if iscell(dsEventsTarget.(tempCompTarget))
%                     if iComb == 1
%                         curr_comp = num2str(curr_comp);
%                     end
%                     matchIndicator = matchIndicator & ( strcmp(curr_comp,dsEventsTarget.(tempCompTarget)) );
%                 else
%                     matchIndicator = matchIndicator & ( curr_comp == dsEventsTarget.(tempCompTarget) );
%                 end
%             end
%             
%             curr_channel_dsEventsTarget = dsEventsTarget(matchIndicator,:);
%             
%             if strcmp(UseSecondColumnAndBothOffsets,'yes')
%                 curr_begins = (curr_channel_dsEventsTarget.(EventsTarget1TimePointColumn) + EventTarget1TimeWindowOffsetTime);
%                 curr_ends = (curr_channel_dsEventsTarget.(EventsTarget1TimePointColumn2) + EventTarget1TimeWindowOffsetTime2);
%             else
%                 curr_begins = ((curr_channel_dsEventsTarget.(EventsTarget1TimePointColumn) + EventTarget1TimeWindowOffsetTime) - EventTarget1TimeWindowPreOffsetTime);
%                 
%                 curr_ends = ((curr_channel_dsEventsTarget.(EventsTarget1TimePointColumn) + EventTarget1TimeWindowOffsetTime) + EventTarget1TimeWindowPostOffsetTime);
%             end
%             
%             
%             cfg.begin_end_events{iCompChan} = round([curr_begins curr_ends]*data.fsample);
%             
%         end
%         
%         
%     end
%     
%     
%     if strcmp(ApplyEventsSelection2,'yes')
%         
%         eventsTargetPath = listOfEventsTarget2Paths{iData};
%         dsEventsTarget = dataset('File',eventsTargetPath,'Delimiter',',');
%         
%         if ~isempty(EventsTarget2FilterForColumn)
%             matchIndicator = zeros(size(dsEventsTarget,1),1);
%             for iComb = 1:length(EventsTarget2FilterValues)
%                 %iComp = 1
%                 tempCompTarget = EventsTarget2FilterValues{iComb};
%                 
%                 if iscell(dsEventsTarget.(EventsTarget2FilterForColumn))
%                     matchIndicator = matchIndicator | ( strcmp(dsEventsTarget.(EventsTarget2FilterForColumn), tempCompTarget) );
%                 else
%                     matchIndicator = matchIndicator | ( dsEventsTarget.(EventsTarget2FilterForColumn) ==  tempCompTarget);
%                 end
%             end
%             
%             dsEventsTarget = dsEventsTarget(matchIndicator,:);
%         end
%         
%         nEventsTarget = size(dsEventsTarget,1);
%         
%         cfg.begin_end_events2 = {};
%         for iCompChan = 1:length(data.label)
%             cfg.begin_end_events2{iCompChan} = [];
%             curr_dat_label = data.label(iCompChan);
%             
%             curr_comp_channel = curr_dat_label;
%             
%             if strcmp(ApplyEventmappingSettings,'yes')
%                 curr_chan_map_number = find(strcmp(curr_dat_label,curr_channel_eventmapping_settings_table.channel_label));
%                 if ~isempty(curr_chan_map_number)
%                     curr_mapped_label = curr_channel_eventmapping_settings_table.event_channel_label(curr_chan_map_number);
%                     curr_comp_channel = curr_mapped_label;
%                 end
%             end
%             
%             
%             matchIndicator = ones(nEventsTarget,1);
%             
%             
%             for iComb = 1:length(EventsTarget2CompareColumns)
%                 %iComp = 1 % datasetnum
%                 %iComp = 2 % channel
%                 tempCompTarget = EventsTarget2CompareColumns{iComb};
%                 if iComb == 1
%                     curr_comp = iData;
%                 elseif iComb == 2
%                     curr_comp = curr_comp_channel;
%                 end
%                 if iscell(dsEventsTarget.(tempCompTarget))
%                     if iComb == 1
%                         curr_comp = num2str(curr_comp);
%                     end
%                     matchIndicator = matchIndicator & ( strcmp(curr_comp,dsEventsTarget.(tempCompTarget)) );
%                 else
%                     matchIndicator = matchIndicator & ( curr_comp == dsEventsTarget.(tempCompTarget) );
%                 end
%             end
%             
%             curr_channel_dsEventsTarget = dsEventsTarget(matchIndicator,:);
%             
%             if strcmp(UseSecondColumnAndBothOffsets2,'yes')
%                 curr_begins = (curr_channel_dsEventsTarget.(EventsTarget2TimePointColumn) + EventTarget2TimeWindowOffsetTime);
%                 curr_ends = (curr_channel_dsEventsTarget.(EventsTarget2TimePointColumn2) + EventTarget2TimeWindowOffsetTime2);
%             else
%                 curr_begins = ((curr_channel_dsEventsTarget.(EventsTarget2TimePointColumn) + EventTarget2TimeWindowOffsetTime) - EventTarget2TimeWindowPreOffsetTime);
%                 
%                 curr_ends = ((curr_channel_dsEventsTarget.(EventsTarget2TimePointColumn) + EventTarget2TimeWindowOffsetTime) + EventTarget2TimeWindowPostOffsetTime);
%             end
%             
%             
%             cfg.begin_end_events2{iCompChan} = round([curr_begins curr_ends]*data.fsample);
%             
%         end
%         
%         
%     end
    
    epochLengthSamples = round(cfg.epochlength*data.fsample);
    nEpochs = floor(size(data.trial{1},2)/epochLengthSamples);
    
    if cfg.startepoch > nEpochs
        cfg.startepoch = nEpochs;
    end
    if cfg.startepoch < 1
        cfg.startepoch = 1;
    end
    
    if isfield(cfg,'scoring')
        scoring_tmp = cfg.scoring;
        cfg_tmp = [];
        cfg_tmp.to = cfg.standard;
        scoring_tmp = st_scoringconvert(cfg_tmp,scoring_tmp);
        
        cfg_tmp = [];
        cfg_tmp.to = 'number';
        scoring_tmp = st_scoringconvert(cfg_tmp,scoring_tmp);
        hypn = [cellfun(@str2num,scoring_tmp.epochs,'UniformOutput',true)' scoring_tmp.excluded']; 
            if size(hypn,1) < nEpochs
                missingEpochs = nEpochs - size(hypn,1);
                %hypn(end+1:end+missingEpochs,:) = [ones(1,missingEpochs,1)*-1 zeros(0,missingEpochs,1)];
                hypn_missing = zeros(missingEpochs,size(hypn,2));
                hypn_missing(:,1) = -1;
                hypn = [hypn ; hypn_missing];
            end
    else
        hypn = [ones(nEpochs,1)*-1 zeros(nEpochs,1)];
    end
    
    
    cfg.hyp_fample = 1;
    cfg.hyp_epochLengthSamples = cfg.hyp_fample*cfg.epochlength;
    
    %if strcmp(ReadInHypnogram,'yes') || strcmp(DoSleepScoring,'yes')
        plot_MA_offset = -5.5;
        plot_confidence_offset = 0.25;
        
        [hypn_plot_interpol hypn_plot_interpol_MA] = interpolate_hypn_for_plot(hypn,cfg.hyp_epochLengthSamples,plot_MA_offset,istrue(cfg.yaxdisteqi));

        %[hypn_plot_interpol hypn_plot_interpol_MA] = interpolate_hypn_for_plot(hypn,cfg.hyp_epochLengthSamples,plot_MA_offset);
        
        %         if (signalOffsetSamples ~= 0)
        %             signalOffsetSamples_downsampled = floor(signalOffsetSeconds*data.fsample);
        %             hypn_plot_interpol = [repmat(0,signalOffsetSamples_downsampled,1); hypn_plot_interpol];
        %             hypn_plot_interpol_MA = [repmat(plot_MA_offset,signalOffsetSamples_downsampled,1); hypn_plot_interpol_MA];
        %         end
        
        cfg.hypn_plot_interpol = hypn_plot_interpol;
        cfg.hypn_plot_interpol_MA = hypn_plot_interpol_MA;
        cfg.plot_MA_offset = plot_MA_offset;
        cfg.plot_confidence_offset = plot_confidence_offset;
        cfg.plotHyp = 'yes';
        cfg.hypn = hypn;
        
        cfg.hypn_plot_interpol_confidence = [];
        if size(cfg.hypn,2) > 2
            if (max(cfg.hypn(:,3)) <= 1) && (min(cfg.hypn(:,3)) >=0)
                [dummy_temp curr_ep_hypn_plot_interpol_confidence] = interpolate_hypn_for_plot(cfg.hypn(:,2:3),cfg.hyp_epochLengthSamples,cfg.plot_confidence_offset,istrue(cfg.yaxdisteqi));
                cfg.hypn_plot_interpol_confidence = curr_ep_hypn_plot_interpol_confidence;
            end
        end
        
    %end
    
        
        cfg.has_more_than_3_channels = true;
        cfg.has_ECG = false;
        
        if numel(data.label) < 3
            cfg.has_more_than_3_channels = false;
            cfg.highlight_scoring_channels = 'no';
            ft_warning('Channel highlighting of EEG, EOG, EMG and ECG disabled');

            %ft_error('Data must conatain at least 3 channels one EOG one EEG and one EMG, please make sure you selected sufficient channels for scoring!');
        end
        
        [numberEEG numberEEG_frontal numberEEG_occipital numberEOG numberEMG numberECG] = getScoringChannelNumbers(data.label);
        
        cfg.score_channel_eeg_number = numberEEG;
        cfg.score_channel_eeg_frontal_number = numberEEG_frontal;
        cfg.score_channel_eeg_occipital_number = numberEEG_occipital;
        cfg.score_channel_eog_number = numberEOG;
        cfg.score_channel_emg_number = numberEMG;
        
        cfg.score_channel_eeg_color           = [1 1 0.8];
        cfg.score_channel_eeg_frontal_color   = [1 0.9 1];
        cfg.score_channel_eeg_occipital_color = [1 0.88 0.72];
        cfg.score_channel_eog_color           = [0.8 1 1];
        cfg.score_channel_emg_color           = [1 0.8 0.8];
        
        %[numberECG] = getECGChannelNumbers(data.label);
        
        if numberECG > 0 % has an ECG channel
            cfg.has_ECG = true;
            cfg.score_channel_ecg_number = numberECG;
            cfg.score_channel_ecg_color = [0.8 1 0.8];
        else
            cfg.has_ECG = false;
        end
    
        
% set core parameters
load_core_cfg
% core_cfg

if (strcmp(core_cfg.bpfilttype,'FIRdesigned') || strcmp(core_cfg.bpfilttype,'IIRdesigned')) && (~strcmp(core_cfg.bpinstabilityfix,'no'))
    error(['band pass filter instability fix not supported for FIRdesigned or IIRdesigned'])
end

if (strcmp(core_cfg.lpfilttype,'FIRdesigned') || strcmp(core_cfg.lpfilttype,'IIRdesigned')) && (~strcmp(core_cfg.lpinstabilityfix,'no'))
    error(['low pass filter instability fix not supported for FIRdesigned or IIRdesigned'])
end

if (strcmp(core_cfg.hpfilttype,'FIRdesigned') || strcmp(core_cfg.hpfilttype,'IIRdesigned')) && (~strcmp(core_cfg.hpinstabilityfix,'no'))
    error(['high pass filter instability fix not supported for FIRdesigned or IIRdesigned'])
end


if ~strcmp(core_cfg.bpfilttype,'FIRdesigned')
    error(['filter type for band pass not supported, only FIRdesigned allowed'])
end

if ~strcmp(core_cfg.lpfilttype,'FIRdesigned')
    error(['filter type for low pass not supported, only FIRdesigned allowed'])
end

if ~(strcmp(core_cfg.hpfilttype,'FIRdesigned') || strcmp(core_cfg.hpfilttype,'IIRdesigned') )
    error(['filter type for band pass not supported, only FIRdesigned or IIRdesigned allowed'])
end



fsample = data.fsample;
use_hp = true;
use_lp = true;
use_bp = true;
adapt_filter_settings_to_toolbox
        
    cfg.core_cfg = core_cfg;
    
    cfg.outputfilespath = [pwd filesep];%[pathOutputFolder filesep];
    %cfg.ouputFilesPrefixString = ouputFilesPrefixString;
    
    
    %cfg.datasetnum = 1;
    %cfg.datasetsPath = datasetsPath;
    
    cfg.StartWithOpenSession = false;
    %cfg.StartWithOpenSession = true;

        
                                        
fprintf([functionname ' function initialized\n']);

if ~isfield(cfg, 'viewmode')
    % butterfly, vertical, component
    if hascomp
        cfg.viewmode = 'component';
    else
        cfg.viewmode = 'butterfly';
    end
end

if ~isempty(cfg.chanscale)
    if ~isfield(cfg,'channel')
        warning('ignoring cfg.chanscale; this should only be used when an explicit channel selection is being made');
        cfg.chanscale = [];
    elseif numel(cfg.channel) ~= numel(cfg.chanscale)
        error('cfg.chanscale should have the same number of elements as cfg.channel');
    end
    
    % make sure chanscale is a column vector, not a row vector
    if size(cfg.chanscale,2) > size(cfg.chanscale,1)
        cfg.chanscale = cfg.chanscale';
    end
end

if ~isempty(cfg.mychanscale) && ~isfield(cfg,'mychan')
    warning('ignoring cfg.mychanscale; no channels specified in cfg.mychan');
    cfg.mychanscale = [];
end

if ~isfield(cfg, 'channel'),
    if hascomp
        if size(data.topo,2)>9
            cfg.channel = 1:10;
        else
            cfg.channel = 1:size(data.topo,2);
        end
    else
        cfg.channel = 'all';
    end
end




%%% check if scoring excedes the data length

if isfield(cfg, 'scoring')
    
    epochLengthSamples = round(cfg.scoring.epochlength*data.fsample);
    offsetSamples = round(cfg.scoring.dataoffset*data.fsample);
    
    if iscell(data.time)
        nSamplesInData = numel(data.time{1});
    else
        nSamplesInData = numel(data.time);
    end
    nPossibleEpochsInData = floor((nSamplesInData-offsetSamples)/epochLengthSamples);
    nEpochsInScoring = numel(cfg.scoring.epochs);
    if (nEpochsInScoring > nPossibleEpochsInData)
        cfg.scoring.epochs = cfg.scoring.epochs(1:nPossibleEpochsInData);
        cfg.scoring.excluded = cfg.scoring.excluded(1:nPossibleEpochsInData);
        ft_warning('only %d epochs possible in data with epoch length of %d and scoring offset of %f seconds \nbut scoring has %d epochs.\nPlease check if scoring matches to the data.\nScoring was shortened to fit data, and some epochs were discarded.',nPossibleEpochsInData,cfg.scoring.epochlength,cfg.scoring.dataoffset,nEpochsInScoring);
    end
end

%%% cut the scoring if necessary
if strcmp(cfg.cuttoscoring,'yes')
    if ~isfield(cfg, 'scoring')
        ft_error('if cfg.cuttoscoring = ''yes'' then also cfg.scoring needs to be set.')
    end
end

if strcmp(cfg.cuttoscoring,'yes') && isfield(cfg, 'scoring')
    epochLengthSamples = round(cfg.scoring.epoch_length * data.fsample);
    data_cut_samples = epochLengthSamples*numel(scorng.epochs);
    
    dataoffset = 0;
    if isfield(cfg.scoring,'dataoffset')
        dataoffset = cfg.scoring.dataoffset;
    end
    dataoffset_samples = round(dataoffset*data.fsample);
    
    
    if (dataoffset_samples ~= 0)
        cfg_tmp = [];
        if (dataoffset_samples > 0) && (dataoffset_samples+1+epochLengthSamples < size(data.trial{1},2))
            cfg_tmp.begsample = dataoffset_samples+1;
            cfg_tmp.endsample = numel(data.time{1});
            
            data = ft_redefinetrial(cfg_tmp,data);
            cfg_tmp = [];
            cfg_tmp.offset = dataoffset_samples+1;
            data = ft_redefinetrial(cfg_tmp,data);
            data.sampleinfo = [1 numel(data.time{1})];
            ft_warning(['IMPORTANT: JUST CUT ' num2str(dataoffset_samples) ' samples at the beginning corresponding to : ' num2str(dataoffset_samples/data.fsample) ' seconds because of scoring data offset!']);
        elseif dataoffset_samples < 0
            for iTrTr = 1:numel(data.trial)
                cfg_tmp.padtype = 'zero';
                data.trial{iTrTr} = ft_preproc_padding(data.trial{iTrTr}, cfg_tmp.padtype, -dataoffset_samples, 0);
            end
            data.time{1} = (0:(size(data.trial{1},2)-1))/data.fsample;
            data.sampleinfo = [1 numel(data.time{1})];
            warning(['IMPORTANT: JUST ADDED ' num2str(dataoffset_samples) ' samples at the beginning corresponding to : ' num2str(dataoffset_samples/data.fsample) ' seconds because of scoring data offset!']);
        end
    end
    
%     if dataoffset_samples > 0
%         data.trial{1}(:,1:dataoffset_samples) = [];
%         data.time{1}(1:dataoffset_samples) = [];
%         data.sampleinfo = [dataoffset_samples+1 data.sampleinfo(2)];
%     end
    
    if size(data.trial{1},2) > data_cut_samples
        data.trial{1}(:,(data_cut_samples+1):end) = [];
        data.time{1}((data_cut_samples+1):end) = [];
        data.sampleinfo = [data.sampleinfo(1) data_cut_samples];
    end
end



if strcmp(cfg.viewmode, 'component')
    % read or create the layout that will be used for the topoplots
    tmpcfg = [];
    tmpcfg.layout = cfg.layout;
    if isempty(cfg.layout)
        warning('No layout specified - will try to construct one using sensor positions');
        if ft_datatype(data, 'meg')
            tmpcfg.grad = ft_fetch_sens(cfg, data);
        elseif ft_datatype(data, 'eeg')
            tmpcfg.elec = ft_fetch_sens(cfg, data);
        else
            error('cannot infer sensor type');
        end
    end
    cfg.layout = ft_prepare_layout(tmpcfg);
elseif ~isempty(cfg.layout)
    tmpcfg = [];
    tmpcfg.layout = cfg.layout;
    cfg.layout = ft_prepare_layout(tmpcfg);
end

if ~isfield(cfg,'drawgrid')
    cfg.drawgrid = 'no';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the defaults and do some preparation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if hasdata
    % save whether data came from a timelock structure
    istimelock = strcmp(ft_datatype(data),'timelock');
    
    % check if the input data is valid for this function
    data = ft_checkdata(data, 'datatype', {'raw', 'comp'}, 'feedback', 'yes', 'hassampleinfo', 'yes');
    % fetch the header from the data structure in memory
    hdr = ft_fetch_header(data);
    
    if isfield(data, 'cfg') && ~isempty(ft_findcfg(data.cfg, 'origfs'))
        % don't use the events in case the data has been resampled
        warning('the data has been resampled, not showing the events');
        event = [];
    elseif ~isempty(cfg.event)
        % use the events that the user passed in the configuration
        event = cfg.event;
    else
        % fetch the events from the data structure in memory
        event = ft_fetch_event(data);
    end
    
    cfg.channel = ft_channelselection(cfg.channel, hdr.label);
    chansel = match_str(data.label, cfg.channel);
    Nchans  = length(chansel);
    
    if isempty(cfg.continuous)
        if numel(data.trial) == 1 && ~istimelock
            cfg.continuous = 'yes';
        else
            cfg.continuous = 'no';
        end
    else
        if strcmp(cfg.continuous, 'yes') && (numel(data.trial) > 1)
            warning('interpreting trial-based data as continous, time-axis is no longer appropriate. t(0) now corresponds to the first sample of the first trial, and t(end) to the last sample of the last trial')
        end
    end
    
    % this is how the input data is segmented
    trlorg = zeros(numel(data.trial), 3);
    trlorg(:,[1 2]) = data.sampleinfo;
    
    % recreate offset vector (databrowser depends on this for visualisation)
    for ntrl = 1:numel(data.trial)
        trlorg(ntrl,3) = time2offset(data.time{ntrl}, data.fsample);
    end
    Ntrials = size(trlorg, 1);
    
else
    % check if the input cfg is valid for this function
    cfg = ft_checkconfig(cfg, 'dataset2files', {'yes'});
    cfg = ft_checkconfig(cfg, 'required', {'headerfile', 'datafile'});
    cfg = ft_checkconfig(cfg, 'renamed',    {'datatype', 'continuous'});
    cfg = ft_checkconfig(cfg, 'renamedval', {'continuous', 'continuous', 'yes'});
    % read the header from file
    hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat);
    
    if isempty(cfg.continuous)
        if hdr.nTrials==1
            cfg.continuous = 'yes';
        else
            cfg.continuous = 'no';
        end
    end
    
    if ~isempty(cfg.event)
        % use the events that the user passed in the configuration
        event = cfg.event;
    else
        % read the events from file
        event = ft_read_event(cfg.dataset);
    end
    
    cfg.channel = ft_channelselection(cfg.channel, hdr.label);
    chansel = match_str(hdr.label, cfg.channel);
    Nchans  = length(chansel);
    
    if strcmp(cfg.continuous, 'yes')
        Ntrials = 1;
    else
        Ntrials = hdr.nTrials;
    end
    
    % FIXME in case of continuous=yes the trl should be [1 hdr.nSamples*nTrials 0]
    % and a scrollbar should be used
    
    % construct trl-matrix for data from file on disk
    trlorg = zeros(Ntrials,3);
    if strcmp(cfg.continuous, 'yes')
        trlorg(1, [1 2]) = [1 hdr.nSamples*hdr.nTrials];
    else
        for k = 1:Ntrials
            trlorg(k,[1 2]) = [1 hdr.nSamples] + [hdr.nSamples hdr.nSamples] .* (k-1);
        end
    end
end % if hasdata
if strcmp(cfg.continuous,'no') && isempty(cfg.blocksize)
    cfg.blocksize = (trlorg(1,2) - trlorg(1,1)+1) ./ hdr.Fs;
elseif strcmp(cfg.continuous,'yes') && isempty(cfg.blocksize)
    cfg.blocksize = 1;
end



% FIXME make a check for the consistency of cfg.continous, cfg.blocksize, cfg.trl and the data header

if Nchans == 0
    error('no channels to display');
end

if Ntrials == 0
    error('no trials to display');
end

if ischar(cfg.selectfeature)
    % ensure that it is a cell array
    cfg.selectfeature = {cfg.selectfeature};
end
if ~isempty(cfg.selectfeature)
    for i=1:length(cfg.selectfeature)
        if ~isfield(cfg.artfctdef, cfg.selectfeature{i})
            cfg.artfctdef.(cfg.selectfeature{i})          = [];
            cfg.artfctdef.(cfg.selectfeature{i}).artifact = zeros(0,2);
        end
    end
end

% determine the vertical scaling
if ischar(cfg.ylim)
    if hasdata
        % the first trial is used to determine the vertical scaling
        dat = data.trial{1}(chansel,:);
    else
        % one second of data is read from file to determine the vertical scaling
        dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', 1, 'endsample', round(hdr.Fs), 'chanindx', chansel, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat, 'headerformat', cfg.headerformat);
    end % if hasdata
    minval = min(dat(:));
    maxval = max(dat(:));
    switch cfg.ylim
        case 'maxabs'
            maxabs   = max(abs([minval maxval]));
            scalefac = 10^(fix(log10(maxabs)));
            maxabs   = (round(maxabs / scalefac * 100) / 100) * scalefac;
            cfg.ylim = [-maxabs maxabs];
        case 'maxmin'
            cfg.ylim = [minval maxval];
        otherwise
            error('unsupported value for cfg.ylim');
    end % switch ylim
    % zoom in a bit when viemode is vertical
    if strcmp(cfg.viewmode,'vertical')
        cfg.ylim = cfg.ylim/10;
    end
else
    if (numel(cfg.ylim) ~= 2) || ~isnumeric(cfg.ylim)
        error('cfg.ylim needs to be a 1x2 vector [ymin ymax], describing the upper and lower limits')
    end
end

% determine coloring of channels
if hasdata
    labels_all = data.label;
else
    labels_all= hdr.label;
end
if size(cfg.channelcolormap,2) ~= 3
    error('cfg.channelcolormap is not valid, size should be Nx3')
end

if isnumeric(cfg.colorgroups)
    % groups defined by user
    if length(labels_all) ~= length(cfg.colorgroups)
        error('length(cfg.colorgroups) should be length(data/hdr.label)')
    end
    R = cfg.channelcolormap(:,1);
    G = cfg.channelcolormap(:,2);
    B = cfg.channelcolormap(:,3);
    chancolors = [R(cfg.colorgroups(:)) G(cfg.colorgroups(:)) B(cfg.colorgroups(:))];
    
elseif strcmp(cfg.colorgroups, 'allblack')
    chancolors = zeros(length(labels_all),3);
    
elseif strcmp(cfg.colorgroups, 'jet')
    chancolors = jet(length(labels_all));
    
elseif strcmp(cfg.colorgroups, 'hsv')
    chancolors = hsv(length(labels_all));
    
elseif strcmp(cfg.colorgroups, 'chantype')
    type = ft_chantype(labels_all);
    [tmp1 tmp2 cfg.colorgroups] = unique(type);
    fprintf('%3d colorgroups were identified\n',length(tmp1))
    R = cfg.channelcolormap(:,1);
    G = cfg.channelcolormap(:,2);
    B = cfg.channelcolormap(:,3);
    chancolors = [R(cfg.colorgroups(:)) G(cfg.colorgroups(:)) B(cfg.colorgroups(:))];
    
elseif strcmp(cfg.colorgroups(1:9), 'labelchar')
    % groups determined by xth letter of label
    labelchar_num = str2double(cfg.colorgroups(10));
    vec_letters = num2str(zeros(length(labels_all),1));
    for iChan = 1:length(labels_all)
        vec_letters(iChan) = labels_all{iChan}(labelchar_num);
    end
    [tmp1 tmp2 cfg.colorgroups] = unique(vec_letters);
    fprintf('%3d colorgroups were identified\n',length(tmp1))
    R = cfg.channelcolormap(:,1);
    G = cfg.channelcolormap(:,2);
    B = cfg.channelcolormap(:,3);
    chancolors = [R(cfg.colorgroups(:)) G(cfg.colorgroups(:)) B(cfg.colorgroups(:))];
    
elseif strcmp(cfg.colorgroups, 'sequential')
    % no grouping
    chancolors = lines(length(labels_all));
    
else
    error('do not understand cfg.colorgroups')
end

% collect the artifacts that have been detected from cfg.artfctdef.xxx.artifact
artlabel = fieldnames(cfg.artfctdef);
sel      = zeros(size(artlabel));
artifact = cell(size(artlabel));
for i=1:length(artlabel)
    sel(i) = isfield(cfg.artfctdef.(artlabel{i}), 'artifact');
    if sel(i)
        artifact{i} = cfg.artfctdef.(artlabel{i}).artifact;
        fprintf('detected %3d %s artifacts\n', size(artifact{i}, 1), artlabel{i});
    end
end

% get the subset of the artfctdef fields that seem to contain artifacts
artifact = artifact(sel==1);
artlabel = artlabel(sel==1);

if length(artlabel) > 9
    error('only up to 9 artifacts groups supported')
end

% make artdata representing all artifacts in a "raw data" format
datendsample = max(trlorg(:,2));

artdata = [];
artdata.trial{1}       = convert_event(artifact, 'boolvec', 'endsample', datendsample); % every artifact is a "channel"
artdata.time{1}        = offset2time(0, hdr.Fs, datendsample);
artdata.label          = artlabel;
artdata.fsample        = hdr.Fs;
artdata.cfg.trl        = [1 datendsample 0];

% determine amount of unique event types (for cfg.ploteventlabels)
if ~isempty(event) && isstruct(event)
    eventtypes = unique({event.type});
else
    eventtypes = [];
end


%FW begin
if isfield(cfg,'begin_end_events')
    begin_end_events_per_channel = cfg.begin_end_events;
    cfg.times_ind_per_channel = {};
    data_samples_length = length(data.time{1,1});
    for channelIndex = 1:length(cfg.channel)
        curr_begins_ends = begin_end_events_per_channel{channelIndex};
        curr_times_ind = zeros(1,data_samples_length);
        if ~isempty(curr_begins_ends)
            for iEv = 1:size(curr_begins_ends,1)
                temp_beg = curr_begins_ends(iEv,1);
                temp_end = curr_begins_ends(iEv,2);
                curr_times_ind(temp_beg:temp_end) = 1;
            end
        end
        cfg.times_ind_per_channel{channelIndex} = curr_times_ind;
    end % for each of the limEvents each channels
end

if isfield(cfg,'begin_end_events2')
    begin_end_events_per_channel = cfg.begin_end_events2;
    cfg.times_ind_per_channel2 = {};
    data_samples_length = length(data.time{1,1});
    for channelIndex = 1:length(cfg.channel)
        curr_begins_ends = begin_end_events_per_channel{channelIndex};
        curr_times_ind = zeros(1,data_samples_length);
        if ~isempty(curr_begins_ends)
            for iEv = 1:size(curr_begins_ends,1)
                temp_beg = curr_begins_ends(iEv,1);
                temp_end = curr_begins_ends(iEv,2);
                curr_times_ind(temp_beg:temp_end) = 1;
            end
        end
        cfg.times_ind_per_channel2{channelIndex} = curr_times_ind;
    end % for each of the limEvents each channels
end
%FW end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up default defselfuns
% use two cfg sections
% cfg.selfun - labels that are presented in rightclick menu, and is appended using ft_getuserfun(..., 'browse') later on to create a function handle
% cfg.selcfg - cfgs for functions to be executed
defselfun = [];
defselcfg = [];
% simplefft
defselfun{1} = 'simpleFFT';
tmpcfg = [];
tmpcfg.chancolors = chancolors;
defselcfg{1}  = tmpcfg;
if ~strcmp(cfg.doSleepScoring,'yes')
    % multiplotER
    defselfun{2}  = 'multiplotER';
    tmpcfg = [];
    tmpcfg.layout = cfg.layout;
    defselcfg{2}  = tmpcfg;
    % topoplotER
    defselfun{3}  = 'topoplotER';
    tmpcfg = [];
    tmpcfg.layout = cfg.layout;
    defselcfg{3}  = tmpcfg;
    % topoplotVAR
    defselfun{4}  = 'topoplotVAR';
    tmpcfg = [];
    tmpcfg.layout = cfg.layout;
    defselcfg{4}  = tmpcfg;
    % movieplotER
    defselfun{5}  = 'movieplotER';
    tmpcfg = [];
    tmpcfg.layout = cfg.layout;
    tmpcfg.interactive = 'yes';
    defselcfg{5}  = tmpcfg;
end

% add defselfuns to user-specified defselfuns
if ~iscell(cfg.selfun) && ~isempty(cfg.selfun)
    cfg.selfun = {cfg.selfun};
    cfg.selfun = [cfg.selfun defselfun];
    % do the same for the cfgs
    cfg.selcfg = {cfg.selcfg}; % assume the cfg isnt a cell
    cfg.selcfg = [cfg.selcfg defselcfg];
else
    cfg.selfun = defselfun;
    cfg.selcfg = defselcfg;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up the data structures used in the GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% opt represents the global data/settings, it should contain
% - the original data, epoched or continuous
% - the artifacts represented as continuous data
% - the redraw_cb settings
% - the preproc   settings
% - the select_range_cb settings (also used in keyboard_cb)

% these elements are stored inside the figure so that the callback routines can modify them
opt = [];
if hasdata
    opt.orgdata   = data;
    % FW begin
    data = [];
    % FW end
else
    opt.orgdata   = [];      % this means that it will look in cfg.dataset
end
if strcmp(cfg.continuous, 'yes')
    opt.trialviewtype = 'epoch';
else
    opt.trialviewtype = 'trial';
end
opt.artdata     = artdata;
opt.hdr         = hdr;
opt.event       = event;
opt.trlop       = cfg.startepoch;          % the active trial being displayed
opt.ftsel       = find(strcmp(artlabel,cfg.selectfeature)); % current artifact/feature being selected
opt.trlorg      = trlorg;
opt.fsample     = hdr.Fs;
opt.artcolors   = [0.9686 0.7608 0.7686; 0.7529 0.7098 0.9647; 0.7373 0.9725 0.6824;0.8118 0.8118 0.8118; 0.9725 0.6745 0.4784; 0.9765 0.9176 0.5686; 0.6863 1 1; 1 0.6863 1; 0 1 0.6000];
opt.chancolors  = chancolors;
opt.cleanup     = false;      % this is needed for a corrent handling if the figure is closed (either in the corner or by "q")
opt.chanindx    = [];         % this is used to check whether the component topographies need to be redrawn
opt.eventtypes  = eventtypes;
opt.eventtypescolors = [0 0 0; 1 0 0; 0 0 1; 0 1 0; 1 0 1; 0.5 0.5 0.5; 0 1 1; 1 1 0];
opt.eventtypecolorlabels = {'black', 'red', 'blue', 'green', 'cyan', 'grey', 'light blue', 'yellow'};
opt.nanpaddata  = []; % this is used to allow horizontal scaling to be constant (when looking at last segment continous data, or when looking at segmented/zoomed-out non-continous data)
opt.trllock     = []; % this is used when zooming into trial based data

% save original layout when viewmode = component
if strcmp(cfg.viewmode,'component')
    opt.layorg    = cfg.layout;
end

% determine labelling of channels
if strcmp(cfg.plotlabels, 'yes')
    opt.plotLabelFlag = 1;
elseif strcmp(cfg.plotlabels, 'some')
    opt.plotLabelFlag = 2;
else
    opt.plotLabelFlag = 0;
end




h = figure;

if ~isempty(cfg.renderer)
    try
        if ~(strcmp(cfg.renderer,'painters') || strcmp(cfg.renderer,'zbuffer') || strcmp(cfg.renderer,'oppengl'))
            error(['cfg.renderer parameter must either be ''painters'' or ''zbuffer'' or ''opengl'' but given was ' cfg.renderer ''])
        end
        set(h, 'renderer', cfg.renderer);
    catch e
    end
end

if strcmp(cfg.bgcolor,'dark')
    whitebg(h,'k');
    opt.chancolors = 1-opt.chancolors;
end
set(h,'MenuBar','none');
axes('Parent',h,'Position',[0.035 0.09 0.965 0.91]);
set(h,'color',[0 0 0]);
set(gca, 'YColor', [0.5 0.5 0.5]);

set(gca,'TickDir','out');
set(gca,'TickLength',[0.005 0.01])
% a = gca;
% b = copyobj(a, h);
% set(b,  'YColor', [0.3 0.3 0.3], 'XTickLabel', [], 'YTickLabel', [])
set(gca,'Fontsize',5,'FontUnits','normalized');

if isfield(cfg,'begin_end_events') || isfield(cfg,'begin_end_events')
    cfg.displayEvents = 'yes';
else
    cfg.displayEvents = 'no';
end

cfg.browserversion = '3.0.0';

if strcmp(cfg.doSleepScoring,'yes')
    
    
    cfg.score_channel_eeg_scale = 50;
    temp_eeg_span = (1/min(cfg.chanscale))*cfg.score_channel_eeg_scale*cfg.chanscale(cfg.score_channel_eeg_number);
    cfg.ylim = [-temp_eeg_span temp_eeg_span];
    
    
    cfg.markSpindles = 'no';
    cfg.markSO = 'no';
    
    cfg.underlaySpindleSignal = 'no';
    cfg.underlayAlphaSignal = 'no';
    cfg.underlaySOSignal = 'no';

    
    
    
    cfg.spindle_mark_color = [0 1 0];
    cfg.slowoscillation_mark_color = [1 0 0];
    
    cfg.underlaySpindleSignal_color = [150/255 192/255 150/255];
    cfg.underlayAlphaSignal_color = [0/255 192/255 0/255];
    cfg.underlaySOSignal_color = [255/255 165/255 0/255];
    
    cfg.color_text_on_bg = [0.8 0.8 0.8];
    
    cfg.curr_displayed_detected_slowosci_perc_display_ind = 0;
    cfg.curr_displayed_detected_slowosci_perc = -1;
    cfg.curr_displayed_detected_slowosci_number = -1;
    cfg.curr_displayed_detected_slowosci_perc_cumulative = 0;
    cfg.curr_displayed_detected_slowosci_perc_substractive = 0;
    cfg.curr_displayed_detected_spindels_number = -1;
    cfg.curr_displayed_marking = -1;
    
    cfg.SOdetection_orientation = 1;
    
    cfg.sp_thresholdForDetectionBeginEnd = 10;
    cfg.sp_thresholdForDetectionCriterion = 15;
    cfg.sp_minSec = 0.5;
    cfg.sp_maxSec = 2;
    cfg.sp_minFreq = 12;
    cfg.sp_maxFreq = 14;
    cfg.al_minFreq = 8;
    cfg.al_maxFreq = 11;
    
    cfg.markSO_filter = 'yes';
    
    
    cfg.so_minFreq = 0.5;
    cfg.so_maxFreq = 2;
    
    cfg.so_filter_minFreq = cfg.so_minFreq;
    cfg.so_filter_maxFreq = cfg.so_maxFreq;
    
    cfg.so_thresholdAmplitudeForDetection = 75;
    
    cfg.emg_thresholdAmplitudeForDetection = 45;
    
    
    
    cfg.nEpochsBuffer = 0.34;
    
    cfg.FrqOfSmpl = opt.orgdata.fsample;
    
    cfg.UseFixedFilterOrder_bp = 'no';
    cfg.FilterOrder_bp = 100;
    cfg.AstopLeft_bp = 100;
    cfg.Apass_bp = 0.001;
    cfg.AstopRight_bp = 100;
    cfg.StopToPassTransitionWidth_bp = 1.25;
    cfg.PassToStopTransitionWidth_bp = 1.25;
    
    cfg.UseFixedFilterOrder_hp = 'yes';
    cfg.FilterOrder_hp = 4;
    cfg.AstopLeft_hp = 100;
    cfg.Apass_hp = 0.001;
    cfg.StopToPassTransitionWidth_hp = 0.2;
    
    
    cfg.UseFixedFilterOrder_lp = 'no';
    cfg.FilterOrder_lp = 100;
    cfg.AstopRight_lp = 100;
    cfg.Apass_lp = 0.001;
    cfg.PassToStopTransitionWidth_lp = 1.25;
    
    
    
    
    cfg = update_filters(cfg);
    
    cfg.markECG = 'no';
    cfg.ECG_signalMultiplicator = 1;
    cfg.ecg_peak_filter_minFreq = 20;
    cfg.ecg_peak_filter_maxFreq = 35;
    
    if cfg.has_ECG
        FpassLeft = cfg.ecg_peak_filter_minFreq; %left pass frequency in Hz
        FpassRight = cfg.ecg_peak_filter_maxFreq; %right pass frequency in Hz
        
        FstopLeft = FpassLeft - cfg.StopToPassTransitionWidth_bp; %left stop frequency in Hz
        FstopRight = FpassRight + cfg.PassToStopTransitionWidth_bp; %left stop frequency in Hz
        
        %usedFilterOrder_bp = NaN;
        %bp_hdm = NaN;
        if strcmp(cfg.core_cfg.bpfilttype,'IIRdesigned') || strcmp(cfg.core_cfg.bpfilttype,'FIRdesigned')
            bp_d = [];
            bp_hd = [];
            if strcmp(cfg.UseFixedFilterOrder_bp,'yes')
                bp_d = fdesign.bandpass('N,Fst1,Fp1,Fp2,Fst2',cfg.FilterOrder_bp,FstopLeft,FpassLeft,FpassRight,FstopRight,cfg.FrqOfSmpl);
                bp_hd = design(bp_d,'equiripple');
            else
                bp_d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',FstopLeft,FpassLeft,FpassRight,FstopRight,cfg.AstopLeft_bp,cfg.Apass_bp,cfg.AstopRight_bp,cfg.FrqOfSmpl);
                bp_hd = design(bp_d,'equiripple','MinOrder', 'even');
            end
            %usedFilterOrder_bp = bp_hd.order;
            cfg.hr_bpfilterdesign = bp_hd;
            %bp_hdm = measure(bp_hd);
        end
    end
    
    opt.markingstatus = 'off';
    opt.zoomstatus = 'off';
    cfg.use_ruler = 'no';
    
    %cfg.autosave_hypnogram = 'yes';
    cfg.autosave_hypnogram_every_number_change = 1;
    opt.autosave_hypnogram_change_interator = 0;
    
    cfg.skip_to_next = 'always'; %'always' 'firstscore' 'unknown' 'stay' 
    cfg.confidence_skip_to_lower_than_threshold = 0;
    
    cfg.display_power_spectrum = 'no';
    cfg.freq_borders = [0.5 4; 4 8; 8 12; 12 15; 20 30;];
    cfg.freq_colors = jet(size(cfg.freq_borders,1));
    
    cfg.display_time_frequency = 'no';
    
    
    cfg.artifact_export_delimiter = ',';
    
    cfg.toggle_epoch_marker = 0;
    
    
    if ~isfield(cfg,'hypn_mult')
        cfg.hypn_mult = {};
        cfg.hypn_plot_interpol_mult = {};
        cfg.hypn_plot_interpol_MA_mult = {};
        cfg.hypn_plot_interpol_confidence_mult = {};
        cfg.hypn_mult_idx = 1;
    end
    
end

setappdata(h, 'opt', opt);
setappdata(h, 'cfg', cfg);

% set the figure window title
funcname = mfilename();
if nargin < 2
    if isfield(cfg, 'dataset')
        dataname = cfg.dataset;
    elseif isfield(cfg, 'datafile')
        dataname = cfg.datafile;
    else
        dataname = [];
    end
else
    dataname = inputname(2);
end
hfig = gcf;
if ft_platform_supports('matlabversion',-Inf, '2014a')
    handlenum = hfig;
else
    handlenum = hfig.Number;
end
set(gcf, 'Name', sprintf('%d: %s: %s', handlenum, funcname, join_str(', ',dataname)));
set(gcf, 'NumberTitle', 'off');

% set zoom option to on
% zoom(h,'on')
% set(zoom(h),'actionPostCallback',@zoom_drawlabels_cb)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up the figure and callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(h, 'KeyPressFcn',           @keyboard_cb);
set(h, 'WindowButtonDownFcn',   {@ft_select_range, 'multiple', false, 'xrange', true, 'yrange', false, 'clear', true, 'contextmenu', cfg.selfun, 'callback', {@select_range_cb, h}, 'event', 'WindowButtonDownFcn'});
set(h, 'WindowButtonUpFcn',     {@ft_select_range, 'multiple', false, 'xrange', true, 'yrange', false, 'clear', true, 'contextmenu', cfg.selfun, 'callback', {@select_range_cb, h}, 'event', 'WindowButtonUpFcn'});
set(h, 'WindowButtonMotionFcn', {@ft_select_range, 'multiple', false, 'xrange', true, 'yrange', false, 'clear', true, 'contextmenu', cfg.selfun, 'callback', {@select_range_cb, h}, 'event', 'WindowButtonMotionFcn'});

% make the user interface elements for the data view
% uicontrol('tag', 'labels',  'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', opt.trialviewtype, 'userdata', 't')
% uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', 'leftarrow')
% uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', 'rightarrow')
%
% uicontrol('tag', 'labels',  'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'channel','userdata', 'c')
% uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', 'uparrow')
% uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', 'downarrow')
%
% uicontrol('tag', 'labels',  'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'horizontal', 'userdata', 'h')
% uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '-', 'userdata', 'shift+leftarrow')
% uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '+', 'userdata', 'shift+rightarrow')
%
% uicontrol('tag', 'labels',  'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'vertical', 'userdata', 'v')
% uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '-', 'userdata', 'shift+downarrow')
% uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '+', 'userdata', 'shift+uparrow')

temp_lower_line_y = 0.00;
temp_lower_line_y2 = temp_lower_line_y+0.04;

uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'x', 'position', [0.0, 0.99 , 0.01, 0.01],'backgroundcolor',[0.5 0.5 0.5],'foregroundcolor',[1 1 1], 'userdata', 'x')

uicontrol('tag', 'labels',  'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', opt.trialviewtype, 'position', [0.05, temp_lower_line_y2 , 0.04, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 't')
uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'position', [0.05, temp_lower_line_y , 0.02, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'leftarrow')
uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'position', [0.07, temp_lower_line_y , 0.02, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'rightarrow')

uicontrol('tag', 'labels',  'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'channel', 'position', [0.09, temp_lower_line_y2 , 0.04, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1],'userdata', 'shift+c')
uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'position', [0.09, temp_lower_line_y , 0.02, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'shift+uparrow')
uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'position', [0.11, temp_lower_line_y , 0.02, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'shift+downarrow')


uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '+', 'position', [0.13, temp_lower_line_y + 0.08/3 + 0.08/3, 0.04, 0.08/3],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'uparrow')
uicontrol('tag', 'labels',  'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'scale', 'position', [0.13, temp_lower_line_y + 0.08/3, 0.04, 0.08/3],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'y')
uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '-', 'position', [0.13, temp_lower_line_y, 0.04, 0.03],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'downarrow')

if ~strcmp(cfg.doSleepScoring,'yes')
    uicontrol('tag', 'labels',  'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'horizontal', 'position', [0.17,temp_lower_line_y2 , 0.04, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'h')
    uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '-', 'position', [0.17, temp_lower_line_y , 0.02, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'shift+leftarrow')
    uicontrol('tag', 'buttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '+', 'position', [0.17, temp_lower_line_y , 0.02, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'shift+rightarrow')
else
    uicontrol('tag', 'scindlabelSOperc', 'parent', h, 'units', 'normalized', 'style', 'text', 'string','SO%: ?','Fontsize',10,'FontUnits','normalized','FontWeight','bold', 'ForegroundColor',[1 0 0],'HorizontalAlignment','left','position', [0.2, temp_lower_line_y+0.08/3+0.08/3 , 0.15, 0.08/3],'backgroundcolor',[0 0 0])
    %uicontrol('tag', 'scindlabelSOnum', 'parent', h, 'units', 'normalized', 'style', 'text', 'string', '#SO: ?','Fontsize',10,'FontUnits','normalized','FontWeight','bold', 'ForegroundColor',[1 0 0],'HorizontalAlignment','left','position', [0.2, temp_lower_line_y+0.04 , 0.15, 0.02],'backgroundcolor',[0 0 0])
    uicontrol('tag', 'scindlabelSpnum', 'parent', h, 'units', 'normalized', 'style', 'text', 'string', '#Sp: ?','Fontsize',10,'FontUnits','normalized','FontWeight','bold', 'ForegroundColor',[0 1 0],'HorizontalAlignment','left','position', [0.2, temp_lower_line_y+0.08/3 , 0.15, 0.08/3],'backgroundcolor',[0 0 0])
    uicontrol('tag', 'scindlabelMarkperc', 'parent', h, 'units', 'normalized', 'style', 'text', 'string', 'Mark%: ?','Fontsize',10,'FontUnits','normalized','FontWeight','bold', 'ForegroundColor',[1 0 0],'HorizontalAlignment','left','position', [0.2, temp_lower_line_y , 0.15, 0.08/3],'backgroundcolor',[0 0 0])
end

if strcmp(cfg.doSleepScoring,'yes')
    
    uicontrol('tag', 'scbuttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'W(0/W)','position', [0.36, temp_lower_line_y2 , 0.04, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', '0')
    uicontrol('tag', 'scbuttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'REM(5/R)','position', [0.40, temp_lower_line_y2 , 0.04, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', '5')
    if strcmp(cfg.standard,'rk')
        uicontrol('tag', 'scbuttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'MT(8)','position', [0.44, temp_lower_line_y2, 0.04, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', '8')
    end
    uicontrol('tag', 'scbuttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'MA(9)','position', [0.48, temp_lower_line_y2 , 0.04, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', '9')
    
    uicontrol('tag', 'scbuttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'S1(1)','position', [0.36, temp_lower_line_y , 0.04, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', '1')
    uicontrol('tag', 'scbuttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'S2(2)','position', [0.40 temp_lower_line_y , 0.04, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', '2')
    uicontrol('tag', 'scbuttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'S3(3)','position', [0.44, temp_lower_line_y , 0.04, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', '3')
    if strcmp(cfg.standard,'rk')
        uicontrol('tag', 'scbuttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'S4(4)','position', [0.48, temp_lower_line_y , 0.04, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', '4')
    end
    uicontrol('tag', 'scoptbuttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '?(7/D)','position', [0.52, temp_lower_line_y2 , 0.04, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'd')
    uicontrol('tag', 'scoptbuttons_nextunk', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '[>]','position', [0.52, temp_lower_line_y , 0.04, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'u')

    uicontrol('tag', 'scoptbuttons_SOdet', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '(+)_.�\_.�\_.�','Fontsize',6,'FontUnits','normalized','position', [0.56, temp_lower_line_y+0.08/3+0.08/3 , 0.05, 0.08/3],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'k')
    uicontrol('tag', 'scoptbuttons_SOdisp', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '--(_.�\)--','Fontsize',6,'FontUnits','normalized','position', [0.56, temp_lower_line_y+0.08/3 , 0.05, 0.08/3],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'o')
    uicontrol('tag', 'scoptbuttons_SPdet', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '--~~-~~--','Fontsize',6,'FontUnits','normalized', 'position', [0.61, temp_lower_line_y+0.08/3+0.08/3 , 0.05, 0.08/3],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'j')
    uicontrol('tag', 'scoptbuttons_SPdisp', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '--(~~~)--','Fontsize',6,'FontUnits','normalized', 'position', [0.61, temp_lower_line_y+0.08/3 , 0.05, 0.08/3],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'n')
    if cfg.has_ECG
        uicontrol('tag', 'scoptbuttons_HRdisp', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<3 <3 <3','Fontsize',6,'FontUnits','normalized', 'position', [0.56, temp_lower_line_y , 0.05, 0.08/3],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'b')
    end
    uicontrol('tag', 'scoptbuttons_ALdisp', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '--wWwWw--','Fontsize',6,'FontUnits','normalized', 'position', [0.61, temp_lower_line_y , 0.05, 0.08/3],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'a')
    
    uicontrol('tag', 'scoptbuttons_pow', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'pow','position',  [0.675, temp_lower_line_y+0.08/3+0.08/3, 0.025, 0.08/3],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'p')
    uicontrol('tag', 'scoptbuttons_tfr', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'tfr','position', [0.675, temp_lower_line_y+0.08/3, 0.025, 0.08/3],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'f')
    uicontrol('tag', 'scoptbuttons_min', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'min','position', [0.675, temp_lower_line_y, 0.025, 0.08/3],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'shift+m')

    uicontrol('tag', 'scoptbuttons_mark', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'mark','position', [0.7, temp_lower_line_y+0.08/3+0.08/3, 0.025, 0.08/3],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'm')
    uicontrol('tag', 'scoptbuttons_zoom', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'zoom','position', [0.7, temp_lower_line_y+0.08/3, 0.025, 0.08/3],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'z')
    uicontrol('tag', 'scoptbuttons_grid', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'grid','position', [0.725, temp_lower_line_y+0.08/3+0.08/3, 0.025, 0.08/3],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'g')
    uicontrol('tag', 'scoptbuttons_ruler', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'ruler','position', [0.725, temp_lower_line_y+0.08/3, 0.025, 0.08/3],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'e')
    if strcmp(cfg.displayEvents,'yes')
        uicontrol('tag', 'scoptbuttons_EvDisp', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'events','position', [0.7, temp_lower_line_y , 0.05, 0.08/3],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'v')
    else
        uicontrol('tag', 'scoptbuttons_EvDisp', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'EVENTS','position', [0.7, temp_lower_line_y , 0.05, 0.08/3],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'v')
    end
    uicontrol('tag', 'scoptbuttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'save','position', [0.75, temp_lower_line_y2 , 0.05, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'shift+s')
    uicontrol('tag', 'scoptbuttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'open','position', [0.75, temp_lower_line_y , 0.05, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'shift+o')
    
    uicontrol('tag', 'scoptbuttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'export hyp','position', [0.8, temp_lower_line_y2 , 0.05, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'shift+e')
    uicontrol('tag', 'scoptbuttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'import hyp','position', [0.8, temp_lower_line_y , 0.05, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'shift+i')
    
    %uicontrol('tag', 'scoptbuttons_focusEEG', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', ['Focus EEG: ' opt.hdr.label{cfg.score_channel_eeg_number}],'position', [0.85, temp_lower_line_y+0.08/3+0.08/3 , 0.08, 0.08/3], 'userdata', 'alt+e','backgroundcolor', cfg.score_channel_eeg_color)
    %uicontrol('tag', 'scoptbuttons_focusEOG', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', ['Focus EOG: ' opt.hdr.label{cfg.score_channel_eog_number}],'position', [0.85, temp_lower_line_y+0.08/3 , 0.08, 0.08/3], 'userdata', 'alt+o','backgroundcolor', cfg.score_channel_eog_color)
    %uicontrol('tag', 'scoptbuttons_focusEMG', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', ['Focus EMG: ' opt.hdr.label{cfg.score_channel_emg_number}],'position', [0.85, temp_lower_line_y , 0.08, 0.08/3], 'userdata', 'alt+m','backgroundcolor', cfg.score_channel_emg_color)
    
    
    
    
    
    uicontrol('tag', 'scoptbuttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'shortcuts','position', [0.9, temp_lower_line_y2 , 0.05, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'shift+h')
    
    
    uicontrol('tag', 'scoptbuttons', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'thresholds','position', [0.95, temp_lower_line_y2 , 0.04, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1], 'userdata', 'l')
end

% legend artifacts/features
for iArt = 1:length(artlabel)
    %   uicontrol('tag', 'artifactui', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', artlabel{iArt}, 'userdata', num2str(iArt), 'position', [0.91, 0.9 - ((iArt-1)*0.09), 0.08, 0.04], 'backgroundcolor', opt.artcolors(iArt,:))
    %   uicontrol('tag', 'artifactui', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', ['shift+' num2str(iArt)], 'position', [0.91, 0.855 - ((iArt-1)*0.09), 0.03, 0.04], 'backgroundcolor', opt.artcolors(iArt,:))
    %   uicontrol('tag', 'artifactui', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', ['control+' num2str(iArt)], 'position', [0.96, 0.855 - ((iArt-1)*0.09), 0.03, 0.04], 'backgroundcolor', opt.artcolors(iArt,:))
    uicontrol('tag', 'artifactui_button', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', ['artifact(' opt.artdata.label{opt.ftsel} ')'], 'userdata', 'shift+a', 'position', [0.01, temp_lower_line_y2 - ((iArt-1)*0.09), 0.04, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1])
    uicontrol('tag', 'artifactui', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', ['shift+' num2str(iArt)], 'position', [0.01, temp_lower_line_y - ((iArt-1)*0.09), 0.02, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1])
    uicontrol('tag', 'artifactui', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', ['control+' num2str(iArt)], 'position', [0.03, temp_lower_line_y - ((iArt-1)*0.09), 0.02, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1])
    
end


if strcmp(cfg.viewmode, 'butterfly')
    % button to find label of nearest channel to datapoint
    uicontrol('tag', 'artifactui_button', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'identify', 'userdata', 'i', 'position', [0.91, temp_lower_line_y2, 0.08, 0.05], 'backgroundcolor', [1 1 1])
end

% 'edit preproc'-button
%uicontrol('tag', 'preproccfg', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string','preproc cfg','position', [0.91, 0.55 - ((iArt-1)*0.09), 0.08, 0.04],'callback',@preproc_cfg1_cb)
uicontrol('tag', 'preproccfg', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string','cfg','position', [0.95, temp_lower_line_y , 0.04, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1],'callback',@preproc_cfg1_cb)



%ft_uilayout(h, 'tag', 'labels',  'width', 0.10, 'height', 0.05);
%ft_uilayout(h, 'tag', 'buttons', 'width', 0.05, 'height', 0.05);


ft_uilayout(h, 'tag', 'labels',     'style', 'pushbutton', 'callback', @keyboard_cb);
ft_uilayout(h, 'tag', 'buttons',    'style', 'pushbutton', 'callback', @keyboard_cb);
ft_uilayout(h, 'tag', 'artifactui', 'style', 'pushbutton', 'callback', @keyboard_cb);
ft_uilayout(h, 'tag', 'artifactui_button', 'style', 'pushbutton', 'callback', @keyboard_cb);

if strcmp(cfg.doSleepScoring,'yes');
    ft_uilayout(h, 'tag', 'scbuttons', 'style', 'pushbutton', 'callback', @keyboard_cb);
    ft_uilayout(h, 'tag', 'scoptbuttons', 'style', 'pushbutton', 'callback', @keyboard_cb);
    
    ft_uilayout(h, 'tag', 'scoptbuttons_SOdet', 'style', 'pushbutton', 'callback', @keyboard_cb);
    ft_uilayout(h, 'tag', 'scoptbuttons_SOdisp', 'style', 'pushbutton', 'callback', @keyboard_cb);
    ft_uilayout(h, 'tag', 'scoptbuttons_SPdet', 'style', 'pushbutton', 'callback', @keyboard_cb);
    ft_uilayout(h, 'tag', 'scoptbuttons_SPdisp', 'style', 'pushbutton', 'callback', @keyboard_cb);
    if cfg.has_ECG
        ft_uilayout(h, 'tag', 'scoptbuttons_HRdisp', 'style', 'pushbutton', 'callback', @keyboard_cb);
    end
    ft_uilayout(h, 'tag', 'scoptbuttons_ALdisp', 'style', 'pushbutton', 'callback', @keyboard_cb);
    ft_uilayout(h, 'tag', 'scoptbuttons_EvDisp', 'style', 'pushbutton', 'callback', @keyboard_cb);
    ft_uilayout(h, 'tag', 'scoptbuttons_mark', 'style', 'pushbutton', 'callback', @keyboard_cb);
    ft_uilayout(h, 'tag', 'scoptbuttons_zoom', 'style', 'pushbutton', 'callback', @keyboard_cb);
    ft_uilayout(h, 'tag', 'scoptbuttons_grid', 'style', 'pushbutton', 'callback', @keyboard_cb);
    ft_uilayout(h, 'tag', 'scoptbuttons_nextunk', 'style', 'pushbutton', 'callback', @keyboard_cb);
    ft_uilayout(h, 'tag', 'scoptbuttons_min', 'style', 'pushbutton', 'callback', @keyboard_cb);
    ft_uilayout(h, 'tag', 'scoptbuttons_ruler', 'style', 'pushbutton', 'callback', @keyboard_cb);
    ft_uilayout(h, 'tag', 'scoptbuttons_pow', 'style', 'pushbutton', 'callback', @keyboard_cb);
    ft_uilayout(h, 'tag', 'scoptbuttons_tfr', 'style', 'pushbutton', 'callback', @keyboard_cb);
    %ft_uilayout(h, 'tag', 'scoptbuttons_focusEEG', 'style', 'pushbutton', 'callback', @keyboard_cb);
    %ft_uilayout(h, 'tag', 'scoptbuttons_focusEOG', 'style', 'pushbutton', 'callback', @keyboard_cb);
    %ft_uilayout(h, 'tag', 'scoptbuttons_focusEMG', 'style', 'pushbutton', 'callback', @keyboard_cb);
end

% ft_uilayout(h, 'tag', 'labels',  'retag', 'viewui');
% ft_uilayout(h, 'tag', 'buttons', 'retag', 'viewui2');
% if strcmp(cfg.doSleepScoring,'yes');
%     ft_uilayout(h, 'tag', 'scbuttons', 'retag', 'viewui');
%     ft_uilayout(h, 'tag', 'scbuttons2', 'retag', 'viewui');
%
% end
% ft_uilayout(h, 'tag', 'viewui', 'BackgroundColor', [0.8 0.8 0.8], 'hpos', 'auto', 'vpos', 0.1);
% ft_uilayout(h, 'tag', 'viewui2', 'BackgroundColor', [0.8 0.8 0.8], 'hpos', 'auto', 'vpos', 0);




definetrial_cb(h);

if cfg.StartWithOpenSession
    try
        [h,opt,cfg,hyp_file_filterindex] = loadSession(h,opt,cfg);
        setappdata(h, 'opt', opt);
        setappdata(h, 'cfg', cfg);
        %redraw_cb(h);
        if hyp_file_filterindex ~= 0
            msgbox('Opening session successful!' ,'Opening successful','modal');
        end
    catch err
        msgbox('Open the session failed!' ,'Open failed','error','modal');
    end
end

redraw_cb(h);




% %% Scrollbar
%
% % set initial scrollbar value
% dx = maxtime;
%
% % set scrollbar position
% fig_pos=get(gca,'position');
% scroll_pos=[fig_pos(1) fig_pos(2) fig_pos(3) 0.02];
%
% % define callback
% S=['set(gca,''xlim'',get(gcbo,''value'')+[ ' num2str(mintime) ',' num2str(maxtime) '])'];
%
% % Creating Uicontrol
% s=uicontrol('style','slider',...
%     'units','normalized','position',scroll_pos,...
%     'callback',S,'min',0,'max',0, ...
%     'visible', 'off'); %'value', xmin

% set initial scrollbar value
% dx = maxtime;
%
% % set scrollbar position
% fig_pos=get(gca,'position');
% scroll_pos=[fig_pos(1) fig_pos(2) fig_pos(3) 0.02];
%
% % define callback
% S=['set(gca,''xlim'',get(gcbo,''value'')+[ ' num2str(mintime) ',' num2str(maxtime) '])'];
%
% % Creating Uicontrol
% s=uicontrol('style','slider',...
%     'units','normalized','position',scroll_pos,...
%     'callback',S,'min',0,'max',0, ...
%     'visible', 'off'); %'value', xmin
%initialize postion of plot
% set(gca,'xlim',[xmin xmin+dx]);





if nargout
    % wait until the user interface is closed, get the user data with the updated artifact details
    set(h, 'CloseRequestFcn', @cleanup_cb);
    
    while ishandle(h)
        uiwait(h);
        opt = getappdata(h, 'opt');
        if opt.cleanup
            delete(h);
        end
    end
    
    % add the updated artifact definitions to the output cfg
    for i=1:length(opt.artdata.label)
        cfg.artfctdef.(opt.artdata.label{i}).artifact = convert_event(opt.artdata.trial{1}(i,:), 'artifact');
    end
    
    % add the updated preproc to the output
    try
        browsecfg = getappdata(h, 'cfg');
        cfg.preproc = browsecfg.preproc;
    end
    
    % add the update event to the output cfg
    cfg.event = opt.event;
    
    % do the general cleanup and bookkeeping at the end of the function
    %if (~isdeployed)
    ft_postamble debug
    ft_postamble trackconfig
    ft_postamble provenance
    ft_postamble previous data
    %end
    
end % if nargout








% 
% 
% epochLength = str2num(getParam('epochLength',listOfCoreParameters)); % in seconds
% 
% 
% 
% ApplyEventmappingSettings = getParam('ApplyEventmappingSettings',listOfParameters);%either yes or no
% if (useDummyDataset)
%     ApplyEventmappingSettings = 'no';
% end
% EventmappingSettingsDefinitionsFileName = getParam('EventmappingSettingsDefinitionsFileName',listOfParameters);
% listOfEventmappingSettingsFiles = {};
% if strcmp(ApplyEventmappingSettings,'yes')
%     
%     if strcmp(ReadSingleDataset,'yes')
%         if exist([EventmappingSettingsDefinitionsFileName],'file') ~= 2
%             error(['EventmappingSettingsDefinitionsFileName file ' [EventmappingSettingsDefinitionsFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
%         end
%         listOfEventmappingSettingsFiles = {EventmappingSettingsDefinitionsFileName};
%     else
%         if exist([pathInputFolder filesep EventmappingSettingsDefinitionsFileName],'file') ~= 2
%             error(['EventmappingSettingsDefinitionsFileName file ' [pathInputFolder filesep EventmappingSettingsDefinitionsFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
%         end
%         listOfEventmappingSettingsFiles = read_mixed_csv([pathInputFolder filesep EventmappingSettingsDefinitionsFileName],',');
%     end
%     if ~(all(size(listOfDatasetsPaths) == size(listOfEventmappingSettingsFiles)))
%         error('files or number of Datasetspaths listOfEventmappingSettingsFiles are invalid or do not aggree')
%     end
%     
%     for iDefFiles = 1:size(listOfEventmappingSettingsFiles,1)
%         if exist([pathInputFolder filesep listOfEventmappingSettingsFiles{iDefFiles}],'file') ~= 2
%             error(['The channel event mapping settings file listed in EventmappingSettingsDefinitionsFileName for dataset number ' num2str(iDefFiles)  ' does not exist'])
%         end
%         
%         try
%             readtable([pathInputFolder filesep listOfEventmappingSettingsFiles{iDefFiles}],'Delimiter',',');
%         catch err
%             error(['The channel event mapping settings file listed in EventmappingSettingsDefinitionsFileName for dataset number ' num2str(iDefFiles)  ' is not readable'])
%         end
%     end
%     
% end
% 
% ApplyEventsSelection = getParam('ApplyEventsSelection',listOfParameters);%either yes or no
% if (useDummyDataset)
%     ApplyEventsSelection = 'no';
% end
% if strcmp(ApplyEventsSelection,'yes')
%     
%     EventsTarget1FilePathsFileName = getParam('EventsTarget1FilePathsFileName',listOfParameters);
%     
%     
%     
%     
%     EventsTarget1TimePointColumn = getParam('EventsTarget1TimePointColumn',listOfParameters);
%     
%     EventsTarget1CompareColumns = strsplit(getParam('EventsTarget1CompareColumns',listOfParameters),' ');
%     
%     EventTarget1TimeWindowOffsetTime = str2num(getParam('EventTarget1TimeWindowOffsetTime',listOfParameters)); % in units of EventsTarget1TimePointColumn
%     UseSecondColumnAndBothOffsets = getParam('UseSecondColumnAndBothOffsets',listOfParameters);
%     EventsTarget1TimePointColumn2 = getParam('EventsTarget1TimePointColumn2',listOfParameters);
%     EventTarget1TimeWindowOffsetTime2 = str2num(getParam('EventTarget1TimeWindowOffsetTime2',listOfParameters)); % in units of EventsTarget1TimePointColumn
%     EventTarget1TimeWindowPreOffsetTime = str2num(getParam('EventTarget1TimeWindowPreOffsetTime',listOfParameters)); % in units of EventsTarget1TimePointColumn
%     EventTarget1TimeWindowPostOffsetTime = str2num(getParam('EventTarget1TimeWindowPostOffsetTime',listOfParameters)); % in units of EventsTarget1TimePointColumn
%     
%     
%     EventsFilesWhich = getParam('EventsFilesWhich',listOfParameters);%Event files to be processed either all or subset if subset then DataSetsNumbers is used for selection default all
%     EventsFilesNumbers = str2num(getParam('EventsFilesNumbers',listOfParameters));%The line numbers of the events file to be processed if EventsFilesWhich parameter is set to subset
%     
%     if strcmp(ReadSingleDataset,'yes')
%         if exist([EventsTarget1FilePathsFileName],'file') ~= 2
%             error(['EventsTarget1FilePathsFileName file ' [EventsTarget1FilePathsFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
%         end
%         listOfEventsTarget1Paths = {EventsTarget1FilePathsFileName};
%     else
%         if exist([pathInputFolder filesep EventsTarget1FilePathsFileName],'file') ~= 2
%             error(['EventsTarget1FilePathsFileName file ' [pathInputFolder filesep EventsTarget1FilePathsFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
%         end
%         listOfEventsTarget1Paths = read_mixed_csv([pathInputFolder filesep EventsTarget1FilePathsFileName],',');
%     end
%     iDatas_Events = 1:(length(listOfEventsTarget1Paths));
%     
%     if strcmp(EventsFilesWhich,'subset')
%         if ~(ismember(min(EventsFilesNumbers),iDatas) && ismember(max(EventsFilesNumbers),iDatas))
%             error('Parameter EventsFilesNumbers contains numbers not matching to any line number, e.g. too less EventsTestFilePaths in EventsTestFilePathsFileName!')
%         end
%         iDatas_Events = EventsFilesNumbers;
%     end
%     
%     
%     FilterValuesSplitString = ' ';
%     EventsTarget1FilterForColumn = getParam('EventsTarget1FilterForColumn',listOfParameters);%variable name for event value in test events files to apply text filter to if nothing is entered it is not filtered e.g. channel default is no value entered
%     EventsTarget1FilterValues = strsplit(getParam('EventsTarget1FilterValues',listOfParameters),FilterValuesSplitString);%variable values for EventsTarget1FilterForColum in test events files to apply text filter e.g. Cz
%     
% end
% 
% 
% ApplyEventsSelection2 = getParam('ApplyEventsSelection2',listOfParameters);%either yes or no
% if (useDummyDataset)
%     ApplyEventsSelection2 = 'no';
% end
% if strcmp(ApplyEventsSelection2,'yes')
%     
%     EventsTarget2FilePathsFileName = getParam('EventsTarget2FilePathsFileName',listOfParameters);
%     
%     
%     
%     
%     EventsTarget2TimePointColumn = getParam('EventsTarget2TimePointColumn',listOfParameters);
%     
%     EventsTarget2CompareColumns = strsplit(getParam('EventsTarget2CompareColumns',listOfParameters),' ');
%     
%     EventTarget2TimeWindowOffsetTime = str2num(getParam('EventTarget2TimeWindowOffsetTime',listOfParameters)); % in units of EventsTarget2TimePointColumn
%     UseSecondColumnAndBothOffsets2 = getParam('UseSecondColumnAndBothOffsets2',listOfParameters);
%     EventsTarget2TimePointColumn2 = getParam('EventsTarget2TimePointColumn2',listOfParameters);
%     EventTarget2TimeWindowOffsetTime2 = str2num(getParam('EventTarget2TimeWindowOffsetTime2',listOfParameters)); % in units of EventsTarget2TimePointColumn
%     EventTarget2TimeWindowPreOffsetTime = str2num(getParam('EventTarget2TimeWindowPreOffsetTime',listOfParameters)); % in units of EventsTarget2TimePointColumn
%     EventTarget2TimeWindowPostOffsetTime = str2num(getParam('EventTarget2TimeWindowPostOffsetTime',listOfParameters)); % in units of EventsTarget2TimePointColumn
%     
%     
%     EventsFilesWhich2 = getParam('EventsFilesWhich2',listOfParameters);%Event files to be processed either all or subset if subset then DataSetsNumbers is used for selection default all
%     EventsFilesNumbers2 = str2num(getParam('EventsFilesNumbers2',listOfParameters));%The line numbers of the events file to be processed if EventsFilesWhich parameter is set to subset
%     
%     if strcmp(ReadSingleDataset,'yes')
%         if exist([EventsTarget2FilePathsFileName],'file') ~= 2
%             error(['EventsTarget2FilePathsFileName file ' [EventsTarget2FilePathsFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
%         end
%         listOfEventsTarget2Paths = {EventsTarget2FilePathsFileName};
%     else
%         if exist([pathInputFolder filesep EventsTarget2FilePathsFileName],'file') ~= 2
%             error(['EventsTarget2FilePathsFileName file ' [pathInputFolder filesep EventsTarget2FilePathsFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
%         end
%         listOfEventsTarget2Paths = read_mixed_csv([pathInputFolder filesep EventsTarget2FilePathsFileName],',');
%     end
%     iDatas_Events = 1:(length(listOfEventsTarget2Paths));
%     
%     if strcmp(EventsFilesWhich,'subset')
%         if ~(ismember(min(EventsFilesNumbers),iDatas) && ismember(max(EventsFilesNumbers),iDatas))
%             error('Parameter EventsFilesNumbers contains numbers not matching to any line number, e.g. too less EventsTestFilePaths in EventsTestFilePathsFileName!')
%         end
%         iDatas_Events = EventsFilesNumbers;
%     end
%     
%     
%     FilterValuesSplitString = ' ';
%     EventsTarget2FilterForColumn = getParam('EventsTarget2FilterForColumn',listOfParameters);%variable name for event value in test events files to apply text filter to if nothing is entered it is not filtered e.g. channel default is no value entered
%     EventsTarget2FilterValues = strsplit(getParam('EventsTarget2FilterValues',listOfParameters),FilterValuesSplitString);%variable values for EventsTarget2FilterForColum in test events files to apply text filter e.g. Cz
%     
% end


    %     for second read in of hypnogramm after downsampling
    %     if (signalOffsetSamples ~= 0)
    %         signalOffsetSamples_downsampled = floor(signalOffsetSeconds*FrqOfSmpl);
    %         roiBegins = roiBegins + signalOffsetSamples_downsampled;
    %         roiEnds = roiEnds + signalOffsetSamples_downsampled;
    %     end
    
    
    
    
    
    
fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cleanup_cb(h, eventdata)
selection = questdlg('Close Browser?',...
    'Close Request',...
    'Yes','No','No');
switch selection,
    case 'Yes',
        
        cfg = getappdata(h, 'cfg');
        if isfield(cfg,'hhyp')
            if ishandle(cfg.hhyp)
                close(cfg.hhyp);
            end
        end
        if isfield(cfg,'f_ps')
            if ishandle(cfg.f_ps)
                close(cfg.f_ps);
            end
        end
        if isfield(cfg,'f_tfr')
            if ishandle(cfg.f_tfr)
                close(cfg.f_tfr);
            end
        end
        opt = getappdata(h, 'opt');
        opt.cleanup = true;
        setappdata(h, 'opt', opt);
        uiresume
        opt = [];
        cfg = [];
        h = [];
    case 'No'
        return
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function definetrial_cb(h, eventdata)
opt = getappdata(h, 'opt');
cfg = getappdata(h, 'cfg');
if strcmp(cfg.continuous, 'no')
    
    % when zooming in, lock the trial! one can only go to the next trial when horizontal scaling doesn't segment the data - from ft-meeting: this might be relaxed later on - roevdmei
    if isempty(opt.trllock)
        opt.trllock = opt.trlop;
    end
    locktrllen = ((opt.trlorg(opt.trllock,2)-opt.trlorg(opt.trllock,1)+1) ./ opt.fsample);
    % if cfg.blocksize is close to the length of the locked trial, set it to that
    if (abs(locktrllen-cfg.blocksize) / locktrllen) < 0.1
        cfg.blocksize = locktrllen;
    end
    
    %%%%%%%%%
    % trial is locked, change subdivision of trial
    if cfg.blocksize < locktrllen
        % lock the trial if it wasn't locked (and thus trlop refers to the actual trial)
        if isempty(opt.trllock)
            opt.trllock = trlop;
        end
        % save current position if already
        if isfield(opt, 'trlvis')
            thissegbeg = opt.trlvis(opt.trlop,1);
        end
        datbegsample = min(opt.trlorg(opt.trllock,1));
        datendsample = max(opt.trlorg(opt.trllock,2));
        smpperseg  = round(opt.fsample * cfg.blocksize);
        begsamples = datbegsample:smpperseg:datendsample;
        endsamples = datbegsample+smpperseg-1:smpperseg:datendsample;
        offset     = (((1:numel(begsamples))-1)*smpperseg) + opt.trlorg(opt.trllock,3);
        if numel(endsamples)<numel(begsamples)
            endsamples(end+1) = datendsample;
        end
        trlvis = [];
        trlvis(:,1) = begsamples';
        trlvis(:,2) = endsamples';
        trlvis(:,3) = offset;
        % determine length of each trial, and determine the offset with the current requested zoom-level
        trllen   = (trlvis(:,2) - trlvis(:,1)+1);
        sizediff = smpperseg - trllen;
        opt.nanpaddata = sizediff;
        
        if isfield(opt, 'trlvis')
            % update the current trial counter and try to keep the current sample the same
            opt.trlop   = nearest(begsamples, thissegbeg);
        end
        % update trialviewtype
        opt.trialviewtype = 'trialsegment';
        % update button
        set(findobj(get(h,'children'),'string','trial'),'string','segment');
        %%%%%%%%%
        
        
        %%%%%%%%%
        % trial is not locked, go to original trial division and zoom out
    elseif cfg.blocksize >= locktrllen
        trlvis = opt.trlorg;
        % set current trlop to locked trial if it was locked before
        if ~isempty(opt.trllock)
            opt.trlop = opt.trllock;
        end
        smpperseg  = round(opt.fsample * cfg.blocksize);
        % determine length of each trial, and determine the offset with the current requested zoom-level
        trllen   = (trlvis(:,2) - trlvis(:,1)+1);
        sizediff = smpperseg - trllen;
        opt.nanpaddata = sizediff;
        
        % update trialviewtype
        opt.trialviewtype = 'trial';
        % update button
        set(findobj(get(h,'children'),'string','trialepoch'),'string',opt.trialviewtype);
        
        % release trial lock
        opt.trllock = [];
        %%%%%%%%%
    end
    
    % save trlvis
    opt.trlvis  = trlvis;
else
    % construct a trial definition for visualisation
    if isfield(opt, 'trlvis') % if present, remember where we were
        thistrlbeg = opt.trlvis(opt.trlop,1);
    end
    % look at cfg.blocksize and make opt.trl accordingly
    datbegsample = min(opt.trlorg(:,1));
    datendsample = max(opt.trlorg(:,2));
    smpperseg  = round(opt.fsample * cfg.blocksize);
    begsamples = datbegsample:smpperseg:datendsample;
    endsamples = datbegsample+smpperseg-1:smpperseg:datendsample;
    if numel(endsamples)<numel(begsamples)
        endsamples(end+1) = datendsample;
    end
    trlvis = [];
    trlvis(:,1) = begsamples';
    trlvis(:,2) = endsamples';
    % compute the offset. In case if opt.trlorg has multiple trials, the first sample is t=0, otherwise use the offset in opt.trlorg
    if size(opt.trlorg,1)==1
        offset = begsamples - repmat(begsamples(1),[1 numel(begsamples)]); % offset for all epochs compared to the first
        offset = offset + opt.trlorg(1,3);
        trlvis(:,3) = offset;
    else
        offset = begsamples - repmat(begsamples(1),[1 numel(begsamples)]);
        trlvis(:,3) = offset;
    end
    
    if isfield(opt, 'trlvis')
        % update the current trial counter and try to keep the current sample the same
        % opt.trlop   = nearest(round((begsamples+endsamples)/2), thissample);
        opt.trlop   = nearest(begsamples, thistrlbeg);
    end
    opt.trlvis  = trlvis;
    
    % NaN-padding when horizontal scaling is bigger than the data
    % two possible situations, 1) zoomed out so far that all data is one epoch, or 2) multiple epochs but last epoch is smaller than the rest
    sizediff = smpperseg-(endsamples-begsamples+1);
    opt.nanpaddata = sizediff;
end % if continuous
setappdata(h, 'opt', opt);
setappdata(h, 'cfg', cfg);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function helptext = help_text()
helptext = [ ...
    '\nKEYBOARD SHORTCUTS ONLY WILL WORK \n' ...
    'IF MOUSE WAS CLICKED IN DISPLAY AREA!\n\n' ...
    'Keyboard shortcuts for... \n' ...
    ' Scoring:\n'...
    '  0 or W: Wake (W) \n'...
    '  1: Stage 1 (S1/N1) \n'...
    '  2: Stage 2 (S2/N2) \n'...
    '  3: Stage 3 (S3/N3) \n'...
    '  4: Stage 4 (S4) \n'...
    '  R: REM (R)\n'...
    '  7 or D: Delete Stage \n'...
    '  8: Movement Time (MT) \n'...
    '  9 or W: Movement Arousal/Epoch with artifact(s) (MA, additional) \n'...
    ' Control:\n'...
    '  Left-arrow: go to previous epoch \n'...
    '  Right-arrow: go to next epoch \n'...
    '  Up-arrow: decrease scaling (zoom in all channels) \n'...
    '  Down-arrow: increase scaling (zoom out all channels) \n'...
    '  Y: Vertical Y-scaling basis (equivalent to scaling factor = 1) (dialog) \n'...
    '  Shift + C: Channel settings \n' ...
    '             select, color, scaling, order, foci (EEG, EMG, EOG) \n'...
    '  Shift + Up-arrow: skip up through undisplayed channels \n'...
    '  Shift + Down-arrow: skip down through undisplayed channels \n'...
    '  Shift + Q: Quit Program \n'...
    '  U: (automatic) epoch skipping modus change \n'...
    '  P: Power spectrum of the current epoch \n'...
    '  F: Time-frequency plot of the current EEG scoring channel \n'...
    '  X: Hide/Display menue bar \n'...
    '  T: Select epoch to jump to (dialog) \n'...
    '  Shift + T and Ctrl + T: toggle and un-toggle epoch start-marker\n'...
    '  L: Set Thresholds (dialog) \n'...
    '  Shift + H: Shortcut help \n'...
    ' Scoring aids:\n'...
    '  M: enable/disable marking, cummulative time \n'...
    '  N: enable/disable spindle signal display (EEG filtered in spindle band) \n'...
    '  V: enable/disable display of pre-readin Events, if processed already \n'...
    '  Q: delete previously made marking \n'...
    '  E: enable/disable Ruler (aka "the Score-ship") \n'...
    '  J: enable/disable spindle marking \n'...
    '  K: enable/[polartity-switch]/disable K-Complex \n'...
    '     or Slow oscillation marking \n'...
    '  B: enable/[polartity-switch]/disable heart rate \n'...
    '     indication and coloring \n'...
    '  G: enable/disable display of the time grid \n'...
    '  Z: enable/disable zooming (double-click = return to full view) \n'...
    '     left-click: zoom in \n'...
    '     Alt + left-click zoom out \n'...
    '     double-click: fully zoom out \n'...
    '  Shift + A: switch to next artifact marking type \n'...
    '  Ctrl + 1..9 OR Alt + 1..9: skip to next artifact type 1..9 \n'...
    ' Saving/Export:\n'...
    '  Shift + S: Save Session (dialog) \n'...
    '  Shift + O: Open Session (dialog) \n'...
    '  Shift + I: Import Hypnogram (dialog) \n'...
    '  Shift + E: Export Hypnogram (dialog) \n'...
    ' Non-Scoring:\n'...
    '  Shift + Left-arrow: decrease epoch length\n'...
    '  Shift + Right-arrow: increase epoch length\n'...
    '  H: Horizontal scaling \n'...
    '  (S: switch highlighting/marking style, deprecated)\n'...
    ];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function help_cb(h, eventdata)
fprintf('------------------------------------------------------------------------------------\n')
fprintf(help_text());
fprintf('------------------------------------------------------------------------------------\n')
fprintf('\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_range_cb(h, range, cmenulab) %range 1X4 in sec relative to current trial
opt = getappdata(h, 'opt');
cfg = getappdata(h, 'cfg');

% the range should be in the displayed box
range(1) = max(opt.hpos-opt.width/2, range(1));
range(2) = max(opt.hpos-opt.width/2, range(2));
range(1) = min(opt.hpos+opt.width/2, range(1));
range(2) = min(opt.hpos+opt.width/2, range(2));
range = (range-(opt.hpos-opt.width/2)) / opt.width; % left side of the box becomes 0, right side becomes 1
range = range * (opt.hlim(2) - opt.hlim(1)) + opt.hlim(1);   % 0 becomes hlim(1), 1 becomes hlim(2)

begsample = opt.trlvis(opt.trlop,1);
endsample = opt.trlvis(opt.trlop,2);
offset    = opt.trlvis(opt.trlop,3);

% determine the selection
begsel = round(range(1)*opt.fsample+begsample-offset-1);
endsel = round(range(2)*opt.fsample+begsample-offset);
% artifact selection is now always based on begsample/endsample/offset
% -roevdmei

% the selection should always be confined to the current trial
begsel = max(begsample, begsel);
endsel = min(endsample, endsel);

% mark or execute selfun
if isempty(cmenulab)
    % the left button was clicked INSIDE a selected range, update the artifact definition or event
    
    if strcmp(cfg.selectmode, 'markartifact')
        % mark or unmark artifacts
        artval = opt.artdata.trial{1}(opt.ftsel, begsel:endsel);
        artval = any(artval,1);
        if any(artval)
            fprintf('there is overlap with the active artifact (%s), disabling this artifact\n',opt.artdata.label{opt.ftsel});
            opt.artdata.trial{1}(opt.ftsel, begsel:endsel) = 0;
        else
            fprintf('there is no overlap with the active artifact (%s), marking this as a new artifact\n',opt.artdata.label{opt.ftsel});
            opt.artdata.trial{1}(opt.ftsel, begsel:endsel) = 1;
        end
        
        % redraw only when marking (so the focus doesn't go back to the databrowser after calling selfuns
        setappdata(h, 'opt', opt);
        setappdata(h, 'cfg', cfg);
        redraw_cb(h);
    elseif strcmp(cfg.selectmode, 'markpeakevent') || strcmp(cfg.selectmode, 'marktroughevent')
        %mark or unmark events, marking at peak/trough of window
        if any(intersect(begsel:endsel, [opt.event.sample]))
            fprintf('there is overlap with one or more event(s), disabling this/these event(s)\n');
            ind_rem = intersect(begsel:endsel, [opt.event.sample]);
            for iRemove = 1:length(ind_rem)
                opt.event([opt.event.sample]==ind_rem(iRemove)) = [];
            end
        else
            fprintf('there is no overlap with any event, adding an event to the peak/trough value\n');
            % check if only 1 chan, other wise not clear max in which channel. %
            % ingnie: would be cool to add the option to select the channel when multiple channels
            if size(opt.curdat.trial{1},1) > 1
                error('cfg.selectmode = ''markpeakevent'' and ''marktroughevent'' only supported with 1 channel in the data')
            end
            if strcmp(cfg.selectmode, 'markpeakevent')
                [dum ind_minmax] = max(opt.curdat.trial{1}(begsel-begsample+1:endsel-begsample+1));
                val = 'peak';
            elseif strcmp(cfg.selectmode, 'marktroughevent')
                [dum ind_minmax] = min(opt.curdat.trial{1}(begsel-begsample+1:endsel-begsample+1));
                val = 'trough';
            end
            samp_minmax = begsel + ind_minmax - 1;
            event_new.type     = 'ft_databrowser_manual';
            event_new.sample   = samp_minmax;
            event_new.value    = val;
            event_new.duration = 1;
            event_new.offset   = 0;
            % add new event to end opt.event
            % check if events are in order now
            if  min(diff([opt.event.sample]))>0
                % add new event in line with old ones
                nearest_event = nearest([opt.event.sample], samp_minmax);
                if opt.event(nearest_event).sample > samp_minmax
                    %place new event before nearest
                    ind_event_new = nearest_event;
                else
                    %place new event after nearest
                    ind_event_new = nearest_event +1;
                end
                event_lastpart = opt.event(ind_event_new:end);
                opt.event(ind_event_new) = event_new;
                opt.event(ind_event_new+1:end+1) = event_lastpart;
            else
                %just add to end
                opt.event(end+1) = event_new;
            end
            clear event_new ind_event_new event_lastpart val dum ind_minmax
        end
        % redraw only when marking (so the focus doesn't go back to the databrowser after calling selfuns
        setappdata(h, 'opt', opt);
        setappdata(h, 'cfg', cfg);
        redraw_cb(h);
    end
    
else
    % the right button was used to activate the context menu and the user made a selection from that menu
    % execute the corresponding function
    
    % get index into cfgs
    selfunind = strcmp(cfg.selfun, cmenulab);
    
    % cut out the requested data segment
    seldata.label    = opt.curdat.label;
    seldata.time{1}  = offset2time(offset+begsel-begsample, opt.fsample, endsel-begsel+1);
    seldata.trial{1} = ft_fetch_data(opt.curdat, 'begsample', begsel, 'endsample', endsel);
    seldata.fsample  = opt.fsample;
    seldata.cfg.trl  = [begsel endsel offset];
    
    % prepare input
    funhandle = ft_getuserfun(cmenulab,'browse');
    funcfg = cfg.selcfg{selfunind};
    % get windowname and give as input (can be used for the other functions as well, not implemented yet)
    if ~strcmp(opt.trialviewtype,'trialsegment')
        str = sprintf('%s %d/%d, time from %g to %g min', opt.trialviewtype, opt.trlop, size(opt.trlvis,1), seldata.time{1}(1)./60, seldata.time{1}(end)./60);
    else
        str = sprintf('trial %d/%d: epoch: %d/%d , time from %g to %g min', opt.trllock, size(opt.trlorg,1), opt.trlop, size(opt.trlvis,1), seldata.time{1}(1)./60, seldata.time{1}(end)./60);
    end
    funcfg.figurename = [cmenulab ': ' str];
    feval(funhandle, funcfg, seldata);
end

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function preproc_cfg1_cb(h,eventdata)
parent = get(h,'parent');
cfg = getappdata(parent, 'cfg');

% parse cfg.preproc
if ~isempty(cfg.preproc)
    tmpcfg = cfg.preproc;
    cfg = [];
    cfg.preproc = tmpcfg;
    code = printstruct('cfg', cfg);
else
    code = [];
end

% add descriptive lines
sep     = sprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
descrip = sprintf('%% Add/change cfg options for on-the-fly preprocessing\n%% Use as cfg.preproc.xxx\n');
code    = [sep descrip sep code];


% make figure displaying the edit box
pph = figure;
axis off
% add save button
uicontrol('tag', 'preproccfg_l2', 'parent', pph, 'units', 'normalized', 'style', 'pushbutton', 'string','save and close','position', [0.81, 0.6 , 0.18, 0.10],'callback',@preproc_cfg2_cb);

% add edit box
ppeh = uicontrol('style', 'edit');
set(pph, 'toolBar', 'none')
set(pph, 'menuBar', 'none')
set(pph, 'Name', 'cfg.preproc editor')
set(pph, 'NumberTitle', 'off')
set(ppeh, 'Units', 'normalized');
set(ppeh, 'Position', [0 0 .8 1]);
set(ppeh, 'backgroundColor', [1 1 1]);
set(ppeh, 'horizontalAlign', 'left');
set(ppeh, 'max', 2);
set(ppeh, 'min', 0);
set(ppeh, 'FontName', 'Courier');
set(ppeh, 'FontSize', 12);
set(ppeh, 'string', code);


% add handle for the edit style to figure
setappdata(pph,'superparent', parent); % superparent is the main ft_databrowser window
setappdata(pph,'ppeh', ppeh);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function preproc_cfg2_cb(h,eventdata)
parent = get(h,'parent');
superparent = getappdata(parent,'superparent');
ppeh   = getappdata(parent,'ppeh');
code = get(ppeh, 'string');

% remove descriptive lines (so they don't display on command line)
code = cellstr(code(5:end,:));
% get rid of empty lines and white space
remind = [];
for iline = 1:numel(code)
    code{iline} = strtrim(code{iline});
    if isempty(code{iline})
        remind = [remind iline];
    end
end
code(remind) = [];

if ~isempty(code)
    ispreproccfg = strncmp(code,'cfg.preproc.',12);
    if ~all(ispreproccfg)
        errordlg('cfg-options must be specified as cfg.preproc.xxx','cfg.preproc editor','modal')
    end
    % eval the code
    for icomm = 1:numel(code)
        eval([code{icomm} ';']);
    end
    
    % check for cfg and output into the original appdata-window
    if ~exist('cfg','var')
        cfg = [];
        cfg.preproc = [];
    end
    maincfg = getappdata(superparent,'cfg');
    maincfg.preproc = cfg.preproc;
    setappdata(superparent,'cfg',maincfg)
end

close(parent)
redraw_cb(superparent)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function keyboard_cb(h, eventdata)

if (isempty(eventdata) && ft_platform_supports('matlabversion',-Inf, '2014a')) || isa(eventdata, 'matlab.ui.eventdata.ActionData')
  % determine the key that corresponds to the uicontrol element that was activated
  key = get(h, 'userdata');
else
  % determine the key that was pressed on the keyboard
  key = parsekeyboardevent(eventdata);
end

% get focus back to figure
if ~strcmp(get(h, 'type'), 'figure')
    set(h, 'enable', 'off');
    drawnow;
    set(h, 'enable', 'on');
end

h = getparent(h);
opt = getappdata(h, 'opt');
cfg = getappdata(h, 'cfg');

if strcmp(cfg.standard,'rk')
    keySet = {'0' '1' '2' '3' '4' '5' '7' '8' '9' 'd' 'r' 'w'};
else
    keySet = {'0' '1' '2' '3' '5' '7' '9' 'd' 'r' 'w'};
end
                        
switch key
    case 'g'
        if strcmp(cfg.drawgrid,'yes')
            cfg.drawgrid = 'no';
        else
            cfg.drawgrid = 'yes';
        end
        setappdata(h, 'opt', opt);
        setappdata(h, 'cfg', cfg);
        redraw_cb(h, eventdata);
    case 'shift+m'
        if strcmp(cfg.drawgrid,'yes')
            cfg.drawgrid = 'no';
            %         else
            %             cfg.drawgrid = 'yes';
        end
        if strcmp(cfg.plot_stage_signatures,'yes')
            cfg.plot_stage_signatures = 'no';
        else
            cfg.plot_stage_signatures = 'yes';
        end
        if strcmp(cfg.highlight_scoring_channels,'yes')
            cfg.highlight_scoring_channels = 'no';
        else
            cfg.highlight_scoring_channels = 'yes';
        end

        setappdata(h, 'opt', opt);
        setappdata(h, 'cfg', cfg);
        redraw_cb(h, eventdata);
    case 'z'
        if ~isfield(opt,'zoomstatus')
            opt.zoomstatus = 'on';
        else
            switch opt.zoomstatus
                case 'on'
                    opt.zoomstatus = 'off';
                    %case 'out'
                    %    opt.zoomstatus = 'off';
                case 'off'
                    opt.zoomstatus = 'on';
                    if isfield(opt,'markingstatus')
                        opt.markingstatus = 'off';
                    end
                otherwise
                    opt.zoomstatus = 'off';
            end
        end
        
        %zoom(h,'off')
        % oldkeyhook = get(h, 'KeyPressFcn');
        %     newkeyhook = @keyboard_cb;
        zoom(h,opt.zoomstatus)
        %     oldWindowButtonDownFcnhook = get(h, 'WindowButtonDownFcn');
        %     oldWindowButtonUpFcnhook = get(h, 'WindowButtonUpFcn');
        %     hManager = uigetmodemanager(h);
        %     set(hManager.WindowListenerHandles,'Enable','off');
        %     set(h, 'KeyPressFcn', newkeyhook);
        %     set(h, 'WindowButtonDownFcn', oldWindowButtonDownFcnhook);
        %     set(h, 'WindowButtonUpFcn', oldWindowButtonUpFcnhook);
        setappdata(h, 'opt', opt);
    case 'm'
        %marking
        if ~isfield(opt,'markingstatus')
            opt.markingstatus = 'on';
            opt.marks = [];
            opt.markSecondClick = 0;
            if isfield(opt,'zoomstatus')
                opt.zoomstatus = 'off';
                zoom(h,opt.zoomstatus)
            end
        else
            switch opt.markingstatus
                case 'on'
                    opt.markingstatus = 'off';
                    %cfg.curr_displayed_marking = -1;
                    %opt.marks = [];
                case 'off'
                    opt.markingstatus = 'on';
                    %opt.marks = [];
                    %opt.markSecondClick = 0;
                    if isfield(opt,'zoomstatus')
                        opt.zoomstatus = 'off';
                        zoom(h,opt.zoomstatus)
                    end
                otherwise
                    opt.markingstatus = 'off';
                    cfg.curr_displayed_marking = -1;
                    opt.marks = [];
            end
            if ~isfield(opt,'marks')
                opt.marks = [];
            end
            opt.markSecondClick = 0;
        end
        
        %zoom(h,'off')
        setappdata(h, 'opt', opt);
        setappdata(h, 'cfg', cfg);
        updateLabels(h);
    case 'q'
        %delete last marking
        if isfield(opt,'marks')
            if ~isempty(opt.marks)
                opt.marks(end,:) = [];
            end
            setappdata(h, 'opt', opt);
            redraw_marks_cb(h);
        end
        if isfield(opt,'markSecondClick')
            if opt.markSecondClick == 1
                opt.markSecondClick = 0;
            end
        end
        setappdata(h, 'opt', opt);
    case 'shift+a'
        % switch to another artifact type
        opt.ftsel = opt.ftsel + 1;
        if opt.ftsel > numel(opt.artdata.label)
            opt.ftsel = 1;
        end
        %opt.ftsel = str2double(strrep(key, 'shift+', ''));
        numart = size(opt.artdata.trial{1}, 1);
        if opt.ftsel > numart
            fprintf('data has no artifact type %i \n', opt.ftsel)
        else
            setappdata(h, 'opt', opt);
            setappdata(h, 'cfg', cfg);
            fprintf('switching to the "%s" artifact\n', opt.artdata.label{opt.ftsel});
            redraw_cb(h, eventdata);
        end
        %FW begin
    case {'shift+1' 'shift+2' 'shift+3' 'shift+4' 'shift+5' 'shift+6' 'shift+7' 'shift+8' 'shift+9'}
        % go to previous artifact
        %opt.ftsel = str2double(key(end));
        numart = size(opt.artdata.trial{1}, 1);
        if opt.ftsel > numart
            fprintf('data has no artifact type %i \n', opt.ftsel)
        else
            cursam = opt.trlvis(opt.trlop,1);
            artsam = find(opt.artdata.trial{1}(opt.ftsel,1:cursam-1), 1, 'last');
            if isempty(artsam)
                fprintf('no earlier "%s" artifact found\n', opt.artdata.label{opt.ftsel});
            else
                fprintf('going to previous "%s" artifact\n', opt.artdata.label{opt.ftsel});
                if opt.trlvis(nearest(opt.trlvis(:,1),artsam),1) < artsam
                    arttrl = nearest(opt.trlvis(:,1),artsam);
                else
                    arttrl = nearest(opt.trlvis(:,1),artsam)-1;
                end
                opt.trlop = arttrl;
                setappdata(h, 'opt', opt);
                setappdata(h, 'cfg', cfg);
                redraw_cb(h, eventdata);
            end
        end
    case keySet
        %profile on
        if strcmp(cfg.doSleepScoring,'yes')
            switch key
	        case 'r'
	            key = '5';
	        case 'w'
	            key = '0';
            end
            
            curr_epoch = opt.trlop;
            
            temp_epochLengthSamples = opt.trlvis(1, 2) - opt.trlvis(1, 1) + 1;
            nEpochs = floor(size(opt.orgdata.trial{1},2)/temp_epochLengthSamples);
            if curr_epoch > nEpochs
                curr_epoch = curr_epoch - 1;
            end
            curr_hypn = cfg.hypn(curr_epoch,:);
            
            switch key
                case {'0' '1' '2' '3' '4' '5' '8'}
                    h1 = str2num(key);
                    h2 = curr_hypn(:,2);
                case '9'
                    h1 = curr_hypn(:,1);
                    h2 = 1;
                case {'d' '7'}
                    h1 = -1;
                    h2 = 0;
            end
            if h1 == 8
                h2 = 2;
            end
            
            switch cfg.skip_to_next
                case 'always'
                    temp_skip_to_next_epoch = (h1 ~= -1);
                case 'firstscore'
                    temp_skip_to_next_epoch = ((cfg.hypn(curr_epoch,1) == -1) && (h1 ~= -1));
                case 'unknown'
                    % see code below
                case 'stay'
                    temp_skip_to_next_epoch = false;
            end
            
            [stagestring h1_str h2_str] = getStageStringByHypnValue(h1,h2);
            opt.curr_stage = stagestring;
            
            next_low_conf_epoch = [];
            if cfg.confidence_skip_to_lower_than_threshold ~= 0
                if size(cfg.hypn,2) > 2
                    if curr_epoch < size(cfg.hypn,1)
                        next_low_conf_epoch = find(cfg.hypn((curr_epoch+1):end,3) < cfg.confidence_skip_to_lower_than_threshold,1,'first');
                        if ~isempty(next_low_conf_epoch)
                            next_low_conf_epoch = curr_epoch + next_low_conf_epoch;
                            temp_skip_to_next_epoch = true;
                        end
                    end
                end
            end
            
            next_unknown_epoch = [];
            if strcmp(cfg.skip_to_next, 'unknown')
                    if curr_epoch < size(cfg.hypn,1)
                        next_unknown_epoch = find(cfg.hypn((curr_epoch+1):end,1) == -1,1,'first');
                        if ~isempty(next_unknown_epoch)
                            next_unknown_epoch = curr_epoch + next_unknown_epoch;
                            temp_skip_to_next_epoch = (h1 ~= -1);
                        end
                    end
            end

            
            
            opt.prev_stages = getPrevStageString_stage(cfg.hypn,curr_epoch,6);
            opt.next_stages = getNextStageString_stage(cfg.hypn,curr_epoch,6);
            
            
            
            if isfield(cfg, 'hypn_mult')
                if ~isempty(cfg.hypn_mult)
                    for iHypMult = 1:numel(cfg.hypn_mult)
                        opt.curr_stage_mult{iHypMult} = getStageStringByHypnValue(cfg.hypn_mult{iHypMult}(curr_epoch,1),cfg.hypn_mult{iHypMult}(curr_epoch,2));
                        %                 opt.prev_stages_mult{iHypMult}
                        %                 opt.next_stages_mult{iHypMult}
                        %                 opt.curr_stage_mult = '?';
                    end
                end
            end
            
            
            
            
            if (cfg.toggle_epoch_marker ~= 0) && (cfg.toggle_epoch_marker <= curr_epoch)
                for iTempEpoch = cfg.toggle_epoch_marker:curr_epoch
                    if size(cfg.hypn,2) > 2
                        cfg.hypn(iTempEpoch,:) = [h1 h2 1];
                    else
                        cfg.hypn(iTempEpoch,:) = [h1 h2];
                    end
                    hyp_begsample = cfg.hyp_epochLengthSamples*(cfg.toggle_epoch_marker-1)+1; % opt.trlvis(opt.trlop,1);
                    
                end
                temp_hyp_part = cfg.hypn(cfg.toggle_epoch_marker:curr_epoch,:);
            else
                hyp_begsample = cfg.hyp_epochLengthSamples*(curr_epoch-1)+1; % opt.trlvis(opt.trlop,1);
                if size(cfg.hypn,2) > 2
                    temp_hyp_part = [h1 h2 1];
                else
                    temp_hyp_part = [h1 h2];
                end
                cfg.hypn(curr_epoch,:) = temp_hyp_part;
                
            end
            
            
            
            
            hyp_begsample_this_epoch = cfg.hyp_epochLengthSamples*(curr_epoch-1)+1; % opt.trlvis(opt.trlop,1);
            hyp_endsample = cfg.hyp_epochLengthSamples*(curr_epoch);% opt.trlvis(opt.trlop,2);
            [curr_ep_hypn_plot_interpol, curr_ep_hypn_plot_interpol_MA] = interpolate_hypn_for_plot(temp_hyp_part,length(hyp_begsample_this_epoch:hyp_endsample),cfg.plot_MA_offset,istrue(cfg.yaxdisteqi));
            
            cfg.toggle_epoch_marker = 0;
            
            
            cfg.hypn_plot_interpol(hyp_begsample:hyp_endsample) = curr_ep_hypn_plot_interpol;
            cfg.hypn_plot_interpol_MA(hyp_begsample:hyp_endsample) = curr_ep_hypn_plot_interpol_MA;
            
            if ~isempty(cfg.hypn_plot_interpol_confidence)
                [temp_dummy, curr_ep_hypn_plot_interpol_confidence] = interpolate_hypn_for_plot(temp_hyp_part(:,2:3),length(hyp_begsample_this_epoch:hyp_endsample),cfg.plot_confidence_offset,istrue(cfg.yaxdisteqi));
                cfg.hypn_plot_interpol_confidence(hyp_begsample:hyp_endsample) = curr_ep_hypn_plot_interpol_confidence;
            end
            
            if temp_skip_to_next_epoch
                next_epoch = opt.trlop + 1;
                if ~isempty(next_low_conf_epoch)
                    next_epoch = next_low_conf_epoch;
                end
                if ~isempty(next_unknown_epoch)
                    if ~isempty(next_low_conf_epoch)
                    	next_epoch = min(next_unknown_epoch,next_low_conf_epoch);
                    else
                        next_epoch = next_unknown_epoch;
                    end
                end
                opt.trlop = min(next_epoch, size(opt.trlvis,1)); % should not be larger than the number of trials
                if isfield(opt,'marks')
                    opt.marks = [];
                end
                cfg.curr_displayed_marking = -1;
            end
            
            %if strcmp(cfg.autosave_hypnogram,'yes')
            if (cfg.autosave_hypnogram_every_number_change ~= 0)
                
                if isfield(cfg,'autosave_hypfilepath')
                    opt.autosave_hypnogram_change_interator = opt.autosave_hypnogram_change_interator + 1;
                    
                    if mod(opt.autosave_hypnogram_change_interator, cfg.autosave_hypnogram_every_number_change) == 0
                        opt.autosave_hypnogram_change_interator = 0;
                        
                        try
                            if ~isfield(cfg,'hypnogram_delimiter_autosave')
                                cfg.hypnogram_delimiter_autosave = '\t';
                            end
                            writeHypnogramFile(cfg.autosave_hypfilepath,cfg.hypn,cfg.hypnogram_delimiter_autosave);
                            
                            temp_ArtifactPath = [cfg.autosave_hypfilepath '.artifacts.csv'];
                            writeArtifactFile(temp_ArtifactPath,opt,cfg.artifact_export_delimiter);
                            
                        catch err
                            
                        end
                    end
                else
                    
                end
            end
            
            
            setappdata(h, 'opt', opt);
            setappdata(h, 'cfg', cfg);
            redraw_cb(h, eventdata);
            
        end
        
        %profile viewer
        % FW end
    case {'control+1' 'control+2' 'control+3' 'control+4' 'control+5' 'control+6' 'control+7' 'control+8' 'control+9' 'alt+1' 'alt+2' 'alt+3' 'alt+4' 'alt+5' 'alt+6' 'alt+7' 'alt+8' 'alt+9'}
        % go to next artifact
        %opt.ftsel = str2double(key(end));
        numart = size(opt.artdata.trial{1}, 1);
        if opt.ftsel > numart
            fprintf('data has no artifact type %i \n', opt.ftsel)
        else
            cursam = opt.trlvis(opt.trlop,2);
            artsam = find(opt.artdata.trial{1}(opt.ftsel,cursam+1:end), 1, 'first') + cursam;
            if isempty(artsam)
                fprintf('no later "%s" artifact found\n', opt.artdata.label{opt.ftsel});
            else
                fprintf('going to next "%s" artifact\n', opt.artdata.label{opt.ftsel});
                if opt.trlvis(nearest(opt.trlvis(:,1),artsam),1) < artsam
                    arttrl = nearest(opt.trlvis(:,1),artsam);
                else
                    arttrl = nearest(opt.trlvis(:,1),artsam)-1;
                end
                opt.trlop = arttrl;
                setappdata(h, 'opt', opt);
                setappdata(h, 'cfg', cfg);
                redraw_cb(h, eventdata);
            end
        end
    case 'leftarrow'
        opt.trlop = max(opt.trlop - 1, 1); % should not be smaller than 1
        if isfield(opt,'marks')
            opt.marks = [];
        end
        cfg.curr_displayed_marking = -1;
        setappdata(h, 'opt', opt);
        setappdata(h, 'cfg', cfg);
        redraw_cb(h, eventdata);
    case 'rightarrow'
        opt.trlop = min(opt.trlop + 1, size(opt.trlvis,1)); % should not be larger than the number of trials
        if isfield(opt,'marks')
            opt.marks = [];
        end
        cfg.curr_displayed_marking = -1;
        setappdata(h, 'opt', opt);
        setappdata(h, 'cfg', cfg);
        redraw_cb(h, eventdata);
    case 'shift+uparrow'
        chansel = match_str(opt.hdr.label, cfg.channel);
        minchan = min(chansel);
        numchan = length(chansel);
        chansel = minchan - numchan : minchan - 1;
        if min(chansel)<1
            chansel = chansel - min(chansel) + 1;
        end
        % convert numeric array into cell-array with channel labels
        cfg.channel = opt.hdr.label(chansel);
        setappdata(h, 'opt', opt);
        setappdata(h, 'cfg', cfg);
        delete(findobj(h,'tag', 'chanlabel'));  % remove channel labels here, and not in redrawing to save significant execution time (see bug 2065)
        redraw_cb(h, eventdata);
    case 'shift+downarrow'
        chansel = match_str(opt.hdr.label, cfg.channel);
        maxchan = max(chansel);
        numchan = length(chansel);
        chansel = maxchan + 1 : maxchan + numchan;
        if max(chansel)>length(opt.hdr.label)
            chansel = chansel - (max(chansel) - length(opt.hdr.label));
        end
        % convert numeric array into cell-array with channel labels
        cfg.channel = opt.hdr.label(chansel);
        setappdata(h, 'opt', opt);
        setappdata(h, 'cfg', cfg);
        delete(findobj(h,'tag', 'chanlabel'));  % remove channel labels here, and not in redrawing to save significant execution time (see bug 2065)
        redraw_cb(h, eventdata);
    case 'shift+leftarrow'
        if ~strcmp(cfg.doSleepScoring,'yes')
            cfg.blocksize = cfg.blocksize*sqrt(2);
            setappdata(h, 'opt', opt);
            setappdata(h, 'cfg', cfg);
            definetrial_cb(h, eventdata);
            redraw_cb(h, eventdata);
        end
    case 'shift+rightarrow'
        if ~strcmp(cfg.doSleepScoring,'yes')
            cfg.blocksize = cfg.blocksize/sqrt(2);
            setappdata(h, 'opt', opt);
            setappdata(h, 'cfg', cfg);
            definetrial_cb(h, eventdata);
            redraw_cb(h, eventdata);
        end
    case 'uparrow'
        cfg.ylim = cfg.ylim/sqrt(2);
        setappdata(h, 'opt', opt);
        setappdata(h, 'cfg', cfg);
        redraw_cb(h, eventdata);
    case 'downarrow'
        cfg.ylim = cfg.ylim*sqrt(2);
        setappdata(h, 'opt', opt);
        setappdata(h, 'cfg', cfg);
        redraw_cb(h, eventdata);
    case 'shift+q'
        setappdata(h, 'opt', opt);
        setappdata(h, 'cfg', cfg);
        cleanup_cb(h);
    case 'x'
        status = get(h,'MenuBar');
        if strcmp(status,'figure')
            set(h,'MenuBar','none');
        else
            set(h,'MenuBar','figure');
        end
        redraw_cb(h, eventdata);
    case 'u'
        if strcmp(cfg.doSleepScoring,'yes')
            if ~isfield(cfg,'skip_to_next')
                cfg.skip_to_next = 'always';
            else
                if strcmp(cfg.skip_to_next, 'always')
                    cfg.skip_to_next= 'firstscore';
                elseif strcmp(cfg.skip_to_next, 'firstscore')
                    cfg.skip_to_next= 'unknown';
                elseif strcmp(cfg.skip_to_next, 'unknown')
                    cfg.skip_to_next= 'stay';
                elseif strcmp(cfg.skip_to_next, 'stay')
                    cfg.skip_to_next= 'always';
                end
            end
            setappdata(h, 'opt', opt);
            setappdata(h, 'cfg', cfg);
        end
    case 'e'
        if strcmp(cfg.doSleepScoring,'yes')
            if ~isfield(cfg,'use_ruler')
                cfg.use_ruler = 'yes';
            else
                if strcmp(cfg.use_ruler,'yes')
                    cfg.use_ruler = 'no';
                else
                    cfg.use_ruler = 'yes';
                end
            end
            setappdata(h, 'opt', opt);
            setappdata(h, 'cfg', cfg);
        end
    case 'k'
        if strcmp(cfg.doSleepScoring,'yes')
            if ~isfield(cfg,'markSO')
                cfg.markSO = 'yes';
                cfg.SOdetection_orientation = 1;
            else
                if strcmp(cfg.markSO,'yes')
                    cfg.markSO = 'no';
                    cfg.curr_displayed_detected_slowosci_perc = -1;
                    cfg.curr_displayed_detected_slowosci_number = -1;
                    cfg.curr_displayed_detected_slowosci_perc_cumulative = 0;
                    cfg.curr_displayed_detected_slowosci_perc_substractive = 0;
                else
                    cfg.markSO = 'yes';
                    if cfg.SOdetection_orientation == 1
                        cfg.SOdetection_orientation = -1;
                    else
                        cfg.SOdetection_orientation = 1;
                    end
                end
            end
            setappdata(h, 'opt', opt);
            setappdata(h, 'cfg', cfg);
            redraw_cb(h, eventdata);
        end
    case 'j'
        if strcmp(cfg.doSleepScoring,'yes')
            if ~isfield(cfg,'markSpindles')
                cfg.markSpindles = 'yes';
            else
                if strcmp(cfg.markSpindles,'yes')
                    cfg.markSpindles = 'no';
                    cfg.curr_displayed_detected_spindels_number = -1;
                else
                    cfg.markSpindles = 'yes';
                end
            end
            setappdata(h, 'opt', opt);
            setappdata(h, 'cfg', cfg);
            redraw_cb(h, eventdata);
        end
    case 'o'
        if strcmp(cfg.doSleepScoring,'yes')
            if ~isfield(cfg,'underlaySOSignal')
                cfg.underlaySOSignal = 'yes';
            else
                if strcmp(cfg.underlaySOSignal,'yes')
                    cfg.underlaySOSignal = 'no';
                else
                    cfg.underlaySOSignal = 'yes';
                end
            end
            setappdata(h, 'opt', opt);
            setappdata(h, 'cfg', cfg);
            redraw_cb(h, eventdata);
        end
    case 'n'
        if strcmp(cfg.doSleepScoring,'yes')
            if ~isfield(cfg,'underlaySpindleSignal')
                cfg.underlaySpindleSignal = 'yes';
            else
                if strcmp(cfg.underlaySpindleSignal,'yes')
                    cfg.underlaySpindleSignal = 'no';
                else
                    cfg.underlaySpindleSignal = 'yes';
                end
            end
            setappdata(h, 'opt', opt);
            setappdata(h, 'cfg', cfg);
            redraw_cb(h, eventdata);
        end
    case 'a'
        if strcmp(cfg.doSleepScoring,'yes')
            if ~isfield(cfg,'underlayAlphaSignal')
                cfg.underlayAlphaSignal = 'yes';
            else
                if strcmp(cfg.underlayAlphaSignal,'yes')
                    cfg.underlayAlphaSignal = 'no';
                else
                    cfg.underlayAlphaSignal = 'yes';
                end
            end
            setappdata(h, 'opt', opt);
            setappdata(h, 'cfg', cfg);
            redraw_cb(h, eventdata);
        end
    case 'v'
        if strcmp(cfg.doSleepScoring,'yes')
            if ~isfield(cfg,'displayEvents')
                cfg.displayEvents = 'yes';
            else
                if strcmp(cfg.displayEvents,'yes')
                    cfg.displayEvents = 'no';
                else
                    cfg.displayEvents = 'yes';
                end
            end
            setappdata(h, 'opt', opt);
            setappdata(h, 'cfg', cfg);
            redraw_cb(h, eventdata);
        end
        
    case 'p'
        if strcmp(cfg.doSleepScoring,'yes')
            if ~isfield(cfg,'display_power_spectrum')
                cfg.display_power_spectrum = 'yes';
            else
                if strcmp(cfg.display_power_spectrum,'yes')
                    cfg.display_power_spectrum = 'no';
                else
                    cfg.display_power_spectrum = 'yes';
                end
            end
            setappdata(h, 'opt', opt);
            setappdata(h, 'cfg', cfg);
            redraw_cb(h, eventdata);
        end
    case 'f'
        if strcmp(cfg.doSleepScoring,'yes')
            if ~isfield(cfg,'display_time_frequency')
                cfg.display_time_frequency = 'yes';
            else
                if strcmp(cfg.display_time_frequency,'yes')
                    cfg.display_time_frequency = 'no';
                else
                    cfg.display_time_frequency = 'yes';
                end
            end
            setappdata(h, 'opt', opt);
            setappdata(h, 'cfg', cfg);
            redraw_cb(h, eventdata);
        end
    case 'b'
        % toggle heart (b)eat detection
        if cfg.has_ECG
            if strcmp(cfg.markECG,'yes') && (cfg.ECG_signalMultiplicator == 1)
                cfg.ECG_signalMultiplicator = -1;
            elseif strcmp(cfg.markECG,'yes') && (cfg.ECG_signalMultiplicator == -1)
                cfg.markECG = 'no';
                cfg.ECG_signalMultiplicator = 1;
            else
                cfg.markECG = 'yes';
                cfg.ECG_signalMultiplicator = 1;
            end
            setappdata(h, 'opt', opt);
            setappdata(h, 'cfg', cfg);
            redraw_cb(h, eventdata);
        end
    case 't'
        % select the trial to display
        if ~strcmp(opt.trialviewtype,'trialsegment')
            str = sprintf('%s to display (current epoch = %d/%d)', opt.trialviewtype, opt.trlop, size(opt.trlvis,1));
        else
            str = sprintf('epoch to display (current epoch = %d/%d)', opt.trlop, size(opt.trlvis,1));
        end
        response = inputdlg(str, 'specify', 1, {num2str(opt.trlop)});
        if ~isempty(response)
            opt.trlop = str2double(response);
            opt.trlop = min(opt.trlop, size(opt.trlvis,1)); % should not be larger than the number of trials
            opt.trlop = max(opt.trlop, 1); % should not be smaller than 1
            if isfield(opt,'marks')
                opt.marks = [];
            end
            cfg.curr_displayed_marking = -1;
            setappdata(h, 'opt', opt);
            setappdata(h, 'cfg', cfg);
            redraw_cb(h, eventdata);
        end
        
    case 'shift+t'
        if strcmp(cfg.doSleepScoring,'yes')
            curr_epoch = opt.trlop;
            if (curr_epoch >= cfg.toggle_epoch_marker)
                cfg.toggle_epoch_marker = curr_epoch;
            else
                cfg.toggle_epoch_marker = 0;
            end
            setappdata(h, 'cfg', cfg);
            redraw_cb(h, eventdata);
        end
    case 'control+t'
        if strcmp(cfg.doSleepScoring,'yes')
            cfg.toggle_epoch_marker = 0;
            setappdata(h, 'cfg', cfg);
            redraw_cb(h, eventdata);
        end
    case 'h'
        if ~strcmp(cfg.doSleepScoring,'yes')
            % select the horizontal scaling
            response = inputdlg('horizontal scale', 'specify', 1, {num2str(cfg.blocksize)});
            if ~isempty(response)
                cfg.blocksize = str2double(response);
                setappdata(h, 'opt', opt);
                setappdata(h, 'cfg', cfg);
                definetrial_cb(h, eventdata);
                redraw_cb(h, eventdata);
            end
        end
    case 'shift+h'
        helpdlg(sprintf(help_text()),'Shortcut help');
    case 'l'
        if strcmp(cfg.doSleepScoring,'yes')
            try
                prompt = {'Spindle min Freq Ruler/Detection/Display [Hz]:','Spindle max Freq Ruler/Detection/Display[Hz]:',...
                    'Spindle threshold borders Ruler/Detection [potential]:','Spindle threshold peak Ruler/Detection [potential]:',...
                    'Spindle min duration Ruler/Detection [s]:','Spindle max duration Ruler/Detection [s]:',...
                    'SO min Freq Ruler/Detection [Hz]:','SO max Freq Ruler/Detection [Hz]:'...
                    'SO min Freq Display [Hz]:','SO max Freq Display [Hz]:'...
                    'SO threshold amplitude Ruler/Detection [potential]:'...
                    'Alpha min Freq Display [Hz]:','Alpha max Freq Display [Hz]:'...
                    'EMG threshold amplitude Ruler [potential]:'...
                    'Autosave hypnogram, #changes needed to autosave (0 = disabled):'...
                    'Confidence skip to lower than [0...1] (0 = disabled):'};
                dlg_title = 'Thresholds';
                num_lines = 1;
                defaultans = {num2str(cfg.sp_minFreq),num2str(cfg.sp_maxFreq),...
                    num2str(cfg.sp_thresholdForDetectionBeginEnd),num2str(cfg.sp_thresholdForDetectionCriterion),...
                    num2str(cfg.sp_minSec),num2str(cfg.sp_maxSec),...
                    num2str(cfg.so_minFreq),num2str(cfg.so_maxFreq),...
                    num2str(cfg.so_filter_minFreq),num2str(cfg.so_filter_maxFreq),...
                    num2str(cfg.so_thresholdAmplitudeForDetection),...
                    num2str(cfg.al_minFreq),num2str(cfg.al_maxFreq),...
                    num2str(cfg.emg_thresholdAmplitudeForDetection),...
                    num2str(cfg.autosave_hypnogram_every_number_change),...
                    num2str(cfg.confidence_skip_to_lower_than_threshold)};
                response = inputdlg(prompt,dlg_title,num_lines,defaultans);
                
                
                if ~isempty(response)
                    if numel(response)==16
                        cfg.sp_minFreq = str2num(['[' response{1} ']']);
                        cfg.sp_maxFreq = str2num(['[' response{2} ']']);
                        cfg.sp_thresholdForDetectionBeginEnd = str2num(['[' response{3} ']']);
                        cfg.sp_thresholdForDetectionCriterion = str2num(['[' response{4} ']']);
                        cfg.sp_minSec = str2num(['[' response{5} ']']);
                        cfg.sp_maxSec = str2num(['[' response{6} ']']);
                        
                        
                        cfg.so_minFreq = str2num(['[' response{7} ']']);
                        cfg.so_maxFreq = str2num(['[' response{8} ']']);
                        cfg.so_filter_minFreq = str2num(['[' response{9} ']']);
                        cfg.so_filter_maxFreq = str2num(['[' response{10} ']']);
                        cfg.so_thresholdAmplitudeForDetection = str2num(['[' response{11} ']']);
                        
                        cfg.al_minFreq = str2num(['[' response{12} ']']);
                        cfg.al_maxFreq = str2num(['[' response{13} ']']);
                        
                        cfg.emg_thresholdAmplitudeForDetection = str2num(['[' response{14} ']']);
                        
                        cfg.autosave_hypnogram_every_number_change = str2num(['[' response{15} ']']);
                        
                        cfg.confidence_skip_to_lower_than_threshold = str2num(['[' response{16} ']']);
                    else
                        error('not 16 elements')
                    end
                    
                    cfg = update_filters(cfg);
                    setappdata(h, 'opt', opt);
                    setappdata(h, 'cfg', cfg);
                    redraw_cb(h, eventdata);
                end
            catch err
                msgbox('Setting the thresholds failed!' ,'Thresholds failed','error','modal');
            end
        end
    case 'y'
        % select the vertical scaling
        response = inputdlg('vertical scale references, [ymin ymax], ''maxabs'' or ''maxmin''', 'specify', 1, {['[ ' num2str(cfg.ylim) ' ]']});
        if ~isempty(response)
            response = ['[' response{1} ']']; % convert to string and add brackets, just to ensure that str2num will work
            if strcmp(response, '[maxmin]')
                minval = min(opt.curdat.trial{1}(:));
                maxval = max(opt.curdat.trial{1}(:));
                cfg.ylim = [minval maxval];
            elseif strcmp(response, '[maxabs]')
                minval = min(opt.curdat.trial{1}(:));
                maxval = max(opt.curdat.trial{1}(:));
                cfg.ylim = [-max(abs([minval maxval])) max(abs([minval maxval]))];
            else
                tmp = str2num(response);
                if numel(tmp)==2
                    cfg.ylim = tmp;
                else
                    warning('incorrect specification of cfg.ylim, not changing the limits for the vertical axes')
                end
            end
            setappdata(h, 'opt', opt);
            setappdata(h, 'cfg', cfg);
            redraw_cb(h, eventdata);
        end
    case 'shift+c'
        %         % select channels
        %         select = match_str(opt.hdr.label, cfg.channel);
        %         select = select_channel_list(opt.hdr.label, select);
        %         cfg.channel = opt.hdr.label(select);
        %         setappdata(h, 'opt', opt);
        %         setappdata(h, 'cfg', cfg);
        %         delete(findobj(h,'tag', 'chanlabel'));  % remove channel labels here, and not in redrawing to save significant execution time (see bug 2065)
        %         redraw_cb(h, eventdata);
        
        temp_all_channels = opt.hdr.label;
        select = match_str(temp_all_channels, cfg.channel);
        [selected_channels, cfg, opt] = channelDialog(temp_all_channels,select,cfg,opt);
        if ~isempty(selected_channels)
            cfg.channel = selected_channels;
            setappdata(h, 'opt', opt);
            setappdata(h, 'cfg', cfg);
            delete(findobj(h,'tag', 'chanlabel'));  % remove channel labels here, and not in redrawing to save significant execution time (see bug 2065)
            redraw_cb(h, eventdata);
        end
        
    case 'i'
        if strcmp(cfg.viewmode, 'butterfly')
            delete(findobj(h,'tag', 'identify'));
            % click in data and get name of nearest channel
            fprintf('click in the figure to identify the name of the closest channel\n');
            val = ginput(1);
            pos = val(1);
            % transform 'val' to match data
            val(1) = val(1) * range(opt.hlim) + opt.hlim(1);
            val(2) = val(2) * range(opt.vlim) + opt.vlim(1);
            channame = val2nearestchan(opt.curdat,val);
            channb = match_str(opt.curdat.label,channame);
            fprintf('channel name: %s\n',channame);
            redraw_cb(h, eventdata);
            ft_plot_text(pos, 0.9, channame, 'FontSize', 16, 'tag', 'identify','interpreter','none');
            if ~ishold
                hold on
                ft_plot_vector(opt.curdat.time{1}, opt.curdat.trial{1}(channb,:), 'box', false, 'tag', 'identify', ...
                    'hpos', opt.laytime.pos(1,1), 'vpos', opt.laytime.pos(1,2), 'width', opt.laytime.width(1), 'height', opt.laytime.height(1), 'hlim', opt.hlim, 'vlim', opt.vlim, ...
                    'color', 'k', 'linewidth', 2);
                hold off
            else
                ft_plot_vector(opt.curdat.time{1}, opt.curdat.trial{1}(channb,:), 'box', false, 'tag', 'identify', ...
                    'hpos', opt.laytime.pos(1,1), 'vpos', opt.laytime.pos(1,2), 'width', opt.laytime.width(1), 'height', opt.laytime.height(1), 'hlim', opt.hlim, 'vlim', opt.vlim, ...
                    'color', 'k', 'linewidth', 2);
            end
        else
            warning('only supported with cfg.viewmode=''butterfly''');
        end
    case 'shift+s'
        if strcmp(cfg.doSleepScoring,'yes')
            
            try
                tempfilepath = [cfg.outputfilespath 'session_saved_' '.mat'];
                
                [hyp_file_name hyp_file_path hyp_file_filterindex] = uiputfile(...
                    {'*.mat','MAT-files (*.mat)';...
                    '*.*',  'All Files (*.*)'},...
                    'Save session as',...
                    tempfilepath);
                if hyp_file_filterindex ~= 0
                    save([hyp_file_path hyp_file_name],'cfg','opt','-v7.3')
                    msgbox('Saving session successful!' ,'Saving successful','modal');
                end
            catch err
                msgbox('Saving session failed!' ,'Saving failed','error','modal');
            end
        end
    case 'shift+o'
        if strcmp(cfg.doSleepScoring,'yes')
            try
                [h,opt,cfg,hyp_file_filterindex] = loadSession(h,opt,cfg);
                setappdata(h, 'opt', opt);
                setappdata(h, 'cfg', cfg);
                redraw_cb(h, eventdata);
                if hyp_file_filterindex ~= 0
                    msgbox('Opening session successful!' ,'Opening successful','modal');
                end
            catch err
                msgbox('Open the session failed!' ,'Open failed','error','modal');
            end
        end
    case 'shift+i'
        if strcmp(cfg.doSleepScoring,'yes')
            
            loaded_hypnogram_from_file = false;
            
            %%% first ask if to remove another additional hypnogram
            if numel(cfg.hypn_mult) > 0
                try
                    hypn_list_request = cellstr(num2str((1:numel(cfg.hypn_mult))'))';
                    [selection,ok] = listdlg('PromptString',['Choose hypnogram to Remove'],...
                        'SelectionMode','single',...
                        'InitialValue',numel(cfg.hypn_mult),...
                        'ListString',hypn_list_request,...
                        'OkString','Delete/Replace',...
                        'CancelString','Add');
                    if ok
                        
                        answer_hyp_replace = questdlg('Do you want to replace the primary hypnogram with the deleted one?', ...
                            'Replace Primary with deleted?', ...
                            'Yes','No','No');
                        switch answer_hyp_replace
                            case 'Yes'
                                
                                cfg.hypn = cfg.hypn_mult{selection};
                                cfg.hypn_plot_interpol = cfg.hypn_plot_interpol_mult{selection};
                                cfg.hypn_plot_interpol_MA = cfg.hypn_plot_interpol_MA_mult{selection};
                                cfg.hypn_plot_interpol_confidence = cfg.hypn_plot_interpol_confidence_mult{selection};

                                if isfield(cfg,'autosave_hypfilepath')
                                    temp_ArtifactPath = [cfg.autosave_hypfilepath '.artifacts.csv'];
                                    if exist(temp_ArtifactPath) == 2
                                        [o, c] = readArtifactFile(temp_ArtifactPath,opt,cfg,cfg.artifact_export_delimiter);
                                        opt = o;
                                        cfg = c;
                                        c = [];
                                        o = [];
                                    end
                                    
                                    [temp_pathstr,temp_name,temp_ext] = fileparts(cfg.autosave_hypfilepath);
                                    
                                    iAutosave = 0;
                                    cfg.autosave_hypfilepath = [temp_pathstr filesep temp_name '_autosave' num2str(iAutosave) temp_ext];
                                    while exist(cfg.autosave_hypfilepath) == 2
                                        iAutosave = iAutosave + 1;
                                        cfg.autosave_hypfilepath = [temp_pathstr filesep temp_name '_autosave' num2str(iAutosave) temp_ext];
                                    end
                                end
                        end
                        
                        [cfg] = removeAdditionalHypnogram(selection,cfg);
                        opt.hyp_figure_reload = true;
                        setappdata(h, 'cfg', cfg);
                        redraw_cb(h, eventdata);
                        return
                    end
                    
                catch err
                    msgbox('Remmoval of a hypnogram failed!' ,'Import failed','error','modal');
                    return
                end
            end
            if numel(cfg.hypn_mult) >= 4
                msgbox('Can only load 4 hypnograms at a time. Please remove one first before adding new' ,'Hypngoram limit reached','error','modal');
                return
            end
            
            
            try
                import_success = false;
                answer_hyp_type = questdlg('File or Automatic?', ...
                    'Read in Hypnogram?', ...
                    'From file','From Z3Score','From file');
                switch answer_hyp_type
                    case 'From file'
                        temppath = [cfg.outputfilespath];
                        
                        
                        list_formats = {'SpiSOP/Schlafaus/sleepin','Zmax','Somnomedics English', 'FASST', 'SleepTrip'};
                        list_formats_st = {'spisop','zmax','somnomedics_english','fasst','sleeptrip',};
                        [indx_file_formats, selected_file_format] = listdlg('ListString',list_formats,'SelectionMode','single','PromptString',{'Select a file format.',['Scoring will be converted to scoring standard = ' cfg.standard '.'],''},'InitialValue',1,'Name','File format');
                        
                        scoringformat = list_formats_st{indx_file_formats};
                        
                        [hyp_file_name hyp_file_path hyp_file_filterindex] = uigetfile(...
                            {'*.txt;*.csv;*.tsv','Import formats (*.txt,*.csv)';...
                            '*.txt','Text - Tab delimited (*.txt)';...
                            '*.tsv','Text - Tab delimited (*.tsv)';...
                            '*.csv','Comma Separated Values (*.csv)';...
                            % '*.m', 'program files (*.m)';...
                            % '*.fig','Figures (*.fig)';...
                            '*.mat','MAT-files from SleepTrip scoring export(*.mat)';...
                            '*.*',  'All Files (*.*)'},...
                            'Import scoring',...
                            temppath);
                        
                        if hyp_file_filterindex ~= 0
                            import_success = selected_file_format;
                        end
                        
                        begsample = opt.trlvis(opt.trlop,1);
                        endsample = opt.trlvis(opt.trlop,2);
                        temp_epochLengthSamples = endsample - begsample + 1;
                        
                        temp_hypnogramPath = [hyp_file_path hyp_file_name];
                        cfg_rs = [];
                        
                        cfg_rs.scoringfile      = temp_hypnogramPath;
                        cfg_rs.scoringformat    = scoringformat;
                        
                        cfg_rs.standard =  cfg.standard;

                        temp_scoring = st_read_scoring(cfg_rs);
                        cfg_tmp = [];
                        cfg_tmp.to = 'number';
                        temp_scoring = st_scoringconvert(cfg_tmp,temp_scoring);
                        temp_hypn = [cellfun(@str2num,temp_scoring.epochs,'UniformOutput',true)' temp_scoring.excluded']; 
                        loaded_hypnogram_from_file = true;
                    case 'From Z3Score'
                        ask_again = true;
                        while ask_again
                            cfg = [];
                            sca = AutoSleepScoringZ3score(cfg,opt.orgdata);
                            answer_temp = questdlg(sca.status(),'Automatic Scoring Status', ...
                                'OK','OK');
                            
                            %%sca = sca.update_channel();
                            answer_sca = questdlg('Perform Scoring?', ...
                                'Sure to try Automatic Sleep Scoring now?', ...
                                'Yes','Change','Cancel','Cancel');
                            switch answer_sca
                                case 'Yes'
                                    sca = sca.score();
                                    temp_hypn = sca.hypnogram;
                                    temp_hypn(:,3) = sca.hypnogram_confidence;
                                    import_success = true;
                                    ask_again = false;
                                case 'Change'
                                    ask_again = true;
                                case 'Cancel'
                                    ask_again = false;
                                    return
                            end
                        end
                end
                
                
                
                
                
                answer_hyp = questdlg('Import as additional hypnogram?', ...
                    'Hypnogram import', ...
                    'No, replace primary', ...
                    'Yes, create additional',...
                    'Cancel','No, replace primary');
                
                
                
                if ~import_success
                    msgbox('Importing the hypnogram failed!' ,'Import failed','error','modal');
                    return
                end
                
                
                
                temp_epochLengthSamples = opt.trlvis(1, 2) - opt.trlvis(1, 1) + 1;
                nEpochs = floor(size(opt.orgdata.trial{1},2)/temp_epochLengthSamples);
                
                if (size(temp_hypn,1) > nEpochs)
                    msgbox(sprintf(['Wrong Hypnogram?\n' 'It is too long!\n' ' data: ' num2str(nEpochs) ' ep\n' ' import: ' num2str(size(temp_hypn,1)) ' ep\n' 'hypnogram will be truncated to data by cutting its tail']) ,'Wrong Hypnogram imported?', 'warn','modal');
                    temp_hypn = temp_hypn(1:nEpochs,:);
                end
                
                
                
                if size(temp_hypn,1) < nEpochs
                    missingEpochs = nEpochs - size(temp_hypn,1);
                    %hypn(end+1:end+missingEpochs,:) = [ones(1,missingEpochs,1)*-1 zeros(0,missingEpochs,1)];
                    hypn_missing = zeros(missingEpochs,size(temp_hypn,2));
                    hypn_missing(:,1) = -1;
                    temp_hypn = [temp_hypn ; hypn_missing];
                end
                
                
                % Handle response
                switch answer_hyp
                    case 'No, replace primary'
                        
                        cfg.hypn = temp_hypn;
                        
                        [curr_ep_hypn_plot_interpol curr_ep_hypn_plot_interpol_MA] = interpolate_hypn_for_plot(cfg.hypn,cfg.hyp_epochLengthSamples,cfg.plot_MA_offset,istrue(cfg.yaxdisteqi));
                        
                        cfg.hypn_plot_interpol = curr_ep_hypn_plot_interpol;
                        cfg.hypn_plot_interpol_MA = curr_ep_hypn_plot_interpol_MA;
                        
                        %has confidence ratings
                        cfg.hypn_plot_interpol_confidence = [];
                        if size(cfg.hypn,2) > 2
                            if (max(cfg.hypn(:,3)) <= 1) && (min(cfg.hypn(:,3)) >=0)
                                [dummy_temp curr_ep_hypn_plot_interpol_confidence] = interpolate_hypn_for_plot(cfg.hypn(:,2:3),cfg.hyp_epochLengthSamples,cfg.plot_confidence_offset,istrue(cfg.yaxdisteqi));
                                cfg.hypn_plot_interpol_confidence = curr_ep_hypn_plot_interpol_confidence;
                            end
                        end
                        
                        if loaded_hypnogram_from_file
                            
                            temp_ArtifactPath = [temp_hypnogramPath '.artifacts.csv'];
                            if exist(temp_ArtifactPath) == 2
                                [o, c] = readArtifactFile(temp_ArtifactPath,opt,cfg,cfg.artifact_export_delimiter);
                                opt = o;
                                cfg = c;
                                c = [];
                                o = [];
                            end
                            
                            [temp_pathstr,temp_name,temp_ext] = fileparts(temp_hypnogramPath);
                            
                            iAutosave = 0;
                            cfg.autosave_hypfilepath = [temp_pathstr filesep temp_name '_autosave' num2str(iAutosave) temp_ext];
                            while exist(cfg.autosave_hypfilepath) == 2
                                iAutosave = iAutosave + 1;
                                cfg.autosave_hypfilepath = [temp_pathstr filesep temp_name '_autosave' num2str(iAutosave) temp_ext];
                            end
                        end
                        
                        
                    case 'Yes, create additional'
                        
                        if ~isfield(cfg,'hypn_mult')
                            cfg.hypn_mult = {};
                            cfg.hypn_plot_interpol_mult = {};
                            cfg.hypn_plot_interpol_MA_mult = {};
                             cfg.hypn_plot_interpol_confidence_mult = {};
                            cfg.hypn_mult_idx = 1;
                        end
                        if ~isfield(cfg,'hypn_mult_idx')
                            cfg.hypn_mult_idx = 1;
                        end
                        
                        if cfg.hypn_mult_idx > 4
                            if numel(cfg.hypn_mult) >= 4
                                msgbox('Importing the hypnogram failed!' ,'Only 4 Hypnograms allowed','error','modal');
                                return
                            end
                        end
                        %cfg.hypn_mult_idx = cfg.hypn_mult_idx+1;
                        
                        
                        
                        cfg.hypn_mult{cfg.hypn_mult_idx} = temp_hypn;
                        
                        [curr_ep_hypn_plot_interpol curr_ep_hypn_plot_interpol_MA] = interpolate_hypn_for_plot(cfg.hypn_mult{cfg.hypn_mult_idx},cfg.hyp_epochLengthSamples,cfg.plot_MA_offset,istrue(cfg.yaxdisteqi));
                        
                        curr_ep_hypn_plot_interpol_confidence = [];
                        if size(cfg.hypn_mult{cfg.hypn_mult_idx},2) > 2
                            if (max(cfg.hypn_mult{cfg.hypn_mult_idx}(:,3)) <= 1) && (min(cfg.hypn_mult{cfg.hypn_mult_idx}(:,3)) >=0)
                                [dummy_temp curr_ep_hypn_plot_interpol_confidence] = interpolate_hypn_for_plot(cfg.hypn_mult{cfg.hypn_mult_idx}(:,2:3),cfg.hyp_epochLengthSamples,cfg.plot_confidence_offset,istrue(cfg.yaxdisteqi));
                            end
                        end
                        
                        cfg.hypn_plot_interpol_mult{cfg.hypn_mult_idx} = curr_ep_hypn_plot_interpol;
                        cfg.hypn_plot_interpol_MA_mult{cfg.hypn_mult_idx} = curr_ep_hypn_plot_interpol_MA;
                        cfg.hypn_plot_interpol_confidence_mult{cfg.hypn_mult_idx} = curr_ep_hypn_plot_interpol_confidence;
                        cfg.hypn_mult_idx = cfg.hypn_mult_idx+1;
                        
                    case 'Cancel'
                        return
                end
                
                opt.hyp_figure_reload = true;
                
                setappdata(h, 'opt', opt);
                setappdata(h, 'cfg', cfg);
                redraw_cb(h, eventdata);
                if import_success
                    msgbox('Importing the hypnogram successful!' ,'Import successful','modal');
                end
            catch err
                msgbox('Importing the hypnogram failed!' ,'Import failed','error','modal');
            end
        end
    case 'shift+e'
        if strcmp(cfg.doSleepScoring,'yes')
            
            try
                
                tempfilepath = [cfg.outputfilespath 'hypn_export' '.txt'];
                
                [hyp_file_name hyp_file_path hyp_file_filterindex] = uiputfile(...
                    {'*.txt;*.csv','Export formats (*.txt,*.csv)';...
                    '*.txt','Text - Tab delimited (*.txt)';...
                    '*.csv','Comma Separated Values (*.csv)';...
                    % '*.m', 'program files (*.m)';...
                    % '*.fig','Figures (*.fig)';...
                    % '*.mat','MAT-files (*.mat)';...
                    '*.*',  'All Files (*.*)'},...
                    'Export Hypnogram as',...
                    tempfilepath);
                if hyp_file_filterindex ~= 0
                    [dummy_pathstr,dummy_name,temp_ext] = fileparts([hyp_file_path hyp_file_name]);
                    if strcmp(temp_ext,'.csv')
                        cfg.hypnogram_delimiter_autosave = ',';
                    else
                        cfg.hypnogram_delimiter_autosave = '\t';
                    end
                    
                    temp_hypnogramPath = [hyp_file_path hyp_file_name];
                    writeHypnogramFile(temp_hypnogramPath,cfg.hypn,cfg.hypnogram_delimiter_autosave);
                    
                    temp_ArtifactPath = [hyp_file_path hyp_file_name '.artifacts.csv'];
                    writeArtifactFile(temp_ArtifactPath,opt,cfg.artifact_export_delimiter);
                    
                    
                    [temp_pathstr,temp_name,temp_ext] = fileparts(temp_hypnogramPath);
                    iAutosave = 0;
                    cfg.autosave_hypfilepath = [temp_pathstr filesep temp_name '_autosave' num2str(iAutosave) temp_ext];
                    while exist(cfg.autosave_hypfilepath) == 2
                        iAutosave = iAutosave + 1;
                        cfg.autosave_hypfilepath = [temp_pathstr filesep temp_name '_autosave' num2str(iAutosave) temp_ext];
                    end
                    
                    
                    
                    
                    
                    setappdata(h, 'cfg', cfg);
                    msgbox('Exporting the hypnogram successful!' ,'Export successful','modal');
                end
            catch err
                msgbox('Exporting the hypnogram failed!' ,'Export failed','error','modal');
            end
        end
        %     case 'alt+e'
        %         if strcmp(cfg.doSleepScoring,'yes')
        %
        %             % select one channel
        %             cfg.score_channel_eeg_number = select_channel_list(opt.hdr.label, cfg.score_channel_eeg_number,'select one channel');
        %             cfg.score_channel_eeg_number = cfg.score_channel_eeg_number(1);
        %             ft_uilayout(h, 'tag', 'scoptbuttons_focusEEG', 'string', ['Focus EEG: ' opt.hdr.label{cfg.score_channel_eeg_number}]);
        %             setappdata(h, 'opt', opt);
        %             setappdata(h, 'cfg', cfg);
        %             redraw_cb(h, eventdata);
        %         end
        %     case 'alt+o'
        %         if strcmp(cfg.doSleepScoring,'yes')
        %
        %             % select one channel
        %             cfg.score_channel_eog_number = select_channel_list(opt.hdr.label, cfg.score_channel_eog_number,'select one channel');
        %             cfg.score_channel_eog_number = cfg.score_channel_eog_number(1);
        %             ft_uilayout(h, 'tag', 'scoptbuttons_focusEOG', 'string', ['Focus EOG: ' opt.hdr.label{cfg.score_channel_eog_number}]);
        %             setappdata(h, 'opt', opt);
        %             setappdata(h, 'cfg', cfg);
        %             redraw_cb(h, eventdata);
        %         end
        %     case 'alt+m'
        %         if strcmp(cfg.doSleepScoring,'yes')
        %
        %             % select one channel
        %             cfg.score_channel_emg_number = select_channel_list(opt.hdr.label, cfg.score_channel_emg_number,'select one channel');
        %             cfg.score_channel_emg_number = cfg.score_channel_emg_number(1);
        %             ft_uilayout(h, 'tag', 'scoptbuttons_focusEMG', 'string', ['Focus EMG: ' opt.hdr.label{cfg.score_channel_emg_number}]);
        %             setappdata(h, 'opt', opt);
        %             setappdata(h, 'cfg', cfg);
        %             redraw_cb(h, eventdata);
        %         end
    case 's'
        % toggle between selectmode options: switch from 'markartifact', to 'markpeakevent' to 'marktroughevent' and back with on screen feedback
        %         curstate = find(strcmp(cfg.selectmode, {'markartifact', 'markpeakevent', 'marktroughevent'}));
        %         if curstate == 1
        %             cfg.selectmode = 'markpeakevent';
        %         elseif curstate == 2
        %             cfg.selectmode = 'marktroughevent';
        %         elseif curstate == 3
        cfg.selectmode = 'markartifact';
        %         end
        
        fprintf('switching to selectmode = %s\n',cfg.selectmode);
        setappdata(h, 'opt', opt);
        setappdata(h, 'cfg', cfg);
        redraw_cb(h, eventdata);
    case 'control+control'
        % do nothing
    case 'shift+shift'
        % do nothing
    case 'alt+alt'
        % do nothing
    otherwise
        setappdata(h, 'opt', opt);
        setappdata(h, 'cfg', cfg);
        help_cb(h);
end

if strcmp(cfg.doSleepScoring,'yes')
    if strcmp(cfg.markSO,'yes')
        ft_uilayout(h, 'tag', 'scoptbuttons_SOdet', 'BackgroundColor',  cfg.slowoscillation_mark_color);
        ft_uilayout(h, 'tag', 'scoptbuttons_SOdet', 'FontWeight', 'bold');
        if cfg.SOdetection_orientation == 1;
            ft_uilayout(h, 'tag', 'scoptbuttons_SOdet', 'string', '(+)_.�\_.�\_.�');
        else
            ft_uilayout(h, 'tag', 'scoptbuttons_SOdet', 'string', '(-)`�/�`�/�`�/');
        end
    else
        ft_uilayout(h, 'tag', 'scoptbuttons_SOdet', 'BackgroundColor', [0.5 0.5 0.5]);
        ft_uilayout(h, 'tag', 'scoptbuttons_SOdet', 'FontWeight', 'normal');
        %ft_uilayout(h, 'tag', 'scoptbuttons_SOdet', 'string', 'off');
    end
    
    if strcmp(cfg.markSpindles,'yes')
        ft_uilayout(h, 'tag', 'scoptbuttons_SPdet', 'BackgroundColor', cfg.spindle_mark_color);
        ft_uilayout(h, 'tag', 'scoptbuttons_SPdet', 'FontWeight', 'bold');
        
    else
        ft_uilayout(h, 'tag', 'scoptbuttons_SPdet', 'BackgroundColor', [0.5 0.5 0.5]);
        ft_uilayout(h, 'tag', 'scoptbuttons_SPdet', 'FontWeight', 'normal');
        
    end
    
    if strcmp(cfg.underlaySpindleSignal,'yes')
        ft_uilayout(h, 'tag', 'scoptbuttons_SPdisp', 'BackgroundColor', cfg.underlaySpindleSignal_color);
        ft_uilayout(h, 'tag', 'scoptbuttons_SPdisp', 'FontWeight', 'bold');
    else
        ft_uilayout(h, 'tag', 'scoptbuttons_SPdisp', 'BackgroundColor', [0.5 0.5 0.5]);
        ft_uilayout(h, 'tag', 'scoptbuttons_SPdisp', 'FontWeight', 'normal');
    end
    
    if strcmp(cfg.underlayAlphaSignal,'yes')
        ft_uilayout(h, 'tag', 'scoptbuttons_ALdisp', 'BackgroundColor', cfg.underlayAlphaSignal_color);
        ft_uilayout(h, 'tag', 'scoptbuttons_ALdisp', 'FontWeight', 'bold');
    else
        ft_uilayout(h, 'tag', 'scoptbuttons_ALdisp', 'BackgroundColor', [0.5 0.5 0.5]);
        ft_uilayout(h, 'tag', 'scoptbuttons_ALdisp', 'FontWeight', 'normal');
    end
        
    if strcmp(cfg.underlaySOSignal,'yes')
        ft_uilayout(h, 'tag', 'scoptbuttons_SOdisp', 'BackgroundColor', cfg.underlaySOSignal_color);
        ft_uilayout(h, 'tag', 'scoptbuttons_SOdisp', 'FontWeight', 'bold');
    else
        ft_uilayout(h, 'tag', 'scoptbuttons_SOdisp', 'BackgroundColor', [0.5 0.5 0.5]);
        ft_uilayout(h, 'tag', 'scoptbuttons_SOdisp', 'FontWeight', 'normal');
    end
    
    if strcmp(cfg.markECG,'yes') && (cfg.ECG_signalMultiplicator == 1)
        ft_uilayout(h, 'tag', 'scoptbuttons_HRdisp', 'string', ['+<3+<3+<3']);
        ft_uilayout(h, 'tag', 'scoptbuttons_HRdisp', 'FontWeight', 'bold')
        ft_uilayout(h, 'tag', 'scoptbuttons_HRdisp', 'BackgroundColor', cfg.score_channel_ecg_color);
        
    elseif strcmp(cfg.markECG,'yes') && (cfg.ECG_signalMultiplicator == -1)
        ft_uilayout(h, 'tag', 'scoptbuttons_HRdisp', 'string', ['-<3-<3-<3']);
        ft_uilayout(h, 'tag', 'scoptbuttons_HRdisp', 'FontWeight', 'bold')
        ft_uilayout(h, 'tag', 'scoptbuttons_HRdisp', 'BackgroundColor', cfg.score_channel_ecg_color);
    else
        ft_uilayout(h, 'tag', 'scoptbuttons_HRdisp', 'string', ['<3 <3 <3']);
        ft_uilayout(h, 'tag', 'scoptbuttons_HRdisp', 'FontWeight', 'normal')
        ft_uilayout(h, 'tag', 'scoptbuttons_HRdisp', 'BackgroundColor', [0.5 0.5 0.5]);
    end
    
    
    if strcmp(cfg.displayEvents,'yes')
        %ft_uilayout(h, 'tag', 'scoptbuttons_SPdisp', 'BackgroundColor', cfg.underlaySpindleSignal_color);
        ft_uilayout(h, 'tag', 'scoptbuttons_EvDisp', 'FontWeight', 'bold');
        ft_uilayout(h, 'tag', 'scoptbuttons_EvDisp', 'string', 'EVENTS');
        
    else
        %ft_uilayout(h, 'tag', 'scoptbuttons_SPdisp', 'BackgroundColor', [0.5 0.5 0.5]);
        ft_uilayout(h, 'tag', 'scoptbuttons_EvDisp', 'FontWeight', 'normal');
        ft_uilayout(h, 'tag', 'scoptbuttons_EvDisp', 'string', 'events');
        
    end
    
    if strcmp(opt.zoomstatus,'on')
        ft_uilayout(h, 'tag', 'scoptbuttons_zoom', 'string', ['zOOm']);
    else
        ft_uilayout(h, 'tag', 'scoptbuttons_zoom', 'string', ['zoom']);
        zoom(h,'out')
    end
    if strcmp(opt.markingstatus,'on')
        set(h, 'WindowButtonMotionFcn',   {@mouse_move_cb, 'h_main',h});
        set(h, 'WindowButtonDownFcn',   {@select_marks_cb, 'h_main',h,'event','WindowButtonDownFcn'});
        set(h, 'WindowButtonUpFcn',   {@select_marks_cb, 'h_main',h,'event','WindowButtonUpFcn'});
        
        ft_uilayout(h, 'tag', 'scoptbuttons_mark', 'string', ['MARK']);
        
    else
        if strcmp(cfg.use_ruler,'no')
            set(h, 'WindowButtonMotionFcn', {@ft_select_range, 'multiple', false, 'xrange', true, 'yrange', false, 'clear', true, 'contextmenu', cfg.selfun, 'callback', {@select_range_cb, h}, 'event', 'WindowButtonMotionFcn'});
        end
        set(h, 'WindowButtonDownFcn',   {@ft_select_range, 'multiple', false, 'xrange', true, 'yrange', false, 'clear', true, 'contextmenu', cfg.selfun, 'callback', {@select_range_cb, h}, 'event', 'WindowButtonDownFcn'});
        set(h, 'WindowButtonUpFcn',     {@ft_select_range, 'multiple', false, 'xrange', true, 'yrange', false, 'clear', true, 'contextmenu', cfg.selfun, 'callback', {@select_range_cb, h}, 'event', 'WindowButtonUpFcn'});
        ft_uilayout(h, 'tag', 'scoptbuttons_mark', 'string', ['mark']);
        
    end
    
    if strcmp(cfg.use_ruler,'yes')
        set(h, 'WindowButtonMotionFcn',   {@mouse_move_cb, 'h_main',h});
        ft_uilayout(h, 'tag', 'scoptbuttons_ruler', 'string', ['RULER']);
    else
        delete(findobj(h, 'tag', 'scale_help'));
        
        if strcmp(opt.markingstatus,'off')
            set(h, 'WindowButtonMotionFcn', {@ft_select_range, 'multiple', false, 'xrange', true, 'yrange', false, 'clear', true, 'contextmenu', cfg.selfun, 'callback', {@select_range_cb, h}, 'event', 'WindowButtonMotionFcn'});
        end
        ft_uilayout(h, 'tag', 'scoptbuttons_ruler', 'string', ['ruler']);
    end
    
    switch cfg.skip_to_next
        case 'always'
            ft_uilayout(h, 'tag', 'scoptbuttons_nextunk', 'string', ['[>]']);
        case 'firstscore'
            ft_uilayout(h, 'tag', 'scoptbuttons_nextunk', 'string', ['[?>]']);
        case 'unknown'
            ft_uilayout(h, 'tag', 'scoptbuttons_nextunk', 'string', ['[>?]']);
        case 'stay'
            ft_uilayout(h, 'tag', 'scoptbuttons_nextunk', 'string', ['[><]']);
    end

    
    if istrue(cfg.drawgrid)
        ft_uilayout(h, 'tag', 'scoptbuttons_grid', 'string', ['|gr:id|']);
    else
        ft_uilayout(h, 'tag', 'scoptbuttons_grid', 'string', ['grid']);
    end
    
     if istrue(cfg.plot_stage_signatures) && (istrue(cfg.highlight_scoring_channels) || ~istrue(cfg.has_more_than_3_channels))
        ft_uilayout(h, 'tag', 'scoptbuttons_min', 'string', ['min']);
     else
        ft_uilayout(h, 'tag', 'scoptbuttons_min', 'string', ['MIN']);
     end
    
    if istrue(cfg.display_power_spectrum)
        ft_uilayout(h, 'tag', 'scoptbuttons_pow', 'string', ['POW']);
    else
        ft_uilayout(h, 'tag', 'scoptbuttons_pow', 'string', ['pow']);
    end
    
    if istrue(cfg.display_time_frequency)
        ft_uilayout(h, 'tag', 'scoptbuttons_tfr', 'string', ['TFR']);
    else
        ft_uilayout(h, 'tag', 'scoptbuttons_tfr', 'string', ['tfr']);
    end
    
    
    
    ft_uilayout(h, 'tag', 'artifactui_button', 'string', ['artfct(' opt.artdata.label{opt.ftsel} ')']);
    
end
uiresume(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function toggle_viewmode_cb(h, eventdata, varargin)
% FIXME should be used
opt = guidata(getparent(h));
if ~isempty(varargin) && ischar(varargin{1})
    cfg.viewmode = varargin{1};
end
guidata(getparent(h), opt);
delete(findobj(h,'tag', 'chanlabel'));  % remove channel labels here, and not in redrawing to save significant execution time (see bug 2065)
redraw_cb(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = getparent(h)
p = h;
while p~=0
    h = p;
    p = get(h, 'parent');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function redraw_cb(h, eventdata)
h = getparent(h);
opt = getappdata(h, 'opt');
cfg = getappdata(h, 'cfg');

%fprintf('redrawing with viewmode %s\n', cfg.viewmode);

begsample = opt.trlvis(opt.trlop, 1);
endsample = opt.trlvis(opt.trlop, 2);
offset    = opt.trlvis(opt.trlop, 3);
chanindx  = match_str(opt.hdr.label, cfg.channel);


% FW begin
if strcmp(cfg.bgcolor,'dark')
    whitebg(h,'k');
end

if strcmp(cfg.doSleepScoring,'yes')
    if isfield(cfg,'plotHyp')
        if isfield(cfg,'hhyp')
            if ishandle(cfg.hhyp)
                %figure(cfg.hhyp)
            else
                cfg.hhyp = figure;
                set(cfg.hhyp, 'WindowButtonDownFcn',   {@select_sleep_stage_cb, 'h_main',h});
                set(cfg.hhyp, 'NumberTitle', 'off');
                cfg.hhypfigax = gca;
            end
            
        else
            cfg.hhyp = figure;
            figure(cfg.hhyp)
            set(cfg.hhyp, 'WindowButtonDownFcn',   {@select_sleep_stage_cb, 'h_main',h});
            set(cfg.hhyp, 'NumberTitle', 'off');
            cfg.hhypfigax = gca;
            %cfg.hhypfig = gcf;
        end


        
        if strcmp(cfg.bgcolor,'dark')
            whitebg(cfg.hhyp,'k');
            set(cfg.hhyp,'color',[0 0 0]);
        else
            set(cfg.hhyp,'color',[1 1 1]);
        end
        
        
        set(cfg.hhypfigax,'Fontsize',8,'FontUnits','normalized');
        %figure(cfg.hhyp);
        %h.hhyp = getparent(cfg.hhyp);
        temp_max_y = 1.25;
        if ~(isfield(cfg,'hyp_x_time') && isfield(cfg,'hyp_x_time_hyp'))
            opt.hyp_figure_reload = true;
        end
        
        
        curr_epoch = opt.trlop;
        
        temp_epochLengthSamples = opt.trlvis(1, 2) - opt.trlvis(1, 1) + 1;
        nEpochs = floor(size(opt.orgdata.trial{1},2)/temp_epochLengthSamples);
        if curr_epoch > nEpochs
            curr_epoch = curr_epoch - 1;
        end
        
        hyp_begsample = cfg.hyp_epochLengthSamples*(curr_epoch-1)+1; % opt.trlvis(opt.trlop,1);
        hyp_endsample = cfg.hyp_epochLengthSamples*(curr_epoch);% opt.trlvis(opt.trlop,2);
        
        if opt.hyp_figure_reload
            beg_time = opt.orgdata.time{1,1}(1);
            end_time = opt.orgdata.time{1,1}(end);
            cfg.hyp_x_time = linspace(beg_time,end_time,ceil((end_time-beg_time)*cfg.hyp_fample));
            cfg.hyp_x_time = cfg.hyp_x_time/60;
            %cfg.hypn_plot_interpol = cfg.hypn_plot_interpol(1:min(length(cfg.hypn_plot_interpol),length(cfg.hyp_x_time)));
            %cfg.hypn_plot_interpol_MA = cfg.hypn_plot_interpol_MA(1:min(length(cfg.hypn_plot_interpol_MA),length(cfg.hyp_x_time)));
            cfg.hyp_x_time_hyp = cfg.hyp_x_time(1:min(length(cfg.hypn_plot_interpol),length(cfg.hyp_x_time)));
            %x_time_hyp = x_time(1:min(length(x_time),length(cfg.hypn_plot_interpol)));
            opt.hyp_figure_reload = false;
        end
        axh = cfg.hhypfigax;
        
        if strcmp(cfg.bgcolor,'dark')
            tempcolor = [1 1 1];
        else
            tempcolor = [0 0 0];
        end
        
        %%% plot primary hypnogram
        plot(axh,cfg.hyp_x_time_hyp(1:min(length(cfg.hyp_x_time_hyp),length(cfg.hypn_plot_interpol))),cfg.hypn_plot_interpol(1:min(length(cfg.hyp_x_time_hyp),length(cfg.hypn_plot_interpol))),'Color',tempcolor)
        hold(axh,'all');
        
        xlim(axh,[0 max(cfg.hyp_x_time)]);
        ylabel(axh,'Stages');
        ylim(axh,[cfg.plot_MA_offset temp_max_y])
        
        plot(axh,cfg.hyp_x_time_hyp(1:min(length(cfg.hyp_x_time_hyp),length(cfg.hypn_plot_interpol_MA))),cfg.hypn_plot_interpol_MA(1:min(length(cfg.hyp_x_time_hyp),length(cfg.hypn_plot_interpol_MA))),'Color',[1 0 0])
        
        temp_x_other_artifact = ~(cfg.hypn_plot_interpol_MA > -4.5);
        if any(~temp_x_other_artifact)
            temp_hypn_plot_interpol_MA = cfg.hypn_plot_interpol_MA;
            temp_hypn_plot_interpol_MA(temp_x_other_artifact) = -4.45;
            %temp_x_other_artifact_index = [temp_x_other_artifact; temp_x_other_artifact+1];
            %temp_x_other_artifact_index = temp_x_other_artifact_index(1:end-1);
            %temp_x_other_artifact_index = sort(temp_x_other_artifact_index);
            plot(axh,cfg.hyp_x_time_hyp(1:min(length(cfg.hyp_x_time_hyp),length(temp_hypn_plot_interpol_MA))),temp_hypn_plot_interpol_MA(1:min(length(cfg.hyp_x_time_hyp),length(temp_hypn_plot_interpol_MA))),'Color',[0.75 0.75 0.75])
        end
        
        temp_do_plot_confidence = false;
        if ~isempty(cfg.hypn_plot_interpol_confidence)
            temp_do_plot_confidence = true;
            plot(axh,cfg.hyp_x_time_hyp(1:min(length(cfg.hyp_x_time_hyp),length(cfg.hypn_plot_interpol_confidence))),cfg.hypn_plot_interpol_confidence(1:min(length(cfg.hyp_x_time_hyp),length(cfg.hypn_plot_interpol_confidence))),'Color',[0 0 0])
        end
        
        %%% end plot primary hypnogram
        
        %%% plot secondary hypnograms
        
        if ~isempty(cfg.hypn_mult)
            
            hyp_mult_colors = jet(4);
            hyp_mult_colors_MA = hyp_mult_colors;
            tempHypCol_ridx = 1;
            y_shift_step = 0.05;
            for iHypMult = 1:numel(cfg.hypn_mult)
                
                tempcolor = hyp_mult_colors(tempHypCol_ridx,:);
                
                plot(axh,cfg.hyp_x_time_hyp(1:min(length(cfg.hyp_x_time_hyp),length(cfg.hypn_plot_interpol_mult{iHypMult}))),cfg.hypn_plot_interpol_mult{iHypMult}(1:min(length(cfg.hyp_x_time_hyp),length(cfg.hypn_plot_interpol_mult{iHypMult})))-iHypMult*y_shift_step,'Color',tempcolor)
                hold(axh,'all');
                
                xlim(axh,[0 max(cfg.hyp_x_time)]);
                ylabel(axh,'Stages');
                ylim(axh,[cfg.plot_MA_offset temp_max_y])
                
                plot(axh,cfg.hyp_x_time_hyp(1:min(length(cfg.hyp_x_time_hyp),length(cfg.hypn_plot_interpol_MA_mult{iHypMult}))),cfg.hypn_plot_interpol_MA_mult{iHypMult}(1:min(length(cfg.hyp_x_time_hyp),length(cfg.hypn_plot_interpol_MA_mult{iHypMult})))-iHypMult*y_shift_step,'Color',tempcolor)
                
                temp_x_other_artifact = ~((cfg.hypn_plot_interpol_MA_mult{iHypMult}-iHypMult*y_shift_step) > -4.5);
                if any(~temp_x_other_artifact)
                    temp_hypn_plot_interpol_MA = cfg.hypn_plot_interpol_MA_mult{iHypMult};
                    temp_hypn_plot_interpol_MA(temp_x_other_artifact) = -4.45-iHypMult*y_shift_step;
                    %temp_x_other_artifact_index = [temp_x_other_artifact; temp_x_other_artifact+1];
                    %temp_x_other_artifact_index = temp_x_other_artifact_index(1:end-1);
                    %temp_x_other_artifact_index = sort(temp_x_other_artifact_index);
                    plot(axh,cfg.hyp_x_time_hyp(1:min(length(cfg.hyp_x_time_hyp),length(temp_hypn_plot_interpol_MA))),temp_hypn_plot_interpol_MA(1:min(length(cfg.hyp_x_time_hyp),length(temp_hypn_plot_interpol_MA)))-iHypMult*y_shift_step,'Color',tempcolor)
                end
                
                %             cfg.hypn_mult = {};
                %             cfg.hypn_plot_interpol_mult = {};
                %             cfg.hypn_plot_interpol_MA_mult = {};
                %             cfg.hypn_mult_idx = 1;
                
                if ~isempty(cfg.hypn_plot_interpol_confidence_mult{iHypMult})
                    temp_do_plot_confidence = true;
                    plot(axh,cfg.hyp_x_time_hyp(1:min(length(cfg.hyp_x_time_hyp),length(cfg.hypn_plot_interpol_confidence_mult{iHypMult}))),cfg.hypn_plot_interpol_confidence_mult{iHypMult}(1:min(length(cfg.hyp_x_time_hyp),length(cfg.hypn_plot_interpol_confidence_mult{iHypMult})))-iHypMult*y_shift_step,'Color',tempcolor)
                end
                
                tempHypCol_ridx = tempHypCol_ridx+1;
                if tempHypCol_ridx > size(hyp_mult_colors,1)
                    tempHypCol_ridx = 1;
                end
            end
            
        end
        %%% end plot secondary hypnograms
        
        
        if temp_do_plot_confidence
            yTick = [0.75 0.25 0 -0.5 -1 -2 -3 -4 cfg.plot_MA_offset+1 cfg.plot_MA_offset+0.5];
            yTickLabel = {'max conf' 'min conf' 'W' 'REM' 'S1' 'S2' 'S3' 'S4' 'MT' 'MA'};
        else
            yTick = [0 -0.5 -1 -2 -3 -4 cfg.plot_MA_offset+1 cfg.plot_MA_offset+0.5];
            yTickLabel = {'W' 'REM' 'S1' 'S2' 'S3' 'S4' 'MT' 'MA'};
        end
        if strcmp(cfg.displayEvents,'yes')
            
            ev1_offset = 0.9;
            ev2_offset = 0.8;
            
            if isfield(cfg,'begin_end_events2')
                yTick = [ev2_offset yTick];
                yTickLabel = {'Ev2' yTickLabel{:}};
                for channelIndex = 1:length(chanindx)
                    curr_begins_ends = cfg.begin_end_events2{chanindx(channelIndex)};
                    if strcmp(cfg.viewmode, 'component')
                        color = 'k';
                    else
                        %color = opt.chancolors(chanindx(i),:);
                        color = opt.chancolors(chanindx(channelIndex),:);
                    end
                    
                    temp_x = (curr_begins_ends(:,1)/opt.fsample)/60;
                    temp_y = repmat(ev2_offset,size(curr_begins_ends,1),1);
                    plot(axh,[temp_x temp_x]',[temp_y-0.04 temp_y+0.04]','Color',color)
                    %scatter(axh,(curr_begins_ends(:,1)/opt.fsample)/60,repmat(ev1_offset,1,size(curr_begins_ends,1)),'MarkerEdgeColor',color)
                end
            end
            
            if isfield(cfg,'begin_end_events')
                yTick = [ev1_offset yTick];
                yTickLabel = {'Ev1' yTickLabel{:}};
                for channelIndex = 1:length(chanindx)
                    curr_begins_ends = cfg.begin_end_events{chanindx(channelIndex)};
                    if strcmp(cfg.viewmode, 'component')
                        color = 'k';
                    else
                        %color = opt.chancolors(chanindx(i),:);
                        color = opt.chancolors(chanindx(channelIndex),:);
                    end
                    
                    temp_x = (curr_begins_ends(:,1)/opt.fsample)/60;
                    temp_y = repmat(ev1_offset,size(curr_begins_ends,1),1);
                    plot(axh,[temp_x temp_x]',[temp_y-0.05 temp_y+0.05]','Color',color)
                    %scatter(axh,(curr_begins_ends(:,1)/opt.fsample)/60,repmat(ev1_offset,1,size(curr_begins_ends,1)),'MarkerEdgeColor',color)
                end
            end
            
            
        end
        
        yTick = [1 yTick];
        yTickLabel = {'?' yTickLabel{:}};
        
        set(axh, 'yTick', flip(yTick));
        set(axh, 'yTickLabel', flip(yTickLabel));
        set(axh,'TickDir','out');
        xTick = [0:20:max(cfg.hyp_x_time)];
        set(axh, 'xTick', xTick);
        x_pos_begin = cfg.hyp_x_time(hyp_begsample);
        x_pos_end = cfg.hyp_x_time(hyp_endsample);
        x_pos = [x_pos_begin x_pos_end x_pos_end x_pos_begin];
        y_pos = [cfg.plot_MA_offset cfg.plot_MA_offset 1 1];
        pos_now = patch(x_pos,y_pos,[0.5 0.25 1],'parent',axh);
        set(pos_now,'FaceAlpha',0.4);
        set(pos_now,'EdgeColor','none');
        
        if cfg.toggle_epoch_marker
            x_pos_begin_toggle_epoch_marker = cfg.hyp_x_time(cfg.hyp_epochLengthSamples*(cfg.toggle_epoch_marker-1)+1); % opt.trlvis(opt.trlop,1);
            line([x_pos_begin_toggle_epoch_marker x_pos_begin_toggle_epoch_marker],[cfg.plot_MA_offset temp_max_y],'color',[0.25 1 0.125],'LineWidth',2,'parent',axh);
        end
        
        line([x_pos_begin x_pos_begin],[cfg.plot_MA_offset temp_max_y],'color',[0.25 0.125 1],'LineWidth',2,'parent',axh);
        
        
        
        set(cfg.hhyp, 'Name', sprintf('Hypnogram'));
        set(axh, 'box', 'off');
        
        hold(axh,'off');
        
        figure(h)
    end
end
% FW end

figure(h); % ensure that the calling figure is in the front
%hold all;
%set(h, 'Name', sprintf('Datasetnum %d: %s',cfg.datasetnum,cfg.datasetsPath));

if ~isempty(opt.event) && isstruct(opt.event)
    % select only the events in the current time window
    event     = opt.event;
    evtsample = [event(:).sample];
    event     = event(evtsample>=begsample & evtsample<=endsample);
else
    event = [];
end

if isempty(opt.orgdata)
    dat = ft_read_data(cfg.datafile, 'header', opt.hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', strcmp(cfg.continuous, 'no'), 'dataformat', cfg.dataformat, 'headerformat', cfg.headerformat);
else
    dat = ft_fetch_data(opt.orgdata, 'header', opt.hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'allowoverlap', true); % ALLOWING OVERLAPPING TRIALS
end
art = ft_fetch_data(opt.artdata, 'begsample', begsample, 'endsample', endsample);

% apply preprocessing and determine the time axis
[dat, lab, tim] = preproc(dat, opt.hdr.label(chanindx), offset2time(offset, opt.fsample, size(dat,2)), cfg.preproc);

% add NaNs to data for plotting purposes. NaNs will be added when requested horizontal scaling is longer than the data.
nsamplepad = opt.nanpaddata(opt.trlop);
if nsamplepad>0
    dat = [dat NaN(numel(lab), opt.nanpaddata(opt.trlop))];
    tim = [tim linspace(tim(end),tim(end)+nsamplepad*mean(diff(tim)),nsamplepad)];  % possible machine precision error here
end
opt.curdat.label      = lab;
opt.curdat.time{1}    = tim;
opt.curdat.trial{1}   = dat;
opt.curdat.fsample    = opt.fsample;
opt.curdat.sampleinfo = [begsample endsample offset];

% apply scaling to selected channels
% using wildcard to support subselection of channels
if ~isempty(cfg.eegscale)
    chansel = match_str(lab, ft_channelselection('EEG', lab));
    dat(chansel,:) = dat(chansel,:) .* cfg.eegscale;
end
if ~isempty(cfg.eogscale)
    chansel = match_str(lab, ft_channelselection('EOG', lab));
    dat(chansel,:) = dat(chansel,:) .* cfg.eogscale;
end
if ~isempty(cfg.ecgscale)
    chansel = match_str(lab, ft_channelselection('ECG', lab));
    dat(chansel,:) = dat(chansel,:) .* cfg.ecgscale;
end
if ~isempty(cfg.emgscale)
    chansel = match_str(lab, ft_channelselection('EMG', lab));
    dat(chansel,:) = dat(chansel,:) .* cfg.emgscale;
end
if ~isempty(cfg.megscale)
    type = opt.hdr.grad.type;
    chansel = match_str(lab, ft_channelselection('MEG', lab, type));
    dat(chansel,:) = dat(chansel,:) .* cfg.megscale;
end
if ~isempty(cfg.magscale)
    chansel = match_str(lab, ft_channelselection('MEGMAG', lab));
    dat(chansel,:) = dat(chansel,:) .* cfg.magscale;
end
if ~isempty(cfg.gradscale)
    chansel = match_str(lab, ft_channelselection('MEGGRAD', lab));
    dat(chansel,:) = dat(chansel,:) .* cfg.gradscale;
end
if ~isempty(cfg.chanscale)
    chansel = match_str(lab, ft_channelselection(cfg.channel, lab));
    dat(chansel,:) = dat(chansel,:) .* repmat(cfg.chanscale(chanindx),1,size(dat,2));
end
if ~isempty(cfg.mychanscale)
    chansel = match_str(lab, ft_channelselection(cfg.mychan, lab));
    dat(chansel,:) = dat(chansel,:) .* cfg.mychanscale;
end

% to assure current feature is plotted on top
ordervec = 1:length(opt.artdata.label);
ordervec(opt.ftsel) = [];
ordervec(end+1) = opt.ftsel;

% FIXME speedup ft_prepare_layout
if strcmp(cfg.viewmode, 'butterfly')
    laytime = [];
    laytime.label = {'dummy'};
    laytime.pos = [0.5 0.5];
    laytime.width = 1;
    laytime.height = 1;
    opt.laytime = laytime;
else
    % this needs to be reconstructed if the channel selection changes
    tmpcfg = [];
    tmpcfg.layout  = 'vertical';
    tmpcfg.channel = cfg.channel;
    tmpcfg.skipcomnt = 'yes';
    tmpcfg.skipscale = 'yes';
    tmpcfg.showcallinfo = 'no';
    opt.laytime = ft_prepare_layout(tmpcfg, opt.orgdata);
end

% determine the position of the channel/component labels relative to the real axes
% FIXME needs a shift to the left for components
labelx = opt.laytime.pos(:,1) - opt.laytime.width/2 - 0.01;
labely = opt.laytime.pos(:,2);

% determine the total extent of all virtual axes relative to the real axes
ax(1) = min(opt.laytime.pos(:,1) - opt.laytime.width/2);
ax(2) = max(opt.laytime.pos(:,1) + opt.laytime.width/2);
ax(3) = min(opt.laytime.pos(:,2) - opt.laytime.height/2);
ax(4) = max(opt.laytime.pos(:,2) + opt.laytime.height/2);
axis(ax)

% determine a single local axis that encompasses all channels
% this is in relative figure units
opt.hpos   = (ax(1)+ax(2))/2;
opt.vpos   = (ax(3)+ax(4))/2;
opt.width  = ax(2)-ax(1);
opt.height = ax(4)-ax(3);

% these determine the scaling inside the virtual axes
% the hlim will be in seconds, the vlim will be in Tesla or Volt
opt.hlim = [tim(1) tim(end)];
opt.vlim = cfg.ylim;

% FW begin
%hold all
if strcmp(cfg.doSleepScoring,'yes')
    
    delete(findobj(h,'tag', 'mark_spind'));
    delete(findobj(h,'tag', 'mark_slowosci'));
    
    if ~strcmp(cfg.display_power_spectrum,'yes')
        if isfield(cfg,'f_ps')
            if ishandle(cfg.f_ps)
                close(cfg.f_ps);
            end
        end
    end
    
    if ~strcmp(cfg.display_time_frequency,'yes')
        if isfield(cfg,'f_tfr')
            if ishandle(cfg.f_tfr)
                close(cfg.f_tfr);
            end
        end
    end
    
    
    nSamples_data = size(opt.orgdata.trial{1},2);
    FrqOfSmpl = opt.orgdata.fsample;
    
    cfg_buffered_signal_redef = [];
    lengthEpochSamples = endsample - begsample+1;
    buff_begsample = begsample-fix(cfg.nEpochsBuffer*lengthEpochSamples);
    buff_endsample = endsample+fix(cfg.nEpochsBuffer*lengthEpochSamples);
    padd_samples_left = [];
    padd_samples_right = [];
    if buff_begsample < 1
        cfg_buffered_signal_redef.begsample = 1;
        padd_samples_left = repmat(0,1,1 - buff_begsample);
    else
        cfg_buffered_signal_redef.begsample = buff_begsample;
    end
    
    if buff_endsample > nSamples_data
        cfg_buffered_signal_redef.endsample = nSamples_data;
        padd_samples_right = repmat(0,1,buff_endsample - nSamples_data);
    else
        cfg_buffered_signal_redef.endsample = buff_endsample;
    end
    
    if strcmp(cfg.markSO,'yes') || strcmp(cfg.markSpindles,'yes') || strcmp(cfg.underlaySpindleSignal,'yes') || strcmp(cfg.underlayAlphaSignal,'yes') || strcmp(cfg.underlaySOSignal,'yes') || strcmp(cfg.display_power_spectrum,'yes') || strcmp(cfg.display_time_frequency,'yes')
        
        data_det_signal_eeg_data = [];
        
        temp_channel_number_in_curr_display = find(chanindx == cfg.score_channel_eeg_frontal_number);
        for iChanDisplayed = temp_channel_number_in_curr_display
            
            
            cfg_eeg_redef_channel = [];
            cfg_eeg_redef_channel.channel = cfg.score_channel_eeg_frontal_number;
            data_det_signal_eeg_data = ft_redefinetrial(cfg_buffered_signal_redef,ft_selectdata(cfg_eeg_redef_channel,opt.orgdata));
            
            if strcmp(cfg.markSO,'yes')
                
                data_det_signal_eeg_so = [ padd_samples_left data_det_signal_eeg_data.trial{1} padd_samples_right];
                
                data_det_signal_eeg_so = data_det_signal_eeg_so*cfg.SOdetection_orientation;
                
                lengthSignal = length(data_det_signal_eeg_so);
                
                minPeakDistanceSamples_peakdet = fix(FrqOfSmpl*(1/(cfg.so_maxFreq)));
                minPeakDistanceSamples = fix(FrqOfSmpl*(1/(cfg.so_maxFreq)));
                maxPeakDistanceSamples = fix(FrqOfSmpl*(1/cfg.so_minFreq));
                
                [tempCandSignalPeaks, tempCandSignalPeaksSamples] = findpeaks(data_det_signal_eeg_so,'MINPEAKHEIGHT',0,'MINPEAKDISTANCE',minPeakDistanceSamples_peakdet);
                
                %     data_det_signal_eeg_so = sin(0:0.01:3)+sin(0:0.1:30)*0.6
                %     [y x] = findpeaks(data_det_signal_eeg_so,'MINPEAKHEIGHT',0,'MINPEAKDISTANCE',150)
                %
                %     plot(data_det_signal_eeg_so)
                %     hold on
                %     scatter(x,y)
                
                so_cand_pairsamples = [];
                for iCand = 1:(length(tempCandSignalPeaksSamples)-1)
                    for iNextCand = (iCand+1):length(tempCandSignalPeaksSamples);
                        if ((tempCandSignalPeaksSamples(iNextCand) - tempCandSignalPeaksSamples(iCand)) <= maxPeakDistanceSamples) %&& ...
                            %((tempCandSignalPeaksSamples(iNextCand) - tempCandSignalPeaksSamples(iCand)) >= minPeakDistanceSamples)
                            so_cand_pairsamples = [so_cand_pairsamples; [tempCandSignalPeaksSamples(iCand) tempCandSignalPeaksSamples(iNextCand)]];
                        else
                            break;
                        end
                    end
                end
                %    cfg.so_thresholdAmplitudeForDetection = 75;
                
                cfg.so_thresholdAmplitudeForDetection_max = 750;
                so_candidateIndex = [];
                for iCand = 1:size(so_cand_pairsamples,1)
                    samp1 = so_cand_pairsamples(iCand,1);
                    samp2 = so_cand_pairsamples(iCand,2);
                    peak1 = data_det_signal_eeg_so(samp1);
                    peak2 = data_det_signal_eeg_so(samp2);
                    mid_sig = data_det_signal_eeg_so(samp1:samp2);
                    [min_trough_mid min_trough_mid_ind] = min(mid_sig);
                    if (range([min_trough_mid,peak1]) > cfg.so_thresholdAmplitudeForDetection) && ...
                            (range([min_trough_mid,peak2]) > cfg.so_thresholdAmplitudeForDetection) && ...
                            (range([min_trough_mid,peak1]) <= cfg.so_thresholdAmplitudeForDetection_max) && ...
                            (range([min_trough_mid,peak2]) <= cfg.so_thresholdAmplitudeForDetection_max)
                        so_candidateIndex = [so_candidateIndex iCand];
                    end
                end
                so_cand_pairsamples = so_cand_pairsamples(so_candidateIndex,:);
                
                %     so_candidateIndex_include1 = [];
                %     so_candidateIndex_include2 = [];
                %     so_candidateIndex_exclude = [];
                %     for iCand = 1:size(so_cand_pairsamples,1)-1
                %         samp1 = so_cand_pairsamples(iCand,1);
                %         samp2 = so_cand_pairsamples(iCand,2);
                %         samp1_next = so_cand_pairsamples(iCand+1,1);
                %         samp2_next = so_cand_pairsamples(iCand+1,2);
                %         peak1 = data_det_signal_eeg_so(samp1);
                %         peak2 = data_det_signal_eeg_so(samp2);
                %         peak1_next = data_det_signal_eeg_so(samp1_next);
                %         peak2_next = data_det_signal_eeg_so(samp2_next);
                %         if  ( (samp1 == samp1_next) && (samp2 ~= samp2_next))%first peak shared second differnt
                %             if (peak2 >= peak1)
                %                 %so_candidateIndex_exclude = [so_candidateIndex_exclude; iCand+1];
                %                 %continue;
                %             else %(peak2 < peak1)
                %                 temp_mid_sig = data_det_signal_eeg_so(samp1:samp2_next);
                %                 temp_min_sig = min(temp_mid_sig);
                %                 if (peak1 - temp_min_sig)*0.5 > (peak1 - peak2) %peak goes higher than 50% of major peak and might split slow wave
                %                     %so_candidateIndex_exclude = [so_candidateIndex_exclude; iCand];
                %                     so_candidateIndex_exclude = [so_candidateIndex_exclude; iCand+1];
                %                     %continue;
                %                 else
                %                     %so_candidateIndex_include1 = [so_candidateIndex_include1; iCand+1];
                %                 end
                %             end
                %         else
                %             %nothing
                %         end
                %         if (samp1 ~= samp1_next) && (samp2 == samp2_next) %first peak different second shared
                %             if (peak1_next >= peak1)
                %                 so_candidateIndex_exclude = [so_candidateIndex_exclude; iCand];
                %                 continue;
                %             else %(peak1_next < peak1)
                %                 temp_mid_sig = data_det_signal_eeg_so(samp1:samp2_next);
                %                 temp_min_sig = min(temp_mid_sig);
                %                 if range([peak2_next temp_min_sig])*0.5 > range([peak1_next peak2_next]) %peak goes higher than 50% of major peak and might split slow wave
                %                     so_candidateIndex_exclude = [so_candidateIndex_exclude; iCand];
                %                     %so_candidateIndex_exclude = [so_candidateIndex_exclude; iCand+1];
                %                     continue;
                %                 else
                %                     %so_candidateIndex_include2 = [so_candidateIndex_include2; iCand];
                %                 end
                %             end
                %         end
                %
                %     end
                %     so_cand_pairsamples = so_cand_pairsamples(unique([so_candidateIndex_include1; so_candidateIndex_include2; setdiff(1:size(so_cand_pairsamples,1),so_candidateIndex_exclude)]),:);
                
                
                
                iCand = 1;
                while iCand < size(so_cand_pairsamples,1)
                    samp1 = so_cand_pairsamples(iCand,1);
                    samp2 = so_cand_pairsamples(iCand,2);
                    peak1 = data_det_signal_eeg_so(samp1);
                    peak2 = data_det_signal_eeg_so(samp2);
                    iCandNext = iCand+1;
                    samp1_next = so_cand_pairsamples(iCandNext,1);
                    samp2_next = so_cand_pairsamples(iCandNext,2);
                    peak1_next = data_det_signal_eeg_so(samp1_next);
                    peak2_next = data_det_signal_eeg_so(samp2_next);
                    
                    while (iCandNext <= size(so_cand_pairsamples,1)) && ( (samp1 == samp1_next) && (samp2 ~= samp2_next))
                        
                        samp1_next = so_cand_pairsamples(iCandNext,1);
                        samp2_next = so_cand_pairsamples(iCandNext,2);
                        peak1_next = data_det_signal_eeg_so(samp1_next);
                        peak2_next = data_det_signal_eeg_so(samp2_next);
                        
                        if  ( (samp1 == samp1_next) && (samp2 ~= samp2_next))%first peak shared second differnt
                            if (peak2 >= peak1)
                                %so_candidateIndex_include1 = [so_candidateIndex_include1 iCand];
                                %so_candidateIndex_exclude = [so_candidateIndex_exclude; iCand+1];
                                %continue;
                                so_cand_pairsamples(iCandNext,:) = [];
                                %iCandNext = iCandNext + 1;
                            else %(peak2 < peak1)
                                temp_mid_sig = data_det_signal_eeg_so(samp1:samp2_next);
                                temp_min_sig = min(temp_mid_sig);
                                if (peak1 - temp_min_sig)*0.5 < (peak1 - peak2) %peak goes higher than 50% of major peak and might split slow wave
                                    %so_candidateIndex_exclude = [so_candidateIndex_exclude; iCand];
                                    %so_candidateIndex_exclude = [so_candidateIndex_exclude; iCandNext];
                                    %iCand = iCand + 2;
                                    so_cand_pairsamples(iCandNext,:) = [];
                                    %continue;
                                else
                                    iCandNext = iCandNext + 1;
                                    %so_candidateIndex_include1 = [so_candidateIndex_include1; iCand+1];
                                end
                            end
                        else
                            break
                        end
                    end
                    
                    if iCand < size(so_cand_pairsamples,1)
                        samp1 = so_cand_pairsamples(iCand,1);
                        samp2 = so_cand_pairsamples(iCand,2);
                        peak1 = data_det_signal_eeg_so(samp1);
                        peak2 = data_det_signal_eeg_so(samp2);
                        iCandNext = iCand + 1;
                        samp1_next = so_cand_pairsamples(iCandNext,1);
                        samp2_next = so_cand_pairsamples(iCandNext,2);
                        peak1_next = data_det_signal_eeg_so(samp1_next);
                        peak2_next = data_det_signal_eeg_so(samp2_next);
                        while (iCandNext <= size(so_cand_pairsamples,1)) && ( (samp1 ~= samp1_next) && (samp2 == samp2_next))
                            
                            samp1 = so_cand_pairsamples(iCand,1);
                            samp2 = so_cand_pairsamples(iCand,2);
                            peak1 = data_det_signal_eeg_so(samp1);
                            peak2 = data_det_signal_eeg_so(samp2);
                            
                            if (samp1 ~= samp1_next) && (samp2 == samp2_next) %first peak different second shared
                                if (peak1_next >= peak1)
                                    %so_candidateIndex_exclude = [so_candidateIndex_exclude; iCand];
                                    so_cand_pairsamples(iCand,:) = [];
                                else %(peak1_next < peak1)
                                    temp_mid_sig = data_det_signal_eeg_so(samp1:samp2);
                                    temp_min_sig = min(temp_mid_sig);
                                    if (peak2 - temp_min_sig)*0.5 < (peak2 - peak1_next) %peak goes higher than 50% of major peak and might split slow wave
                                        %so_candidateIndex_exclude = [so_candidateIndex_exclude; iCand];
                                        %so_candidateIndex_exclude = [so_candidateIndex_exclude; iCand+1];
                                        so_cand_pairsamples(iCand,:) = [];
                                        %continue;
                                    else
                                        iCand = iCand + 1;
                                        %so_candidateIndex_include2 = [so_candidateIndex_include2; iCand];
                                    end
                                end
                            else
                                %iCand = iCand + 1;
                                break;
                            end
                        end
                    else
                        break;
                    end
                    iCand = iCand + 1;
                end
                
                iCand = 1;
                while iCand < size(so_cand_pairsamples,1)
                    iCandNext = iCand + 1;
                    samp1 = so_cand_pairsamples(iCand,1);
                    samp2 = so_cand_pairsamples(iCand,2);
                    peak1 = data_det_signal_eeg_so(samp1);
                    peak2 = data_det_signal_eeg_so(samp2);
                    samp1_next = so_cand_pairsamples(iCandNext,1);
                    samp2_next = so_cand_pairsamples(iCandNext,2);
                    peak1_next = data_det_signal_eeg_so(samp1_next);
                    peak2_next = data_det_signal_eeg_so(samp2_next);
                    while (iCandNext <= size(so_cand_pairsamples,1)) && (samp1_next <= samp2)
                        
                        if (samp1_next <  samp2) && (samp2_next >=  samp2)
                            so_cand_pairsamples(iCand,1) = samp1_next;
                            so_cand_pairsamples(iCand,2) = samp2;
                            so_cand_pairsamples(iCandNext,:) = [];
                        else
                            iCandNext = iCandNext + 1;
                        end
                        if iCandNext < size(so_cand_pairsamples,1)
                            samp1_next = so_cand_pairsamples(iCandNext,1);
                            samp2_next = so_cand_pairsamples(iCandNext,2);
                            peak1_next = data_det_signal_eeg_so(samp1_next);
                            peak2_next = data_det_signal_eeg_so(samp2_next);
                        end
                        
                        
                    end
                    iCand = iCand + 1;
                end
                
                
                so_display_ind = zeros(1,lengthSignal);
                for iCand = 1:size(so_cand_pairsamples,1)
                    so_display_ind(so_cand_pairsamples(iCand,1):so_cand_pairsamples(iCand,2)) = 1;
                end
                
                tmp = diff([0 so_display_ind(fix(cfg.nEpochsBuffer*lengthEpochSamples)+1:fix((cfg.nEpochsBuffer+1)*lengthEpochSamples)) 0]);
                evbeg = find(tmp==+1);
                evend = find(tmp==-1) - 1;
                
                if cfg.SOdetection_orientation == 1
                    temp_so_y_offset1 = 0.9;
                    temp_so_y_offset2 = 1;
                else
                    temp_so_y_offset1 = -1;
                    temp_so_y_offset2 = -0.9;
                end
                
                %init first axis
                temp_ax = [];
                
                h_so_event_begin_end = temp_ax;
                for k=1:numel(evbeg)
                    
                    h_so_event_begin_end = ft_plot_box([tim(evbeg(k)) tim(evend(k)) temp_so_y_offset1 temp_so_y_offset2],'facealpha',0.38, 'facecolor', cfg.slowoscillation_mark_color, 'edgecolor', 'none', 'tag', 'mark_slowosci',  ...
                        'hpos', opt.laytime.pos(iChanDisplayed,1), 'vpos', opt.laytime.pos(iChanDisplayed,2), 'width', opt.width, 'height', opt.laytime.height(iChanDisplayed), 'hlim', opt.hlim, 'vlim', [-1 1],'axis', h_so_event_begin_end);
                    
                end
                temp_ax = h_so_event_begin_end;
                
                so_points_for_display = so_cand_pairsamples(:);
                so_points_for_display = so_points_for_display(((fix(cfg.nEpochsBuffer*lengthEpochSamples)+1) <= so_points_for_display) & ...
                    (so_points_for_display <= fix(cfg.nEpochsBuffer+1)*lengthEpochSamples)) - fix(cfg.nEpochsBuffer*lengthEpochSamples);
                
                h_so_points_begin_end = ft_plot_vector([tim(so_points_for_display)' tim(so_points_for_display)'],[repmat(temp_so_y_offset1,length(so_points_for_display),1) repmat(temp_so_y_offset2,length(so_points_for_display),1)], 'color', cfg.slowoscillation_mark_color, 'linewidth', 1, 'tag', 'mark_slowosci',  ...
                    'hpos', opt.laytime.pos(iChanDisplayed,1), 'vpos', opt.laytime.pos(iChanDisplayed,2), 'width', opt.width, 'height', opt.laytime.height(iChanDisplayed), 'hlim', opt.hlim, 'vlim', [-1 1]);
                
                
                so_points_for_display_begins = [];
                so_points_for_display_ends = [];
                so_points_for_display_pairs = [];
                if ~isempty(so_cand_pairsamples)
                    so_points_for_display_begins = so_cand_pairsamples(((fix(cfg.nEpochsBuffer*lengthEpochSamples+1) <= so_cand_pairsamples(:,1)) & ...
                        (so_cand_pairsamples(:,1) <= fix((cfg.nEpochsBuffer+1)*lengthEpochSamples))),:) - fix(cfg.nEpochsBuffer*lengthEpochSamples);
                    so_points_for_display_ends = so_cand_pairsamples(((fix(cfg.nEpochsBuffer*lengthEpochSamples+1) <= so_cand_pairsamples(:,2)) & ...
                        (so_cand_pairsamples(:,2) <= fix((cfg.nEpochsBuffer+1)*lengthEpochSamples))),:) - fix(cfg.nEpochsBuffer*lengthEpochSamples);
                    so_points_for_display_pairs = so_cand_pairsamples((((cfg.nEpochsBuffer*lengthEpochSamples+1) <= so_cand_pairsamples(:,1)) & ...
                        (so_cand_pairsamples(:,2) <= fix((cfg.nEpochsBuffer+1)*lengthEpochSamples))),:) - fix(cfg.nEpochsBuffer*lengthEpochSamples);
                end
                %                 temp_prevbeg = -1;
                %                 temp_count_overlapp = 0;
                %                 for k=1:size(so_points_for_display_pairs,1)
                %                     if (temp_prevbeg == so_points_for_display_pairs(k,1))
                %                         temp_count_overlapp = temp_count_overlapp + 1;
                %                     else
                %                         temp_count_overlapp = 0;
                %                     end
                %                     h_so_event_begin_end = ft_plot_box([tim(so_points_for_display_pairs(k,1)) tim(so_points_for_display_pairs(k,2)) (0.8-temp_count_overlapp*0.1) (0.9-temp_count_overlapp*0.1)],'facealpha',0.6, 'facecolor', [0 0 1], 'edgecolor', 'none', 'tag', 'mark_slowosci',  ...
                %                         'hpos', opt.laytime.pos(iChanDisplayed,1), 'vpos', opt.laytime.pos(iChanDisplayed,2), 'width', opt.width, 'height', opt.laytime.height(iChanDisplayed), 'hlim', opt.hlim, 'vlim', [-1 1]);
                %                     h_so_event_begin_end = ft_plot_box([tim(so_points_for_display_pairs(k,2))-0.1 tim(so_points_for_display_pairs(k,2)) (0.8-temp_count_overlapp*0.1) (0.9-temp_count_overlapp*0.1)],'facealpha',0.6, 'facecolor', [1 0 0], 'edgecolor', 'none', 'tag', 'mark_slowosci',  ...
                %                         'hpos', opt.laytime.pos(iChanDisplayed,1), 'vpos', opt.laytime.pos(iChanDisplayed,2), 'width', opt.width, 'height', opt.laytime.height(iChanDisplayed), 'hlim', opt.hlim, 'vlim', [-1 1]);
                %                     temp_prevbeg = so_points_for_display_pairs(k,1);
                %                 end
                
                temp_epochLengthSamples = endsample - begsample + 1;
                cfg.curr_displayed_detected_slowosci_perc_display_ind = so_display_ind(fix(cfg.nEpochsBuffer*lengthEpochSamples+1):fix((cfg.nEpochsBuffer+1)*lengthEpochSamples));
                cfg.curr_displayed_detected_slowosci_perc = sum(cfg.curr_displayed_detected_slowosci_perc_display_ind )/temp_epochLengthSamples;
                
                
                cfg.curr_displayed_detected_slowosci_number = ( size(so_points_for_display_pairs,1) + 0.5*abs(size(so_points_for_display_begins,1)-size(so_points_for_display_pairs,1)) + 0.5*abs(size(so_points_for_display_ends,1)-size(so_points_for_display_pairs,1)) );
                
            end
            
            if strcmp(cfg.underlaySOSignal,'yes')
                
                cfg_eeg_so = [];
                cfg_eeg_so = cfg.core_cfg;
                cfg_eeg_so.hpfilter = 'yes';
                %cfg_eeg_so.hpfilterdesign = cfg.so_hpfilterdesign;
                
                if strcmp(cfg.UseFixedFilterOrder_hp,'yes')
                    cfg_eeg_so.hpfiltord     = cfg.FilterOrder_hp;
                end
                cfg_eeg_so.hpfreq        = cfg.so_filter_minFreq;%dummy values are overwritten by low level function
                cfg_eeg_so.feedback = 'no';
                
                data_det_signal_eeg_so_disp1 = st_preprocessing(cfg_eeg_so,data_det_signal_eeg_data);
                
                
                
                
                cfg_eeg_so2 = [];
                cfg_eeg_so2 = cfg.core_cfg;
                cfg_eeg_so2.lpfilter = 'yes';
                %cfg_eeg_so2.lpfilterdesign = cfg.so_lpfilterdesign;
                if strcmp(cfg.UseFixedFilterOrder_lp,'yes')
                    cfg_eeg_so2.lpfiltord     = cfg.FilterOrder_lp;
                end
                cfg_eeg_so2.lpfreq        = cfg.so_filter_maxFreq;%dummy values are overwritten by low level function
                data_det_signal_eeg_so_disp2 = st_preprocessing(cfg_eeg_so2,data_det_signal_eeg_so_disp1);
                
                
                
                
                SOdispsignal = [ padd_samples_left data_det_signal_eeg_so_disp2.trial{1} padd_samples_right];
                
                cfg.so_signal_display = SOdispsignal(fix(cfg.nEpochsBuffer*lengthEpochSamples+1):fix((cfg.nEpochsBuffer+1)*lengthEpochSamples));
                
            end
        end
        
        
        temp_channel_number_in_curr_display = find(chanindx == cfg.score_channel_eeg_number);
        for iChanDisplayed = temp_channel_number_in_curr_display
            
            
            
            if (cfg.score_channel_eeg_frontal_number ~= cfg.score_channel_eeg_number) || isempty(data_det_signal_eeg_data)
                cfg_eeg_redef_channel = [];
                cfg_eeg_redef_channel.channel = cfg.score_channel_eeg_number;
                data_det_signal_eeg_data = ft_redefinetrial(cfg_buffered_signal_redef,ft_selectdata(cfg_eeg_redef_channel,opt.orgdata));
            end
            
            
            %%%%%% Time-frequency #####
            
            % if     isfield(cfg,'display_power_spectrum') && isfield(cfg,'display_time_frequency')
            if strcmp(cfg.display_power_spectrum,'yes') || strcmp(cfg.display_time_frequency,'yes')
                
                
                data_det_signal_eeg_data_tfr = data_det_signal_eeg_data;
                
                minFreq = 0.5;
                maxFreq = 30;
                FreqSteps = 0.5;
                TimeSteps = 0.1;
                Xtick = fix(minFreq):2:fix(maxFreq);
                cfg_tfr = [];
                cfg_tfr.method    = 'wavelet';%single number (in unit of time, typically seconds) of the required snippets
                cfg_tfr.output   = 'pow';%single number (between 0 and 1 (exclusive)) specifying the fraction of overlap between snippets (0 = no overlap)
                cfg_tfr.foi = [minFreq:FreqSteps:maxFreq];%
                cfg_tfr.width = 4;%7
                cfg_tfr.pad = 'maxperlen';
                cfg_tfr.feedback = 'no';
                cfg_tfr.keeptrials = 'no';
                cfg_tfr.toi = [min(cellfun(@min,data_det_signal_eeg_data_tfr.time)):TimeSteps:max(cellfun(@max,data_det_signal_eeg_data_tfr.time))];
                data_tfr = ft_freqanalysis(cfg_tfr,data_det_signal_eeg_data_tfr);
                
                
                if strcmp(cfg.display_power_spectrum,'yes')
                    
                    
                    if isfield(cfg,'f_ps')
                        if ishandle(cfg.f_ps)
                            figure(cfg.f_ps);
                        else
                            cfg.f_ps = figure;
                            figure(cfg.f_ps)
                        end
                    else
                        cfg.f_ps = figure;
                        figure(cfg.f_ps)
                    end
                    
                    
                    if isfield(cfg, 'f_ps_gca')
                        if ishandle(cfg.f_ps_gca)
                            delete(cfg.f_ps_gca);
                        end
                    end
                    powerspectrum = nanmean(squeeze(data_tfr.powspctrm),2);
                    freq = data_tfr.freq;
                    bars_y = 0;
                    for iFreq =  1:size(cfg.freq_borders,1)
                        temp_freq_border_left = cfg.freq_borders(iFreq,1);
                        temp_freq_border_right = cfg.freq_borders(iFreq,2);
                        temp_mean_power = 10*log10(nanmean(powerspectrum((temp_freq_border_left <= freq) & (freq <= temp_freq_border_right))));
                        rectangle('Position',[temp_freq_border_left,bars_y,temp_freq_border_right-temp_freq_border_left,temp_mean_power],'FaceColor',cfg.freq_colors(iFreq,:),'EdgeColor','k','LineWidth',1)
                        if iFreq == 1 && ~ishold
                            hold all
                        end
                    end
                    plot(freq,10*log10(powerspectrum),'color','k','LineWidth',2);
                    ylabel('Power [dB]');
                    xlabel('Frequency [Hz]');
                    chname = data_tfr.label{1};
                    cfg.f_ps_gca = gca;
                    title(cfg.f_ps_gca,['EEG log power Spectrum (' chname ')' ],'interpreter','none');
                    set(cfg.f_ps_gca, 'TickDir', 'out','Xtick', Xtick);
                    set(cfg.f_ps, 'Name', 'Power Spectrum');
                    hold off
                    
                    
                end
                
                if strcmp(cfg.display_time_frequency,'yes')
                    
                    if isfield(cfg,'f_tfr')
                        if ishandle(cfg.f_tfr)
                            figure(cfg.f_tfr);
                        else
                            cfg.f_tfr = figure;
                            figure(cfg.f_tfr)
                            %set(cfg.f_tfr, 'WindowButtonDownFcn',   {@select_sleep_stage_cb, 'h_main',h});
                            %set(cfg.f_tfr, 'NumberTitle', 'off');
                            cfg.f_tfr_gca = gca;
                            %cfg.hhypfig = gcf;
                        end
                    else
                        cfg.f_tfr = figure;
                        figure(cfg.f_tfr)
                        %set(cfg.f_tfr, 'WindowButtonDownFcn',   {@select_sleep_stage_cb, 'h_main',h});
                        %set(cfg.f_tfr, 'NumberTitle', 'off');
                        cfg.f_tfr_gca = gca;
                        %cfg.hhypfig = gcf;
                    end
                    
                    
                    
                    
                    
                    data_tfr.powspctrm(isnan(data_tfr.powspctrm(:))) = 10E-12;
                    
                    time_bins_epoch = (lengthEpochSamples/data_det_signal_eeg_data_tfr.fsample)/TimeSteps;
                    time_bins_to_substract_left = (fix(cfg.nEpochsBuffer*lengthEpochSamples)/data_det_signal_eeg_data_tfr.fsample)/TimeSteps;
                    time_bins_to_substract_right = time_bins_to_substract_left;
                    
                    
                    if ~isempty(padd_samples_left)
                        temp_index_display_tfr_time_points = 1:time_bins_epoch;
                    else%if ~isempty(padd_samples_right)
                        temp_index_display_tfr_time_points = (time_bins_to_substract_left+1):(time_bins_to_substract_left+time_bins_epoch);
                        %             else
                        %                 temp_index_display_tfr_time_points = (time_bins_to_substract_left+1):(time_bins_to_substract_left+time_bins_to_substract_right+time_bins_epoch);
                    end
                    
                    temp_index_display_tfr_time_points = fix(temp_index_display_tfr_time_points);
                    
                    
                    
                    
                    
                    temp_curr_tfr_channel_signal_ylim = cfg.ylim/cfg.chanscale(cfg.score_channel_eeg_number);
                    temp_curr_tfr_channel_signal_ylim = temp_curr_tfr_channel_signal_ylim*2;
                    
                    Ysteps = (max(temp_curr_tfr_channel_signal_ylim)-min(temp_curr_tfr_channel_signal_ylim))/5;
                    
                    Ytick1 = [minFreq 4:4:maxFreq];
                    Ytick2 = [min(temp_curr_tfr_channel_signal_ylim):Ysteps:max(temp_curr_tfr_channel_signal_ylim)];
                    Ztick = -8:4:8;
                    cfg_tfr = [];
                    temp_time_interval_display = [min(data_tfr.time(temp_index_display_tfr_time_points)) max(data_tfr.time(temp_index_display_tfr_time_points))];
                    cfg_tfr.baseline     = temp_time_interval_display;%normalized with reference to average power in +-0.9 s interval
                    cfg_tfr.baselinetype = 'db'; % power in dB = 10*log_10(pwr/mean) 10*log10(data ./ meanVals);
                    cfg_tfr.zlim         = [min(Ztick) max(Ztick)];
                    cfg_tfr.xlim         = temp_time_interval_display;
                    cfg_tfr.ylim         = [minFreq maxFreq];
                    
                    x2 = data_det_signal_eeg_data_tfr.time{:}; %(-spindle_trough_prestim_actual <= timelock.events.time & timelock.events.time <= spindle_trough_poststim_actual);
                    cfg_tfr.x2range = temp_time_interval_display;
                    
                    
                    y2 = [];
                    y2(1,:) = data_det_signal_eeg_data_tfr.trial{1};
                    %y2(2,:) = timelock.nonevents.avg * 1000000;
                    
                    cfg_tfr.y2range = [min(Ytick2) max(Ytick2)];
                    cfg_tfr.y2label = 'signal units';
                    cfg_tfr.y2colors = [[0 0 0]];%color for second y axis default 'b'
                    cfg_tfr.y2linestyles = ['-'];%linestylw for second y axis default '-'
                    cfg_tfr.y2linewidths = [1.5];%lable for second y axis default 1
                    cfg_tfr.y2alphas = [0.55];%lable for second y axis default 1
                    
                    %       cfg_tfr.y2colors = [[1 1 1]; [.8 .8 .8]];%color for second y axis default 'b'
                    %       cfg_tfr.y2linestyles = ['-';'-'];%linestylw for second y axis default '-'
                    %       cfg_tfr.y2linewidths = [3 ; 1];%lable for second y axis default 1
                    %       cfg_tfr.y2alphas = [0.55 ; 0.55];%lable for second y axis default 1
                    
                    cfg_tfr.colormap = jet(128);%individual_color_map_insertion(min(Ztick),max(Ztick),{stat.critval(1), stat.critval(2)},[1 1 1],jet(256)); %excludes a little bit more t-values due to imprecicion errors
                    %cfg.colormap = cfg.colormap(end:-1:1,:);
                    cfg_tfr.interactive = 'no';
                    
                    if isfield(cfg, 'f_tfr_p1gca')
                        if ishandle(cfg.f_tfr_p1gca)
                            delete(cfg.f_tfr_p1gca);
                        end
                        if ishandle(cfg.f_tfr_p2gca)
                            delete(cfg.f_tfr_p2gca);
                        end
                        %delete(cfg.f_tfr_p1gccb);
                    end
                    figure(cfg.f_tfr)
                    
                    %set(0,'DefaultFigureVisible','off');
                    %cfg_tfr.colorbar = 'no';
                    [tempcfg p1gca p2gca p1gccb] = ft_fw_singleplotTFR_yy(cfg_tfr, data_tfr,x2,y2);%ft_singleplotTFR edit ft_freqbaseline
                    %                         tfr_fig = gcf;
                    %                         %set(0,'DefaultFigureVisible','on');
                    %                         tt = figure
                    %                         s1 = copyobj(p1gca,tt)
                    %                         s2 = copyobj(p2gca,tt)
                    %                         close(tfr_fig)
                    set(p1gccb,'location','northoutside')
                    set(p1gca,'xlim',cfg_tfr.x2range)
                    set(p1gca,'ylim',cfg_tfr.ylim)
                    pos = get(p1gca, 'Position');
                    pos(1) = 0.055;
                    pos(3) = 0.9;
                    set(p1gca, 'Position', pos)
                    set(p2gca, 'Position', pos)
                    
                    cfg.f_tfr_p1gca = p1gca;
                    cfg.f_tfr_p2gca = p2gca;
                    cfg.f_tfr_p1gccb = p1gccb;
                    chname = data_tfr.label{1};
                    title(cfg.f_tfr_p1gca,['EEG log power (' chname ') normalized to average in time window' ]);
                    xlabel(cfg.f_tfr_p1gca,'Time [sec]');
                    ylabel(cfg.f_tfr_p1gca,'Frequency [Hz]');
                    nticks = 11;
                    xTick = round(linspace(temp_time_interval_display(1), temp_time_interval_display(2), nticks));
                    temp_fontsize = 10;
                    set(cfg.f_tfr_p1gca,'Fontsize',temp_fontsize,'FontUnits','normalized', 'TickDir', 'out', 'Xtick', xTick, 'Ytick', Ytick1);
                    set(cfg.f_tfr_p1gccb ,'Fontsize',temp_fontsize, 'TickDir','out','Ytick',Ztick);
                    ylabel(cfg.f_tfr_p1gccb ,'Power/mean(power) [dB]');
                    set(cfg.f_tfr_p2gca,'Fontsize',temp_fontsize,'FontUnits','normalized', 'TickDir', 'out','Ytick', Ytick2);
                    set(cfg.f_tfr, 'Name', 'Time-frequency of EEG, normalized');
                    
                    %myaa(4)
                    
                    %figure_width = 12;    % Width in inches
                    %figure_height = 5;    % Height in inches
                    
                    %pos = get(cfg.f_tfr, 'Position');
                    %set(cfg.f_tfr, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
                    
                end
                
                figure(h); % ensure that the calling figure is in the front
            end
            %end
            
            
            %%%%%% Time-frequncy end #####
            
            
            %cfg_emg_redef_channel = [];
            %cfg_emg_redef_channel.channel = cfg.score_channel_emg_number;
            %data_det_signal_emg_data = ft_selectdata(cfg_emg_redef_channel,ft_redefinetrial(cfg_emg_redef,opt.orgdata));
            
            %data_det_signal_emg_artifact =
            
            
            if strcmp(cfg.markSpindles,'yes') || strcmp(cfg.underlaySpindleSignal,'yes')
                
                EnvelopeMethod = 'hilbertEnv';
                smplsMinDetectionLength = fix(opt.orgdata.fsample*cfg.sp_minSec);
                smplsMaxDetectionLength = fix(opt.orgdata.fsample*cfg.sp_maxSec);
                
                cfg_eeg_sp = [];
                cfg_eeg_sp = cfg.core_cfg;
                cfg_eeg_sp.bpfilter = 'yes';
                
                %cfg_eeg_sp.bpfilterdesign = cfg.sp_bpfilterdesign;
                
                if strcmp(cfg.UseFixedFilterOrder_bp,'yes')
                    cfg_eeg_sp.bpfiltord     = cfg.FilterOrder_bp;
                end
                cfg_eeg_sp.bpfreq        = [cfg.sp_minFreq cfg.sp_maxFreq];%dummy values are overwritten by low level function
                cfg_eeg_sp.feedback = 'no';
                
                data_det_signal_eeg_sp = st_preprocessing(cfg_eeg_sp,data_det_signal_eeg_data);
                
                frqBndPssSignal = [ padd_samples_left data_det_signal_eeg_sp.trial{1} padd_samples_right];
                
                cfg.spindsignal_display = frqBndPssSignal(fix(cfg.nEpochsBuffer*lengthEpochSamples+1):fix((cfg.nEpochsBuffer+1)*lengthEpochSamples));
                
                if strcmp(cfg.markSpindles,'yes')
                    
                    lengthSignal = length(frqBndPssSignal);
                    
                    smplsRMSTimeWndw = 0.2*FrqOfSmpl;
                    smplsMovAvgTimeWndw = 0.2*FrqOfSmpl;
                    envelope = [];
                    if strcmp(EnvelopeMethod,'hilbertEnv')
                        envelope = abs(hilbert(frqBndPssSignal))';
                        %     elseif strcmp(EnvelopeMethod,'smoothedRMSwd')
                        %         envelope = smoothRMSwd(frqBndPssSignal,smplsRMSTimeWndw);
                        %         if exist('smooth','file') == 2
                        %             envelope = smooth(envelope,smplsMovAvgTimeWndw);
                        %         else
                        %             envelope = smoothwd(envelope,smplsMovAvgTimeWndw)';
                        %         end
                    end
                    
                    cfg.spindsignal_envelope_display = envelope(fix(cfg.nEpochsBuffer*lengthEpochSamples+1):fix((cfg.nEpochsBuffer+1)*lengthEpochSamples));
                    
                    [begins, ends] = getBeginsAndCorrespondingEndsIndicesAboveThreshold(envelope,cfg.sp_thresholdForDetectionBeginEnd/2);
                    ind_valid_criterion = [];
                    sp_events = [];
                    for iBeg = 1:numel(begins)
                        envmax = max(envelope(begins(iBeg):ends(iBeg)));
                        if ((cfg.sp_thresholdForDetectionCriterion/2) <= envmax)
                            ind_valid_criterion = [ind_valid_criterion;iBeg];
                            sp_events = [sp_events; begins(iBeg) + find(envelope(begins(iBeg):ends(iBeg)) == envmax,1,'first') - 1];
                        end
                    end
                    
                    begins = begins(ind_valid_criterion);
                    ends = ends(ind_valid_criterion);
                    
                    tempCandidatesLengths = ends - begins + 1;
                    
                    indicesCandiates = find((tempCandidatesLengths >= smplsMinDetectionLength) & (tempCandidatesLengths <= smplsMaxDetectionLength));
                    begins = begins(indicesCandiates);
                    ends = ends(indicesCandiates);
                    sp_events = sp_events(indicesCandiates);
                    curr_nDetected_spindels = length(indicesCandiates);
                    
                    
                    
                    sp_events_for_display = sp_events((fix(cfg.nEpochsBuffer*lengthEpochSamples+1) <= sp_events) & (sp_events <= fix((cfg.nEpochsBuffer+1)*lengthEpochSamples))) - fix(cfg.nEpochsBuffer*lengthEpochSamples);
                    
                    cfg.curr_displayed_detected_spindels_number = numel(sp_events_for_display);
                    
                    if cfg.SOdetection_orientation == 1
                        temp_sp_y_offset1 = -1;
                        temp_sp_y_offset2 = -0.9;
                    else
                        temp_sp_y_offset1 = 0.9;
                        temp_sp_y_offset2 = 1;
                    end
                    
                    temp_ax = [];
                    if curr_nDetected_spindels > 0
                        sp_display_ind = zeros(1,lengthSignal);
                        for iBeg = 1:numel(begins)
                            sp_display_ind(begins(iBeg):ends(iBeg)) = 1;
                        end
                        
                        
                        
                        tmp = diff([0 sp_display_ind(fix(cfg.nEpochsBuffer*lengthEpochSamples+1):fix((cfg.nEpochsBuffer+1)*lengthEpochSamples)) 0]);
                        evbeg = find(tmp==+1);
                        evend = find(tmp==-1) - 1;
                        
                        h_sp_event_begin_end = temp_ax;
                        for k=1:numel(evbeg)
                            h_sp_event_begin_end = ft_plot_box([tim(evbeg(k)) tim(evend(k)) temp_sp_y_offset1 temp_sp_y_offset2],'facealpha',0.38, 'facecolor', cfg.spindle_mark_color, 'edgecolor', 'none', 'tag', 'mark_spind',  ...
                                'hpos', opt.laytime.pos(iChanDisplayed,1), 'vpos', opt.laytime.pos(iChanDisplayed,2), 'width', opt.width, 'height', opt.laytime.height(iChanDisplayed), 'hlim', opt.hlim, 'vlim', [-1 1],'axis',h_sp_event_begin_end);
                        end
                        temp_ax = h_sp_event_begin_end;
                        
                    end
                    h_sp_event = temp_ax;
                    for iEv = 1:numel(sp_events_for_display)
                        h_sp_event = ft_plot_line([tim(sp_events_for_display(iEv)) tim(sp_events_for_display(iEv))], [temp_sp_y_offset1 temp_sp_y_offset2],'facealpha',0.8, 'color', cfg.spindle_mark_color, 'edgecolor', 'none', 'tag', 'mark_spind','linewidth',2,  ...
                            'hpos', opt.laytime.pos(iChanDisplayed,1), 'vpos', opt.laytime.pos(iChanDisplayed,2), 'width', opt.width, 'height', opt.laytime.height(iChanDisplayed), 'hlim', opt.hlim, 'vlim', [-1 1], 'axis', h_sp_event);
                        
                    end
                    temp_ax = h_sp_event;
                    
                end
            end
            
            
            
            if strcmp(cfg.underlayAlphaSignal,'yes')
                
                cfg_eeg_al = [];
                cfg_eeg_al = cfg.core_cfg;
                cfg_eeg_al.bpfilter = 'yes';
                
                %cfg_eeg_sp.bpfilterdesign = cfg.sp_bpfilterdesign;
                
                if strcmp(cfg.UseFixedFilterOrder_bp,'yes')
                    cfg_eeg_al.bpfiltord     = cfg.FilterOrder_bp;
                end
                cfg_eeg_al.bpfreq        = [cfg.al_minFreq cfg.al_maxFreq];%dummy values are overwritten by low level function
                cfg_eeg_al.feedback = 'no';
                
                data_det_signal_eeg_al = st_preprocessing(cfg_eeg_al,data_det_signal_eeg_data);
                
                frqBndPssSignal = [ padd_samples_left data_det_signal_eeg_al.trial{1} padd_samples_right];
                
                cfg.alphasignal_display = frqBndPssSignal(fix(cfg.nEpochsBuffer*lengthEpochSamples+1):fix((cfg.nEpochsBuffer+1)*lengthEpochSamples));
                
            end
            
            
            
        end
        
        
    end
    
    if cfg.has_ECG
        if strcmp(cfg.markECG,'yes')
            
            temp_channel_number_in_curr_display = find(chanindx == cfg.score_channel_ecg_number);
            for iChanDisplayed = temp_channel_number_in_curr_display
                
                cfg_ecg_redef_channel = [];
                cfg_ecg_redef_channel.channel = cfg.score_channel_ecg_number;
                data_det_signal_ecg_data = ft_redefinetrial(cfg_buffered_signal_redef,ft_selectdata(cfg_ecg_redef_channel,opt.orgdata));
                
                
                cfg_ecg_peak = [];
                cfg_ecg_peak = cfg.core_cfg;
                cfg_ecg_peak.bpfilter = 'yes';
                %cfg_ecg_peak.bpfilterdesign = cfg.hr_bpfilterdesign;
                
                if strcmp(cfg.UseFixedFilterOrder_bp,'yes')
                    cfg_ecg_peak.bpfiltord     = cfg.FilterOrder_bp;
                end
                cfg_ecg_peak.bpfreq        = [cfg.ecg_peak_filter_minFreq cfg.ecg_peak_filter_maxFreq];%dummy values are overwritten by low level function
                cfg_ecg_peak.feedback = 'no';
                
                %cfg_ecg_peak.hilbert = 'abs';
                data_det_signal_ecg_peak = st_preprocessing(cfg_ecg_peak,data_det_signal_ecg_data);
                
                if (cfg.ECG_signalMultiplicator ~= 1)
                    data_det_signal_ecg_peak = ft_fw_factorMultiplicationOnSignal(data_det_signal_ecg_peak,'trial',cfg.ECG_signalMultiplicator);
                end
                
                ECGpeakfrqBndPssSignal_hilbert = abs(hilbert([ padd_samples_left data_det_signal_ecg_peak.trial{1} padd_samples_right]));
                
                candSignal = ECGpeakfrqBndPssSignal_hilbert;
                threshold_ECG = 2*std(ECGpeakfrqBndPssSignal_hilbert);
                
                
                minPeakDistanceSamples = FrqOfSmpl * 0.2;%200ms the refractory time of a heart beat
                [tempCandSignalPeaks, tempCandSignalPeaksSamples] = findpeaks(candSignal,'MINPEAKHEIGHT',threshold_ECG,'MINPEAKDISTANCE',minPeakDistanceSamples);
                candSignalPeaks = candSignal(tempCandSignalPeaksSamples);
                candSignalPeaksSamples = tempCandSignalPeaksSamples;%    candSignalPeaksSamples = currentRawDataSampleOffset + tempCandSignalPeaksSamples;
                
                if numel(candSignalPeaksSamples) > 2
                    samples_diff = diff(candSignalPeaksSamples);
                    %tempinstHRmin_pre = cat(2,0,60./(samples_diff(1:end-1)/FrqOfSmpl));
                    %tempinstHRmin_post = cat(2,0,60./(samples_diff(2:end)/FrqOfSmpl));
                    %tempinstHRmin_change = [0 , tempinstHRmin_post-tempinstHRmin_pre];
                    tempinstHRmin = cat(2,0,60./(samples_diff/FrqOfSmpl));
                    
                    tempinstHRmin = interp1(candSignalPeaksSamples,tempinstHRmin,1:length(candSignal),'linear');
                    
                    %tempinstHRmin_change = interp1(candSignalPeaksSamples,tempinstHRmin_change,1:length(candSignal),'linear');
                    
                    cfg.ECG_instHRsmin = tempinstHRmin(fix(cfg.nEpochsBuffer*lengthEpochSamples+1):fix((cfg.nEpochsBuffer+1)*lengthEpochSamples));
                    %cfg.ECG_instHRsmin_change = tempinstHRmin_change(fix(cfg.nEpochsBuffer*lengthEpochSamples+1):fix((cfg.nEpochsBuffer+1)*lengthEpochSamples));
                    cfg.ECG_instHRsmin(isnan(cfg.ECG_instHRsmin)) = 0;
                    %cfg.ECG_instHRsmin_change(isnan(cfg.ECG_instHRsmin_change)) = 0;
                    cfg.ECG_instHRsmin_SignalPeaksSamples = candSignalPeaksSamples((fix(cfg.nEpochsBuffer*lengthEpochSamples+1) <= candSignalPeaksSamples) & (candSignalPeaksSamples <= fix((cfg.nEpochsBuffer+1)*lengthEpochSamples))) - fix(cfg.nEpochsBuffer*lengthEpochSamples);
                else
                    cfg.ECG_instHRsmin = zeros(1,lengthEpochSamples);
                    %cfg.ECG_instHRsmin_change = zeros(1,lengthEpochSamples);
                    cfg.ECG_instHRsmin_SignalPeaksSamples = [];
                    
                end
            end
            
        end
    end
    
    setappdata(h, 'cfg', cfg);
    updateLabels(h);
    
end


if strcmp(cfg.doSleepScoring,'yes')
    curr_epoch = opt.trlop;
    
    temp_epochLengthSamples = opt.trlvis(1, 2) - opt.trlvis(1, 1) + 1;
    nEpochs = floor(size(opt.orgdata.trial{1},2)/temp_epochLengthSamples);
    if curr_epoch > nEpochs
        curr_epoch = curr_epoch - 1;
    end
    curr_hypn = cfg.hypn(curr_epoch,:);
    h1 = curr_hypn(:,1);
    h2 = curr_hypn(:,2);
    [stagestring h1_str h2_str] = getStageStringByHypnValue(h1,h2);
    opt.curr_stage = stagestring;
    
    opt.prev_stages = getPrevStageString_stage(cfg.hypn,curr_epoch,6);
    opt.next_stages = getNextStageString_stage(cfg.hypn,curr_epoch,6);
    
    if isfield(cfg, 'hypn_mult')
        if ~isempty(cfg.hypn_mult)
            for iHypMult = 1:numel(cfg.hypn_mult)
                opt.curr_stage_mult{iHypMult} = getStageStringByHypnValue(cfg.hypn_mult{iHypMult}(curr_epoch,1),cfg.hypn_mult{iHypMult}(curr_epoch,2));
                %                 opt.prev_stages_mult{iHypMult}
                %                 opt.next_stages_mult{iHypMult}
                %                 opt.curr_stage_mult = '?';
            end
        end
    end
end


if strcmp(cfg.doSleepScoring,'yes')
    
    temp_channel_number_in_curr_display = find(chanindx == cfg.score_channel_eeg_number);
    
    temp_ax = [];
    
    delete(findobj(h,'tag', 'scorechan_eeg'));
    
    for iChanDisplayed = temp_channel_number_in_curr_display;
        
        if cfg.has_more_than_3_channels && istrue(cfg.highlight_scoring_channels)
            h_scorechan_eeg = ft_plot_box([tim(1) tim(end) -1 1],'facealpha',0.5, 'facecolor', cfg.score_channel_eeg_color , 'edgecolor', 'none', 'tag', 'scorechan_eeg',  ...
                'hpos', opt.laytime.pos(iChanDisplayed,1), 'vpos', opt.laytime.pos(iChanDisplayed,2), 'width', opt.width, 'height', opt.laytime.height(iChanDisplayed), 'hlim', opt.hlim, 'vlim', [-1 1], 'axis', temp_ax);
            temp_ax = h_scorechan_eeg;
            
            h_scorechan_eeg_middle_line = ft_plot_line([tim(1) tim(end)],[0 0], 'color', cfg.score_channel_eeg_color , 'linewidth', 1, 'tag', 'scorechan_eeg',  ...
                'hpos', opt.laytime.pos(iChanDisplayed,1), 'vpos', opt.laytime.pos(iChanDisplayed,2), 'width', opt.width, 'height', opt.laytime.height(iChanDisplayed), 'hlim', opt.hlim, 'vlim', [-1 1], 'axis', temp_ax);
            temp_ax = h_scorechan_eeg_middle_line;
            %        ft_plot_matrix([tim(1) tim(end)],[-1 1],datamatrix, 'clim',[-1,1], 'tag', 'TFR',  ...
            %             'hpos', opt.laytime.pos(iChanDisplayed,1), 'vpos', opt.laytime.pos(iChanDisplayed,2), 'width', opt.width, 'height', opt.laytime.height(iChanDisplayed), 'hlim', opt.hlim, 'vlim', [-1 1]);
            
            
            
            %         if ~isfield(opt,'curr_stage')
            %             opt.curr_stage = '?';
            %         end
            %
            %
            %         temp_curr_channels_displayed = numel(opt.laytime.label);
            %
            %         h_curr_stage = ft_plot_text(tim(floor(end/2)), 0, opt.curr_stage, 'tag', 'curr_stage', 'Color', [0.9 0.9 0.9], 'FontSize', 64, 'FontUnits',  'normalized', ...
            %             'hpos', opt.laytime.pos(iChanDisplayed,1), 'vpos', opt.laytime.pos(iChanDisplayed,2), 'width', opt.width, 'height', opt.laytime.height(iChanDisplayed), 'hlim', opt.hlim, 'vlim', [-1 1],'interpreter','none', 'axis', temp_ax);
            %         temp_ax = h_curr_stage;
            
        end
    end
    
    %half channel hight
    temp_channel_number_in_curr_display = ceil(numel(chanindx)/2);
    
    %temp_ax = [];
    
    delete(findobj(h,'tag', 'curr_stage'));
    
    if istrue(cfg.plot_stage_signatures)
        for iChanDisplayed = temp_channel_number_in_curr_display;
        
        if numel(cfg.hypn_mult) > 0
            if isfield(opt,'curr_stage_mult')
                temp_col = jet(4);
                for iHypMult = 1:numel(opt.curr_stage_mult)
                    
                    h_curr_stage = ft_plot_text(tim(floor(end/2)), 0, opt.curr_stage_mult{iHypMult}, 'tag', 'curr_stage', 'Color', temp_col(iHypMult,:), 'FontSize', 0.15-iHypMult*0.02, 'FontUnits',  'normalized', 'FontName', 'FixedWidth', ...
                        'hpos', opt.laytime.pos(iChanDisplayed,1), 'vpos', opt.laytime.pos(iChanDisplayed,2), 'width', opt.width, 'height', opt.laytime.height(iChanDisplayed), 'hlim', opt.hlim, 'vlim', [-1 1],'interpreter','none', 'axis', temp_ax);
                    temp_ax = h_curr_stage;
                    
                    
                    %                 temp_iChanDisplayed = iChanDisplayed;%iChanDisplayed-1;
                    %                 if temp_iChanDisplayed < 1
                    %                     temp_iChanDisplayed = 1;
                    %                 end
                    %                 h_prev_stage = ft_plot_text(tim(floor(end/2)), 0, [opt.prev_stages_mult{iHypMult} '                    ' opt.next_stages_mult{iHypMult}], 'tag', 'curr_stage', 'Color', [0.9 0.9 0.9]-iHypMult*0.2, 'FontSize', 24-iHypMult*2, 'FontUnits',  'normalized', ...
                    %                     'hpos', opt.laytime.pos(temp_iChanDisplayed,1), 'vpos', opt.laytime.pos(temp_iChanDisplayed,2), 'width', opt.width, 'height', opt.laytime.height(temp_iChanDisplayed), 'hlim', opt.hlim, 'vlim', [-1 1],'interpreter','none', 'axis', temp_ax);
                    %                 temp_ax = h_prev_stage;
                    %
                end
            end
        end
        
        if ~isfield(opt,'curr_stage')
            opt.curr_stage = '?';
        end
        
        if ~isfield(opt,'prev_stages')
            opt.prev_stages = '';
        end
        
        if ~isfield(opt,'next_stages')
            opt.next_stages = '';
        end
        
        
        %temp_curr_channels_displayed = numel(opt.laytime.label);
        
        h_curr_stage = ft_plot_text(tim(floor(end/2)), 0, opt.curr_stage, 'tag', 'curr_stage', 'Color', cfg.color_text_on_bg , 'FontSize', 0.15, 'FontUnits',  'normalized', 'FontName', 'FixedWidth', ...
            'hpos', opt.laytime.pos(iChanDisplayed,1), 'vpos', opt.laytime.pos(iChanDisplayed,2), 'width', opt.width, 'height', opt.laytime.height(iChanDisplayed), 'hlim', opt.hlim, 'vlim', [-1 1],'interpreter','none', 'axis', temp_ax);
        temp_ax = h_curr_stage;
        
        
        temp_iChanDisplayed = iChanDisplayed;%iChanDisplayed-1;
        if temp_iChanDisplayed < 1
            temp_iChanDisplayed = 1;
        end
        temp_prev_stages = opt.prev_stages;
        temp_next_stages = opt.next_stages;
        if length(temp_prev_stages) < 6
            temp_prev_stages = [temp_prev_stages [repmat(' ',1,6-length(temp_prev_stages))]];
        end
        if length(temp_next_stages) < 6
            temp_next_stages = [[repmat(' ',1,6-length(temp_next_stages))] temp_next_stages];
        end
        
        h_prev_stage = ft_plot_text(tim(floor(end/2)), 0, [temp_prev_stages '                     ' temp_next_stages], 'tag', 'curr_stage', 'Color', cfg.color_text_on_bg, 'FontSize', 0.05, 'FontUnits',  'normalized', 'FontName', 'FixedWidth', ...
            'hpos', opt.laytime.pos(temp_iChanDisplayed,1), 'vpos', opt.laytime.pos(temp_iChanDisplayed,2), 'width', opt.width, 'height', opt.laytime.height(temp_iChanDisplayed), 'hlim', opt.hlim, 'vlim', [-1 1],'interpreter','none', 'axis', temp_ax);
        temp_ax = h_prev_stage;
        
        
        
        
        
    end
    end
    
    
    delete(findobj(h, 'tag', 'scorechan_eeg_frontal'));
    
    temp_channel_number_in_curr_display = find(chanindx == cfg.score_channel_eeg_frontal_number);
    
    for iChanDisplayed = temp_channel_number_in_curr_display;
        if cfg.has_more_than_3_channels && istrue(cfg.highlight_scoring_channels) && (cfg.score_channel_eeg_frontal_number ~= cfg.score_channel_eeg_number)
            h_scorechan_eeg_frontal = ft_plot_box([tim(1) tim(end) -1 1],'facealpha',0.5, 'facecolor', cfg.score_channel_eeg_frontal_color , 'edgecolor', 'none', 'tag', 'scorechan_eeg_frontal',  ...
                'hpos', opt.laytime.pos(iChanDisplayed,1), 'vpos', opt.laytime.pos(iChanDisplayed,2), 'width', opt.width, 'height', opt.laytime.height(iChanDisplayed), 'hlim', opt.hlim, 'vlim', [-1 1], 'axis', temp_ax);
            temp_ax = h_scorechan_eeg_frontal;
        end
    end
    
    delete(findobj(h, 'tag', 'scorechan_eeg_occipital'));
    
    temp_channel_number_in_curr_display = find(chanindx == cfg.score_channel_eeg_occipital_number);
    
    for iChanDisplayed = temp_channel_number_in_curr_display;
        if cfg.has_more_than_3_channels && istrue(cfg.highlight_scoring_channels) && (cfg.score_channel_eeg_occipital_number ~= cfg.score_channel_eeg_number)
            h_scorechan_eeg_occipital = ft_plot_box([tim(1) tim(end) -1 1],'facealpha',0.5, 'facecolor', cfg.score_channel_eeg_occipital_color , 'edgecolor', 'none', 'tag', 'scorechan_eeg_occipital',  ...
                'hpos', opt.laytime.pos(iChanDisplayed,1), 'vpos', opt.laytime.pos(iChanDisplayed,2), 'width', opt.width, 'height', opt.laytime.height(iChanDisplayed), 'hlim', opt.hlim, 'vlim', [-1 1], 'axis', temp_ax);
            temp_ax = h_scorechan_eeg_occipital;
        end
    end
    
    delete(findobj(h, 'tag', 'scorechan_eog'));
    
    temp_channel_number_in_curr_display = find(chanindx == cfg.score_channel_eog_number);
    
    for iChanDisplayed = temp_channel_number_in_curr_display;
        if cfg.has_more_than_3_channels && istrue(cfg.highlight_scoring_channels)
            h_scorechan_eog = ft_plot_box([tim(1) tim(end) -1 1],'facealpha',0.5, 'facecolor', cfg.score_channel_eog_color , 'edgecolor', 'none', 'tag', 'scorechan_eog',  ...
                'hpos', opt.laytime.pos(iChanDisplayed,1), 'vpos', opt.laytime.pos(iChanDisplayed,2), 'width', opt.width, 'height', opt.laytime.height(iChanDisplayed), 'hlim', opt.hlim, 'vlim', [-1 1], 'axis', temp_ax);
            temp_ax = h_scorechan_eog;
        end
    end
    
    
    delete(findobj(h, 'tag', 'scorechan_emg'));
    
    temp_channel_number_in_curr_display = find(chanindx == cfg.score_channel_emg_number);
    
    for iChanDisplayed = temp_channel_number_in_curr_display;
        if cfg.has_more_than_3_channels && istrue(cfg.highlight_scoring_channels)
            h_scorechan_emg = ft_plot_box([tim(1) tim(end) -1 1],'facealpha',0.5, 'facecolor', cfg.score_channel_emg_color , 'edgecolor', 'none', 'tag', 'scorechan_emg',  ...
                'hpos', opt.laytime.pos(iChanDisplayed,1), 'vpos', opt.laytime.pos(iChanDisplayed,2), 'width', opt.width, 'height', opt.laytime.height(iChanDisplayed), 'hlim', opt.hlim, 'vlim', [-1 1], 'axis', temp_ax);
            temp_ax = h_scorechan_emg;
        end
    end
    
    
    delete(findobj(h, 'tag', 'scorechan_ecg'));
    
    if cfg.has_ECG
        
        temp_channel_number_in_curr_display = find(chanindx == cfg.score_channel_ecg_number);
        
        for iChanDisplayed = temp_channel_number_in_curr_display;
            if cfg.has_more_than_3_channels && istrue(cfg.highlight_scoring_channels)
                
                h_scorechan_ecg = ft_plot_box([tim(1) tim(end) -1 1],'facealpha',0.3, 'facecolor', cfg.score_channel_ecg_color , 'edgecolor', 'none', 'tag', 'scorechan_ecg',  ...
                    'hpos', opt.laytime.pos(iChanDisplayed,1), 'vpos', opt.laytime.pos(iChanDisplayed,2), 'width', opt.width, 'height', opt.laytime.height(iChanDisplayed), 'hlim', opt.hlim, 'vlim', [-1 1], 'axis', temp_ax);
                temp_ax = h_scorechan_ecg;
            end
        end
        
    end
end

delete(findobj(h,'tag', 'grid'));

if strcmp(cfg.drawgrid,'yes')
    %hold all
    h_curr_gridline = [];
    for iSec = 1:numel(cfg.drawgrid_seconds)
        for iSecLine = tim(1):cfg.drawgrid_seconds(iSec):tim(end)
            h_curr_gridline = ft_plot_line([iSecLine iSecLine], [-1 1], 'tag', 'grid', 'Color', cfg.drawgrid_colors{iSec}, 'linestyle',cfg.drawgrid_LineStyle{iSec}, ...
                'hpos', opt.hpos, 'vpos', opt.vpos, 'width', opt.width, 'height', opt.height, 'hlim', opt.hlim, 'vlim', [-1 1],'axis',h_curr_gridline);
        end
    end
    % hold off
end



% FW end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('plotting artifacts...\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delete(findobj(h,'tag', 'artifact'));

for j = ordervec
    tmp = diff([0 art(j,:) 0]);
    artbeg = find(tmp==+1);
    artend = find(tmp==-1) - 1;
    
    for k=1:numel(artbeg)
        
        h_artifact = ft_plot_box([tim(artbeg(k)) tim(artend(k)) -1 1], 'facecolor', opt.artcolors(j,:), 'edgecolor', 'none', 'tag', 'artifact',  ...
            'hpos', opt.hpos, 'vpos', opt.vpos, 'width', opt.width, 'height', opt.height, 'hlim', opt.hlim, 'vlim', [-1 1]);
        
        
    end
end % for each of the artifact channels

%FW begin

delete(findobj(h,'tag', 'event_begin_end2'));
delete(findobj(h,'tag', 'event_begin_end1'));

if strcmp(cfg.displayEvents,'yes')
    if isfield(cfg,'begin_end_events')
        
        for channelIndex = 1:length(chanindx)
            tmp = diff([0 cfg.times_ind_per_channel{chanindx(channelIndex)}(begsample:endsample) 0]);
            evbeg = find(tmp==+1);
            evend = find(tmp==-1) - 1;
            
            for k=1:numel(evbeg)
                h_event_begin_end1 = ft_plot_box([tim(evbeg(k)) tim(evend(k)) -0.8 0.8],'facealpha',0.38, 'facecolor', cfg.event_begin_end_color, 'edgecolor', 'none', 'tag', 'event_begin_end1',  ...
                    'hpos', opt.laytime.pos(channelIndex,1), 'vpos', opt.laytime.pos(channelIndex,2), 'width', opt.width, 'height', opt.laytime.height(channelIndex), 'hlim', opt.hlim, 'vlim', [-1 1]);
            end
            
            evbegend_ind = cfg.begin_end_events{chanindx(channelIndex)}(((cfg.begin_end_events{chanindx(channelIndex)}(:,1) >= begsample) & (cfg.begin_end_events{chanindx(channelIndex)}(:,1) <= endsample)),:);
            evbeg_ind = evbegend_ind(:,1);
            evend_ind = floor(evbeg_ind+ 0.1*(evbegend_ind(:,2) - evbeg_ind));
            
            evend_ind((evend_ind-evbeg_ind) < 1) = evend_ind((evend_ind-evbeg_ind) < 1)+1;
            evend_ind(evend_ind > endsample) = endsample;
            
            evbeg_ind = evbeg_ind-(begsample-1);
            evend_ind = evend_ind-(begsample-1);
            
            for k=1:numel(evbeg_ind)
                h_event_begin_end1_ind = ft_plot_box([tim(evbeg_ind(k)) tim(evend_ind(k)) 0.8 1],'facealpha',0.38, 'facecolor', cfg.event_begin_end_color, 'edgecolor', 'none', 'tag', 'event_begin_end1',  ...
                    'hpos', opt.laytime.pos(channelIndex,1), 'vpos', opt.laytime.pos(channelIndex,2), 'width', opt.width, 'height', opt.laytime.height(channelIndex), 'hlim', opt.hlim, 'vlim', [-1 1]);
                
            end
            
        end
    end
    
    if isfield(cfg,'begin_end_events2')
        
        for channelIndex = 1:length(chanindx)
            tmp = diff([0 cfg.times_ind_per_channel2{chanindx(channelIndex)}(begsample:endsample) 0]);
            evbeg = find(tmp==+1);
            evend = find(tmp==-1) - 1;
            
            for k=1:numel(evbeg)
                h_event_begin_end2 = ft_plot_box([tim(evbeg(k)) tim(evend(k)) -0.8 0.8],'facealpha',0.38, 'facecolor', cfg.event_begin_end_color2, 'edgecolor', 'none', 'tag', 'event_begin_end2',  ...
                    'hpos', opt.laytime.pos(channelIndex,1), 'vpos', opt.laytime.pos(channelIndex,2), 'width', opt.width, 'height', opt.laytime.height(channelIndex), 'hlim', opt.hlim, 'vlim', [-1 1]);
                
            end
            
            evbegend_ind = cfg.begin_end_events2{chanindx(channelIndex)}(((cfg.begin_end_events2{chanindx(channelIndex)}(:,1) >= begsample) & (cfg.begin_end_events2{chanindx(channelIndex)}(:,1) <= endsample)),:);
            evbeg_ind = evbegend_ind(:,1);
            evend_ind = floor(evbeg_ind+ 0.1*(evbegend_ind(:,2) - evbeg_ind));
            
            evend_ind((evend_ind-evbeg_ind) < 1) = evend_ind((evend_ind-evbeg_ind) < 1)+1;
            evend_ind(evend_ind > endsample) = endsample;
            
            evbeg_ind = evbeg_ind-(begsample-1);
            evend_ind = evend_ind-(begsample-1);
            
            for k=1:numel(evbeg_ind)
                h_event_begin_end2_ind = ft_plot_box([tim(evbeg_ind(k)) tim(evend_ind(k)) 0.8 1],'facealpha',0.38, 'facecolor', cfg.event_begin_end_color2, 'edgecolor', 'none', 'tag', 'event_begin_end2',  ...
                    'hpos', opt.laytime.pos(channelIndex,1), 'vpos', opt.laytime.pos(channelIndex,2), 'width', opt.width, 'height', opt.laytime.height(channelIndex), 'hlim', opt.hlim, 'vlim', [-1 1]);
                
            end
            
        end
    end
end
% begin_end_events_per_channel = cfg.begin_end_events;
% for channelIndex = 1:length(chanindx)
%     curr_begins_ends = begin_end_events_per_channel{chanindx(channelIndex)};
%     if ~isempty(curr_begins_ends)
%         curr_begins_ends = curr_begins_ends((curr_begins_ends(:,1) <= endsample) | (curr_begins_ends(:,2) >= begsample),:);
%         curr_begins_ends_seconds = curr_begins_ends ./opt.fsample;
%         full_time = opt.orgdata.time{1,1};
%         begsec = begsample./opt.fsample;
%         endsec = endsample./opt.fsample;
%
%         full_time_ind_segment = ((full_time >= begsec) & (full_time <= endsec));
%         %full_time_ind_events = zeros(1,length(full_time_ind_segment));
%         for iEv = 1:size(curr_begins_ends_seconds,1)
%             temp_beg = curr_begins_ends_seconds(iEv,1);
%             temp_end = curr_begins_ends_seconds(iEv,2);
%             full_time_ind_events = ((full_time >= temp_beg) & (full_time <= temp_end));
%             ind_highlight = full_time_ind_events(begsample:endsample);
%             h_event_begin_end = ft_plot_box([temp_beg temp_end -1 1], 'facecolor', cfg.event_begin_end_color, 'edgecolor', 'none', 'tag', 'artifact',  ...
%             'hpos', opt.laytime.pos(channelIndex,1), 'vpos', opt.laytime.pos(channelIndex,2), 'width', opt.width, 'height', opt.laytime.height(channelIndex), 'hlim', opt.hlim, 'vlim', [-1 1]);
%         end
%     end
% end % for each of the limEvents each channels


% %FW end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fprintf('plotting events...\n');
if strcmp(cfg.ploteventlabels , 'colorvalue') && ~isempty(opt.event)
    eventlabellegend = [];
    for iType = 1:length(opt.eventtypes)
        eventlabellegend = [eventlabellegend sprintf('%s = %s\n',opt.eventtypes{iType},opt.eventtypecolorlabels{iType})];
    end
    fprintf(eventlabellegend);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delete(findobj(h,'tag', 'event'));
% save stuff to able to shift event labels downwards when they occur at the same time-point
eventcol = cell(1,numel(event));
eventstr = cell(1,numel(event));
eventtim = NaN(1,numel(event));
% gather event info and plot lines
for ievent = 1:numel(event)
    try
        if strcmp(cfg.ploteventlabels , 'type=value')
            if isempty(event(ievent).value)
                eventstr{ievent} = '';
            else
                eventstr{ievent} = sprintf('%s = %s', event(ievent).type, num2str(event(ievent).value)); % value can be both number and string
            end
            eventcol{ievent} = 'k';
        elseif strcmp(cfg.ploteventlabels , 'colorvalue')
            eventcol{ievent} = opt.eventtypescolors(match_str(opt.eventtypes, event(ievent).type),:);
            eventstr{ievent} = sprintf('%s', num2str(event(ievent).value)); % value can be both number and string
        end
    catch
        eventstr{ievent} = 'unknown';
        eventcol{ievent} = 'k';
    end
    eventtim(ievent) = (event(ievent).sample-begsample)/opt.fsample + opt.hlim(1);
    ft_plot_line([eventtim(ievent) eventtim(ievent)], [-1 1], 'tag', 'event', 'color', eventcol{ievent}, ...
        'hpos', opt.hpos, 'vpos', opt.vpos, 'width', opt.width, 'height', opt.height, 'hlim', opt.hlim, 'vlim', [-1 1]);
end
% count the consecutive occurrence of each time point
concount = NaN(1,numel(event));
for ievent = 1:numel(event)
    concount(ievent) = sum(eventtim(ievent)==eventtim(1:ievent-1));
end
% plot labels
for ievent = 1:numel(event)
    ft_plot_text(eventtim(ievent), 0.9-concount(ievent)*.06, eventstr{ievent}, 'tag', 'event', 'Color', eventcol{ievent}, ...
        'hpos', opt.hpos, 'vpos', opt.vpos, 'width', opt.width, 'height', opt.height, 'hlim', opt.hlim, 'vlim', [-1 1],'interpreter','none');
end


delete(findobj(h,'tag', 'toggle_marker'));

if strcmp(cfg.doSleepScoring,'yes')
    if cfg.toggle_epoch_marker == opt.trlop
        h_top = ft_plot_line([tim(1) tim(end)],[1 1], 'color', [0.25 1 0.125] , 'linewidth', 8, 'tag', 'toggle_marker',  ...
            'hpos', opt.laytime.pos(1,1), 'vpos', opt.laytime.pos(1,2), 'width', opt.width, 'height', opt.laytime.height(1), 'hlim', opt.hlim, 'vlim', [-1 1]);
    end
end

delete(findobj(h,'tag', 'title_time'));

startim = tim(1);
if nsamplepad>0
    endtim = tim(end-nsamplepad);
else
    endtim = tim(end);
end

if ~strcmp(opt.trialviewtype,'trialsegment')
    str = sprintf('%s %d/%d, time from %0.2f to %0.2f min', opt.trialviewtype, opt.trlop, size(opt.trlvis,1), startim/60, endtim/60);
else
    str = sprintf('trial %d/%d: epoch: %d/%d , time from %g to %g min', opt.trllock, size(opt.trlorg,1), opt.trlop, size(opt.trlvis,1), startim/60, endtim/60);
end


%title(str,'interpreter','none');
temp_ylim = get(gca, 'ylim');
temp_xlim = get(gca, 'xlim');
delete(findobj(h,'tag', 'title_time'));
text(temp_xlim(end)-0.01*range(temp_xlim),temp_ylim(end)-0.005*range(temp_ylim),str,'HorizontalAlignment','right','VerticalAlignment','top','color',cfg.color_text_on_bg,'tag','title_time','FontSize', 12, 'FontUnits',  'points', 'interpreter','none')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fprintf('plotting data...\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delete(findobj(h,'tag', 'timecourse'));
delete(findobj(h,'tag', 'identify'));
delete(findobj(h,'tag', 'spindle_timecourse'));
delete(findobj(h,'tag', 'alpha_timecourse'));
delete(findobj(h,'tag', 'so_timecourse'));
delete(findobj(h,'tag', 'ecg_HR_timecourse'));
delete(findobj(h,'tag', 'ecg_HR_peaks_markers'));




% not removing channel labels, they cause the bulk of redrawing time for the slow text function (note, interpreter = none hardly helps)
% warning, when deleting the labels using the line below, one can easily tripple the excution time of redrawing in viewmode = vertical (see bug 2065)
%delete(findobj(h,'tag', 'chanlabel'));




if strcmp(cfg.viewmode, 'butterfly')
    set(gca,'ColorOrder',opt.chancolors(chanindx,:)) % plot vector does not clear axis, therefore this is possible
    ft_plot_vector(tim, dat, 'box', false, 'tag', 'timecourse', ...
        'hpos', opt.laytime.pos(1,1), 'vpos', opt.laytime.pos(1,2), 'width', opt.laytime.width(1), 'height', opt.laytime.height(1), 'hlim', opt.hlim, 'vlim', opt.vlim);
    
    
    % two ticks per channel
    yTick = sort([opt.laytime.pos(:,2)+(opt.laytime.height/2); ...
        opt.laytime.pos(:,2)+(opt.laytime.height/4); ...
        opt.laytime.pos(:,2);                        ...
        opt.laytime.pos(:,2)-(opt.laytime.height/4); ...
        opt.laytime.pos(:,2)-(opt.laytime.height/2)]);
    
    yTickLabel = {num2str(yTick.*range(opt.vlim) + opt.vlim(1))};
    
    set(gca, 'yTick', yTick);
    set(gca, 'yTickLabel', yTickLabel)
    
elseif any(strcmp(cfg.viewmode, {'vertical' 'component'}))
    
    % determine channel indices into data outside of loop
    laysels = match_str(opt.laytime.label, opt.hdr.label);
    
    for i = 1:length(chanindx)
        if strcmp(cfg.viewmode, 'component')
            color = 'k';
        else
            color = opt.chancolors(chanindx(i),:);
        end
        datsel = i;
        laysel = laysels(i);
        if ~isempty(datsel) && ~isempty(laysel)
            
            % only plot labels when current chanlabel objects are less then the total number of channels (see bug 2065)
            % this is a cheap quick fix. If it causes error in plotting components, do this conditional on viewmode
            if numel(findobj(h,'tag', 'chanlabel'))<numel(chanindx)
                if opt.plotLabelFlag == 1 || (opt.plotLabelFlag == 2 && mod(i,10)==0)

                    h_chanlabel = ft_plot_text(tim(1)+range(tim)*0.0125, 0.5, opt.hdr.label(chanindx(i)), 'tag', 'chanlabel', 'Color', cfg.color_text_on_bg, 'FontSize', 0.025, 'FontUnits',  'normalized','HorizontalAlignment', 'left', ...
                        'hpos', opt.laytime.pos(laysel,1), 'vpos', opt.laytime.pos(laysel,2), 'width', opt.width, 'height', opt.laytime.height(laysel), 'hlim', opt.hlim, 'vlim', [-1 1],'interpreter','none');
                    
%                     h_chanlabel = ft_plot_text(tim(1)+range(tim)*0.0125, 0.5, opt.hdr.label(chanindx(i)), 'tag', 'chanlabel', 'Color', cfg.color_text_on_bg, 'FontSize', 0.9/2/numel(chanindx), 'FontUnits',  'normalized','HorizontalAlignment', 'left', ...
%                         'hpos', opt.laytime.pos(laysel,1), 'vpos', opt.laytime.pos(laysel,2), 'width', opt.width, 'height', opt.laytime.height(laysel), 'hlim', opt.hlim, 'vlim', [-1 1],'interpreter','none');
                    
                    %ft_plot_text(labelx(laysel), labely(laysel), opt.hdr.label(chanindx(i)), 'tag', 'chanlabel', 'HorizontalAlignment', 'right','interpreter','none','FontUnits','normalized','FontSize',0.9/2/numel(chanindx));
                end
            end
            
            if strcmp(cfg.doSleepScoring,'yes')
                
                if laysel == find(chanindx == cfg.score_channel_eeg_number);
                    
                    if strcmp(cfg.underlaySpindleSignal,'yes')
                        if isfield(cfg,'spindsignal_envelope_display')
                            ft_plot_vector(tim(1:numel(cfg.spindsignal_envelope_display)), cfg.spindsignal_envelope_display, 'box', false, 'color', cfg.underlaySpindleSignal_color, 'tag', 'spindle_timecourse', ...
                                'hpos', opt.laytime.pos(laysel,1), 'vpos', opt.laytime.pos(laysel,2), 'width', opt.laytime.width(laysel), 'height', opt.laytime.height(laysel), 'hlim', opt.hlim, 'vlim', opt.vlim);
                            ft_plot_vector(tim(1:numel(cfg.spindsignal_envelope_display)), -cfg.spindsignal_envelope_display, 'box', false, 'color', cfg.underlaySpindleSignal_color, 'tag', 'spindle_timecourse', ...
                                'hpos', opt.laytime.pos(laysel,1), 'vpos', opt.laytime.pos(laysel,2), 'width', opt.laytime.width(laysel), 'height', opt.laytime.height(laysel), 'hlim', opt.hlim, 'vlim', opt.vlim);
                            
                        end
                        
                        if isfield(cfg,'spindsignal_display')
                            ft_plot_vector(tim(1:numel(cfg.spindsignal_display)), cfg.spindsignal_display , 'box', false, 'color', [0.75 0.75 0.75], 'tag', 'spindle_timecourse', ...
                                'hpos', opt.laytime.pos(laysel,1), 'vpos', opt.laytime.pos(laysel,2), 'width', opt.laytime.width(laysel), 'height', opt.laytime.height(laysel), 'hlim', opt.hlim, 'vlim', opt.vlim);
                        end
                    end
                end
                
                
                if laysel == find(chanindx == cfg.score_channel_eeg_occipital_number);
                    
                    if strcmp(cfg.underlayAlphaSignal,'yes')
                        
                        if isfield(cfg,'alphasignal_display')
                            ft_plot_vector(tim(1:numel(cfg.alphasignal_display)), cfg.alphasignal_display , 'box', false, 'color', [0.75 0.75 0.75], 'tag', 'alpha_timecourse', ...
                                'hpos', opt.laytime.pos(laysel,1), 'vpos', opt.laytime.pos(laysel,2), 'width', opt.laytime.width(laysel), 'height', opt.laytime.height(laysel), 'hlim', opt.hlim, 'vlim', opt.vlim);
                        end
                    end
                end
                
                if laysel == find(chanindx == cfg.score_channel_eeg_frontal_number);
                  if strcmp(cfg.underlaySOSignal,'yes')
                        if isfield(cfg,'so_signal_display')
                            ft_plot_vector(tim(1:numel(cfg.so_signal_display)), cfg.so_signal_display , 'box', false, 'color', cfg.underlaySOSignal_color, 'tag', 'so_timecourse', ...
                                'hpos', opt.laytime.pos(laysel,1), 'vpos', opt.laytime.pos(laysel,2), 'width', opt.laytime.width(laysel), 'height', opt.laytime.height(laysel), 'hlim', opt.hlim, 'vlim', opt.vlim);
                        end
                        
                    end
                end
                
            end
            
            %the time course of channels
            ft_plot_vector(tim, dat(datsel, :), 'box', false, 'color', color, 'tag', 'timecourse', ...
                'hpos', opt.laytime.pos(laysel,1), 'vpos', opt.laytime.pos(laysel,2), 'width', opt.laytime.width(laysel), 'height', opt.laytime.height(laysel), 'hlim', opt.hlim, 'vlim', opt.vlim);
            
            
            if strcmp(cfg.doSleepScoring,'yes')
                if cfg.has_ECG
                    
                    if laysel == find(chanindx == cfg.score_channel_ecg_number);
                        
                        if strcmp(cfg.markECG,'yes')
                            
                            %                             %rgb_ind = 1+fix(fw_normalize(cfg.ECG_instHRsmin, min(cfg.ECG_instHRsmin), max(cfg.ECG_instHRsmin), 0, 127));
                            %                             %rgb_ind(isnan(rgb_ind)) = 1;
                            %                             %rgb = flip(autumn(128),1);
                            %                             %ft_plot_vector(tim, dat(datsel, :), 'box', false, 'color', rgb(rgb_ind,:), 'tag', 'ecg_HR_timecourse', ...
                            %                             ft_plot_vector(tim, dat(datsel, :), 'box', false, 'color', 'k', 'tag', 'ecg_HR_timecourse', ...
                            %                             'linewidth', 1.5, ...
                            %                                 'hpos', opt.laytime.pos(laysel,1), 'vpos', opt.laytime.pos(laysel,2), 'width', opt.laytime.width(laysel), 'height', opt.laytime.height(laysel), 'hlim', opt.hlim, 'vlim', opt.vlim);
                            %
                            if ~isempty(cfg.ECG_instHRsmin_SignalPeaksSamples)
                                rgb_ind = 1+fix(fw_normalize(cfg.ECG_instHRsmin(cfg.ECG_instHRsmin_SignalPeaksSamples), min(cfg.ECG_instHRsmin), max(cfg.ECG_instHRsmin), 0, 127));
                                rgb_ind(isnan(rgb_ind)) = 1;
                                rgb = flip(autumn(128),1);
                                for iECGpeak = 1:numel(cfg.ECG_instHRsmin_SignalPeaksSamples)
                                    %ft_plot_text(tim(cfg.ECG_instHRsmin_SignalPeaksSamples(iECGpeak)), dat(datsel, cfg.ECG_instHRsmin_SignalPeaksSamples(iECGpeak)), num2str(round(cfg.ECG_instHRsmin(cfg.ECG_instHRsmin_SignalPeaksSamples(iECGpeak)))),'FontSize', 8, 'tag', 'ecg_HR_peaks_markers','interpreter','none','color', [0.2 0.2 0.2], ...
                                    ft_plot_text(tim(cfg.ECG_instHRsmin_SignalPeaksSamples(iECGpeak)), dat(datsel, cfg.ECG_instHRsmin_SignalPeaksSamples(iECGpeak)), num2str(round(cfg.ECG_instHRsmin(cfg.ECG_instHRsmin_SignalPeaksSamples(iECGpeak)))),'FontSize', 8, 'tag', 'ecg_HR_peaks_markers','interpreter','none','color', rgb(rgb_ind(iECGpeak),:), ...
                                        'hpos', opt.laytime.pos(laysel,1), 'vpos', opt.laytime.pos(laysel,2), 'width', opt.laytime.width(laysel), 'height', opt.laytime.height(laysel), 'hlim', opt.hlim, 'vlim', opt.vlim);
                                end
                                %                             ft_plot_vector(tim(cfg.ECG_instHRsmin_SignalPeaksSamples), dat(datsel, cfg.ECG_instHRsmin_SignalPeaksSamples), 'box', false, 'color', 'none', 'marker', 'o', 'markerfacecolor', [1 90/255 0], 'tag', 'ecg_HR_peaks_markers', ...
                                %                                 'hpos', opt.laytime.pos(laysel,1), 'vpos', opt.laytime.pos(laysel,2), 'width', opt.laytime.width(laysel), 'height', opt.laytime.height(laysel), 'hlim', opt.hlim, 'vlim', opt.vlim);
                            end
                        end
                    end
                    
                end
                
            end
            
        end
    end
    
    if length(chanindx)>19
        % no space for yticks
        yTick = [];
        yTickLabel = [];
    else
        %FW begin
        yTickLabel = [];
        for i = 1:length(chanindx)
            curr_scale = cfg.chanscale(chanindx);
            if length(chanindx) > 6
                % one tick per channel
                yTick = sort([opt.laytime.pos(:,2)+(opt.laytime.height(laysel)/4); ...
                    opt.laytime.pos(:,2)-(opt.laytime.height(laysel)/4)]);
                %FW begin
                temp_tick = [.25 .75] .* range(opt.vlim./curr_scale(i)) + opt.vlim(1)./abs(curr_scale(i));
                yTickLabel_temp = cellfun(@str2num,cellfun(@(x) num2str(x,'%1.2f'),{temp_tick(1) temp_tick(2)},'UniformOutput',false));
                yTickLabel = [yTickLabel {yTickLabel_temp}];
                %yTickLabel = [yTickLabel {[.25 .75] .* temp_factor}];
                %FW end
            else
                %FW begin
                % two ticks per channel
                temp_tick = [0 .25 .75 1] .* range(opt.vlim./curr_scale(i)) + opt.vlim(1)./abs(curr_scale(i));
                %annotation('line',[-0.01 -0.01],[0 1],'color',[0.5 0.5 0.5]) .* range(opt.vlim./curr_scale(i)) + opt.vlim(1)./abs(curr_scale(i)));
                yTick = sort([opt.laytime.pos(:,2)+(opt.laytime.height(laysel)/2); ...
                    opt.laytime.pos(:,2)+(opt.laytime.height(laysel)/4); ...
                    opt.laytime.pos(:,2)-(opt.laytime.height(laysel)/4); ...
                    opt.laytime.pos(:,2)-(opt.laytime.height(laysel)/2)]);
                yTickLabel_temp = cellfun(@str2num,cellfun(@(x) num2str(x,'%1.2f'),{temp_tick(1) temp_tick(2) temp_tick(3) temp_tick(4)},'UniformOutput',false));
                if i==1
                    yTickLabel = [{num2str(yTickLabel_temp(1)) num2str(yTickLabel_temp(2)) num2str(yTickLabel_temp(3)) ' '} yTickLabel];
                elseif i == length(chanindx)
                    yTickLabel = [{' ' num2str(yTickLabel_temp(2)) num2str(yTickLabel_temp(3)) num2str(yTickLabel_temp(4))} yTickLabel];
                else
                    yTickLabel = [{num2str(yTickLabel_temp(1)) num2str(yTickLabel_temp(2)) num2str(yTickLabel_temp(3)) num2str(yTickLabel_temp(4))} yTickLabel];
                end
                %yTickLabel = [yTickLabel {[.0 .25 .75 1] .* temp_factor}];
                %FW end
            end
        end
    end
    %FW end
    %FW begin
    %yTickLabel = repmat(yTickLabel, 1, length(chanindx));
    set(gca, 'yTick', yTick);
    if length(chanindx) > 6
        set(gca, 'yTickLabel', flip(yTickLabel));
    else
        set(gca, 'yTickLabel', yTickLabel);
    end
    %FW end
    
else
    error('unknown viewmode "%s"', cfg.viewmode);
end % if strcmp viewmode


nticks = 11;
if strcmp(cfg.doSleepScoring,'yes')
    xTickLabel = cellstr(num2str( round(linspace(tim(1), tim(end), nticks)') , '%1.1f'))';
    if nsamplepad>0
        nlabindat = sum(linspace(tim(1), tim(end), nticks) < tim(end-nsamplepad));
        xTickLabel(nlabindat+1:end) = repmat({' '},[1 nticks-nlabindat]);
    end
    xTick = linspace(ax(1), ax(2), nticks);
    set(gca, 'xTick',xTick)
    %yTick = get(gca, 'yTick');
    temp_ylim = get(gca, 'ylim');
    temp_xlim = get(gca, 'xlim');
    set(gca,'xticklabel',[])
    delete(findobj(h,'tag', 'xticks'));
    for k = 1:length(xTick)
        if k == 1
            text(xTick(k)+0.001*range(temp_xlim),temp_ylim(1)+0.0025*range(temp_ylim),xTickLabel{k},'HorizontalAlignment','left','VerticalAlignment','bottom','tag','xticks','FontSize', 8, 'FontUnits',  'points')
        else
            text(xTick(k)-0.001*range(temp_xlim),temp_ylim(1)+0.0025*range(temp_ylim),xTickLabel{k},'HorizontalAlignment','right','VerticalAlignment','bottom','tag','xticks','FontSize', 8, 'FontUnits',  'points')
        end
    end
    %set(gca, 'xTickLabel', xTickLabel)
else
    xTickLabel = cellstr(num2str( linspace(tim(1), tim(end), nticks)' , '%1.2f'))';
    if nsamplepad>0
        nlabindat = sum(linspace(tim(1), tim(end), nticks) < tim(end-nsamplepad));
        xTickLabel(nlabindat+1:end) = repmat({' '},[1 nticks-nlabindat]);
    end
    set(gca, 'xTick', linspace(ax(1), ax(2), nticks))
    set(gca, 'xTickLabel', xTickLabel)
end


if strcmp(cfg.viewmode, 'component')
    
    % determine the position of each of the original channels for the topgraphy
    laychan = opt.layorg;
    
    % determine the position of each of the topographies
    laytopo.pos(:,1)  = opt.laytime.pos(:,1) - opt.laytime.width/2 - opt.laytime.height;
    laytopo.pos(:,2)  = opt.laytime.pos(:,2) + opt.laytime.height/2;
    laytopo.width     = opt.laytime.height;
    laytopo.height    = opt.laytime.height;
    laytopo.label     = opt.laytime.label;
    
    if ~isequal(opt.chanindx, chanindx)
        opt.chanindx = chanindx;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fprintf('plotting component topographies...\n');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        delete(findobj(h,'tag', 'topography'));
        
        [sel1, sel2] = match_str(opt.orgdata.topolabel, laychan.label);
        chanx = laychan.pos(sel2,1);
        chany = laychan.pos(sel2,2);
        
        if strcmp(cfg.compscale, 'global')
            for i=1:length(chanindx) % loop through all components to get max and min
                zmin(i) = min(opt.orgdata.topo(sel1,chanindx(i)));
                zmax(i) = max(opt.orgdata.topo(sel1,chanindx(i)));
            end
            
            if strcmp(cfg.zlim, 'maxmin')
                zmin = min(zmin);
                zmax = max(zmax);
            elseif strcmp(cfg.zlim, 'maxabs')
                zmax = max([abs(zmin) abs(zmax)]);
                zmin = -zmax;
            else
                error('configuration option for component scaling could not be recognized');
            end
        end
        
        for i=1:length(chanindx)
            % plot the topography of this component
            laysel = match_str(opt.laytime.label, opt.hdr.label(chanindx(i)));
            chanz = opt.orgdata.topo(sel1,chanindx(i));
            
            if strcmp(cfg.compscale, 'local')
                % compute scaling factors here
                if strcmp(cfg.zlim, 'maxmin')
                    zmin = min(chanz);
                    zmax = max(chanz);
                elseif strcmp(cfg.zlim, 'maxabs')
                    zmax = max(abs(chanz));
                    zmin = -zmax;
                end
            end
            
            % scaling
            chanz = (chanz - zmin) ./  (zmax- zmin);
            
            % laychan is the actual topo layout, in pixel units for .mat files
            % laytopo is a vertical layout determining where to plot each topo,
            %   with one entry per component
            
            
            ft_plot_topo(chanx, chany, chanz, 'mask', ...
                laychan.mask, 'interplim', 'mask', 'outline', ...
                laychan.outline, 'tag', 'topography', ...
                'hpos', laytopo.pos(laysel,1)-laytopo.width(laysel)/2,...
                'vpos', laytopo.pos(laysel,2)-laytopo.height(laysel)/2,...
                'width', laytopo.width(laysel), 'height', laytopo.height(laysel), 'gridscale', 45);
            
            %axis equal
            %drawnow
        end
        
        caxis([0 1]);
        
    end % if redraw_topo
    
    set(gca, 'yTick', [])
    
    ax(1) = min(laytopo.pos(:,1) - laytopo.width);
    ax(2) = max(opt.laytime.pos(:,1) + opt.laytime.width/2);
    ax(3) = min(opt.laytime.pos(:,2) - opt.laytime.height/2);
    ax(4) = max(opt.laytime.pos(:,2) + opt.laytime.height/2);
    axis(ax)
end % plotting topographies





if strcmp(cfg.doSleepScoring,'yes')
    xlabel([],'interpreter','none');
else
    xlabel('time','interpreter','none');
end

if isfield(opt,'marks')
    redraw_marks_cb(h)
end

%hold off;

% possibly adds some responsiveness if the 'thing' is clogged
%drawnow
drawnow expose
%drawnow update
%drawnow expose update




hax = get(h, 'CurrentAxes');
opt.xlim_init = get(hax, 'xlim');
opt.ylim_init = get(hax, 'ylim');

setappdata(h, 'opt', opt);
setappdata(h, 'cfg', cfg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function key = parseKeyboardEvent(eventdata)
% 
% key = eventdata.Key;
% 
% % handle possible numpad events (different for Windows and UNIX systems)
% % NOTE: shift+numpad number does not work on UNIX, since the shift
% % modifier is always sent for numpad events
% if isunix()
%     shiftInd = match_str(eventdata.Modifier, 'shift');
%     if ~isnan(str2double(eventdata.Character)) && ~isempty(shiftInd)
%         % now we now it was a numpad keystroke (numeric character sent AND
%         % shift modifier present)
%         key = eventdata.Character;
%         eventdata.Modifier(shiftInd) = []; % strip the shift modifier
%     end
% elseif ispc()
%     if strfind(eventdata.Key, 'numpad')
%         key = eventdata.Character;
%     end
% end
% 
% if ~isempty(eventdata.Modifier)
%     key = [eventdata.Modifier{1} '+' key];
% end
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_sleep_stage_cb(hObject, eventdata, varargin) %
h_main = ft_getopt(varargin, 'h_main');
h = getparent(h_main);
opt = getappdata(h, 'opt');
cfg = getappdata(h, 'cfg');

%[press_x press_y] = ft_select_point([0 0],'nearest',0);
point = get(cfg.hhypfigax,'CurrentPoint');     % button down detected
press_x = point(1,1);
xlime_hypn = get(cfg.hhypfigax, 'xlim');
if (press_x > xlime_hypn(1)) && (press_x  < xlime_hypn(2))
    press_x_sec = press_x*60;
    
    new_epoch = max(1,fix(press_x_sec/cfg.blocksize));
    if new_epoch ~= opt.trlop
        opt.trlop = max(1,fix(press_x_sec/cfg.blocksize));
        
        setappdata(h, 'opt',opt);
        setappdata(h, 'cfg',cfg);
        redraw_cb(h);
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_marks_cb(hObject, eventdata, varargin) %
h_main = ft_getopt(varargin, 'h_main');
temp_event = ft_getopt(varargin, 'event');
h = getparent(h_main);
opt = getappdata(h, 'opt');
%cfg = getappdata(h, 'cfg');
if strcmp(opt.markingstatus,'on')
    hax = get(h, 'CurrentAxes');
    
    %[press_x press_y] = ft_select_point([0 0],'nearest',0);
    point = get(hax,'CurrentPoint');     % button down detected
    click_x = point(1,1);
    click_y = point(1,2);
    
    
    
    %xlim_h_init = get(hax, 'xlim');
    %ylim_h_init = get(hax, 'ylim');
    xlim_h_init = opt.xlim_init;
    ylim_h_init = opt.ylim_init;
    xlim_h = get(hax, 'xlim');
    ylim_h = get(hax, 'ylim');
    
    
    if (click_x > xlim_h(1)) && (click_x  < xlim_h(2))% && (over_y > ylim_h(1)) && (over_y  < ylim_h(2))
        
        
        if opt.markSecondClick
            if strcmp(temp_event,'WindowButtonUpFcn')
                %draw line marking
                begsample = opt.trlvis(opt.trlop, 1);
                endsample = opt.trlvis(opt.trlop, 2);
                temp_epochLengthSamples = endsample - begsample + 1;
                
                tim_1 =     opt.hlim(1) + ( ( (opt.curr_first_click_x - xlim_h(1)) / range(xlim_h) ) * range(opt.hlim) );
                tim_2 =     opt.hlim(1) + ( ( (click_x - xlim_h(1)) / range(xlim_h) ) * range(opt.hlim) );
                pos_1 =     opt.curr_first_click_x;
                pos_2 =     click_x;
                ind_1 =     max(1, min(temp_epochLengthSamples,round( ((opt.curr_first_click_x - xlim_h(1)) / range(xlim_h) ) * temp_epochLengthSamples )));
                ind_2 =     max(1, min(temp_epochLengthSamples,round( ((click_x - xlim_h(1))                / range(xlim_h) ) * temp_epochLengthSamples )));
                opt.marks = [opt.marks; [tim_1 tim_2 pos_1 pos_2 ind_1 ind_2] ];
                opt.markSecondClick = 0;
                setappdata(h, 'opt',opt);
                redraw_marks_cb(h);
            end
        else
            if strcmp(temp_event,'WindowButtonDownFcn')
                %safe and wait for second click
                opt.curr_first_click_x = click_x;
                opt.curr_first_click_y = click_y;
                opt.markSecondClick = 1;
                setappdata(h, 'opt',opt);
            end
        end
    end
    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function redraw_marks_cb(h)
cfg = getappdata(h, 'cfg');
opt = getappdata(h, 'opt');

chanindx  = match_str(opt.hdr.label, cfg.channel);

temp_channel_number_in_curr_display = find(chanindx == cfg.score_channel_eeg_number);

cfg.marking_color = [1 0 0];
for iChanDisplayed = temp_channel_number_in_curr_display;
    delete(findobj(h, 'tag', 'markering'));
    for iMark = 1:size(opt.marks,1)
        h_scorechan_eeg_mark_middle_line = ft_plot_line([opt.marks(iMark,1) opt.marks(iMark,2)],[0 0], 'color', cfg.marking_color , 'linewidth', 2, 'tag', 'markering',  ...
            'hpos', opt.laytime.pos(iChanDisplayed,1), 'vpos', opt.laytime.pos(iChanDisplayed,2), 'width', opt.width, 'height', opt.laytime.height(iChanDisplayed), 'hlim', opt.hlim, 'vlim', [-1 1]);
    end
    
    
    
    begsample = opt.trlvis(opt.trlop, 1);
    endsample = opt.trlvis(opt.trlop, 2);
    temp_epochLengthSamples = endsample - begsample + 1;
    
    mark_display_ind = zeros(1,temp_epochLengthSamples);
    for iMark = 1:size(opt.marks,1)
        ind_1 = opt.marks(iMark,5);
        ind_2 = opt.marks(iMark,6);
        if ind_1 > ind_2
            mark_display_ind(ind_2:ind_1) = 1;
        else
            mark_display_ind(ind_1:ind_2) = 1;
        end
    end
    
    cfg.curr_displayed_detected_slowosci_perc_cumulative = num2str(100*sum(cfg.curr_displayed_detected_slowosci_perc_display_ind | mark_display_ind)/temp_epochLengthSamples,3);
    cfg.curr_displayed_detected_slowosci_perc_substractive = num2str(100*sum(cfg.curr_displayed_detected_slowosci_perc_display_ind & ~mark_display_ind)/temp_epochLengthSamples,3);
    cfg.curr_displayed_marking = num2str(100*sum(mark_display_ind)/temp_epochLengthSamples,3);
    
    %     ft_uilayout(h, 'tag', 'scindlabelSOperc', 'string', ['SO%: ' num2str(100*cfg.curr_displayed_detected_slowosci_perc,3) ...
    %         ' /+ ' cfg.curr_displayed_detected_slowosci_perc_cumulative  ...
    %         ' /- ' cfg.curr_displayed_detected_slowosci_perc_substractive ...
    %         ]);
    %
    %     ft_uilayout(h, 'tag', 'scindlabelMarkperc', 'string', ['Mark%: ' cfg.curr_displayed_marking ]);
    
    setappdata(h, 'cfg',cfg);
    updateLabels(h)
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouse_move_cb(hObject, eventdata, varargin) %
h_main = ft_getopt(varargin, 'h_main');
h = getparent(h_main);
opt = getappdata(h, 'opt');
cfg = getappdata(h, 'cfg');

if strcmp(cfg.use_ruler,'yes') || strcmp(opt.markingstatus,'on')
    
    hax = get(h, 'CurrentAxes');
    
    %[press_x press_y] = ft_select_point([0 0],'nearest',0);
    point = get(hax,'CurrentPoint');     % button down detected
    over_x = point(1,1);
    over_y = point(1,2);
    
    
    
    %xlim_h_init = get(hax, 'xlim');
    %ylim_h_init = get(hax, 'ylim');
    xlim_h_init = opt.xlim_init;
    ylim_h_init = opt.ylim_init;
    xlim_h = get(hax, 'xlim');
    ylim_h = get(hax, 'ylim');
    
    delete(findobj(h, 'tag', 'scale_help'));
    delete(findobj(h, 'tag', 'marking_move_line'));
    if (over_x > xlim_h(1)) && (over_x  < xlim_h(2)) && (over_y > ylim_h(1)) && (over_y  < ylim_h(2))
        
        if strcmp(cfg.use_ruler,'yes')
            
            ruler_outer_amplitude_units = cfg.so_thresholdAmplitudeForDetection ;
            ruler_outer_duration_seconds_max = 1/cfg.so_minFreq;
            ruler_outer_duration_seconds_min = 1/cfg.so_maxFreq;
            
            ruler_inner2_amplitude_units = cfg.emg_thresholdAmplitudeForDetection;
            
            ruler_inner_amplitude_units = cfg.sp_thresholdForDetectionCriterion;
            ruler_inner_amplitude_units_min = cfg.sp_thresholdForDetectionBeginEnd;
            ruler_inner_duration_seconds_max = cfg.sp_maxSec;
            ruler_inner_duration_seconds_min = cfg.sp_minSec;
            chanindx  = match_str(opt.hdr.label, cfg.channel);
            
            ncurr_Channels = length(chanindx);
            y_range_channel = range(ylim_h_init)/ncurr_Channels;
            channelIndex = 1;
            for curr_lay = 1:length(chanindx)
                %if (opt.laytime.pos(channelIndex,2) < over_y ) && (opt.laytime.height(channelIndex) < over_y)
                if ((min(ylim_h_init)+((curr_lay-1)*y_range_channel)) < over_y ) && ((min(ylim_h_init)+(curr_lay*y_range_channel)) > over_y)
                    channelIndex = length(chanindx)-curr_lay+1;
                    break;
                end
            end
            
            %channelIndex = 3;
            
            curr_channel = chanindx(channelIndex);
            curr_channel_scaling = cfg.chanscale(curr_channel);
            
            y_range_channel_plot = opt.laytime.height(channelIndex);
            
            x_scale_outer = (ruler_outer_duration_seconds_max/cfg.blocksize)*range(xlim_h_init);
            x_scale_outer_min = (ruler_outer_duration_seconds_min/cfg.blocksize)*range(xlim_h_init);
            
            y_scale_outer = ((ruler_outer_amplitude_units/range(cfg.ylim))*curr_channel_scaling)*y_range_channel_plot;
            
            scale_help_outerlow = ft_plot_line([max(xlim_h(1),over_x) min(xlim_h(2),over_x+x_scale_outer)], [max(ylim_h(1),over_y+(y_scale_outer/2)) min(ylim_h(2),over_y+y_scale_outer/2)],'facealpha',0.9, 'color', [1 0 0] , 'linewidth', 1, 'tag', 'scale_help');
            scale_help_outerhigh = ft_plot_line([max(xlim_h(1),over_x) min(xlim_h(2),over_x+x_scale_outer)],[max(ylim_h(1),over_y-y_scale_outer/2) min(ylim_h(2),over_y-(y_scale_outer/2))],'facealpha',0.9, 'color', [1 0 0] , 'linewidth', 1, 'tag', 'scale_help');
            
            %scale_help_outerlow = ft_plot_box([max(xlim_h(1),over_x) min(xlim_h(2),over_x+x_scale_outer) max(ylim_h(1),over_y+(y_scale_outer/2)-(y_scale_outer/8)) min(ylim_h(2),over_y+y_scale_outer/2)],'facealpha',0.5, 'facecolor', [1 0 0] , 'edgecolor', 'none', 'tag', 'scale_help');
            %scale_help_outerhigh = ft_plot_box([max(xlim_h(1),over_x) min(xlim_h(2),over_x+x_scale_outer) max(ylim_h(1),over_y-y_scale_outer/2) min(ylim_h(2),over_y-(y_scale_outer/2)+(y_scale_outer/8))],'facealpha',0.5, 'facecolor', [1 0 0] , 'edgecolor', 'none', 'tag', 'scale_help');
            
            scale_help_outerlow = ft_plot_box([max(xlim_h(1),over_x) min(xlim_h(2),over_x+x_scale_outer_min) max(ylim_h(1),over_y+(y_scale_outer/2)-(y_scale_outer/8)) min(ylim_h(2),over_y+y_scale_outer/2)],'facealpha',0.5, 'facecolor', [1 0 0] , 'edgecolor', 'none', 'tag', 'scale_help');
            scale_help_outerhigh = ft_plot_box([max(xlim_h(1),over_x) min(xlim_h(2),over_x+x_scale_outer_min) max(ylim_h(1),over_y-y_scale_outer/2) min(ylim_h(2),over_y-(y_scale_outer/2)+(y_scale_outer/8))],'facealpha',0.5, 'facecolor', [1 0 0] , 'edgecolor', 'none', 'tag', 'scale_help');
            
            x_scale_inner = (ruler_inner_duration_seconds_max/cfg.blocksize)*range(xlim_h_init);
            x_scale_inner_min = (ruler_inner_duration_seconds_min/cfg.blocksize)*range(xlim_h_init);
            
            y_scale_inner = ((ruler_inner_amplitude_units/range(cfg.ylim))*curr_channel_scaling)*y_range_channel_plot;
            y_scale_inner_min = ((ruler_inner_amplitude_units_min/range(cfg.ylim))*curr_channel_scaling)*y_range_channel_plot;
            
            y_scale_inner2 = ((ruler_inner2_amplitude_units/range(cfg.ylim))*curr_channel_scaling)*y_range_channel_plot;
            x_scale_inner2 = (1/cfg.blocksize)*range(xlim_h_init);
            
            scale_help_inner_emg_max_top = ft_plot_line([max(xlim_h(1),over_x-x_scale_inner2/2) min(xlim_h(2),over_x+x_scale_inner2/2)], [max(ylim_h(1),over_y+y_scale_inner2/2) min(ylim_h(2),over_y+y_scale_inner2/2)],'facealpha',0.5, 'color', [1 0.7 0] , 'linewidth', 1, 'tag', 'scale_help');
            scale_help_inner_emg_max_bottom = ft_plot_line([max(xlim_h(1),over_x-x_scale_inner2/2) min(xlim_h(2),over_x+x_scale_inner2/2)], [max(ylim_h(1),over_y-y_scale_inner2/2) min(ylim_h(2),over_y-y_scale_inner2/2)],'facealpha',0.5, 'color', [1 0.7 0] , 'linewidth', 1, 'tag', 'scale_help');
            
            
            
            scale_help_inner_max_top = ft_plot_line([max(xlim_h(1),over_x-x_scale_inner/2) min(xlim_h(2),over_x+x_scale_inner/2)], [max(ylim_h(1),over_y+y_scale_inner/2) min(ylim_h(2),over_y+y_scale_inner/2)],'facealpha',0.5, 'color', [0 1 0] , 'linewidth', 1, 'tag', 'scale_help');
            scale_help_inner_max_bottom = ft_plot_line([max(xlim_h(1),over_x-x_scale_inner/2) min(xlim_h(2),over_x+x_scale_inner/2)], [max(ylim_h(1),over_y-y_scale_inner/2) min(ylim_h(2),over_y-y_scale_inner/2)],'facealpha',0.5, 'color', [0 1 0] , 'linewidth', 1, 'tag', 'scale_help');
            
            scale_help_inner_min_top = ft_plot_line([max(xlim_h(1),over_x-x_scale_inner/2) min(xlim_h(2),over_x+x_scale_inner/2)], [max(ylim_h(1),over_y+y_scale_inner_min/2) min(ylim_h(2),over_y+y_scale_inner_min/2)],'facealpha',0.5, 'color', [0 1 0] , 'linewidth', 1, 'tag', 'scale_help');
            scale_help_inner_min_bottom = ft_plot_line([max(xlim_h(1),over_x-x_scale_inner/2) min(xlim_h(2),over_x+x_scale_inner/2)], [max(ylim_h(1),over_y-y_scale_inner_min/2) min(ylim_h(2),over_y-y_scale_inner_min/2)],'facealpha',0.5, 'color', [0 1 0] , 'linewidth', 1, 'tag', 'scale_help');
            
            scale_help_inner = ft_plot_box([max(xlim_h(1),over_x-x_scale_inner_min/2) min(xlim_h(2),over_x+x_scale_inner_min/2) max(ylim_h(1),over_y-y_scale_inner_min/2) min(ylim_h(2),over_y+y_scale_inner_min/2)],'facealpha',0.5, 'facecolor', [0 1 0] , 'edgecolor', 'none', 'tag', 'scale_help');
            %scale_help_inner = ft_plot_box([max(xlim_h(1),over_x-x_scale_inner/2) min(xlim_h(2),over_x+x_scale_inner/2) max(ylim_h(1),over_y-y_scale_inner_min/2) min(ylim_h(2),over_y+y_scale_inner_min/2)],'facealpha',0.5, 'facecolor', [0 1 0] , 'edgecolor', 'none', 'tag', 'scale_help');
            
            
        end
        
        if strcmp(opt.markingstatus,'on')
            
            if opt.markSecondClick
                %draw curr line marking
                pos_x =  opt.curr_first_click_x;
                pos_y =  opt.curr_first_click_y;
                h_marking_move_line = ft_plot_line([pos_x over_x],[pos_y over_y], 'color', [1 0 0] , 'linewidth', 1, 'tag', 'marking_move_line');
                
                %else
            end
            
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MAstr] = getStageStringByHypnValue_MA(h2)
MAstr = '';
switch h2
    case 0
        MAstr = '';
    case {1 2}
        MAstr = ' A';
    case 3
        MAstr = ' a';
    otherwise
        MAstr = ' ?A?';
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stagestring] = getStageStringByHypnValue_stage(h1)
stagestring = '';
switch h1
    case 0
        stagestring = 'W';
    case {1 2 3 4}
        stagestring = ['' num2str(h1)];
    case 5
        stagestring = 'R';
    case 8
        stagestring = 'M';
    otherwise
        stagestring = '?';
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stagestring h1_str h2_str] = getStageStringByHypnValue(h1,h2)
h1_str = getStageStringByHypnValue_stage(h1);
h2_str = getStageStringByHypnValue_MA(h2);
stagestring = [h1_str h2_str];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stagesstring] = getPrevStageString_stage(hypn,curr_epoch,n)
idx = max(1,curr_epoch-n):max(0,(curr_epoch-1));
hv = hypn(idx,1);
if isempty(hv)
    stagesstring = '';
else
    stagesstring = [arrayfun(@(h)getStageStringByHypnValue_stage(h),hv')];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stagesstring] = getNextStageString_stage(hypn,curr_epoch,n)
idx = min(size(hypn,1),curr_epoch+1):min(size(hypn,1),(curr_epoch+n));
hv = hypn(idx,1);
if isempty(hv)
    stagesstring = '';
else
    stagesstring = [arrayfun(@(h)getStageStringByHypnValue_stage(h),hv')];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cfg = update_filters(cfg)

FpassLeft = cfg.sp_minFreq; %left pass frequency in Hz
FpassRight = cfg.sp_maxFreq; %right pass frequency in Hz

FstopLeft = FpassLeft - cfg.StopToPassTransitionWidth_bp; %left stop frequency in Hz
FstopRight = FpassRight + cfg.PassToStopTransitionWidth_bp; %left stop frequency in Hz

%usedFilterOrder_bp = NaN;
%bp_hdm = NaN;
if strcmp(cfg.core_cfg.bpfilttype,'IIRdesigned') || strcmp(cfg.core_cfg.bpfilttype,'FIRdesigned')
    bp_d = [];
    bp_hd = [];
    if strcmp(cfg.UseFixedFilterOrder_bp,'yes')
        bp_d = fdesign.bandpass('N,Fst1,Fp1,Fp2,Fst2',cfg.FilterOrder_bp,FstopLeft,FpassLeft,FpassRight,FstopRight,cfg.FrqOfSmpl);
        bp_hd = design(bp_d,'equiripple');
    else
        bp_d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',FstopLeft,FpassLeft,FpassRight,FstopRight,cfg.AstopLeft_bp,cfg.Apass_bp,cfg.AstopRight_bp,cfg.FrqOfSmpl);
        bp_hd = design(bp_d,'equiripple','MinOrder', 'even');
    end
    %usedFilterOrder_bp = bp_hd.order;
    cfg.sp_bpfilterdesign = bp_hd;
    %bp_hdm = measure(bp_hd);
end


FpassLeft = cfg.al_minFreq; %left pass frequency in Hz
FpassRight = cfg.al_maxFreq; %right pass frequency in Hz

FstopLeft = FpassLeft - cfg.StopToPassTransitionWidth_bp; %left stop frequency in Hz
FstopRight = FpassRight + cfg.PassToStopTransitionWidth_bp; %left stop frequency in Hz

%usedFilterOrder_bp = NaN;
%bp_hdm = NaN;
if strcmp(cfg.core_cfg.bpfilttype,'IIRdesigned') || strcmp(cfg.core_cfg.bpfilttype,'FIRdesigned')
    bp_d = [];
    bp_hd = [];
    if strcmp(cfg.UseFixedFilterOrder_bp,'yes')
        bp_d = fdesign.bandpass('N,Fst1,Fp1,Fp2,Fst2',cfg.FilterOrder_bp,FstopLeft,FpassLeft,FpassRight,FstopRight,cfg.FrqOfSmpl);
        bp_hd = design(bp_d,'equiripple');
    else
        bp_d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',FstopLeft,FpassLeft,FpassRight,FstopRight,cfg.AstopLeft_bp,cfg.Apass_bp,cfg.AstopRight_bp,cfg.FrqOfSmpl);
        bp_hd = design(bp_d,'equiripple','MinOrder', 'even');
    end
    %usedFilterOrder_bp = bp_hd.order;
    cfg.al_bpfilterdesign = bp_hd;
    %bp_hdm = measure(bp_hd);
end




FpassLeft = cfg.so_filter_minFreq; %left pass frequency in Hz
FstopLeft = FpassLeft - cfg.StopToPassTransitionWidth_hp; %left stop frequency in Hz

%     usedFilterOrder_hp = NaN;
%     hp_hdm = NaN;
if strcmp(cfg.core_cfg.hpfilttype,'IIRdesigned') || strcmp(cfg.core_cfg.hpfilttype,'FIRdesigned')
    hp_d = [];
    hp_hd = [];
    if strcmp(cfg.UseFixedFilterOrder_hp,'yes')
        hp_d = fdesign.highpass('N,F3db',cfg.FilterOrder_hp,FpassLeft,cfg.FrqOfSmpl);%designmethods(hp_preDS_d); help(hp_preDS_d,'equiripple'); help(hp_preDS_d,'butter')
    else
        hp_d = fdesign.highpass('Fst,Fp,Ast,Ap',FstopLeft,FpassLeft,cfg.AstopLeft_hp,cfg.Apass_hp,cfg.FrqOfSmpl);%designmethods(hp_preDS_d); help(hp_preDS_d,'equiripple'); help(hp_preDS_d,'butter')
    end
    if strcmp(cfg.core_cfg.hpfilttype,'IIRdesigned')
        hp_hd = design(hp_d,'butter'); %isstable(hp_hd)
    elseif strcmp(cfg.core_cfg.hpfilttype,'FIRdesigned')
        hp_hd = design(hp_d,'equiripple','MinOrder', 'even');
    else
        error(['highpass filter type of ' cfg.core_cfg.hpfilttype ' unknown or not allowed'])
    end
    %         usedFilterOrder_hp = hp_hd.order;
    cfg.so_hpfilterdesign = hp_hd;
    %         hp_hdm = measure(hp_hd);
end

FpassRight = cfg.so_filter_maxFreq; %right pass frequency in Hz
FstopRight = FpassRight + cfg.PassToStopTransitionWidth_lp; %right stop frequency in Hz
%     usedFilterOrder_lp = NaN;
%     lp_hdm = NaN;
if strcmp(cfg.core_cfg.lpfilttype,'IIRdesigned') || strcmp(cfg.core_cfg.lpfilttype,'FIRdesigned')
    lp_d = [];
    lp_hd = [];
    if strcmp(cfg.UseFixedFilterOrder_lp,'yes')
        lp_d = fdesign.lowpass('N,Fp,Fst',cfg.FilterOrder_lp,FpassRight,FstopRight,cfg.FrqOfSmpl);
        lp_hd = design(lp_d,'equiripple');
    else
        lp_d = fdesign.lowpass('Fp,Fst,Ap,Ast',FpassRight,FstopRight,cfg.Apass_lp,cfg.AstopRight_lp,cfg.FrqOfSmpl);
        lp_hd = design(lp_d,'equiripple','MinOrder', 'even');
    end
    %         usedFilterOrder_lp = lp_hd.order;
    cfg.so_lpfilterdesign = lp_hd;
    %         lp_hdm = measure(lp_hd);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateLabels(h)
cfg = getappdata(h, 'cfg');
if (cfg.curr_displayed_detected_slowosci_perc ~= -1)
    %ft_uilayout(h, 'tag', 'scindlabelSOnum', 'string', ['#SO: ' num2str(cfg.curr_displayed_detected_slowosci_number,3)]);
    if (cfg.curr_displayed_marking ~= -1)
        ft_uilayout(h, 'tag', 'scindlabelSOperc', 'string', ['SO%: ' num2str(100*cfg.curr_displayed_detected_slowosci_perc,3) ...
            ' /+ ' cfg.curr_displayed_detected_slowosci_perc_cumulative  ...
            ' /- ' cfg.curr_displayed_detected_slowosci_perc_substractive ...
            ' (#' num2str(cfg.curr_displayed_detected_slowosci_number,3) ')' ...
            ]);
    else
        ft_uilayout(h, 'tag', 'scindlabelSOperc', 'string', ['SO%: ' num2str(100*cfg.curr_displayed_detected_slowosci_perc,3) ' (#' num2str(cfg.curr_displayed_detected_slowosci_number,3) ')']);
    end
else
    ft_uilayout(h, 'tag', 'scindlabelSOperc', 'string', ['SO%: ' '?' ' (#' '?' ')']);
    %ft_uilayout(h, 'tag', 'scindlabelSOnum', 'string', ['#SO: ' '?']);
end
if (cfg.curr_displayed_detected_spindels_number ~= -1)
    ft_uilayout(h, 'tag', 'scindlabelSpnum', 'string', ['#Sp: ' num2str(cfg.curr_displayed_detected_spindels_number)]);
else
    ft_uilayout(h, 'tag', 'scindlabelSpnum', 'string', ['#Sp: ' '?']);
end
if (cfg.curr_displayed_marking ~= -1)
    ft_uilayout(h, 'tag', 'scindlabelMarkperc', 'string', ['Mark%: ' cfg.curr_displayed_marking ]);
else
    ft_uilayout(h, 'tag', 'scindlabelMarkperc', 'string', ['Mark%: ' '?' ]);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [channels, scalings, order, enabled, focus_eeg, focus_emg, focus_eog] = channelDialog(channels)
function [selected_channels,cfg,opt] = channelDialog(channels,indices_selected,cfg,opt)
%selected_channels = channels;
selected_channels = channels(indices_selected);


Nchannels = numel(channels);
titlestr = 'Channel Settings';
pos      = get(0,'DefaultFigurePosition');
temp_screensize = get(0,'screensize');
temp_dlgSize = floor(temp_screensize(3:4)/2);
pos(3:4) = temp_dlgSize;

vertical_size_row_pre = 0.075;
vertical_size_total = max(1,(Nchannels+1)*vertical_size_row_pre);

vertical_size_row = vertical_size_row_pre* (1/vertical_size_total);

%dlg = dialog('Name', titlestr, 'Position', pos);
dlg = figure('Name',titlestr,'NumberTitle','off','MenuBar','none');
axis off % explicitly turn of axis, sometimes an axis system appears

panel1 = uipanel('Parent',dlg);
panel2 = uipanel('Parent',panel1);
set(panel1,'Position',[0 0 0.95 1]);
%set(panel2,'Position',[0 1-vertical_size_total 1 vertical_size_total]);
set(panel2,'Position',[0 1-vertical_size_total 1 vertical_size_total]);
set(gca,'Parent',panel2);
if vertical_size_total > 1
    slider = uicontrol('Style','Slider','Parent',dlg,...
        'Units','normalized','Position',[0.95 0 0.05 1],...
        'Value',1,'Callback',{@(src,eventdata,nch,vsizech,vsize,panel) slider_callback1(src,eventdata,nch,vsizech,vsize,panel), Nchannels, vertical_size_row_pre, vertical_size_total, panel2});
       %'Value',1,'Callback',{@(src,eventdata,arg1,arg2) set(arg1,'Position',[0 -get(src,'Value')*Nchannels/(1/vertical_size_row_pre) 1 arg2]), panel2, vertical_size_total});
    
end
indices_selected = indices_selected(:)';     % ensure that it is a row array

% userdata.label    = channels;
% userdata.select   = index_selected;
index_unselected = setdiff(1:length(channels), indices_selected);
% set(dlg, 'userdata', userdata);
% uicontrol(dlg, 'style', 'text',       'position', [ 10 240+20 80  20], 'string', 'unselected');
% uicontrol(dlg, 'style', 'text',       'position', [200 240+20 80  20], 'string', 'selected  ');
% uicontrol(dlg, 'style', 'listbox',    'position', [ 10  40+20 80 200], 'min', 0, 'max', 2, 'tag', 'lbunsel')
% uicontrol(dlg, 'style', 'listbox',    'position', [200  40+20 80 200], 'min', 0, 'max', 2, 'tag', 'lbsel')
% uicontrol(dlg, 'style', 'pushbutton', 'position', [105 175+20 80  20], 'string', 'add all >'   , 'callback', @label_addall);
% uicontrol(dlg, 'style', 'pushbutton', 'position', [105 145+20 80  20], 'string', 'add >'       , 'callback', @label_add);
% uicontrol(dlg, 'style', 'pushbutton', 'position', [105 115+20 80  20], 'string', '< remove'    , 'callback', @label_remove);
% uicontrol(dlg, 'style', 'pushbutton', 'position', [105  85+20 80  20], 'string', '< remove all', 'callback', @label_removeall);
% uicontrol(dlg, 'style', 'pushbutton', 'position', [ 55  10    80  20], 'string', 'Cancel',       'callback', 'close');
% uicontrol(dlg, 'style', 'pushbutton', 'position', [155  10    80  20], 'string', 'OK',           'callback', 'uiresume');
% line_x_offsets = [0   0.2 0.3 0.4 0.5 0.6  0.65 0.7  0.75 0.83];
% line_x_width =   [0.2 0.1 0.1 0.1 0.1 0.05 0.05 0.05 0.05 0.17];

line_x_offsets = [0.00 0.15 0.250 0.325 0.40 0.50 0.55 0.60 0.65 0.70 0.75 0.83];
line_x_width =   [0.15 0.10 0.075 0.075 0.10 0.05 0.05 0.05 0.05 0.05 0.05 0.17];
uicontrol(panel2, 'style', 'text', 'units', 'normalized','HorizontalAlignment','left', 'position', [line_x_offsets(1) 1-1*vertical_size_row line_x_width(1), vertical_size_row], 'string', 'channel','tag',['heading1']);
uicontrol(panel2, 'style', 'text', 'units', 'normalized','HorizontalAlignment','left', 'position', [line_x_offsets(2) 1-1*vertical_size_row line_x_width(2), vertical_size_row], 'string', 'visible','tag',['heading2']);
uicontrol(panel2, 'style', 'text', 'units', 'normalized','HorizontalAlignment','left', 'position', [line_x_offsets(3) 1-1*vertical_size_row line_x_width(3), vertical_size_row], 'string', 'scaling','tag',['heading3']);
uicontrol(panel2, 'style', 'text', 'units', 'normalized','HorizontalAlignment','left', 'position', [line_x_offsets(4) 1-1*vertical_size_row line_x_width(4), vertical_size_row], 'string', 'color','tag',['heading4']);
uicontrol(panel2, 'style', 'text', 'units', 'normalized','HorizontalAlignment','left', 'position', [line_x_offsets(5) 1-1*vertical_size_row line_x_width(5), vertical_size_row], 'string', 'chOrder','tag',['heading5']);

uicontrol(panel2, 'style', 'text', 'units', 'normalized','HorizontalAlignment','left', 'position', [line_x_offsets(6) 1-1*vertical_size_row line_x_width(6), vertical_size_row], 'string', 'EEG','tag',['heading6'],'backgroundcolor', cfg.score_channel_eeg_color);
uicontrol(panel2, 'style', 'text', 'units', 'normalized','HorizontalAlignment','left', 'position', [line_x_offsets(7) 1-1*vertical_size_row line_x_width(7), vertical_size_row], 'string', 'EEG_F','tag',['heading7'],'backgroundcolor', cfg.score_channel_eeg_frontal_color);
uicontrol(panel2, 'style', 'text', 'units', 'normalized','HorizontalAlignment','left', 'position', [line_x_offsets(8) 1-1*vertical_size_row line_x_width(8), vertical_size_row], 'string', 'EEG_O','tag',['heading8'],'backgroundcolor', cfg.score_channel_eeg_occipital_color);
uicontrol(panel2, 'style', 'text', 'units', 'normalized','HorizontalAlignment','left', 'position', [line_x_offsets(9) 1-1*vertical_size_row line_x_width(9), vertical_size_row], 'string', 'EOG','tag',['heading9'],'backgroundcolor', cfg.score_channel_eog_color);
uicontrol(panel2, 'style', 'text', 'units', 'normalized','HorizontalAlignment','left', 'position', [line_x_offsets(10) 1-1*vertical_size_row line_x_width(10), vertical_size_row], 'string', 'EMG','tag',['heading10'],'backgroundcolor', cfg.score_channel_emg_color);
bg_EEG           = uibuttongroup('Visible','off', 'units', 'normalized','Position',[line_x_offsets(6) 1-(numel(channels)+1)*vertical_size_row line_x_width(6), vertical_size_row*numel(channels)],'tag',['rdg_EEG'],'backgroundcolor', cfg.score_channel_eeg_color,'parent',panel2);
bg_EEG_frontal   = uibuttongroup('Visible','off', 'units', 'normalized','Position',[line_x_offsets(7) 1-(numel(channels)+1)*vertical_size_row line_x_width(7), vertical_size_row*numel(channels)],'tag',['rdg_EEG_frontal'],'backgroundcolor', cfg.score_channel_eeg_frontal_color,'parent',panel2);
bg_EEG_occipital = uibuttongroup('Visible','off', 'units', 'normalized','Position',[line_x_offsets(8) 1-(numel(channels)+1)*vertical_size_row line_x_width(8), vertical_size_row*numel(channels)],'tag',['rdg_EEG_occipital'],'backgroundcolor', cfg.score_channel_eeg_occipital_color,'parent',panel2);
bg_EOG           = uibuttongroup('Visible','off', 'units', 'normalized','Position',[line_x_offsets(9) 1-(numel(channels)+1)*vertical_size_row line_x_width(9), vertical_size_row*numel(channels)],'tag',['rdg_EOG'],'backgroundcolor', cfg.score_channel_eog_color,'parent',panel2);
bg_EMG           = uibuttongroup('Visible','off', 'units', 'normalized','Position',[line_x_offsets(10) 1-(numel(channels)+1)*vertical_size_row line_x_width(10), vertical_size_row*numel(channels)],'tag',['rdg_EMG'],'backgroundcolor', cfg.score_channel_emg_color,'parent',panel2);

if cfg.has_ECG
    uicontrol(panel2, 'style', 'text', 'units', 'normalized','HorizontalAlignment','left', 'position', [line_x_offsets(11) 1-1*vertical_size_row line_x_width(11), vertical_size_row], 'string', 'ECG','tag',['heading11'],'backgroundcolor', cfg.score_channel_ecg_color);
    bg_ECG = uibuttongroup('Visible','off', 'units', 'normalized','Position',[line_x_offsets(11) 1-(numel(channels)+1)*vertical_size_row line_x_width(11), vertical_size_row*numel(channels)],'tag',['rdg_ECG'],'backgroundcolor', cfg.score_channel_ecg_color,'parent',panel2);
end

uicontrol(panel2, 'style', 'pushbutton', 'units', 'normalized', 'position', [line_x_offsets(12) 1-1*vertical_size_row line_x_width(12), vertical_size_row], 'string', 'OK','tag',['button_OK'],'Callback','uiresume');
uicontrol(panel2, 'style', 'pushbutton', 'units', 'normalized', 'position', [line_x_offsets(12) 1-2*vertical_size_row line_x_width(12), vertical_size_row], 'string', 'Cancel','tag',['button_Cancel'],'Callback','close');
%uicontrol(panel2, 'style', 'pushbutton', 'units', 'normalized', 'position', [line_x_offsets(12) 1-2*vertical_size_row line_x_width(12), vertical_size_row], 'string', 'Cancel','tag',['button_Cancel'],'Callback','close');
uicontrol(panel2, 'style', 'pushbutton', 'units', 'normalized', 'position', [line_x_offsets(12) 1-3*vertical_size_row line_x_width(12), vertical_size_row], 'string', '(ch colors)','tag',['button_Color'],'Callback',{@cb_channelDialog_color_channels, numel(channels)});

uicontrol(panel2, 'style', 'text', 'units', 'normalized','HorizontalAlignment','left', 'position', [line_x_offsets(12) 1-4*vertical_size_row line_x_width(12), vertical_size_row], 'string', 'all chan.?','tag',['text_all_chan']);
uicontrol(panel2, 'style', 'checkbox', 'units', 'normalized', 'position', [line_x_offsets(12) 1-5*vertical_size_row+vertical_size_row/2 line_x_width(12), vertical_size_row],'tag',['enabled_all_chan'],'Callback',{@cb_channelDialog_all_channels});



for iCh = 1:numel(channels)
    uicontrol(panel2, 'style', 'text', 'units', 'normalized','HorizontalAlignment','left', 'position', [line_x_offsets(1) 1-(iCh+1)*vertical_size_row line_x_width(1), vertical_size_row], 'string', channels(iCh),'tag',['text_chan' num2str(iCh)]);
    uicontrol(panel2, 'style', 'checkbox', 'units', 'normalized', 'position', [line_x_offsets(2) 1-(iCh+1)*vertical_size_row line_x_width(2), vertical_size_row],'tag',['enabled_chan' num2str(iCh)]);
    uicontrol(panel2, 'style', 'edit', 'units', 'normalized', 'position', [line_x_offsets(3) 1-(iCh+1)*vertical_size_row line_x_width(3), vertical_size_row], 'string', num2str(cfg.chanscale(iCh)),'tag',['scale_chan' num2str(iCh)],'Min',0);
    %uicontrol(panel2, 'style', 'edit', 'units', 'normalized', 'position', [line_x_offsets(4) 1-(iCh+1)*vertical_size_row line_x_width(4), vertical_size_row], 'string', '','tag',['order_chan' num2str(iCh)],'Min',1);
    uicontrol(panel2, 'style', 'pushbutton', 'units', 'normalized', 'position', [line_x_offsets(4) 1-(iCh+1)*vertical_size_row line_x_width(4), vertical_size_row], 'string', 'color','ForegroundColor',opt.chancolors(iCh,:),'tag',['color_chan' num2str(iCh)],'Callback',{@cb_channelDialog_Colorchooser,opt.chancolors(iCh,:)});
    uicontrol(panel2, 'style', 'edit', 'units', 'normalized', 'position', [line_x_offsets(5) 1-(iCh+1)*vertical_size_row line_x_width(5), vertical_size_row], 'string', num2str(iCh),'tag',['order_chan' num2str(iCh)],'Min',0);
    
    uicontrol(panel2, 'style', 'radiobutton', 'units', 'normalized', 'position', [0.1 1-(iCh)*vertical_size_row/(vertical_size_row*numel(channels))+vertical_size_row*0.25 1-0.1 vertical_size_row ], 'string', '','tag',['radiobutton_focusEEG' num2str(iCh)],'parent',bg_EEG,'backgroundcolor', cfg.score_channel_eeg_color);
    uicontrol(panel2, 'style', 'radiobutton', 'units', 'normalized', 'position', [0.1 1-(iCh)*vertical_size_row/(vertical_size_row*numel(channels))+vertical_size_row*0.25 1-0.1 vertical_size_row ], 'string', '','tag',['radiobutton_focusEEG_frontal' num2str(iCh)],'parent',bg_EEG_frontal,'backgroundcolor', cfg.score_channel_eeg_frontal_color);
    uicontrol(panel2, 'style', 'radiobutton', 'units', 'normalized', 'position', [0.1 1-(iCh)*vertical_size_row/(vertical_size_row*numel(channels))+vertical_size_row*0.25 1-0.1 vertical_size_row ], 'string', '','tag',['radiobutton_focusEEG_occipital' num2str(iCh)],'parent',bg_EEG_occipital,'backgroundcolor', cfg.score_channel_eeg_occipital_color);
    uicontrol(panel2, 'style', 'radiobutton', 'units', 'normalized', 'position', [0.1 1-(iCh)*vertical_size_row/(vertical_size_row*numel(channels))+vertical_size_row*0.25 1-0.1 vertical_size_row ], 'string', '','tag',['radiobutton_focusEOG' num2str(iCh)],'parent',bg_EOG,'backgroundcolor', cfg.score_channel_eog_color);
    uicontrol(panel2, 'style', 'radiobutton', 'units', 'normalized', 'position', [0.1 1-(iCh)*vertical_size_row/(vertical_size_row*numel(channels))+vertical_size_row*0.25 1-0.1 vertical_size_row ], 'string', '','tag',['radiobutton_focusEMG' num2str(iCh)],'parent',bg_EMG,'backgroundcolor', cfg.score_channel_emg_color);
    if cfg.has_ECG
        uicontrol(panel2, 'style', 'radiobutton', 'units', 'normalized', 'position', [0.1 1-(iCh)*vertical_size_row/(vertical_size_row*numel(channels))+vertical_size_row*0.25 1-0.1 vertical_size_row ], 'string', '','tag',['radiobutton_focusECG' num2str(iCh)],'parent',bg_ECG,'backgroundcolor', cfg.score_channel_ecg_color);
    end
    
    if ismember(iCh,indices_selected)
        ft_uilayout(panel2, 'tag', ['enabled_chan' num2str(iCh)], 'value', 1);
    end
end


set(bg_EEG,'SelectedObject',findobj(get(bg_EEG,'children'),'tag',['radiobutton_focusEEG' num2str(cfg.score_channel_eeg_number)]))
set(bg_EEG_frontal,'SelectedObject',findobj(get(bg_EEG_frontal,'children'),'tag',['radiobutton_focusEEG_frontal' num2str(cfg.score_channel_eeg_frontal_number)]))
set(bg_EEG_occipital,'SelectedObject',findobj(get(bg_EEG_occipital,'children'),'tag',['radiobutton_focusEEG_occipital' num2str(cfg.score_channel_eeg_occipital_number)]))
set(bg_EOG,'SelectedObject',findobj(get(bg_EOG,'children'),'tag',['radiobutton_focusEOG' num2str(cfg.score_channel_eog_number)]))
set(bg_EMG,'SelectedObject',findobj(get(bg_EMG,'children'),'tag',['radiobutton_focusEMG' num2str(cfg.score_channel_emg_number)]))

if cfg.has_ECG
    set(bg_ECG,'SelectedObject',findobj(get(bg_ECG,'children'),'tag',['radiobutton_focusECG' num2str(cfg.score_channel_ecg_number)]))
end




set(bg_EEG,'Visible','on')
set(bg_EEG_frontal,'Visible','on')
set(bg_EEG_occipital,'Visible','on')
set(bg_EOG,'Visible','on')
set(bg_EMG,'Visible','on')
if cfg.has_ECG
    set(bg_ECG,'Visible','on')
end


drawnow
% label_redraw(dlg);
% wait untill the dialog is closed or the user presses OK/Cancel
uiwait(dlg);
if ishandle(dlg)
    % the user pressed OK, return the selection from the dialog
    %   userdata = get(dlg, 'userdata');
    %   select = userdata.select;
    
    
    channel_order = (1:numel(channels))';
    channel_order_compare = channel_order;
    for iCh = 1:numel(channels)
        value = str2num( get(findobj(get(dlg,'children'),'tag',['order_chan' num2str(iCh)]),'String') );
        if ~isempty(value)
            channel_order(iCh) = value;
        end
    end
    
    if ~all(channel_order_compare == channel_order)
        %channel_order = (1:numel(opt.hdr.label))';
        %channel_order(2) = 0.5;
        opt = changeDataChannelOrder(channel_order,opt);
        %cfg.channel = ft_channelselection(cfg.channel, opt.hdr.label);
    end
    
    
    [dummy_chanOrder curr_chanIndexOrder] = sort(channel_order);
    
    chanNum_focusEEG = 1;
    for iCh = 1:numel(channels)
        if strcmp(get(get(bg_EEG,'SelectedObject'),'Tag'),['radiobutton_focusEEG' num2str(iCh)])
            chanNum_focusEEG = iCh;
            break;
        end
    end
    
    chanNum_focusEEG_frontal = 1;
    for iCh = 1:numel(channels)
        if strcmp(get(get(bg_EEG_frontal,'SelectedObject'),'Tag'),['radiobutton_focusEEG_frontal' num2str(iCh)])
            chanNum_focusEEG_frontal = iCh;
            break;
        end
    end
    
    chanNum_focusEEG_occipital = 1;
    for iCh = 1:numel(channels)
        if strcmp(get(get(bg_EEG_occipital,'SelectedObject'),'Tag'),['radiobutton_focusEEG_occipital' num2str(iCh)])
            chanNum_focusEEG_occipital = iCh;
            break;
        end
    end
    
    chanNum_focusEOG = 1;
    for iCh = 1:numel(channels)
        if strcmp(get(get(bg_EOG,'SelectedObject'),'Tag'),['radiobutton_focusEOG' num2str(iCh)])
            chanNum_focusEOG = iCh;
            break;
        end
    end
    
    chanNum_focusEMG = 1;
    for iCh = 1:numel(channels)
        if strcmp(get(get(bg_EMG,'SelectedObject'),'Tag'),['radiobutton_focusEMG' num2str(iCh)])
            chanNum_focusEMG = iCh;
            break;
        end
    end
    
    if cfg.has_ECG
        chanNum_focusECG = 1;
        for iCh = 1:numel(channels)
            if strcmp(get(get(bg_ECG,'SelectedObject'),'Tag'),['radiobutton_focusECG' num2str(iCh)])
                chanNum_focusECG = iCh;
                break;
            end
        end
    end
    
    
    
    cfg.score_channel_eeg_number = find(curr_chanIndexOrder == chanNum_focusEEG,1,'first');
    %ft_uilayout(h, 'tag', 'scoptbuttons_focusEEG', 'string', ['Focus EEG: ' opt.hdr.label{cfg.score_channel_eeg_number}]);
    
    cfg.score_channel_eeg_frontal_number = find(curr_chanIndexOrder == chanNum_focusEEG_frontal,1,'first');
    %ft_uilayout(h, 'tag', 'scoptbuttons_focusEEG', 'string', ['Focus EEG: ' opt.hdr.label{cfg.score_channel_eeg_number}]);
    
    cfg.score_channel_eeg_occipital_number = find(curr_chanIndexOrder == chanNum_focusEEG_occipital ,1,'first');
    %ft_uilayout(h, 'tag', 'scoptbuttons_focusEEG', 'string', ['Focus EEG: ' opt.hdr.label{cfg.score_channel_eeg_number}]);
    
    cfg.score_channel_eog_number = find(curr_chanIndexOrder == chanNum_focusEOG,1,'first');
    %ft_uilayout(h, 'tag', 'scoptbuttons_focusEOG', 'string', ['Focus EOG: ' opt.hdr.label{cfg.score_channel_eog_number}]);
    
    cfg.score_channel_emg_number = find(curr_chanIndexOrder == chanNum_focusEMG,1,'first');
    %ft_uilayout(h, 'tag', 'scoptbuttons_focusEMG', 'string', ['Focus EMG: ' opt.hdr.label{cfg.score_channel_emg_number}]);
    
    if cfg.has_ECG
        cfg.score_channel_ecg_number = find(curr_chanIndexOrder == chanNum_focusECG,1,'first');
    end
    
    if isfield(cfg,'begin_end_events')
        cfg.times_ind_per_channel = cfg.times_ind_per_channel(curr_chanIndexOrder);
        cfg.begin_end_events = cfg.begin_end_events(curr_chanIndexOrder);
    end
    
    if isfield(cfg,'begin_end_events2')
        cfg.times_ind_per_channel2 = cfg.times_ind_per_channel2(curr_chanIndexOrder);
        cfg.begin_end_events2 = cfg.begin_end_events2(curr_chanIndexOrder);
    end
    
    index_selected_channels = [];
    for iCh = 1:numel(channels)
        if get(findobj(get(dlg,'children'),'tag',['enabled_chan' num2str(iCh)]),'value')
            index_selected_channels = [index_selected_channels;iCh];
        end
    end
    
    for iCh = 1:numel(channels)
        value = str2num( get(findobj(get(dlg,'children'),'tag',['scale_chan' num2str(iCh)]),'String') );
        if ~isempty(value)
            cfg.chanscale(find(curr_chanIndexOrder == iCh,1,'first')) = value;
        end
    end
    
    selected_channels = channels(index_selected_channels);
    
    for iCh = 1:numel(channels)
        color = get(findobj(get(dlg,'children'),'tag',['color_chan' num2str(iCh)]),'ForegroundColor');
        if ~isempty(value)
            opt.chancolors(find(curr_chanIndexOrder == iCh,1,'first'),:) = color;
        end
    end
    
    close(dlg);
    return
else
    % the user pressed Cancel or closed the dialog, return the initial selection
    return
end




end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_channelDialog_all_channels(src,eventdata,handles)
val = get(src,'Value');
for obj = findobj('-regexp','tag',['enabled_chan.*'])
    set(obj,'Value',val);
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_channelDialog_color_channels(src,eventdata,nChan)
colSchemes = {'black' 'jet' '1-jet' 'hsv' '1-hsv' 'copper' 'gray' '1-gray' 'bone' 'lines' 'prism' 'white' 'mono'};
nColSchemes = numel(colSchemes);
for obj = findobj('tag',['button_Color'])
    colscheme = get(obj,'string');
    colschemeInd = find(strcmp(colscheme,colSchemes));
    if isempty(colschemeInd) || (colschemeInd == nColSchemes)
        colschemeInd = 1;
    else
        colschemeInd = colschemeInd + 1;
    end
    newColScheme = colSchemes(colschemeInd);
    cols = zeros(nChan,3);
    switch newColScheme{:}
        case 'jet'
            cols = jet(nChan);
        case '1-jet'
            cols = 1-jet(nChan);
        case 'hsv'
            cols = hsv(nChan);
        case '1-hsv'
            cols = 1-hsv(nChan);
        case 'copper'
            cols = copper(nChan);
        case 'gray'
            cols = gray(nChan);
        case '1-gray'
            cols = 1-gray(nChan);
        case 'bone'
            cols = bone(nChan);
        case 'lines'
            cols = lines(nChan);
        case 'prism'
            cols = prism(nChan);
        case 'white'
            cols = white(nChan);
        case 'mono'
            color_new = uisetcolor([0.5 0.5 0.5],'Color for all channels');
            cols = repmat(color_new,nChan,1);
            newColScheme = num2str(round(color_new*255));
    end
    chanobjs = findobj('-regexp','tag',['color_chan.*']);
    for iChan = 1:numel(chanobjs)
        set(chanobjs(iChan),'ForegroundColor',cols(iChan,:));
    end
    set(obj,'string',newColScheme);
end
end


% function  cb_channelDialog_enable_checkbox(source, eventdata, handles)
% %new_value = setxor(get(source,'Value'),[0 1]);
% set(source,'Value', get(source,'Value'))
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % SUBFUNCTION
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function cb_channelDialog_OK(src,eventdata,arg1,arg2)
% uiresume
% end
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % SUBFUNCTION
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function cb_channelDialog_Cancel(src,eventdata,arg1,arg2)
%
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_channelDialog_Colorchooser(src,eventdata,arg1)
color = arg1;
color_new = uisetcolor(color,'new channel color');
set(src,'ForegroundColor',color_new);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider_callback1(src,eventdata,nch,vsizech,vsize,panel)
val = get(src,'Value');
%set(arg1,'Position',[0 -val 1 arg2])
%set(panel,'Position',[0 -val*nch/(1/vsizech) 1 vsize]);
change_range = vsize-1;
set(panel,'Position',[0 (1-vsize)+(1-val)*change_range 1 vsize]);
end

function opt = changeDataChannelOrder(channel_order,opt)
%channel_order(2) = 3.5
[dummy_chanOrder curr_chanIndexOrder] = sort(channel_order);

data_new = {};
iChanCount = 1;
for iChannelbyOrder = curr_chanIndexOrder'
    
    curr_channel_label = opt.orgdata.label(iChannelbyOrder);
    if ~all(ismember(curr_channel_label,opt.orgdata.label))
        continue
    end
    cfg = [];
    cfg.channel = curr_channel_label;
    data_new{iChanCount} = ft_selectdata(cfg,opt.orgdata);
    iChanCount = iChanCount + 1;
end

if numel(data_new) > 1
    opt.orgdata = ft_appenddata([],data_new{:});
else
    opt.orgdata = data_new{:};
end

opt.hdr.label = opt.orgdata.label;
data_new = {};
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % SUBFUNCTION
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeHypnogramFile(filepath,hypn,delimiter)
hyp_export = fopen(filepath, 'wt');
if size(hypn,2) > 2
    for iRow = 1:size(hypn,1)
        fprintf(hyp_export, ['%i' delimiter '%i' [repmat([delimiter '%f'],1,size(hypn,2)-2)]], hypn(iRow,:));
        if iRow ~= size(hypn,1)
            fprintf(hyp_export,['\n']);
        end
    end
else
    for iRow = 1:size(hypn,1)
        fprintf(hyp_export, ['%i' delimiter '%i'], hypn(iRow,:));
        if iRow ~= size(hypn,1)
            fprintf(hyp_export,['\n']);
        end
    end
end
fclose(hyp_export);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % SUBFUNCTION
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeArtifactFile(filepath,opt,delimiter)
art_export = fopen(filepath, 'wt');
fprintf(art_export, ['%s' delimiter '%s' delimiter '%s\n'], 'artifact', 'begin_seconds', 'end_seconds');
temp_cfg = [];
for i=1:length(opt.artdata.label)
    temp_cfg.artfctdef.(opt.artdata.label{i}).artifact = convert_event(opt.artdata.trial{1}(i,:), 'artifact');
end
artTypes = fieldnames(temp_cfg.artfctdef);
for iArtType = 1:numel(artTypes)
    artType = artTypes{iArtType};
    
    artifactSamples = temp_cfg.artfctdef.(artType).artifact;
    for iArt = 1:size(artifactSamples,1)
        beginsample = artifactSamples(iArt,1);
        endsample = artifactSamples(iArt,2);
        fprintf(art_export, ['%s' delimiter '%f' delimiter '%f\n'], artType, beginsample/opt.fsample, endsample/opt.fsample);
    end
end
fclose(art_export);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % SUBFUNCTION
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cfg] = removeAdditionalHypnogram(hypn_mult_idx,cfg)
if ~isempty(cfg.hypn_mult) && (hypn_mult_idx <= numel(cfg.hypn_mult))
    idx = 1:numel(cfg.hypn_mult);
    cfg.hypn_mult = cfg.hypn_mult(idx ~= hypn_mult_idx);
    cfg.hypn_plot_interpol_mult = cfg.hypn_plot_interpol_mult(idx ~= hypn_mult_idx);
    cfg.hypn_plot_interpol_MA_mult = cfg.hypn_plot_interpol_MA_mult(idx ~= hypn_mult_idx);
    cfg.hypn_mult_idx = numel(cfg.hypn_mult)+1;
end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % SUBFUNCTION
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opt, cfg] = readArtifactFile(filepath,opt,cfg,delimiter)

af = dataset('File',[filepath],'Delimiter',delimiter);
if ~isempty(af)
    %unique(af.artifact)
    af.begin_seconds = round(af.begin_seconds*opt.fsample);
    af.end_seconds = round(af.end_seconds*opt.fsample);
    
    
    
    
    
    % collect the artifacts that have been detected from cfg.artfctdef.xxx.artifact
    artlabel = unique(af.artifact);
    sel      = zeros(size(artlabel));
    artifact = cell(size(artlabel));
    
    for i=1:length(artlabel)
        temp_idx = strcmp(artlabel(i),af.artifact);
        artifact{i} = [af.begin_seconds(temp_idx) af.end_seconds(temp_idx)];
        fprintf('detected %3d %s artifacts\n', size(artifact{i}, 1), artlabel{i});
    end
    
    cfg.selectfeature = artlabel(1);
    opt.ftsel       = find(strcmp(artlabel,cfg.selectfeature)); % current artifact/feature being selected
    
    if length(artlabel) > 9
        error(['only up to 9 artifacts groups supported, but ' num2str(length(artlabel)) ' found'])
    end
    
    % make artdata representing all artifacts in a "raw data" format
    datendsample = max(opt.trlorg(:,2));
    
    artdata = [];
    artdata.trial{1}       = convert_event(artifact, 'boolvec', 'endsample', datendsample); % every artifact is a "channel"
    artdata.time{1}        = offset2time(0, opt.fsample, datendsample);
    artdata.label          = artlabel;
    artdata.fsample        = opt.fsample;
    artdata.cfg.trl        = [1 datendsample 0];
    
    
    % % legend artifacts/features
    % for iArt = 1:length(artlabel)
    %     %   uicontrol('tag', 'artifactui', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', artlabel{iArt}, 'userdata', num2str(iArt), 'position', [0.91, 0.9 - ((iArt-1)*0.09), 0.08, 0.04], 'backgroundcolor', opt.artcolors(iArt,:))
    %     %   uicontrol('tag', 'artifactui', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', ['shift+' num2str(iArt)], 'position', [0.91, 0.855 - ((iArt-1)*0.09), 0.03, 0.04], 'backgroundcolor', opt.artcolors(iArt,:))
    %     %   uicontrol('tag', 'artifactui', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', ['control+' num2str(iArt)], 'position', [0.96, 0.855 - ((iArt-1)*0.09), 0.03, 0.04], 'backgroundcolor', opt.artcolors(iArt,:))
    %     uicontrol('tag', 'artifactui_button', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', ['artifact(' opt.artdata.label{opt.ftsel} ')'], 'userdata', 'a', 'position', [0.01, temp_lower_line_y2 - ((iArt-1)*0.09), 0.04, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1])
    %     uicontrol('tag', 'artifactui', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', ['shift+' num2str(iArt)], 'position', [0.01, temp_lower_line_y - ((iArt-1)*0.09), 0.02, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1])
    %     uicontrol('tag', 'artifactui', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', ['control+' num2str(iArt)], 'position', [0.03, temp_lower_line_y - ((iArt-1)*0.09), 0.02, 0.04],'backgroundcolor',[0 0 0],'foregroundcolor',[1 1 1])
    %
    % end
    
    
    
    opt.artdata = artdata;
end

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % SUBFUNCTION
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h, opt, cfg, hyp_file_filterindex] = loadSession(h,opt,cfg)
temppath = [cfg.outputfilespath];

[hyp_file_name hyp_file_path hyp_file_filterindex] = uigetfile(...
    {'*.mat','MAT-files (*.mat)';...
    '*.*',  'All Files (*.*)'},...
    'Open session',...
    temppath);

if hyp_file_filterindex ~= 0
    if ~isfield(cfg,'hhyp')
        cfg.hhyp = figure;
    end
    if  ~ishandle(cfg.hhyp)
        cfg.hhyp = figure;
    end
    figure(cfg.hhyp)
    set(cfg.hhyp, 'WindowButtonDownFcn',   {@select_sleep_stage_cb, 'h_main',h});
    set(cfg.hhyp, 'NumberTitle', 'off');
    cfg.hhypfigax = gca;
    
    hhyp = cfg.hhyp;
    hhypfigax = cfg.hhypfigax;
    
    figure(h)
    
    tempfilepath = [hyp_file_path hyp_file_name];
    tempcfg = load(tempfilepath,'cfg');
    
    if isfield(tempcfg.cfg,'browserversion')
        temp_saved_version = tempcfg.cfg.browserversion;
    else
        temp_saved_version = 'lower than 2.3.8';
        tempcfg.cfg.browserversion = 'lower than 2.3.8';
    end
    if strcmp(temp_saved_version,cfg.browserversion)
    else
        answer_open = questdlg(['Opening this session might fail!, This version ' cfg.browserversion ' does not match the one of the saved session ' temp_saved_version], ...
            'Open might fail', ...
            'Still try', ...
            'Cancel','Cancel');
        switch answer_open
            case 'Still try'
                load(tempfilepath,'opt')
                cfg = tempcfg.cfg;
                tempcfg = [];
                
                cfg.hhyp = hhyp;
                cfg.hhypfigax = hhypfigax;
            otherwise
                tempcfg = [];
                return
        end
        
    end
    
    
else
    error('');
end

if ~isfield(cfg,'color_text_on_bg')
    cfg.color_text_on_bg = [0.8 0.8 0.8];
end

end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % SUBFUNCTION
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hypn_plot_interpol hypn_plot_interpol_MA] = interpolate_hypn_for_plot(hypn,epochLengthSamples,plot_MA_offset,plot_yaxequidist)

if plot_yaxequidist
    remY = -1;
else
    remY  = -0.5;
end

        %plot_MA_offset = -5.5;
        hypn_plot = hypn;
        hypn_plot(hypn_plot(:,1) == 5,1) = 0.5;
        hypn_plot(hypn_plot(:,1) == 8,1) = 0;
        hypn_plot(:,1) = hypn_plot(:,1)*-1;
        hypn_plot_MA = hypn_plot(:,2) ;
        hypn_plot_MA = hypn_plot_MA*0.5;
        hypn_plot_MA(hypn_plot_MA > 1) = 1.35;
        hypn_plot = hypn_plot(:,1) ;
        hypn_plot_interpol = [];
        hypn_plot_interpol_MA = [];
        for iEp = 1:length(hypn_plot)
            temp_samples = repmat(hypn_plot(iEp),epochLengthSamples,1);
            if (hypn_plot(iEp) == remY) %REM
                if plot_yaxequidist
                    temp_samples(1:2:end) = -0.5;
                    temp_samples(2:2:end) = -1.5;
                else
                    temp_samples(1:2:end) = -0.3;
                    temp_samples(2:2:end) = -0.7;
                end
%                 for iSamp = 1:length(temp_samples)
%                     if mod(iSamp,2) == 0
%                         temp_samples(iSamp) = -0.20;
%                     else
%                         temp_samples(iSamp) = -0.70;
%                     end
%                 end
            end
            
            hypn_plot_interpol = [hypn_plot_interpol; temp_samples];
            
            temp_samples_MA = repmat(plot_MA_offset+hypn_plot_MA(iEp),epochLengthSamples,1);
            if (hypn_plot_MA(iEp) > 0) %REM
                for iSamp = 1:length(temp_samples_MA)
                    if mod(iSamp,2) == 1
                        temp_samples_MA(iSamp) = plot_MA_offset;
                    end
                end
            end
            hypn_plot_interpol_MA = [hypn_plot_interpol_MA; temp_samples_MA];
            
        end

end

