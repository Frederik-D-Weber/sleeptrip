function [fh] = st_hypnoplot(cfg, scoring)

% ST_HYPNOPLOT plots a hypnogram from the scoring
%
% Use as
%   [fh] = st_hypnoplot(cfg,scoring)
%
%   scoring is a structure provided by ST_READ_SCORING
%   it returns the figure handle
%
%   config file can be empty, e.g. cfg = []
%
% Optional configuration parameters are
%   cfg.plottype               = string, the type of plot 'classic' plots the line graph as typical or 'colorblocks' plots the colorbocks  or 'colorbar' for only one bar of colors (default = 'classic')
%   cfg.colorscheme            = srting, indicating the color schemes:
%                                       'bright' or 'dark' or 'restless'
%                                       (default = 'dark')
%   cfg.colorblocksconnect     = string, either 'yes' or 'no' if lines between colorblocks should be shown
%                                only has effect for cfg.plottype = 'colorblocks'(default = 'no')
%   cfg.plotlegend             = string, if the legend should be plotted either 'yes' or 'no' (default = 'yes')
%   cfg.plotsleeponset         = string, plot an indicator of sleep onset either 'yes' or 'no' (default = 'yes')
%   cfg.plotsleepoffset        = string, plot an indicator of sleep offset either 'yes' or 'no' (default = 'yes')
%   cfg.plotunknown            = string, plot unscored/unkown epochs or not either 'yes' or 'no' (default = 'yes')
%   cfg.plotexcluded           = string, plot excluded epochs 'yes' or 'no' (default = 'yes')
%   cfg.yaxdisteqi             = string, plot the y-axis ticks in equal distanve from each other 'yes' or 'no' (default = 'no')
%   cfg.sleeponsetdef          = string, sleep onset either 'N1' or 'N1_NR' or 'N1_XR' or
%                                'NR' or 'N2R' or 'XR' or 'AASM' or 'X2R' or 
%                                'N2' or 'N3' or 'SWS' or 'S4' or 'R',
%                                see ST_SLEEPONSET for details (default = 'N1_XR')
%   cfg.title                  = string, title of the figure to export the figure
%   cfg.timeticksdiff          = scalar, time difference in minutes the ticks are places from each other (default = 30);
%   cfg.timemin                = scalar, minimal time in minutes the ticks 
%                                have, e.g. 480 min, will plot tick at least to 480 min (default = 0);
%   cfg.timerange              = vector, [mintime maxtime] of the time axis
%                                limits in minutes, overwrites all the
%                                other contraints
%                                have, e.g. 480 min, will plot tick at least to 480 min (default = display all);
%   cfg.timeunitdisplay        = string, time unit for display in the x-axis labels
%                                either 'minutes' or 'seconds' or 'hours' or 'days'
%                                Note: this will not affect the other parameters like cfg.timerange to be given in minutes
%                                (default = 'minutes')
%   cfg.considerdataoffset     = string, 'yes' or 'no' if dataoffset is represented in time axis (default = 'yes');
%
%  Events can be plotted using the following options
%
%   cfg.eventtimes             = a Nx1 cell containing 1x? vectors of event
%                                  time points (in seconds) representing N
%                                  event types of ? instances to be
%                                  plotted. e.g.
%                                 {[1.5, 233.2, 455.6]; ...
%                                  [98, 3545.9]; ...
%                                  [393.4, 425.8, 900.0, 4001.01]}
%   cfg.eventdurations          = optional, but if set events are given a duration with 
%                                 cfg.eventtimes being the start of the
%                                 event and the respctive duration in
%                                 seconds added to that. The event
%                                 durations then needs to match the number
%                                 of event times for each event type or be
%                                 empty, i.e. having no duration for all
%                                 events of that event type. e.g.
%                                {[30, 120, 600]; ...
%                                             []; ...
%                                    [0, 0, 0, 1]}
%   cfg.eventlabels            = Nx1 cellstr with the labels to the events corresponding to the rows in cfg.eventstimes
%   cfg.eventvalues            = a Nx1 cell containing 1x? vectors of event
%                                values (e.g. amplitude)
%                                 {[20.3, 23.2, 45.6]; ...
%                                  [18, 35.9]; ...
%                                  [39.1, 42.5, 80.0, 42.1]}
%   cfg.eventvalueranges       = a Nx1 cell containing 1x2 vectors of event
%                                values ranges (e.g. min and max of amplitude)
%                                 {[20 40]; ...
%                                  [18, 36]; ...
%                                  [39, 80.0]}
%   cfg.eventvaluerangesrnddec       = round event ranges to that amount of decimal (default = 2)
%   cfg.eventcolors           = a Nx3 color matrix with 3 RGB values from 0 to 1 color for each of
%                               the N event types, (default = lines(N))
%
%
% If you wish to export the figure then define also the following
%   cfg.figureoutputfile       = string, file to export the figure
%   cfg.figureoutputformat     = string, either 'png' or 'epsc' or 'svg' or 'tiff' or
%                                'pdf' or 'bmp' or 'fig' (default = 'png')
%   cfg.figureoutputunit       = string, dimension unit (1 in = 2.54 cm) of hypnograms.
%                                either 'points' or 'normalized' or 'inches'
%                                or 'centimeters' or 'pixels' (default =
%                                'inches')
%   cfg.figureoutputwidth      = scalar, choose format dimensions in inches
%                                (1 in = 2.54 cm) of hypnograms. (default = 9)
%   cfg.figureoutputheight     = scalar, format dimensions in inches (1 in = 2.54 cm) of hypnograms. (default = 3)
%   cfg.figureoutputresolution = scalar, choose resolution in pixesl per inches (1 in = 2.54 cm) of hypnograms. (default = 300)
%   cfg.figureoutputfontsize   = scalar, Font size in units stated in
%                                parameter cfg.figureoutputunit (default = 0.1)
%   cfg.timestamp              = either 'yes' or 'no' if a time stamp should be
%                                added to filename (default = 'yes')
%   cfg.folderstructure        = either 'yes' or 'no' if a folder structure should
%                                be created with the result origin and type 
%                                all results will be stored in "/res/..." (default = 'yes')
%
%
% See also ST_READ_SCORING

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
%dt = now;

ttic = tic;
mtic = memtic;
functionname = getfunctionname();
fprintf([functionname ' function started\n']);


% set the defaults
cfg.plottype                = ft_getopt(cfg, 'plottype', 'classic');
cfg.plotlegend              = ft_getopt(cfg, 'plotlegend', 'yes');
cfg.title                   = ft_getopt(cfg, 'title', '');
cfg.timeticksdiff           = ft_getopt(cfg, 'timeticksdiff', 30);
cfg.timemin                 = ft_getopt(cfg, 'timemin', 0);
cfg.considerdataoffset      = ft_getopt(cfg, 'considerdataoffset', 'yes');
cfg.plotsleeponset          = ft_getopt(cfg, 'plotsleeponset', 'yes');
cfg.plotsleepoffset         = ft_getopt(cfg, 'plotsleepoffset', 'yes');
cfg.plotlightsoff           = ft_getopt(cfg, 'plotlightsoff', 'yes');
cfg.plotlightson            = ft_getopt(cfg, 'plotlightson', 'yes');
cfg.plotunknown             = ft_getopt(cfg, 'plotunknown', 'yes');
cfg.plotexcluded            = ft_getopt(cfg, 'plotexcluded', 'yes');
cfg.sleeponsetdef           = ft_getopt(cfg, 'sleeponsetdef', 'N1_XR');
cfg.colorscheme             = ft_getopt(cfg, 'colorscheme', 'dark');
cfg.colorblocksconnect      = ft_getopt(cfg, 'colorblocksconnect', 'no');

cfg.eventvaluerangesrnddec  = ft_getopt(cfg, 'eventrangernddec', 2);
cfg.timerange               = ft_getopt(cfg, 'timerange', [], true);

cfg.timeunitdisplay         = ft_getopt(cfg, 'timeunitdisplay', 'minutes');

cfg.figureoutputformat      = ft_getopt(cfg, 'figureoutputformat', 'png');
cfg.figureoutputunit        = ft_getopt(cfg, 'figureoutputunit', 'inches');
cfg.figureoutputwidth       = ft_getopt(cfg, 'figureoutputwidth', 9);
cfg.figureoutputheight      = ft_getopt(cfg, 'figureoutputheight', 3);
cfg.figureoutputresolution  = ft_getopt(cfg, 'figureoutputresolution', 300);
cfg.figureoutputfontsize    = ft_getopt(cfg, 'figureoutputfontsize', 0.1);
cfg.timestamp               = ft_getopt(cfg, 'timestamp', 'yes');
cfg.folderstructure         = ft_getopt(cfg, 'folderstructure', 'yes');


if istrue(cfg.colorblocksconnect) && istrue(cfg.plotlegend) && strcmp(cfg.plottype,'colorblocks')
    ft_warning('cfg.colorblocksconnect = ''yes'' currently does not support the legend with cfg.plottype = ''colorblocks'' and thus legend is DISABLED!')
    cfg.plotlegend = 'no';
end
if isfield(cfg, 'eventranges'), ft_error('Changed naming convention, use cfg.eventvalueranges instead of cfg.eventvalueranges'); end


if strcmp(cfg.plottype,'colorbar') || strcmp(cfg.plottype,'colorblocks')
    if isfield(cfg,'yaxdisteqi')
        if ~istrue(cfg.yaxdisteqi)
            ft_warning('cfg.yaxdisteqi is set to ''yes'' because of the cfg.plottype = %s',cfg.plottype)
            cfg.yaxdisteqi = 'yes';
        end
    else
         cfg.yaxdisteqi = 'yes';
    end
else
    cfg.yaxdisteqi = ft_getopt(cfg, 'yaxdisteqi', 'no');
end


if (isfield(cfg, 'eventtimes') && ~isfield(cfg, 'eventlabels')) || (~isfield(cfg, 'eventtimes') && isfield(cfg, 'eventlabels'))  
    ft_error('both cfg.eventtimes and cfg.eventlabels have to be defined togehter.');
end

if isfield(cfg, 'eventtimes')
    if size(cfg.eventtimes,1) ~=  numel(cfg.eventlabels)
        ft_error('dimensions of cfg.eventtimes and cfg.eventlabels do not match.');
    end
end

if (isfield(cfg, 'eventvalues') && ~isfield(cfg, 'eventtimes')) 
    ft_error('both cfg.eventvalues needs a cfg.eventtimes to be defined.');
end

if (isfield(cfg, 'eventvalues') && ~isfield(cfg, 'eventvalueranges')) || (~isfield(cfg, 'eventvalues') && isfield(cfg, 'eventvalueranges'))  
    ft_error('both cfg.eventvalues and cfg.eventvalueranges have to be defined togehter.');
end

if isfield(cfg, 'eventvalues')
    if size(cfg.eventtimes,1) ~=  numel(cfg.eventvalues)
        ft_error('dimensions of cfg.eventtimes and cfg.eventvalues do not match.');
    end
end

if isfield(cfg, 'eventdurations')
    if size(cfg.eventtimes,1) ~=  numel(cfg.eventdurations)
        ft_error('dimensions of cfg.eventtimes and cfg.eventdurations do not match.');
    end
    
    for iEvent = 1:numel(cfg.eventtimes)
        err = false;
        size_evt = size(cfg.eventtimes{iEvent});
        size_evd = size(cfg.eventdurations{iEvent});
        if ~isempty(cfg.eventdurations{iEvent})
            if ~all(size_evt == size_evd)
                err = true;
                ft_warning('some event types in cfg.eventtimes do not match with the ones in cfg.eventdurations, for the %d event type and times (%d %d) not matching dimension of duratoins (%d %d).',iEvent,size_evt(1),size_evt(2),size_evd(1),size_evd(2))
            end
            if err
                ft_error('some event types in cfg.eventtimes do not match with the ones in cfg.eventdurations, see prior warning to find out which.');
            end
        end
    end
end



if isfield(cfg, 'eventvalueranges')
    if size(cfg.eventtimes,1) ~=  numel(cfg.eventvalueranges)
        ft_error('dimensions of cfg.eventtimes and cfg.eventvalueranges do not match.');
    end
end

if isfield(cfg, 'eventtimes')
    nEvents = numel(cfg.eventtimes);
    if isfield(cfg, 'eventcolors')
        nEventColors = size(cfg.eventcolors,1);
        if(nEventColors ~= nEvents)
        	ft_error('number of rows in cfg.eventcolors %d does not match with number of event types %d.',nEvents,nEventColors);
        end
    else %set default colors
        cfg.eventcolors = lines(nEvents);
    end
end

if strcmp(cfg.considerdataoffset, 'yes')
    offsetseconds = scoring.dataoffset;
else
    offsetseconds = 0;
end

saveFigure   = false;
if isfield(cfg, 'figureoutputfile')
    saveFigure = true;
end

hasLightsOff = false;
if isfield(scoring, 'lightsoff')
    hasLightsOff = true;
end

hasLightsOn = false;
if isfield(scoring, 'lightson')
    hasLightsOn = true;
end

fprintf([functionname ' function initialized\n']);

dummySampleRate = 100;
epochLengthSamples = scoring.epochlength * dummySampleRate;
nEpochs = numel(scoring.epochs);

% if hasLightsOff
%     lightsOffSample = scoring.lightsoff*dummySampleRate;
% else
%     lightsOffSample = 0;
% end

%convert the sleep stages to hypnogram numbers
hypn = [cellfun(@(st) sleepStage2hypnNum(st,~istrue(cfg.plotunknown),istrue(cfg.yaxdisteqi)),scoring.epochs','UniformOutput',1) ...
    scoring.excluded'];


hypnStages = [cellfun(@sleepStage2str,scoring.epochs','UniformOutput',0) ...
    cellfun(@sleepStage2str_alt,scoring.epochs','UniformOutput',0) ...
    cellfun(@sleepStage2str_alt2,scoring.epochs','UniformOutput',0)...
    cellfun(@sleepStage2str_alt3,scoring.epochs','UniformOutput',0)];


hypnEpochs = 1:numel(scoring.epochs);
hypnEpochsBeginsSamples = (((hypnEpochs - 1) * epochLengthSamples) + 1)';

%onsetCandidateIndex = getSleepOnsetEpoch(hypnStages,hypnEpochsBeginsSamples,lightsOffSample,cfg.sleeponsetdef);

[onsetCandidateIndex lastsleepstagenumber onsetepoch] = st_sleeponset(cfg,scoring);

if isempty(lastsleepstagenumber)
    lastsleepstagenumber = nEpochs;
end

hasA = false;
if any(ismember('A',scoring.epochs))
    hasA = true;
end

%%% plot hypnogram figure

switch scoring.standard
    case 'aasm'
        if hasA
            if istrue(cfg.yaxdisteqi)
                plot_exclude_offset = -6;
                yTick      = [2  1        0     -1  -2   -3   -4 ];
            else
                
                plot_exclude_offset = -5;
                yTick      = [1.5  1        0     -0.5  -1   -2   -3 ];
            end
            
            yTickLabel = {'?' 'A'      'W'    'R'  'N1' 'N2' 'N3'};
        else
            if istrue(cfg.yaxdisteqi)
                plot_exclude_offset = -6;
                yTick      = [1        0     -1  -2   -3   -4 ];
            else
                
                plot_exclude_offset = -5;
                yTick      = [1        0     -0.5  -1   -2   -3 ];
            end
            
            yTickLabel = {'?'      'W'    'R'  'N1' 'N2' 'N3'};
        end
    case 'rk'
        if hasA
            if istrue(cfg.yaxdisteqi)
                plot_exclude_offset = -8;
                yTick      = [3  2   1  0     -1  -2   -3   -4   -5 ];
            else
                plot_exclude_offset = -7;
                yTick      = [1.5  1   0.5  0     -0.5  -1   -2   -3   -4 ];
            end
            
            yTickLabel = {'?' 'A' 'MT' 'W' 'R' 'S1' 'S2' 'S3' 'S4'};
        else
            if istrue(cfg.yaxdisteqi)
                plot_exclude_offset = -8;
                yTick      = [2   1  0     -1  -2   -3   -4   -5 ];
            else
                plot_exclude_offset = -7;
                yTick      = [ 1   0.5  0     -0.5  -1   -2   -3   -4 ];
            end
            
            yTickLabel = {'?' 'MT' 'W' 'R' 'S1' 'S2' 'S3' 'S4'};
        end
        
    otherwise
        ft_error('scoring standard ''%s'' not supported for ploting.\n Maybe use ST_SCORINGCONVERT to convert the scoring first.', scoring.standard);
end

switch cfg.plottype
    case 'colorbar'
        plot_exclude_offset = 1;
        yTick = [3];
        yTickLabel = {'Stage'};
end

if istrue(cfg.plotexcluded)
    yTickLabel{end+1} = 'Excl';
    yTick(end+1) = plot_exclude_offset;
end

if ~istrue(cfg.plotunknown)
    tempremind = strcmp(yTickLabel,'?');
    yTickLabel(tempremind) = [];
    yTick(tempremind) = [];
end


hhyp = figure;
axh = gca;
set(hhyp,'color',[1 1 1]);
set(axh,'FontUnits',cfg.figureoutputunit)
set(axh,'Fontsize',cfg.figureoutputfontsize);

switch cfg.plottype
    case 'classic'
        [hypn_plot_interpol hypn_plot_interpol_exclude] = interpolate_hypn_for_plot(hypn,epochLengthSamples,plot_exclude_offset,istrue(cfg.yaxdisteqi));
        x_time = (1:length(hypn_plot_interpol))/(dummySampleRate)  - 1/dummySampleRate;
        x_time = x_time + offsetseconds;
        x_time = x_time/60; % minutes
        x_time_hyp = x_time(1:length(hypn_plot_interpol));
        plot(axh,x_time_hyp,hypn_plot_interpol,'Color',[0 0 0])
        hold(axh,'on');
        
    case {'colorblocks', 'colorbar'}
        x_time = (0:numel(scoring.epochs)) * scoring.epochlength;
        x_time = x_time + offsetseconds;
        x_time = x_time/60; % minutes
        x_time_hyp = x_time;
        
        hp = [];
        
        labels = scoring.label;
        [lables_colors_topdown labels_ordered] = st_epoch_colors(labels,cfg.colorscheme);
        idxUsedLabels = [];
        
        incLabel = 1;
        
        [epoch_colors labels_ordered] = st_epoch_colors(scoring.epochs,cfg.colorscheme);
        
        offset_y = -0.5;%(iScoring-0.5);
        height = 1;
        y_hyp_pos_prev = [];
        for iEpoch = 1:numel(scoring.epochs)
            x1 = x_time(iEpoch);
            x2 = x_time(iEpoch+1);
            epoch = scoring.epochs(iEpoch);
            
            switch cfg.plottype
                case 'colorblocks'
                    y_hyp_pos = yTick(ismember(yTickLabel,epoch));
                case 'colorbar'
                    y_hyp_pos = yTick(1);
            end
            
            
            %if isfield(cfg,'plotunknown')
            if ~(~istrue(cfg.plotunknown) && strcmp(epoch,'?'))
                %h = ft_plot_patch([x1 x2 x2 x1], [offset_y offset_y offset_y+height offset_y+height], 'facecolor',epoch_colors(iEpoch,:));
                if istrue(cfg.colorblocksconnect) && ~isempty(y_hyp_pos_prev) && (y_hyp_pos_prev ~= y_hyp_pos)
                    htmp = plot(axh,[x1 x1],[y_hyp_pos+offset_y y_hyp_pos_prev+offset_y],'Color',[0.8 0.8 0.8]);
                    hold(axh,'on');
                end
                y_hyp_pos_prev = y_hyp_pos;
                h = patch([x1 x2 x2 x1], [y_hyp_pos+offset_y y_hyp_pos+offset_y y_hyp_pos+offset_y+height y_hyp_pos+offset_y+height],epoch_colors(iEpoch,:),'edgecolor','none');
                member = find(ismember(labels,epoch),1,'first');
                if ~ismember(member,idxUsedLabels)
                    hp(incLabel) = h;
                    incLabel = incLabel + 1;
                    idxUsedLabels = [idxUsedLabels member];
                end
            end
            %end
            
            if isfield(cfg,'plotexcluded')
                if istrue(cfg.plotexcluded) && scoring.excluded(iEpoch)
                    y_hyp_pos = yTick(end);
                    he = patch([x1 x2 x2 x1], [y_hyp_pos+offset_y y_hyp_pos+offset_y y_hyp_pos+offset_y+height y_hyp_pos+offset_y+height],[1 0 0],'edgecolor','none');
                end
            end
            

        end
        
        
        collabels = labels;
        for iLabel = 1:numel(labels)
            collabels{iLabel} = sprintf(['\\color[rgb]{%.4f,%.4f,%.4f}' labels{iLabel}],lables_colors_topdown(iLabel,1),lables_colors_topdown(iLabel,2),lables_colors_topdown(iLabel,3));
        end
        collabels = collabels(idxUsedLabels);
        [b, idx_ori_labels] = sort(idxUsedLabels);
        if istrue(cfg.plotlegend)
            hLegend = legend(hp(idx_ori_labels),collabels(idx_ori_labels),'Location','northoutside','Orientation','horizontal','Box','off');
        end
        hold(axh,'on');
    otherwise
       ft_error('cfg.plottype = %s is unknown, please see the help for available options.', cfg.plottype)
end


if strcmp(cfg.plotlightsoff, 'yes')
    if hasLightsOff
        lightsoff_time = (scoring.lightsoff/60);%in minutes
        switch cfg.plottype
            case 'classic'
                onset_y_coord_offset = 0.2;
                onset_y_coord = 0+onset_y_coord_offset;
                
            case 'colorblocks'
                onset_y_coord_offset = 0.5;
                onset_y_coord =  yTick(1)+onset_y_coord_offset;
                
            case 'colorbar'
                onset_y_coord_offset = 0.5;
                onset_y_coord =  yTick(1)+onset_y_coord_offset;
        end
        hold(axh,'on');
        scatter(axh,lightsoff_time,onset_y_coord,'filled','>','MarkerFaceColor',[0.1 0.1 0.1])
    end
end

if strcmp(cfg.plotlightson, 'yes')
    if hasLightsOn
        lightson_time = (scoring.lightson/60);%in minutes
        switch cfg.plottype
            case 'classic'
                onset_y_coord_offset = 0.2;
                onset_y_coord = 0+onset_y_coord_offset;
                
            case 'colorblocks'
                onset_y_coord_offset = 0.5;
                onset_y_coord =  yTick(1)+onset_y_coord_offset;
                
            case 'colorbar'
                onset_y_coord_offset = 0.5;
                onset_y_coord =  yTick(1)+onset_y_coord_offset;
        end
        hold(axh,'on');
        scatter(axh,lightson_time,onset_y_coord,'filled','>','MarkerFaceColor',[1 1 0])
    end
end

if strcmp(cfg.plotsleeponset, 'yes')
    if onsetCandidateIndex ~= -1
        onset_time = (onsetCandidateIndex-0.5)*(scoring.epochlength/60) + (offsetseconds/60);%in minutes
        switch cfg.plottype
            case 'classic'
                onset_y_coord_offset = 0.2;
                onset_y_coord = hypn_plot_interpol(find(x_time >=onset_time,1,'first'))+onset_y_coord_offset;
                
            case 'colorblocks'
                onset_y_coord_offset = 0.5;
                onset_y_coord =  yTick(ismember(yTickLabel,scoring.epochs{onsetCandidateIndex}))+onset_y_coord_offset;
                
            case 'colorbar'
                onset_y_coord_offset = 0.5;
                onset_y_coord =  yTick(1)+onset_y_coord_offset;
        end
        hold(axh,'on');
        scatter(axh,onset_time,onset_y_coord,'filled','v','MarkerFaceColor',[0 1 0]);
    end
end


offset_time = max(x_time);
if strcmp(cfg.plotsleepoffset, 'yes')
    if onsetCandidateIndex ~= -1
        offset_time = (lastsleepstagenumber+0.5)*(scoring.epochlength/60)+(offsetseconds/60);%in minutes
        switch cfg.plottype
            case 'classic'
                offset_y_coord_offset = 0.2;
                offset_y_coord = hypn_plot_interpol(find(x_time <=offset_time,1,'last'))+offset_y_coord_offset;
            case 'colorblocks'
                onset_y_coord_offset = 0.5;
                offset_y_coord =  yTick(ismember(yTickLabel,scoring.epochs{lastsleepstagenumber}))+onset_y_coord_offset;
            case 'colorbar'
                onset_y_coord_offset = 0.5;
                offset_y_coord =  yTick(1)+onset_y_coord_offset;
        end
        hold(axh,'on');
        scatter(axh,offset_time,offset_y_coord,'filled','^','MarkerFaceColor',[0 0 1]);
    end
end





eventTimeMaxSeconds = cfg.timemin*60;
offset_step = 0.5;
eventHeight = 0.4;
offset_event_y = max(yTick);

switch cfg.plottype
    case {'colorbar', 'colorblocks'}
        offset_event_y = offset_event_y - offset_y;
end


%find the maximal time of all events
max_temp_x_all = 0;
if isfield(cfg, 'eventtimes')
    for iEvent = 1:numel(cfg.eventtimes)
        if ~isempty(cfg.eventtimes{iEvent})
            if isfield(cfg, 'eventdurations')
                if ~isempty(cfg.eventdurations{iEvent})
                    max_temp_x_all = max(max_temp_x_all,max(max(cfg.eventtimes{iEvent}+cfg.eventdurations{iEvent})));
                else
                    max_temp_x_all = max(max_temp_x_all,max(max(cfg.eventtimes{iEvent})));
                end
            else
                max_temp_x_all = max(max_temp_x_all,max(max(cfg.eventtimes{iEvent})));
            end
        end
    end
end
max_temp_x_all = max_temp_x_all/60;

if isfield(cfg, 'eventtimes')
    
    nEvents = numel(cfg.eventtimes);
    tempcolors = cfg.eventcolors;
    for iEventTypes = 1:nEvents
        currEvents = cfg.eventtimes{iEventTypes};
        currEventsDurations = [];
        if isfield(cfg,'eventdurations')
            currEventsDurations = cfg.eventdurations{iEventTypes};
        end
        if ~isempty(currEvents)
            offset_event_y = offset_event_y + offset_step;
            currEventLabel = cfg.eventlabels{iEventTypes};
            
            yTick = [offset_event_y yTick];
            yTickLabel = {currEventLabel yTickLabel{:}};
            
            color = tempcolors(iEventTypes,:);
            eventTimeMaxSeconds = max([eventTimeMaxSeconds currEvents]);
            temp_x1 = (currEvents/60)';
            if ~isempty(currEventsDurations)
                temp_x2 = ((currEvents+currEventsDurations)/60)';
            end
            temp_y = repmat(offset_event_y,numel(currEvents),1);
            if isfield(cfg, 'eventvalues')
                currEventValues = cfg.eventvalues{iEventTypes};
                currEventValueRanges = cfg.eventvalueranges{iEventTypes};
                currEventValueRanges = round(currEventValueRanges,cfg.eventvaluerangesrnddec);
                event_scale = fw_normalize(currEventValues, min(currEventValueRanges),  max(currEventValueRanges), 0.1, 1)';
                text(max_temp_x_all+1,temp_y(1),['[' num2str(min(currEventValueRanges)) ' ' num2str(max(currEventValueRanges)) ']']);
            else
                event_scale = 1;
            end
            if ~isempty(currEventsDurations)
                temp_plot_y = [temp_y-(eventHeight*event_scale)/2 temp_y+(eventHeight*event_scale)/2];
                %plot(axh,[temp_x1 temp_x2]',temp_plot_y,'Color',color)
                
%                 for iDurEv = 1:numel(temp_x1)
%                     hev = patch([temp_x1(iDurEv) temp_x2(iDurEv) temp_x2(iDurEv) temp_x1(iDurEv)], [temp_plot_y(1,iDurEv) temp_plot_y(1,iDurEv) temp_plot_y(2,iDurEv) temp_plot_y(2,iDurEv)],color,'edgecolor','none');
%                 end
                 hev = patch([temp_x1 temp_x2 temp_x2 temp_x1]', [temp_plot_y(:,1) temp_plot_y(:,1) temp_plot_y(:,2) temp_plot_y(:,2)]',color,'edgecolor','none');
            else
                temp_plot_y = [temp_y-(eventHeight*event_scale)/2 temp_y+(eventHeight*event_scale)/2]';
                plot(axh,[temp_x1 temp_x1]',temp_plot_y,'Color',color)
            end
        end
    end
end


switch cfg.plottype
    case 'classic'
        temp_max_y = max(yTick);
        
        if istrue(cfg.plotexcluded)
            temp_min_y = plot_exclude_offset;
        else
            temp_min_y = min(yTick) - 1;
        end
    case {'colorblocks', 'colorbar'}
        temp_max_y = max(yTick)+0.5;
        temp_min_y = min(yTick)-0.5;
        
end




if isfield(cfg, 'eventtimes')
    temp_max_y = temp_max_y + eventHeight;
end


if isfield(cfg,'plotexcluded')
    if istrue(cfg.plotexcluded)
        if strcmp(cfg.plottype,'classic')
            plot(axh,x_time_hyp,hypn_plot_interpol_exclude,'Color',[1 0 0])
        end
    end
end

if ~isempty(cfg.timerange)
    xlim(axh,[min(cfg.timerange) max(cfg.timerange)]);
else
    xlim(axh,[0 (max([max(x_time), cfg.timemin, eventTimeMaxSeconds/60, offset_time]))]);
end

ylabel(axh,'Stages');
ylim(axh,[temp_min_y temp_max_y])

set(axh, 'yTick', flip(yTick));
set(axh, 'yTickLabel', flip(yTickLabel));
set(axh,'TickDir','out');
xTick = [0:cfg.timeticksdiff:(max([max(x_time),cfg.timemin,eventTimeMaxSeconds/60]))];
set(axh, 'xTick', xTick);
timeunit = 'min';
switch cfg.timeunitdisplay
    case {'m' 'min' 'minute' 'minutes'}
        timeunit = 'min';
    case {'s' 'sec' 'seconds'}
        set(axh, 'xTickLabel', arrayfun(@num2str,round(xTick*60),'UniformOutput',false)); 
        timeunit = 's';
    case {'h' 'hour' 'hours'}
        set(axh, 'xTickLabel', arrayfun(@num2str,round(xTick/60,2),'UniformOutput',false)); 
        timeunit = 'h';
    case {'d' 'day' 'days'}
        set(axh, 'xTickLabel', arrayfun(@num2str,round(xTick/(60*24),3),'UniformOutput',false)); 
        timeunit = 'd';
end
    
set(axh, 'box', 'off')

%     begsample = 0;
%     endsample = 0;
%     x_pos_begin = x_time(begsample);
%     x_pos_end = x_time(endsample);
%     x_pos = [x_pos_begin x_pos_end x_pos_end x_pos_begin];
%     y_pos = [plot_exclude_offset plot_exclude_offset 1 1];
%     pos_now = patch(x_pos,y_pos,[0.5 0.25 1],'parent',axh);
%     set(pos_now,'FaceAlpha',0.4);
%     set(pos_now,'EdgeColor','none');

%     line([x_pos_begin x_pos_begin],[plot_exclude_offset temp_max_y],'color',[0.25 0.125 1],'parent',axh);

%titleName = sprintf('Hypnogram_datasetnum_%d_file_%d',iData,iHyp);
xlabel(['Time [' timeunit ']']);
ylabel('Sleep stage');


cfg = st_adjustfigure(cfg,hhyp);

hold(axh,'off')

if saveFigure
    cfg.functionname = functionname;
    cfg.subfolder = 'hypnograms';
    cfg = st_savefigure(cfg,hhyp);
end
fh = hhyp;

%%% plot hypnogram figure end

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % SUBFUNCTION
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hypn_plot_interpol hypn_plot_interpol_exclude] = interpolate_hypn_for_plot(hypn,epochLengthSamples,plot_exclude_offset, plot_yaxequidist)

if plot_yaxequidist
    remY = -1;
else
    remY  = -0.5;
end


hypn_plot = hypn;
hypn_plot_exclude = hypn_plot(:,2) ;
%hypn_plot_exclude = hypn_plot_exclude*0.5;
%hypn_plot_exclude(hypn_plot_exclude > 1) = 1.35;
hypn_plot = hypn_plot(:,1) ;
hypn_plot_interpol = [];
hypn_plot_interpol_exclude = [];
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
    
    temp_samples_exclude = repmat(plot_exclude_offset+hypn_plot_exclude(iEp),epochLengthSamples,1);
    if (hypn_plot_exclude(iEp) > 0) %excluded
        for iSamp = 1:length(temp_samples_exclude)
            if mod(iSamp,2) == 1
                temp_samples_exclude(iSamp) = plot_exclude_offset;
            end
        end
    end
    hypn_plot_interpol_exclude = [hypn_plot_interpol_exclude; temp_samples_exclude];
end

end
