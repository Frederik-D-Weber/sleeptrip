function [cfg, fh, ah, ch, res] = st_topoplotres(cfg, res)

% ST_TOPOPLOTRES plots the topographic distribution over the head
% of a 2-dimensional property (column of a table in the) result structure as
% retruned from SleepTrip functions (e.g. ST_SPINDLES, ST_POWER) as if they were
% plotted with FT_TOPOPLOTER
%
% Use as
%   [cfg] = st_topoplotRes(cfg, res)
%   [cfg, fh, ah, ch, res] = st_topoplotres(cfg, res)
%
% Required configuration parameters are:
%   cfg.property      = string, with the columname/property in the result
%                       table which to use for the values
%   cfg.layout        = specify the channel layout for plotting using one of
%                       the supported ways (see below).
%
% Optional configuration parameters are:
%   cfg.channel       = the channels you want to select (default = 'all');
%   cfg.renderer      = string, the Renderer of the figure,
%                        either 'painters' use for vector graphic export (recommended)
%                         'zbuffer'  fast and accurate
%                         'OpenGL'   when you have a good opengl function
%   cfg.filtercolumns = either a string or a cellstr of the columns you want to filter for, e.g.
%                        {'resnum', 'freq'}
%   cfg.filtervalues  = either a string/value or cell of cells/cellstr, with the corresponding values to
%                       the columns defined cfg.filtercolumns, e.g. 
%                        {{1}, {4:5}}
%   cfg.average       = if (after filtering) the property values should be
%                       avearaged for each channel (NaNs are ignored)
%                       either 'yes' or 'no', NOT supported with parallel use of cfg.maskproperty (default = 'no')
%   cfg.maskproperty  = column in the data to be used for masking of values. It should have values between 0 and 1, where 0 corresponds to transparent.
%
% Additional parameters from FT_TOPOPLOTER can also be included that are passed to that function internally:
%   cfg.xlim               = limit for 1st dimension in data (e.g., time), can be 'maxmin' or [xmin xmax] (default = 'maxmin')
%   cfg.zlim               = limits for color dimension, 'maxmin', 'maxabs', 'zeromax', 'minzero', or [zmin zmax] (default = 'maxmin')
%   cfg.colormap           = any sized colormap, see COLORMAP
%   cfg.marker             = 'on', 'labels', 'numbers', 'off'
%   cfg.markersymbol       = channel marker symbol (default = 'o')
%   cfg.markercolor        = channel marker color (default = [0 0 0] (black))
%   cfg.markersize         = channel marker size (default = 2)
%   cfg.markerfontsize     = font size of channel labels (default = 8 pt)
%   cfg.highlight          = 'off', 'on', 'labels', 'numbers'
%   cfg.highlightchannel   =  Nx1 cell-array with selection of channels, or vector containing channel indices see FT_CHANNELSELECTION
%   cfg.highlightsymbol    = highlight marker symbol (default = 'o')
%   cfg.highlightcolor     = highlight marker color (default = [0 0 0] (black))
%   cfg.highlightsize      = highlight marker size (default = 6)
%   cfg.highlightfontsize  = highlight marker size (default = 8)
%   cfg.hotkeys            = enables hotkeys (pageup/pagedown/m) for dynamic zoom and translation (ctrl+) of the color limits
%   cfg.colorbar           = 'yes'
%                            'no' (default)
%                            'North'              inside plot box near top
%                            'South'              inside bottom
%                            'East'               inside right
%                            'West'               inside left
%                            'NorthOutside'       outside plot box near top
%                            'SouthOutside'       outside bottom
%                            'EastOutside'        outside right
%                            'WestOutside'        outside left
%   cfg.colorbartext       = string indicating the text next to colorbar
%   cfg.interplimits       = limits for interpolation (default = 'head')
%                            'electrodes' to furthest electrode
%                            'head' to edge of head
%   cfg.interpolation      = 'linear','cubic','nearest','v4' (default = 'v4') see GRIDDATA
%   cfg.style              = plot style (default = 'both')
%                            'straight' colormap only
%                            'contour' contour lines only
%                            'both' (default) both colormap and contour lines
%                            'fill' constant color between lines
%                            'blank' only the head shape
%   cfg.gridscale          = scaling grid size (default = 67)
%                            determines resolution of figure
%   cfg.shading            = 'flat' or 'interp' (default = 'flat')
%   cfg.comment            = 'no', 'auto' or 'xlim' (default = 'auto')
%                            'auto': date, xparam and zparam limits are printed
%                            'xlim': only xparam limits are printed
%   cfg.commentpos         = string or two numbers, position of the comment (default = 'leftbottom')
%                            'lefttop' 'leftbottom' 'middletop' 'middlebottom' 'righttop' 'rightbottom'
%                            'title' to place comment as title
%                            'layout' to place comment as specified for COMNT in layout
%                            [x y] coordinates
%   cfg.interactive        = Interactive plot 'yes' or 'no' (default = 'yes')
%                            In a interactive plot you can select areas and produce a new
%                            interactive plot when a selected area is clicked. Multiple areas
%                            can be selected by holding down the SHIFT key.
%   cfg.interpolatenan     = string 'yes', 'no' (default = 'yes')
%                            interpolate over channels containing NaNs
%
%
% The layout defines how the channels are arranged. You can specify the
% layout in a variety of ways:
%  - you can provide a pre-computed layout structure, see FT_PREPARE_LAYOUT
%  - you can give the name of an ascii layout file with extension *.lay
%  - you can give the name of an electrode file
%  - you can give an electrode definition, i.e. "elec" structure
%  - you can give a gradiometer definition, i.e. "grad" structure
% If you do not specify any of these and the data structure contains an
% electrode or gradiometer structure, that will be used for creating a
% layout. If you want to have more fine-grained control over the layout
% of the subplots, you should create your own layout file.
%
% See also FT_TOPOPLOTER, FT_PREPARE_LAYOUT, ST_POWER, ST_SPINDLES

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

prev_renderer = get(0, 'DefaultFigureRenderer');

% set defaults
cfg.channel  = ft_getopt(cfg, 'channel', 'all', 1);
cfg.renderer  = ft_getopt(cfg, 'renderer', prev_renderer);
cfg.average  = ft_getopt(cfg, 'average', 'no');

set(0, 'DefaultFigureRenderer', cfg.renderer);

fprintf([functionname ' function initialized\n']);

res = st_res_filter(cfg, res);

if ~any(strcmp(res.table.Properties.VariableNames,cfg.property))
    ft_error(['result structure table needs to have a column called %s as defined in cfg.property'], cfg.property)
end

if istrue(cfg.average)
    if isfield(cfg,'maskproperty')
        ft_error(['cfg.maskproperty not supported for cfg.average = ''yes'''])
    end
    res.table = res.table{ismember(fv,cfg.filtervalues{iCol}),:};
    res.table = varfun(@nanmean,res.table,'InputVariables',cfg.property,'GroupingVariables','channel');
    res.table.Properties.VariableNames{end} = cfg.property;
end

cfg2 = cfg;

data = [];
data.dimord = 'chan_time';
data.time = [0];
data.label = res.table.channel;
data.avg = res.table.(cfg.property);
cfg2.parameter  = 'avg';
if isfield(cfg,'maskproperty')
    cfg2.maskparameter  = 'mask';
    data.mask = res.table.(cfg.maskproperty);
end

cfg2 = ft_topoplotER(cfg2,data);

ch = [];
if ~strcmp(cfg2.colorbar, 'no')
    ch = colorbar;
end
fh = gcf;
ah = gca;
%cfg = cfg2;

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end

