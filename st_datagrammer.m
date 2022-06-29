function data = st_datagrammer(cfg, data)

% ST_DATAGRAMMER applies data grammer like filters to the data.
%
% Use as
%   data = st_datagrammer(cfg, data)
%
% Available configuration parameters are:
%   cfg.grammer  = string with data operations, (default = '')
%   cfg.channel  = Nx1 cell-array with selection of channels the grammer
%                  will be applied to.
%                  The rest of channels is left unaffected (default = 'all'),
%                  see FT_CHANNELSELECTION
%
%   GRAMMERS TO USE ARE:
%     'normtomean' norms to a mean ignoring nans
%     '< [threshold] [replace_with_constant]' replaces all values with a contstant below a certain threshold
%     '> [threshold] [replace_with_constant]' replaces all values with a contstant above a certain threshold
%     '<= [threshold] [replace_with_constant]' replaces all values with a contstant below a certain threshold including that threshold
%     '>= [threshold] [replace_with_constant]' replaces all values with a contstant above a certain threshold including that threshold
%     '= [constant] [replace_with_constant]' replaces all values with a contstant that equal to the constant
%     'between [value1] [value2] [replace_with_constant]' replaces all values between value1 and value2 with a contstant 
%     'between_including [value1] [value2] [replace_with_constant]' replaces all values between value1 and value2 with a contstant including values equal to value1 or value2 
%     'min' sets the channel signal to a constant equivalent of the minimum of the signal
%     'max' sets the channel signal to a constant equivalent of the maximum of the signal
%     'ztransform' z-transforms the data
%     'smooth [window_seconds] [[type]]' smoothes the signal with a (carbox) filter with window size of [window_seconds] using the default smoothing type called 'moving' (carbox filter), other types are listed below
%     'medianfilter [window_seconds] [[nan_treatment]] [[paddingType]]' uses a median filter with window size of [window_seconds] and optional nan treatment (either 'includenan' (default) or 'omitnan') and optional padding type ('zeropad' (default) | 'truncate'), see matlab medfilt1 documentation for detail
%     'envpeaks [min_interval_seconds]' get the peaks of the envelope with an interval of at least [min_interval_seconds]
%     'rect' is rectifying the data (i.e. abs(signal))
%     'env' in creating the hilbert envelope of the signal
%     'dev' creating the (first) deviation of the signal, note that the deviation to the first sample is always 0 to keep trial length
%     'mult [factor]' multiplies the data with a constant [factor]
%     'nanreplace [constant]' replaces nans in the signal with a constant number
%     '( ... ) naninsert ( ... )' inserts nans from the first term signal before it into the second term signal after it
%     'bp [low_Hz] [high_Hz] [[filter_order]] [[iir|fir]]' applies a band-pass filter with to preserve signal between [low_Hz] [high_Hz] with optional parameters of [filter_order] and filter type as either fir or iir fitler [[fir|iir]]
%     'hp [low_Hz] [[filter_order]] [[iir|fir]]' applies a high-pass filter with to preserve signal above [low_Hz] with optional parameters of [filter_order] and filter type as either fir or iir fitler [[fir|iir]]
%     'lp [high_Hz] [[filter_order]] [[iir|fir]]' applies a low-pass filter with to preserve signal below [high_Hz] with optional parameters of [filter_order] and filter type as either fir or iir fitler [[fir|iir]]
%     'phase_wrapped_rad' for transformation into wrappted phase in radians [rad, +-pi] from the positive peak being t = 0, please use apply a filter band first 
%     'phase_wrapped_turns' for transformation into wrappted phase in turns [0 to 1 turn]from the positive peak being t = 0, please use apply a filter band first 
%     'phase_wrapped_deg' for transformation into wrappted phase in degree [0 to 360] from the positive peak being t = 0, please use apply a filter band first 
%     'phase_unwrapped_rad' for transformation into unwrappted (non-zigsag but continous) phase in radians [rad, +-pi] from the positive peak being t = 0, please use apply a filter band first 
%     'phase_unwrapped_turns' for transformation into unwrappted (non-zigsag but continous) phase in turns [0 to 1 turn]from the positive peak being t = 0, please use apply a filter band first 
%     'phase_unwrapped_deg' for transformation into unwrappted (non-zigsag but continous) phase in degree [0 to 360] from the positive peak being t = 0, please use apply a filter band first 
%     'instfreq' for transformation into instantaenous frequency, please use apply a filter band first
%     'waveletband [low_Hz] [high_Hz] [[filter_order]]' applies a band-pass filter based on wavelet filtering with to preserve signals activity between [low_Hz] [high_Hz] with optional parameters of [cycle_witdh] fro the width of the wavelet in number of cycles
%     'no' apply nothing
%     '( ... ) + ( ... )' add signals together sample by sample of two terms enclosed in brackets
%     '( ... ) - ( ... )' substract signals from each other sample by sample of two terms  enclosed in brackets
%     '( ... )' enclose operations in brackets
%
%     EXAMPLE:
%     cfg.grammer = 'hp 0.3 lp 35'
%     cfg.grammer = '( ( bp 12 14 mult 4 ) + ( hp 0.5 lp 2 ) ) + ( bp 6 8 mult 2 )'
%
%    smooth types (as decribed in matlab documentation for function smooth):
%         'moving'
%         Moving average (default). A lowpass filter with filter coefficients equal to the reciprocal of the span.
%         'lowess'
%         Local regression using weighted linear least squares and a 1st degree polynomial model
%         'loess'
%         Local regression using weighted linear least squares and a 2nd degree polynomial model
%         'sgolay'
%         Savitzky-Golay filter. A generalized moving average with filter coefficients determined by an unweighted linear least-squares regression and a polynomial model of specified degree (default is 2). The method can accept nonuniform predictor data.
%         'rlowess'
%         A robust version of 'lowess' that assigns lower weight to outliers in the regression. The method assigns zero weight to data outside six mean absolute deviations.
%         'rloess'
%         A robust version of 'loess' that assigns lower weight to outliers in the regression. The method assigns zero weight to data outside six mean absolute deviations.
%
% See also FT_PREPROCESSING, FT_APPLY_MONTAGE

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

% set defaults
cfg.channel  = ft_getopt(cfg, 'channel', 'all', 1);
cfg.grammer  = ft_getopt(cfg, 'grammer', '', 1);


fprintf([functionname ' function initialized\n']);

channel_all = data.label;
channel_to_be_processed = ft_channelselection(cfg.channel, data.label);

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

cfg.grammer = assure_wellformed_grammer(cfg.grammer);

curr_filterdefs = strtrim(cfg.grammer);
curr_filterdefs = strsplit(char(curr_filterdefs));
if (isempty(curr_filterdefs))
    curr_filterdefs = 'no';
end


%init the filter order variables
%usedFilterOrder_lp = -1;
%usedFilterOrder_hp = -1;
%usedFilterOrder_bp = -1;

used_specified_filtering_once_bp = false;
used_specified_filtering_once_hp = false;
used_specified_filtering_once_lp = false;



data_new = {};
iChanCount = 1;
for iChannelbyOrder = 1:numel(channel_all)
    curr_channel_label = channel_all{iChannelbyOrder};
    
    cfg_tmp = [];
    cfg_tmp.channel = curr_channel_label;
    data_new{iChanCount} = ft_selectdata(cfg_tmp, data);
    
    if ~ismember(channel_all{iChannelbyOrder},channel_to_be_processed)
        iChanCount = iChanCount + 1;
        continue
    else
        
        curr_filterdefs_filterPos = 1;
        
        curr_bracket_opened = 0;
        %curr_operator_applied = false;
        curr_operator_wait_for_second_term = {};
        %curr_operator = '';
        data_store = {};%stack for brackets
        data_stack = {};%stack for explict storage
        data_stack_name = {};%stack for explict storage names
        
        used_specified_filtering_once_bp = false;
        used_specified_filtering_once_hp = false;
        used_specified_filtering_once_lp = false;
        
        while curr_filterdefs_filterPos <= length(curr_filterdefs)
            curr_filter = char(curr_filterdefs(curr_filterdefs_filterPos));
            FpassLeft = -1;
            FpassRight = -1;
            MultFactor = 1;
            SmoothSeconds = -1;
            medianfilterSeconds = -1;
            conversion_state_success = true;
            useSpecifiedFilterOrder = false;
            overwriteFilterType = false;
            filterType = '';
            smoothingType = 'moving';
            median_nan_treatment = 'omitnan';
            median_paddingType = 'zeropad';
            nanreplace_constant = 0;
            threshold_value = 0;
            replace_value = 0;
            between_lower = 0;
            between_higher = 0;
            try
                switch curr_filter
                    case 'push'
                        name = char(curr_filterdefs(curr_filterdefs_filterPos+1));
                         data_stack_name{end+1} = name;%stack
                        data_stack{end+1} = data_new{iChanCount};%stack
                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 2;
                        continue;
                    case 'get'
                        name = char(curr_filterdefs(curr_filterdefs_filterPos+1));
                        [~, iStack] = ismember(name, data_stack_name);
                        if iStack > 0
                            data_new{iChanCount} = data_stack{iStack};
                        else 
                            ft_error(['converion of grammer parameter failed for filter settings of channel '  curr_channel_label ', could not retrieve the pushed state called ' name])
                        end
                        data_stack{end+1} = data_new{iChanCount};%stack
                        data_stack_name{end+1} = name;%stack
                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 2;
                        continue;
                    case {'<' '>' '=' '>=' '<='}
                        [threshold_value conversion_state_success1] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+1)));
                        [replace_value   conversion_state_success2] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+2)));
                        conversion_state_success = conversion_state_success1 & conversion_state_success2;
                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 3;
                   case {'between' 'between_including'}
                        [between1 conversion_state_success1] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+1)));
                        [between2 conversion_state_success2] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+2)));
                        [replace_value   conversion_state_success3] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+2)));
                        between_lower = min([between1 between2]);
                        between_higher = max([between1 between2]);
                        conversion_state_success = conversion_state_success1 & conversion_state_success2 & conversion_state_success3;

                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 4;
                    case 'normtomean'
                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 1;
                    case 'min'
                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 1;
                    case 'max'
                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 1;
                    case 'ztransform'
                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 1;
                    case 'smooth'
                        [SmoothSeconds conversion_state_success] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+1)));
                        temp_additional_values = 0;
                        if conversion_state_success
                            if (numel(curr_filterdefs) >= curr_filterdefs_filterPos+2)
                                pot_smoothingType = char(curr_filterdefs(curr_filterdefs_filterPos+3));
                                if any(ismember(pot_smoothingType,{'moving', 'lowess', 'loess', 'sgolay', 'sgolay', 'rlowess', 'rloess'}))
                                    temp_additional_values = temp_additional_values + 1;
                                    smoothingType = pot_smoothingType;
                                end
                            end
                        end
                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 2 + temp_additional_values;
                    case 'medianfilter'
                        % 'medianfilter [window_seconds] [[nan_treatment]] [[paddingType]]' uses a median filter with window size of [window_seconds] and optional nan treatment (either 'includenan' (default) or 'omitnan') and optional padding type ('zeropad' (default) | 'truncate'), see matlab medfilt1 documentation for detail
                        [medianfilterSeconds conversion_state_success] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+1)));
                        temp_additional_values = 0;
                        
                        if conversion_state_success
                            
                            if (numel(curr_filterdefs) >= curr_filterdefs_filterPos+2)
                                pot_param = char(curr_filterdefs(curr_filterdefs_filterPos+3));
                                
                                if any(ismember(pot_smoothingType,{'omitnan', 'includenan'}))
                                    temp_additional_values = temp_additional_values + 1;
                                    median_nan_treatment = pot_param;
                                elseif any(ismember(pot_smoothingType,{'zeropad', 'truncate'}))
                                    temp_additional_values = temp_additional_values + 1;
                                    median_paddingType = pot_param;
                                end
                                if (numel(curr_filterdefs) >= curr_filterdefs_filterPos+3)
                                    pot_param = char(curr_filterdefs(curr_filterdefs_filterPos+4));
                                    if any(ismember(pot_smoothingType,{'omitnan', 'includenan'}))
                                        temp_additional_values = temp_additional_values + 1;
                                        median_nan_treatment = pot_param;
                                    elseif any(ismember(pot_smoothingType,{'zeropad', 'truncate'}))
                                        temp_additional_values = temp_additional_values + 1;
                                        median_paddingType = pot_param;
                                    end
                                end
                                
                            end
                        end
                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 2 + temp_additional_values;
                    case 'envpeaks'
                        [minPeakDistanceSeconds conversion_state_success] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+1)));
                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 2;
                    case 'rect'
                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 1;
                    case 'env'
                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 1;
                    case 'dev'
                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 1;
                    case 'instfreq'
                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 1;
                    case {'phase_wrapped_rad', 'phase_wrapped_deg', 'phase_wrapped_turns', 'phase_unwrapped_rad', 'phase_unwrapped_deg', 'phase_unwrapped_turns'}
                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 1;
                    case 'mult'
                        [MultFactor conversion_state_success] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+1)));
                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 2;
                    case 'nanreplace'
                        [nanreplace_constant conversion_state_success] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+1)));
                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 2;
                    case {'bp', 'waveletband'}
                        [FpassLeft conversion_state_success] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+1)));
                        [FpassRight conversion_state_success] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+2)));
                        temp_additional_values = 0;
                        if (numel(curr_filterdefs) >= curr_filterdefs_filterPos+3)
                            [pot_additional_value conversion_state_success_potaddval] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+3)));
                            if conversion_state_success_potaddval
                                useSpecifiedFilterOrder = true;
                                filterOrder = round(pot_additional_value);
                                temp_additional_values = temp_additional_values + 1;
                                
                                if (numel(curr_filterdefs) >= curr_filterdefs_filterPos+4) && strcmp(curr_filter,'bp')
                                    pot_filterType = char(curr_filterdefs(curr_filterdefs_filterPos+4));
                                    if strcmp(pot_filterType,'iir') || strcmp(pot_filterType,'fir')
                                        overwriteFilterType = true;
                                        used_specified_filtering_once_bp = true;
                                        filterType = pot_filterType;
                                        temp_additional_values = temp_additional_values + 1;
                                    end
                                end
                                if strcmp(curr_filter,'waveletband')
                                    used_specified_filtering_once_bp = true;
                                end
                            end
                        end
                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 3 + temp_additional_values;
                    case 'hp'
                        [FpassLeft conversion_state_success] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+1)));
                        temp_additional_values = 0;
                        if (numel(curr_filterdefs) >= curr_filterdefs_filterPos+2)
                            [pot_additional_value conversion_state_success_potaddval] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+2)));
                            if conversion_state_success_potaddval
                                useSpecifiedFilterOrder = true;
                                filterOrder = round(pot_additional_value);
                                temp_additional_values = temp_additional_values + 1;
                                if (numel(curr_filterdefs) >= curr_filterdefs_filterPos+3)
                                    pot_filterType = char(curr_filterdefs(curr_filterdefs_filterPos+3));
                                    if strcmp(pot_filterType,'iir') || strcmp(pot_filterType,'fir')
                                        overwriteFilterType = true;
                                        used_specified_filtering_once_hp = true;
                                        filterType = pot_filterType;
                                        temp_additional_values = temp_additional_values + 1;
                                    end
                                end
                            end
                        end
                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 2 + temp_additional_values;
                    case 'lp'
                        [FpassRight conversion_state_success] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+1)));
                        temp_additional_values = 0;
                        if (numel(curr_filterdefs) >= curr_filterdefs_filterPos+2)
                            [pot_additional_value conversion_state_success_potaddval] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+2)));
                            if conversion_state_success_potaddval
                                useSpecifiedFilterOrder = true;
                                filterOrder = round(pot_additional_value);
                                temp_additional_values = temp_additional_values + 1;
                                if (numel(curr_filterdefs) >= curr_filterdefs_filterPos+3)
                                    pot_filterType = char(curr_filterdefs(curr_filterdefs_filterPos+3));
                                    if strcmp(pot_filterType,'iir') || strcmp(pot_filterType,'fir')
                                        overwriteFilterType = true;
                                        used_specified_filtering_once_lp = true;
                                        filterType = pot_filterType;
                                        temp_additional_values = temp_additional_values + 1;
                                    end
                                end
                            end
                        end
                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 2 + temp_additional_values;
                    case 'no'
                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 1;
                    case {'+' '-' 'naninsert'}
                        curr_operator_wait_for_second_term{end+1} = curr_filter;
                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 1;
                        continue;
                    case '('
                        curr_bracket_opened = curr_bracket_opened + 1;
                        if numel(curr_operator_wait_for_second_term) > 0
                            data_new{iChanCount} = data_store{end};
                            data_store(end) = [];
                        else
                            if numel(data_store) < 1
                                data_store{end+1} = data_new{iChanCount};%stack
                            else
                                data_store{end+1} = data_store{end};%stack
                            end
                        end
                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 1;
                        continue;
                    case ')'
                        if numel(curr_operator_wait_for_second_term) > 0
                            switch curr_operator_wait_for_second_term{end}
                                case '+'
                                    data_new{iChanCount}.trial{1} = data_term1.trial{1} + data_new{iChanCount}.trial{1};
                                    data_term1 = [];
                                case '-'
                                    data_new{iChanCount}.trial{1} = data_term1.trial{1} - data_new{iChanCount}.trial{1};
                                    data_term1 = [];
                                case 'naninsert'
                                    data_new{iChanCount}.trial{1}(isnan(data_term1.trial{1})) = NaN;
                                    data_term1 = [];
                                otherwise
                                    error(['filter not well defined for channel ' curr_channel_label])
                            end
                            curr_operator_wait_for_second_term(end) = [];
                        else
                            if (curr_bracket_opened > 0) && ~(numel(curr_operator_wait_for_second_term) > 0)
                                data_term1 = data_new{iChanCount};
                            end
                        end
                        curr_bracket_opened = curr_bracket_opened - 1;
                        curr_filterdefs_filterPos = curr_filterdefs_filterPos + 1;
                        continue;
                    otherwise
                        ft_error(['grammer not well defined for channel ' curr_channel_label])
                end
            catch err
                ft_error(['grammer not well defined for channel ' curr_channel_label])
            end
            
            %[numeral conversion_successfull] = str2num(curr_filter);
            %            if conversion_successfull
            %                data_filt{iChanCount}.trial{1} = numeral;
            %            else
            
            
            if ~conversion_state_success
                ft_error(['converion of grammer parameter failed for filter settings of channel '  curr_channel_label])
            end
            
            cfg_tmp = [];
            cfg_tmp = core_cfg;
            if curr_filterdefs_filterPos > 1
                cfg_tmp.dftfilter = 'no';
            end
            cfg_tmp.channel = curr_channel_label;
            
            switch curr_filter
                case '<'
                    data_new{iChanCount}.trial{1}(data_new{iChanCount}.trial{1} < threshold_value) = replace_value;
                case '>'
                    data_new{iChanCount}.trial{1}(data_new{iChanCount}.trial{1} > threshold_value) = replace_value;
                case '='
                    data_new{iChanCount}.trial{1}(data_new{iChanCount}.trial{1} == threshold_value) = replace_value;
                case '>='
                    data_new{iChanCount}.trial{1}(data_new{iChanCount}.trial{1} >= threshold_value) = replace_value;
                case '<='
                    data_new{iChanCount}.trial{1}(data_new{iChanCount}.trial{1} <= threshold_value) = replace_value;
                case 'between'
                    data_new{iChanCount}.trial{1}((between_lower < data_new{iChanCount}.trial{1}) & (data_new{iChanCount}.trial{1} < between_higher)) = replace_value;
                case 'between_including'
                    data_new{iChanCount}.trial{1}((between_lower <= data_new{iChanCount}.trial{1}) & (data_new{iChanCount}.trial{1} <= between_higher)) = replace_value;
                case 'normtomean'
                    data_new{iChanCount}.trial{1} = (data_new{iChanCount}.trial{1} - nanmean(data_new{iChanCount}.trial{1}));
                case 'min'
                    data_new{iChanCount}.trial{1} = repmat(min(data_new{iChanCount}.trial{1}),size(data_new{iChanCount}.trial{1}));
                case 'max'
                    data_new{iChanCount}.trial{1} = repmat(max(data_new{iChanCount}.trial{1}),size(data_new{iChanCount}.trial{1}));
                case 'ztransform'
                    data_new{iChanCount}.trial{1} = (data_new{iChanCount}.trial{1} - nanmean(data_new{iChanCount}.trial{1}))./nanstd(data_new{iChanCount}.trial{1});
                case 'smooth'
                    data_new{iChanCount}.trial{1} = smooth(data_new{iChanCount}.trial{1},max(1,round(SmoothSeconds*data_new{iChanCount}.fsample)), smoothingType)';
                case 'medianfilter'
                    data_new{iChanCount}.trial{1} = medfilt1(data_new{iChanCount}.trial{1},max(1,round(medianfilterSeconds*data_new{iChanCount}.fsample)), median_nan_treatment, median_paddingType)';
                case 'envpeaks'
                    data_new{iChanCount}.trial{1} = envpeaks(data_new{iChanCount}.trial{1},data_new{iChanCount}.fsample, minPeakDistanceSeconds);
                case 'rect'
                    data_new{iChanCount}.trial{1} = abs(data_new{iChanCount}.trial{1});
                case 'env'
                    data_new{iChanCount}.trial{1} = abs(hilbert(data_new{iChanCount}.trial{1}));
                case 'dev'
                    data_new{iChanCount}.trial{1} = [0 diff(data_new{iChanCount}.trial{1})];
                case 'instfreq'
                    instfreq = (diff(unwrap(angle(hilbert(data_new{iChanCount}.trial{1})))) / (2.0 * pi) * fsample); 
                    if numel(instfreq) >= 1
                        data_new{iChanCount}.trial{1} = [instfreq(1) instfreq];
                    end
                    
                    data_new{iChanCount}.trial{1} = angle(hilbert(data_new{iChanCount}.trial{1}));
                    
%                     t = 0:(1/100):10; 
%                     x = cos(pi/0.1*t);
% 
%                     fsample = 100;
%                     phase_raw_rad = angle(hilbert(x));
%                     phase_raw_deg = 360*phase_raw_rad/(2.0 * pi);
%                     phase_raw_turns = phase_raw_rad/(2.0 * pi);
%                     phase_raw_turns_unwarp = unwrap(phase_raw_turns);
%                     instfreq = (diff(phase_raw_rad) / (2.0 * pi) * fsample);
%                     instfreq = [instfreq(1) instfreq];
%                     figure
%                     plot(t,x)
%                     hold on
%                     plot(t,phase_raw_rad)
%                     plot(t,phase_raw_deg)
%                     plot(t,phase_raw_turns)
%                     plot(t,phase_raw_turns_unwarp)
%                     figure
%                     plot(t,instfreq)

                case 'phase_wrapped_rad'
                    data_new{iChanCount}.trial{1} = angle(hilbert(data_new{iChanCount}.trial{1}));
                case 'phase_wrapped_deg'
                    data_new{iChanCount}.trial{1} = 360*angle(hilbert(data_new{iChanCount}.trial{1}))/(2*pi);
                case 'phase_wrapped_turns'
                    data_new{iChanCount}.trial{1} = angle(hilbert(data_new{iChanCount}.trial{1}))/(2*pi);
                case 'phase_unwrapped_rad'
                    data_new{iChanCount}.trial{1} = unwrap(angle(hilbert(data_new{iChanCount}.trial{1})));
                case 'phase_unwrapped_deg'
                    data_new{iChanCount}.trial{1} = 360*unwrap(angle(hilbert(data_new{iChanCount}.trial{1})))/(2*pi);
                case 'phase_unwrapped_turns'
                    data_new{iChanCount}.trial{1} = unwrap(angle(hilbert(data_new{iChanCount}.trial{1})))/(2*pi);
                case 'mult'
                    data_new{iChanCount}.trial{1} = data_new{iChanCount}.trial{1} .* MultFactor;
                case 'nanreplace'
                    data_new{iChanCount}.trial{1}(isnan(data_new{iChanCount}.trial{1})) = nanreplace_constant;
                %case 'naninsert'
                %    data_new{iChanCount}.trial{1}(isnan(data_new{iChanCount}.trial{1})) = nanreplace_constant;
                case 'waveletband'
                    usedFilterOrder_bp = NaN;
                    bp_hdm = NaN;
                    
                    cfg_tmp.method = 'wavelet';
                    cfg_tmp.output = 'pow';
                    FreqSteps = abs(FpassRight-FpassLeft)/10;
                    cfg_tmp.foi = [FpassLeft:FreqSteps:FpassRight];
                    if useSpecifiedFilterOrder
                        cfg_tmp.width = filterOrder;
                    else
                        cfg_tmp.width =  4;
                    end
                    %cfg.pad = 'maxperlen';
                    cfg_tmp.toi = data_new{iChanCount}.time{1};
                    cfg_tmp.pad = ( 2*cfg_tmp.width*(1/min(cfg_tmp.foi)) ) + ( max(cfg_tmp.toi)-min(cfg_tmp.toi) ) +  10 ;
                    cfg_tmp.feedback = core_cfg.feedback;
                    fprintf('reprocess and apply band wavlet filter to %s\n',curr_channel_label);
                    data_freq = ft_freqanalysis(cfg_tmp,data_new{iChanCount});
                    temp_freqdat = squeeze(nanmean(data_freq.powspctrm,2))';
                    temp_freqdat(isnan(temp_freqdat)) = 0;
                    data_new{iChanCount}.trial = {temp_freqdat};
                    temp_freqdat = [];
                    dat_freq = [];
                case 'bp'
                    cfg_tmp.bpfilter = 'yes';
                    FstopLeft = FpassLeft - StopToPassTransitionWidth_bp; %left stop frequency in Hz
                    FstopRight = FpassRight + PassToStopTransitionWidth_bp; %left stop frequency in Hz
                    
                    usedFilterOrder_bp = NaN;
                    bp_hdm = NaN;
                    if strcmp(core_cfg.bpfilttype,'IIRdesigned') || strcmp(core_cfg.bpfilttype,'FIRdesigned') && ~overwriteFilterType
                        bp_d = [];
                        bp_hd = [];
                        fprintf('designing band pass filter\n');
                        if strcmp(UseFixedFilterOrder_bp,'yes')
                            bp_d = fdesign.bandpass('N,Fst1,Fp1,Fp2,Fst2',FilterOrder_bp,FstopLeft,FpassLeft,FpassRight,FstopRight,fsample);
                            bp_hd = design(bp_d,'equiripple');
                        else
                            bp_d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',FstopLeft,FpassLeft,FpassRight,FstopRight,AstopLeft_bp,Apass_bp,AstopRight_bp,fsample);
                            bp_hd = design(bp_d,'equiripple','MinOrder', 'even');
                        end
                        usedFilterOrder_bp = bp_hd.order;
                        cfg_tmp.bpfilterdesign = bp_hd;
                        bp_hdm = measure(bp_hd);
                    else
                        cfg_tmp.bpinstabilityfix = 'split';
                        if isempty(filterType)
                            switch core_cfg.hpfilttype
                                case {'IIRdesigned','but','iir'}
                                    filterType = 'but';
                                case {'FIRdesigned','fir'}
                                    filterType = 'fir';
                                otherwise
                                    ft_warning('filter type chosen as FieldTrip default')
                                    filterType = '';
                            end
                        end
                        switch filterType
                            case 'iir'
                                cfg_tmp.bpfilttype =  'but';
                            case 'fir'
                                cfg_tmp.bpfilttype =  'fir';
                            otherwise
                                %error(['filter type ' filterType ' unknown'])
                        end
                    end
                    
                    if strcmp(UseFixedFilterOrder_bp,'yes') || useSpecifiedFilterOrder
                        if useSpecifiedFilterOrder
                            cfg_tmp.bpfiltord = filterOrder;
                        else
                            cfg_tmp.bpfiltord = FilterOrder_bp;
                        end
                        usedFilterOrder_bp = cfg_tmp.bpfiltord;
                    end
                    cfg_tmp.bpfreq        = [FpassLeft FpassRight];%dummy values are overwritten by low level function
                    cfg_tmp.feedback = core_cfg.feedback;
                    fprintf('reprocess and apply band filter to %s\n',curr_channel_label);
                    data_new{iChanCount} = st_preprocessing(cfg_tmp,data_new{iChanCount});
                    
                case 'hp'
                    cfg_tmp.hpfilter = 'yes';
                    FstopLeft = FpassLeft - StopToPassTransitionWidth_hp; %left stop frequency in Hz
                    
                    usedFilterOrder_hp = NaN;
                    hp_hdm = NaN;
                    if strcmp(core_cfg.hpfilttype,'IIRdesigned') || strcmp(core_cfg.hpfilttype,'FIRdesigned') && ~overwriteFilterType
                        hp_d = [];
                        hp_hd = [];
                        if strcmp(UseFixedFilterOrder_hp,'yes')
                            hp_d = fdesign.highpass('N,F3db',FilterOrder_hp,FpassLeft,fsample);%designmethods(hp_preDS_d); help(hp_preDS_d,'equiripple'); help(hp_preDS_d,'butter')
                        else
                            hp_d = fdesign.highpass('Fst,Fp,Ast,Ap',FstopLeft,FpassLeft,AstopLeft_hp,Apass_hp,fsample);%designmethods(hp_preDS_d); help(hp_preDS_d,'equiripple'); help(hp_preDS_d,'butter')
                        end
                        fprintf('designing high pass filter\n');
                        if strcmp(core_cfg.hpfilttype,'IIRdesigned')
                            hp_hd = design(hp_d,'butter'); %isstable(hp_hd)
                        elseif strcmp(core_cfg.hpfilttype,'FIRdesigned')
                            hp_hd = design(hp_d,'equiripple','MinOrder', 'even');
                        else
                            error(['highpass filter type of ' core_cfg.hpfilttype ' unknown or not allowed'])
                        end
                        usedFilterOrder_hp = hp_hd.order;
                        cfg_tmp.hpfilterdesign = hp_hd;
                        hp_hdm = measure(hp_hd);
                    else
                        cfg_tmp.hpinstabilityfix = 'split';
                        if isempty(filterType)
                            switch core_cfg.hpfilttype
                                case {'IIRdesigned','but','iir'}
                                    filterType = 'but';
                                case {'FIRdesigned','fir'}
                                    filterType = 'fir';
                                otherwise
                                    ft_warning('filter type chosen as FieldTrip default')
                                    filterType = '';
                            end
                        end
                        
                        switch filterType
                            case 'iir'
                                cfg_tmp.hpfilttype =  'but';
                            case 'fir'
                                cfg_tmp.hpfilttype =  'fir';
                            otherwise
                                %error(['filter type ' filterType ' unknown'])
                        end
                    end
                    if strcmp(UseFixedFilterOrder_hp,'yes') || useSpecifiedFilterOrder
                        if useSpecifiedFilterOrder
                            cfg_tmp.hpfiltord = filterOrder;
                        else
                            cfg_tmp.hpfiltord = FilterOrder_hp;
                        end
                        usedFilterOrder_hp = cfg_tmp.hpfiltord;
                    end
                    cfg_tmp.hpfreq        = [FpassLeft];%dummy values are overwritten by low level function
                    cfg_tmp.feedback = core_cfg.feedback;
                    fprintf('reprocess and apply high pass filter to %s\n',curr_channel_label);
                    data_new{iChanCount} = st_preprocessing(cfg_tmp,data_new{iChanCount});
                case 'lp'
                    cfg_tmp.lpfilter = 'yes';
                    FstopRight = FpassRight + PassToStopTransitionWidth_lp; %right stop frequency in Hz
                    usedFilterOrder_lp = NaN;
                    lp_hdm = NaN;
                    if strcmp(core_cfg.lpfilttype,'IIRdesigned') || strcmp(core_cfg.lpfilttype,'FIRdesigned') && ~overwriteFilterType
                        lp_d = [];
                        lp_hd = [];
                        fprintf('dataset: designing low pass filter \n');
                        if strcmp(UseFixedFilterOrder_lp,'yes')
                            lp_d = fdesign.lowpass('N,Fp,Fst',FilterOrder_lp,FpassRight,FstopRight,fsample);
                            lp_hd = design(lp_d,'equiripple');
                        else
                            lp_d = fdesign.lowpass('Fp,Fst,Ap,Ast',FpassRight,FstopRight,Apass_lp,AstopRight_lp,fsample);
                            lp_hd = design(lp_d,'equiripple','MinOrder', 'even');
                        end
                        usedFilterOrder_lp = lp_hd.order;
                        cfg_tmp.lpfilterdesign = lp_hd;
                        lp_hdm = measure(lp_hd);
                    else
                        cfg_tmp.lpinstabilityfix = 'split';
                        if isempty(filterType)
                            switch core_cfg.hpfilttype
                                case {'IIRdesigned','but','iir'}
                                    filterType = 'but';
                                case {'FIRdesigned','fir'}
                                    filterType = 'fir';
                                otherwise
                                    ft_warning('filter type chosen as FieldTrip default')
                                    filterType = '';
                            end
                        end
                        switch filterType
                            case 'iir'
                                cfg_tmp.lpfilttype =  'but';
                            case 'fir'
                                cfg_tmp.lpfilttype =  'fir';
                            otherwise
                                %error(['filter type ' filterType ' unknown'])
                        end
                    end
                    if strcmp(UseFixedFilterOrder_lp,'yes') || useSpecifiedFilterOrder
                        if useSpecifiedFilterOrder
                            cfg_tmp.lpfiltord = filterOrder;
                        else
                            cfg_tmp.lpfiltord = FilterOrder_lp;
                        end
                        usedFilterOrder_lp = cfg_tmp.lpfiltord;
                    end
                    cfg_tmp.lpfreq        = [FpassRight];%dummy values are overwritten by low level function
                    cfg_tmp.feedback = core_cfg.feedback;
                    fprintf('reprocess and apply low pass filter to %s\n',curr_channel_label);
                    data_new{iChanCount} = st_preprocessing(cfg_tmp,data_new{iChanCount});
                case 'no'
                    cfg_tmp.feedback = core_cfg.feedback;
                    fprintf('reprocess without filtering of %s\n',curr_channel_label);
                    data_new{iChanCount} = st_preprocessing(cfg_tmp,data_new{iChanCount});
                otherwise
                    error('filter in the filter definitions are not well defined, please check them')
            end
            
        end
    end
    iChanCount = iChanCount + 1;
end


% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous data

if length(data_new) > 1
    data = ft_appenddata([],data_new{:});
else
    data = data_new{:};
end

ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end

function signal = envpeaks(x,fsample,seconds)
[pks locs] = findpeaks(x,'MINPEAKDISTANCE',max(1,floor(seconds*fsample)));
signal = interp1([1 locs numel(x)],[x(1) pks x(end)],1:numel(x),'linear');
end

function g = assure_wellformed_grammer(g)
g = strrep(g, '(', ' ( ');
g = strrep(g, ')', ' ) ');
g = strrep(g, '+', ' + ');
g = strrep(g, '-', ' - ');

while ~strcmp(g,strrep(g, '  ', ' '))
    g = strrep(g, '  ', ' ');
end

g = strrep(g, ' - - ', ' + ');
g = strrep(g, ' + - ', ' - ');
g = strrep(g, ' - + ', ' - ');
g = strrep(g, ' + + ', ' + ');

g = strtrim(g);
end

