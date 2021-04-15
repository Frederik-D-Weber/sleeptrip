function [res_channel, res_event, res_summary] = st_rems(cfg, data)

% ST_REMS detects rapid eye movements during sleep and their properties.
% results are stored it in a result structure
%
% Use as
%   [res_channel, res_event, res_summary] = st_rems(cfg)
%   [res_channel, res_event, res_summary] = st_rems(cfg, data)
%
% Required configuration parameters are:
%   cfg.scoring         = structure provided by ST_READ_SCORING
%   cfg.channel         = Nx1 cell-array with selection of LEFT and RIGHT
%                         EOG channels in that order
%
%  if no data structure is provided as parameter then you need to define
%   cfg.dataset  = string with the filename
%
% THE DEFAULT options are an adapted and improved version of the Mölle et
%             al. 2002-2013 papers, with improved filters, narrower
%             frequency ranges (i.e. 2 Hz symetrically around center frequency)
%             correct window size of 0.2 s
%             using a stricter criterion to filter some spinldes like 200 uV
%             amplitude (cfg.minamplitude)
%             and minimal threshold factor of 1.75 of the std (instead of cfg.factorthresholdcriterion = 1.5)
%             and the maximal duration set to 2 seconds (instead of 3)
%
% Optional configuration parameters are:
%
%    cfg.ConsiderExcludedEpochs =   string, either 'yes' or 'no' if the excluded epochs in the scoring should not be considered at all
%                            (default = 'yes')
%    cfg.RemoveSignalThatIsNotREM  = string, either 'yes' or 'no' if the signal
%                            that is not REM should be excluded from the detection completely
%                            (default = 'no')
%    cfg.RelaxedThresholdFactorOfBasicThreshold  = number, factor mulitplied to the basic threshold (determined with the cfg.BasicThresholdFactorOfEOGnoise) that results in the relaxed threshold for REMs detection
%                            (default = 0.66)
%    cfg.BasicThresholdFactorOfEOGnoise  = number, factor of EOG noise that results in the basic threshold in REMs detection
%                            (default = 0.7)
%    cfg.amplitudemax       = number, the maximal amplitude for +-
%                             amplitude extension of the REMs in the signal in micro Volts (uV)
%                            (default = 100)
%    cfg.downsamplefs       = downsample the data to this frequency in Hz before doing the anlysis (default = 100/128)
%
% Some additional parameters from FT_PREPROCESSING can also be included
% including the reprocessing options that you can only use for EEG data are:
%
%     cfg.reref         = 'no' or 'yes' (default = 'no')
%     cfg.refchannel    = cell-array with new EEG reference channel(s), this can be 'all' for a common average reference
%     cfg.refmethod     = 'avg', 'median', or 'bipolar' for bipolar derivation of sequential channels (default = 'avg')
%     cfg.implicitref   = 'label' or empty, add the implicit EEG reference as zeros (default = [])
%     cfg.montage       = 'no' or a montage structure, see FT_APPLY_MONTAGE (default = 'no')
%
%
% See also ST_READ_SCORING, FT_PREPROCESSING, FT_APPLY_MONTAGE

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

if ~isfield(cfg, 'scoring')
    cfg.scoring = st_read_scoring(cfg);
end

scoring = cfg.scoring;

downsamplefsNotSet = false;
if ~isfield(cfg, 'downsamplefs')
    downsamplefsNotSet = true;
else
    ft_warning('This algorithm currently only supports only 100 Hz sampling rate and will thus resample at 100 Hz')
end


if isfield(cfg, 'montage') && isfield(cfg, 'channel')
    if ~any(ismember(cfg.channel,cfg.montage.labelnew))
        ft_error('The selected channels with cfg.channel does not match any of the ones defined in cfg.montage.\nPlease make them match or select all channels in the cfg.montage and NOT setting the cfg.channel')
    end
end



% set defaults
%cfg.channel  = ft_getopt(cfg, 'channel', 'all', 1);
%cfg.stages  = ft_getopt(cfg, 'stages', {'N2', 'N3', 'S4'});

cfg.amplitudemax = ft_getopt(cfg, 'amplitudemax', 100);

% EOG_low_pass_filter_freq = 5; % Hz
% EOG_high_pass_filter_freq = 0.5; % Hz
% basic_deviation_per_second_threshold = 261; %?V/s
% relaxed_deviation_per_second_threshold = 165; %?V/s
% 
% cfg.ConsiderExcludedEpochs = 'no';%
% cfg.RemoveSignalThatIsNotREM = 'no';%
% cfg.RelaxedThresholdFactorOfBasicThreshold = 0.66;
% cfg.BasicThresholdFactorOfEOGnoise = 0.7;

cfg.ConsiderExcludedEpochs  = ft_getopt(cfg, 'ConsiderExcludedEpochs', 'yes');
cfg.RemoveSignalThatIsNotREM  = ft_getopt(cfg, 'RemoveSignalThatIsNotREM', 'no');
cfg.RelaxedThresholdFactorOfBasicThreshold  = ft_getopt(cfg, 'RelaxedThresholdFactorOfBasicThreshold', 0.66);
cfg.BasicThresholdFactorOfEOGnoise  = ft_getopt(cfg, 'BasicThresholdFactorOfEOGnoise', 0.7);

% DEL:
% cfg.thresholdstages  = ft_getopt(cfg, 'thresholdstages', cfg.stages);
% cfg.leftofcenterfreq  = ft_getopt(cfg, 'leftofcenterfreq', 1);
% cfg.rightofcenterfreq  = ft_getopt(cfg, 'rightofcenterfreq', 1);
% cfg.thresholdaggmeth  = ft_getopt(cfg, 'thresholdaggmeth', 'respectivechan');
% cfg.envelopemeth  = ft_getopt(cfg, 'envelopemeth','smoothedRMSwd');
% cfg.factorthreshbeginend  = ft_getopt(cfg, 'factorthreshbeginend', 1.5);
% cfg.minduration  = ft_getopt(cfg, 'minduration', 0.5);
% cfg.maxduration  = ft_getopt(cfg, 'maxduration', 2.0);
% cfg.mergewithin  = ft_getopt(cfg, 'mergewithin', 0);
% cfg.minamplitude  = ft_getopt(cfg, 'minamplitude', 0);
% cfg.maxamplitude  = ft_getopt(cfg, 'maxamplitude', 120);
% cfg.filterSDamp  = ft_getopt(cfg, 'filterSDamp', 5);
% cfg.filterSDdur  = ft_getopt(cfg, 'filterSDdur', 5);
% cfg.filterSDfreq  = ft_getopt(cfg, 'filterSDfreq', 5);
cfg.downsamplefs     = ft_getopt(cfg, 'downsamplefs', 100);
% cfg.thresholdformbase  = ft_getopt(cfg, 'thresholdformbase', 'std');
% cfg.thresholdsignal  = ft_getopt(cfg, 'thresholdsignal', 'filtered_signal');
% cfg.rmstimewndw  = ft_getopt(cfg, 'rmstimewndw', 0.2);
% cfg.movavgtimewndw  = ft_getopt(cfg, 'movavgtimewndw', 0.2);
% 
% if ~isfield(cfg, 'factorthresholdcriterion')
%     switch cfg.envelopemeth
%         case 'smoothedRMSwd'
%             cfg.factorthresholdcriterion = 1.5;
%         case 'hilbertEnv'
%             cfg.factorthresholdcriterion = 2.25;
%         otherwise
%             cfg.factorthresholdcriterion = cfg.factorthreshbeginend;
%     end
% end
% 
% if cfg.factorthreshbeginend > cfg.factorthresholdcriterion
%     ft_error('cfg.factorthreshbeginend cannot be bigger than cfg.factorthresholdcriterion')
% end
% 
% UseAbsoluteEnvelopeThreshold = 'no';
% if isfield(cfg, 'envelopethresholds')
%     UseAbsoluteEnvelopeThreshold = 'yes'; % If abosolute positive envelope threshold for all channels should be used. In this case the factorSDbeginEnd and factorSDcriterion is not applied. either yes or no default no
%     AbsoluteEnvelopeThresholdBeginEnd = cfg.envelopethresholds(1); % Abosolute positive envelope potential threshold for all channels for begin and end of event default 4 (for microVolts potential)
%     AbsoluteEnvelopeThresholdCriterion = cfg.envelopethresholds(2); % Abosolute positive envelope potential threshold for all channels for minimum criterion to reach to count as an event. must be greater or equal than AbsoluteEnvelopeThresholdBeginEnd default 4 (for microVolts potential)
%     if AbsoluteEnvelopeThresholdBeginEnd > AbsoluteEnvelopeThresholdCriterion
%         error(['cfg.envelopethresholds must be a vector with two numbers, and the first cannot be greater than the second'])
%     end
% end
% 
% 
% if ~all(ismember(cfg.thresholdstages, cfg.stages))
% 	ft_error(['cfg.thresholdstages must be a subset of cfg.stages'])
% end
% thresholdingInDifferentStages = ~all(strcmp(cfg.stages,cfg.thresholdstages));

hasdata = false;
if nargin > 1
    % data structure provided
    % check if the input data is valid for this function, the input data must be raw
    data = ft_checkdata(data, 'datatype', {'raw+comp', 'raw'}, 'hassampleinfo', 'yes');
    if isfield(data, 'trial') && isfield(data, 'time')
        % check if data structure is likely continous data
        if (numel(data.trial) ~= 1) || ~all(size(data.sampleinfo) == [1 2])
            ft_error('data structure does not look like continous data and has more than one trial')
        end
    end
    hasdata = true;
    nSamplesInData = data.sampleinfo(2);
    fsample = data.fsample;
else
    % no data structure provided
    hdr = ft_read_header(cfg.dataset);
    nSamplesInData = hdr.nSamples;
    fsample = hdr.Fs;
end

% preCenterFreqFilterTo_FpassLeft = cfg.leftofcenterfreq;
% postCenterFreqFilterTo_FpassRight = cfg.rightofcenterfreq;  % in Hz
% 
% 
% 
% if ~(strcmp(cfg.envelopemeth,'hilbertEnv') || strcmp(cfg.envelopemeth,'smoothedRMSwd'))
%     error(['cfg.envelopemeth ' cfg.envelopemeth ' in parameters is unknown, use either hilbertEnv or smoothedRMSwd'])
% end
% 
% if ~(strcmp(cfg.thresholdformbase,'mean') || strcmp(cfg.thresholdformbase,'std'))
%     error(['cfg.thresholdformbase ' cfg.thresholdformbase ' in parameters is unknown, use either mean or std, default std'])
% end
% 
% if ~(strcmp(cfg.thresholdsignal,'filtered_signal') || strcmp(cfg.thresholdsignal,'envelope'))
%     error(['cfg.thresholdsignal ' cfg.thresholdsignal ' in parameters is unknown, use either filtered_signal or envelope, default filtered_signal'])
% end
% 
% factorThresholdCriterion = cfg.factorthresholdcriterion;
% 
% if scoring.epochlength < (1/(cfg.centerfrequency - cfg.leftofcenterfreq))
%     error(['the epoch length ' num2str(scoring.epochlength) 's must not be greater in order to support the requrested maximum of ' num2str((cfg.centerfrequency - cfg.leftofcenterfreq)) ' Hz!'])
% end
% 
% 
% if strcmp(UseAbsoluteEnvelopeThreshold,'yes')
%     cfg.factorthreshbeginend  = 1;
%     AbsoluteEnvelopeThresholdCriterionRatio = AbsoluteEnvelopeThresholdCriterion/AbsoluteEnvelopeThresholdBeginEnd;
%     factorThresholdCriterion = AbsoluteEnvelopeThresholdCriterionRatio;
% end
% 
% 
% 
% 
% set core parameters
load_core_cfg
% core_cfg
% 
% if (strcmp(core_cfg.bpfilttype,'FIRdesigned') || strcmp(core_cfg.bpfilttype,'IIRdesigned')) && (~strcmp(core_cfg.bpinstabilityfix,'no'))
%     error(['band pass filter instability fix not supported for FIRdesigned or IIRdesigned'])
% end
% 
if (strcmp(core_cfg.hpfilttype,'FIRdesigned') || strcmp(core_cfg.hpfilttype,'IIRdesigned')) && (~strcmp(core_cfg.hpinstabilityfix,'no'))
    error(['high pass filter instability fix not supported for FIRdesigned or IIRdesigned'])
end
% 
% Apass_bp = Apass;
% AstopLeft_bp = AstopLeft;
% AstopRight_bp = AstopRight;
% 
Apass_hp = Apass;
AstopLeft_hp = AstopLeft;

PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff = 0.5;

% 
% 
%filtfilt -- The length of the input x must be more than three times the filter order (N) defined as max(length(b)-1,length(a)-1).
%For best results, make sure the sequence you are filtering has length at least three times the filter order and tapers to zero on both edges.i.e.:
minSignalLengthSamples = scoring.epochlength*cfg.downsamplefs;
maxFilterOrder = floor((minSignalLengthSamples) / 3) - 1;
% 
% if strcmp(UseFixedFilterOrder_bp,'yes') && strcmp(useTwoPassFiltering_bp,'yes') && (FilterOrder_bp > maxFilterOrder)
%     error(['filter order for band pass not optimal, use a maximal filter order three times of samples in one epoch, i.e. filter order of ' num2str(maxFilterOrder)])
% elseif (FilterOrder_bp < maxFilterOrder) && strcmp(MaximizeFilterOrderIfFixedFilterOrderIsUsed,'yes')
%     FilterOrder_bp = maxFilterOrder;
% end
% 
if strcmp(UseFixedFilterOrder_hp,'yes') && strcmp(useTwoPassFiltering_hp,'yes') && (FilterOrder_hp > maxFilterOrder)
    error(['filter order for high pass not optimal, use a maximal filter order three times of samples in one epoch, i.e. filter order of ' num2str(maxFilterOrder)])
elseif (FilterOrder_hp < maxFilterOrder)  && strcmp(MaximizeFilterOrderIfFixedFilterOrderIsUsed,'yes')
    FilterOrder_hp = maxFilterOrder;
end
% 
% if strcmp(UseFixedFilterOrder_bp,'yes') && logical(mod(FilterOrder_bp,2))
%     error('band pass filter order must be an even number')
% end
% 
% if strcmp(UseFixedFilterOrder_hp,'yes') && logical(mod(FilterOrder_hp,2))
%     error('high pass order must be an even number')
% end
% 
% 
% % Note that a one- or two-pass filter has consequences for the
% % strength of the filter, i.e. a two-pass filter with the same filter
% % order will attenuate the signal twice as strong.
% if strcmp(useTwoPassFiltering_bp,'yes') && strcmp(UseTwoPassAttenuationCorrection_bp,'yes')
%     Apass_bp = Apass_bp/2;
%     AstopLeft_bp = AstopLeft_bp/2;
%     AstopRight_bp = AstopRight_bp/2;
% end
% 
if strcmp(useTwoPassFiltering_hp,'yes') && strcmp(UseTwoPassAttenuationCorrection_hp,'yes')
    Apass_hp = Apass_hp/2;
    AstopLeft_hp = AstopLeft_hp/2;
end

if strcmp(UseFixedFilterOrder_hp,'yes') && strcmp(core_cfg.hpfilttype,'FIRdesigned')
    error('UseFixedFilterOrder_hp not allowed for high pass filters of type FIRdesigned')
end
% 
% 
% if ~isnan(cfg.centerfrequency)
% 
% minFreq = cfg.centerfrequency - preCenterFreqFilterTo_FpassLeft;
% maxFreq = cfg.centerfrequency + postCenterFreqFilterTo_FpassRight;
% 
% if ~(maxFreq*3 < fsample)
%     error(['sample frequency of ' num2str(fsample) ' Hz must be MORE THAN three-fold (i.e. 3-fold) the maximal band frequency of ' num2str(maxFreq) ' Hz,\n even lower than requested by Nyquist-Shannon sample theorem.\n consider excludding higher frequency bands!'])
% end

fprintf([functionname ' function initalized\n']);



cfg_int = [];

use_hp = true;
use_lp = false;
use_bp = false;
StopToPassTransitionWidth_hp_temp = StopToPassTransitionWidth_hp_predownsample;
adapt_filter_settings_to_toolbox

cfg_int = core_cfg;
if isfield(cfg, 'reref'),       cfg_int.reref = cfg.reref;             end
if isfield(cfg, 'refchannel'),  cfg_int.refchannel = cfg.refchannel;   end
if isfield(cfg, 'refmethod'),   cfg_int.refmethod = cfg.refmethod;     end
if isfield(cfg, 'implicitref'), cfg_int.implicitref = cfg.implicitref; end
if isfield(cfg, 'montage'),     cfg_int.montage = cfg.montage;         end

FpassLeft = PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff; %left pass frequency in Hz
FstopLeft = FpassLeft - StopToPassTransitionWidth_hp_predownsample; %left stop frequency in Hz
usedFilterOrder_hp_preDS = NaN;
hp_preDS_hdm = [];
if PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff ~= 0
    cfg_int.hpfilter = 'yes';
    
    if strcmp(core_cfg.hpfilttype,'IIRdesigned') || strcmp(core_cfg.hpfilttype,'FIRdesigned')
        hp_preDS_d = [];
        hp_preDS_hd = [];
        if strcmp(UseFixedFilterOrder_hp,'yes')
            hp_preDS_d = fdesign.highpass('N,F3db',FilterOrder_hp,FpassLeft,fsample);%designmethods(hp_preDS_d); help(hp_preDS_d,'equiripple'); help(hp_preDS_d,'butter')
        else
            hp_preDS_d = fdesign.highpass('Fst,Fp,Ast,Ap',FstopLeft,FpassLeft,AstopLeft_hp,Apass_hp,fsample);%designmethods(hp_preDS_d); help(hp_preDS_d,'equiripple'); help(hp_preDS_d,'butter')
        end
        fprintf('designing high pass filter for pre downsampling filtering \n');
        if strcmp(core_cfg.hpfilttype,'IIRdesigned')
            hp_preDS_hd = design(hp_preDS_d,'butter'); %isstable(hp_preDS_hd)
        elseif strcmp(core_cfg.hpfilttype,'FIRdesigned')
            hp_preDS_hd = design(hp_preDS_d,'equiripple','MinOrder', 'even');
        else
            error(['highpass filter type of ' core_cfg.hpfilttype ' unknown or not allowed'])
        end
        usedFilterOrder_hp_preDS = hp_preDS_hd.order;
        cfg_int.hpfilterdesign = hp_preDS_hd;
        hp_preDS_hdm = measure(hp_preDS_hd);
    end
else
    cfg_int.hpfilter = 'no';
end
if strcmp(UseFixedFilterOrder_hp,'yes')
    cfg_int.hpfiltord = FilterOrder_hp;
end
cfg_int.hpfreq        = [FpassLeft];%dummy values are overwritten by low level function

if hasdata
    cfg_int.channel = ft_channelselection(cfg.channel, data.label);
else
    cfg_int.channel = ft_channelselection(cfg.channel, hdr.label);
end

cfg_int.feedback = core_cfg.feedback;
fprintf('preprocess and pre filter data\n');


%cfg_int.stages   = cfg.stages;
%cfg_int.scoring  = scoring;

hasROIs = true;
hasROIs_threshold = true;


if hasdata
    data_t = st_preprocessing(cfg_int, data);
    %data_t = st_select_scoring(cfg_int, data_t);
    if isempty(data_t.trial)
        hasROIs = false;
        hasROIs_threshold = false;
        % read in dummy data
        cfg_int.trl = [1 round(fsample*60) 0];
        data = st_preprocessing(cfg_int, data);
        %         data.time = {};
        %         data.trial = {};
        data.sampleinfo = [0 -1];
    else
        data = data_t;
        data_t = [];
    end
else
    cfg_int.dataset  = cfg.dataset;
    %[cfg_int] = st_select_scoring(cfg_int);
    cfg_int.continuous   = 'yes';
%     if isempty(cfg_int.trl)
%         hasROIs = false;
%         hasROIs_threshold = false;
% 
%         % read in dummy data
%         cfg_int.trl = [1 round(fsample*60) 0];
%         data = st_preprocessing(cfg_int);
%         %         data.time = {};
%         %         data.trial = {};
%         data.sampleinfo = [0 -1];
%     else
        data = st_preprocessing(cfg_int);
%     end
    
end

cfg_chan = [];
cfg_chan.channel = ft_channelselection(cfg.channel, data.label);
data = ft_selectdata(cfg_chan, data);

if ~hasROIs
    for iT = 1:numel(data.trial)
        data.trial{iT}(:) = NaN;
    end
end

cfg.downsamplefs = 100;

% if (data.fsample == 128) && downsamplefsNotSet
%     ft_warning('leaving 128 Hz sampling rate as default sampling rate')
%     cfg.downsamplefs = 128;
% end

preDownsampleFreq = data.fsample;
if (cfg.downsamplefs < data.fsample)
    fprintf('resample data from %i to %i Hz\n',data.fsample,cfg.downsamplefs);
    cfg_resample = [];
    cfg_resample.resamplefs = cfg.downsamplefs;
    cfg_resample.detrend = 'no';
    cfg_resample.feedback = core_cfg.feedback;
    data = ft_resampledata(cfg_resample,data);
elseif (cfg.downsamplefs > data.fsample)
    ft_error('upsampling prohibited, choose cfg.downsamplefs (now as %d Hz)  smaller or equal to the sampling frequency of provided data %d Hz', cfg.downsamplefs, data.fsample);
end

fsample = data.fsample;

trlSampleLengths = cellfun(@numel, data.time)';

epochLengthSamples = scoring.epochlength*fsample;
%signalOffsetSamples = 100
    %signalOffsetSeconds = 1;
    if (scoring.dataoffset ~= 0) %&& ~strcmp(DoEpochData,'yes')
        signalOffsetSamples_new = round(scoring.dataoffset*fsample);
        if (signalOffsetSamples_new ~= 0)
            cfg = [];
            if (signalOffsetSamples_new > 0) && (signalOffsetSamples_new+1+epochLengthSamples < size(data.trial{1},2))
                cfg.begsample = signalOffsetSamples_new+1;
                cfg.endsample = numel(data.time{1});
                
                data = ft_redefinetrial(cfg,data);
                cfg = [];
                cfg.offset = signalOffsetSamples_new+1;
                data = ft_redefinetrial(cfg,data);
                data.sampleinfo = [1 numel(data.time{1})];
                warning(['IMPORTANT: JUST CUT ' num2str(signalOffsetSamples_new) ' samples at the beginning corresponding to : ' num2str(signalOffsetSamples_new/fsample) ' seconds because of hypnogram/data offset!']);
            elseif signalOffsetSamples_new < 0
                for iTrTr = 1:numel(data.trial)
                    cfg.padtype = 'zero';
                    data.trial{iTrTr} = ft_preproc_padding(data.trial{iTrTr}, cfg.padtype, -signalOffsetSamples_new, 0);
                end
                data.time{1} = (0:(size(data.trial{1},2)-1))/data.fsample;
                data.sampleinfo = [1 numel(data.time{1})];
                warning(['IMPORTANT: JUST ADDED ' num2str(signalOffsetSamples_new) ' samples at the beginning corresponding to : ' num2str(signalOffsetSamples_new/fsample) ' seconds because of hypnogram/data offset!']);
            end
        end
    end



% use_hp = false;
% use_lp = false;
% use_bp = true;
% %StopToPassTransitionWidth_hp_temp = StopToPassTransitionWidth_hp_predownsample;
% adapt_filter_settings_to_toolbox
% 
% cfg_int = core_cfg;
% 
% cfg_int = [];
% cfg_int = core_cfg;
% cfg_int.bpfilter = 'yes';
% FpassLeft = minFreq; %left pass frequency in Hz
% FpassRight = maxFreq; %right pass frequency in Hz
% FstopLeft = FpassLeft - StopToPassTransitionWidth_bp; %left stop frequency in Hz
% FstopRight = FpassRight + PassToStopTransitionWidth_bp; %left stop frequency in Hz
% 
% usedFilterOrder_bp = NaN;
% bp_hdm = [];
% if strcmp(core_cfg.bpfilttype,'IIRdesigned') || strcmp(core_cfg.bpfilttype,'FIRdesigned')
%     bp_d = [];
%     bp_hd = [];
%     fprintf('designing band pass filter\n');
%     if strcmp(UseFixedFilterOrder_bp,'yes')
%         bp_d = fdesign.bandpass('N,Fst1,Fp1,Fp2,Fst2',FilterOrder_bp,FstopLeft,FpassLeft,FpassRight,FstopRight,fsample);
%         bp_hd = design(bp_d,'equiripple');
%     else
%         bp_d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',FstopLeft,FpassLeft,FpassRight,FstopRight,AstopLeft_bp,Apass_bp,AstopRight_bp,fsample);
%         %bp_hd = design(bp_d,'equiripple','MinOrder', 'even');
%         bp_hd = design(bp_d,'equiripple'); % order is even with 442 for 100 Hz 884 for 200 Hz
%     end
%     usedFilterOrder_bp = bp_hd.order;
%     cfg_int.bpfilterdesign = bp_hd;
%     bp_hdm = measure(bp_hd);
% end
% if strcmp(UseFixedFilterOrder_bp,'yes')
%     cfg_int.bpfiltord     = FilterOrder_bp;
% end
% cfg_int.bpfreq        = [FpassLeft FpassRight];%dummy values are overwritten by low level function
% cfg_int.feedback = core_cfg.feedback;
% fprintf('reprocess and apply band filter\n');
% data = st_preprocessing(cfg_int,data);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %     %ROI2
% %     epochLengthSamples = scoring.epochlength * fsample;
% %     [roiBegins, roiEnds] = getROIsByHypnogram(hypnogramPath,epochLengthSamples,sleepStagesOfInterest);
% %
% %     if (signalOffsetSamples ~= 0)
% %         signalOffsetSamples_downsampled = floor(signalOffsetSeconds*fsample);
% %         roiBegins = roiBegins + signalOffsetSamples_downsampled;
% %         roiEnds = roiEnds + signalOffsetSamples_downsampled;
% %     end
% %
% %     if strcmp(IgnoreDataSetHeader,'no')
% %         roiBegins = roiBegins(1:indexLastIncludedROIinData);
% %         roiEnds = roiEnds(1:indexLastIncludedROIinData);
% %         nSampleLength = floor(nSampleLength*(FrqOfSmplWishedPar/preDownsampleFreq));
% %         if (roiEnds(end) > nSampleLength)
% %             roiEnds(indexLastIncludedROIinData) = nSampleLength;
% %         end;
% %     end
% %
% %
% %     trlSampleBeginsAndEnds = [roiBegins roiEnds];
% 
% sampleLengthsAcrossROIs = sum(trlSampleLengths);
% lengthsAcrossROIsSeconds = sampleLengthsAcrossROIs/fsample; % in seconds
% 
% 
% smplsMinDetectionLength  = round(cfg.minduration*fsample);
% smplsMaxDetectionLength  = round(cfg.maxduration*fsample);
% 
% minFreqPostFreqBorderBufferLength = round(fsample/minFreq);
% 
% %smplsRMSTimeWndw = round(RMSTimeWndw*fsample);
% 
% smplsRMSTimeWndw = 0;
% %assure odd samplesize for moving average time window
% if (mod(floor(cfg.rmstimewndw*fsample),2) == 0)
%     smplsRMSTimeWndw = floor(cfg.rmstimewndw*fsample)+1;
% else
%     smplsRMSTimeWndw = floor(cfg.rmstimewndw*fsample);
% end
% 
% 
% smplsMovAvgTimeWndw = 0;
% %assure odd samplesize for moving average time window
% if (mod(floor(cfg.movavgtimewndw*fsample),2) == 0)
%     smplsMovAvgTimeWndw = floor(cfg.movavgtimewndw*fsample)+1;
% else
%     smplsMovAvgTimeWndw = floor(cfg.movavgtimewndw*fsample);
% end
% 
% 
% nChannels = length(data.label);
% 
% SDfrqBndPssSignal = [];
% if ~strcmp(UseAbsoluteEnvelopeThreshold,'yes')
%     
%     tempAllDataValues = [];
%     
%     for iTr = 1:size(data.trial,2)
%         if iTr == 1
%             tempAllDataValues = data.trial{iTr};
%         else
%             tempAllDataValues = cat(2,tempAllDataValues,data.trial{iTr});
%         end
%     end;
%     
%     if thresholdingInDifferentStages && hasROIs_threshold
%         cfg_int = [];
%         cfg_int.stages   = cfg.thresholdstages;
%         cfg_int.scoring  = scoring;
%         [begins, ends]   = st_select_scoring(cfg_int, data);
%         if isempty(begins)
%             hasROIs_threshold = false;
%             ft_warning('cfg.thresholdstages not present in the data! This will result in no detection')
%         end
%         
%         if ~hasROIs_threshold
%             for iT = 1:numel(data.trial)
%                 data.trial{iT}(:) = NaN;
%             end
%         end
%         
%         tempAllDataValuesTimes = [];
%         for iTr = 1:size(data.trial,2)
%                 tempAllDataValuesTimes = cat(2,tempAllDataValuesTimes,data.time{iTr});
%         end
%         
%         indInStage = repmat(false,1,size(tempAllDataValues,2));
%         for iBegins = 1:numel(begins)
%             indInStage = indInStage | ((begins(iBegins) >= tempAllDataValuesTimes) & (tempAllDataValuesTimes <= ends(iBegins)));
%         end
%         
%         tempAllDataValues = tempAllDataValues(:,indInStage);
%         indInStage = [];
%         tempAllDataValuesTimes = [];
%         
%     end
%     
%     if strcmp(cfg.thresholdsignal,'filtered_signal')
%     elseif strcmp(cfg.thresholdsignal,'envelope')
%         for iChan_temp = 1:size(tempAllDataValues,1)
%             iNaN = isnan(tempAllDataValues(iChan_temp,:));
%             tempAllDataValues(iChan_temp,iNaN) = 0;
%             tempAllDataValues(iChan_temp,:) = abs(hilbert(tempAllDataValues(iChan_temp,:)));
%             tempAllDataValues(iChan_temp,iNaN) = NaN;
%         end
%     end
%     
%     
%     if strcmp('meanoverchan',cfg.thresholdaggmeth)
%         if strcmp(cfg.thresholdformbase,'mean')
%             SDfrqBndPssSignal = nanmean(nanmean(tempAllDataValues,2));
%         elseif strcmp(cfg.thresholdformbase,'std')
%             SDfrqBndPssSignal = nanmean(nanstd(tempAllDataValues,0,2));
%         end
%     elseif strcmp('valuesoverchan',cfg.thresholdaggmeth)
%         if strcmp(cfg.thresholdformbase,'mean')
%             SDfrqBndPssSignal = nanmean(tempAllDataValues(:));
%         elseif strcmp(cfg.thresholdformbase,'std')
%             SDfrqBndPssSignal = nanstd(tempAllDataValues(:));
%         end
%     elseif strcmp('respectivechan',cfg.thresholdaggmeth)
%         if strcmp(cfg.thresholdformbase,'mean')
%             SDfrqBndPssSignal = nanmean(tempAllDataValues,2);
%         elseif strcmp(cfg.thresholdformbase,'std')
%             SDfrqBndPssSignal = nanstd(tempAllDataValues,0,2);
%         end
%         if (1 == numel(SDfrqBndPssSignal)) && any(isnan(SDfrqBndPssSignal))
%            SDfrqBndPssSignal = repmat(NaN,nChannels,1);   
%         end
%     else
%         error('cfg.thresholdaggmeth for calculating mean or standard deviation over channels is unknown');
%     end
%     tempAllDataValues = [];%clear
% end
% 
% 
% dataset_factorThresholdBeginEnd = cfg.factorthreshbeginend ;
% dataset_factorThresholdCriterion = factorThresholdCriterion;
% 
% 
% 
% 
% maxPeaksOrTroughsPerSpindle = ceil(maxFreq*cfg.maxduration+1);




if isfield(cfg, 'channel')
    cfg_rod = [];
    if isnumeric(cfg.channel)
        cfg_rod.order = cfg.channel(:);
    else
        [Lia, Lib ] = ismember(cfg.channel, data.label);
        cfg_rod.order = Lib(:);
    end
    data_ord = st_reorderdata(cfg_rod, data);
end

if numel(data.label) > 2
    ft_error('This funtion expects two channels the left and the right EOG channel that are inversely correlated')
end


hasOnlyOneChannel = false;
if numel(data.label) == 2
    channel_label_EOG_left = data.label{1};
    channel_label_EOG_right = data.label{2};
elseif numel(data.label) == 1
    ft_warning('ONLY ONCE CHANNEL PROVIDED. BUT This funtion expects two channels the left and the right EOG channel that are inversely correlated. Taking the inverse of the first channel as the second one.')
    channel_label_EOG_left = data.label{1};
    channel_label_EOG_right = data.label{1};
    hasOnlyOneChannel = true;
end


epochs = cellfun(@sleepStage2str,scoring.epochs','UniformOutput',0);
hypnStages = [cellfun(@sleepStage2str,epochs,'UniformOutput',0) ...
    cellfun(@sleepStage2str_alt,epochs,'UniformOutput',0) ...
    cellfun(@sleepStage2str_alt2,epochs,'UniformOutput',0)];

hypnEpochsEndsSamples = cumsum(repmat(scoring.epochlength,numel(scoring.epochs),1));
hypnEpochsBeginsSamples = hypnEpochsEndsSamples - scoring.epochlength;

% hypnEpochsEndsSamples = hypnEpochsEndsSamples + scoring.dataoffset;
% hypnEpochsBeginsSamples = hypnEpochsBeginsSamples + scoring.dataoffset;

hypnEpochsEndsSamples = hypnEpochsEndsSamples*fsample;
hypnEpochsBeginsSamples = hypnEpochsBeginsSamples*fsample+1;

%%% Alorythm by Marek_Adamczyk adapted by Frederik D. Weber %%%
    
    
    usedFixedSamplingRate = fsample; %% must be 100 Hz for the filters to work out.

    %usedFixedSamplingRate = 100;
    
    cfg_scc = [];
    cfg_scc.to = 'number';

    %scoring.label = unique(scoring.epochs)';
    scoring_numbers = st_scoringconvert(cfg_scc,scoring);
    scoring.numbers = cellfun(@str2num,scoring_numbers.epochs,'UniformOutput',1);

    if istrue(cfg.ConsiderExcludedEpochs)
        scoring.numbers((scoring.excluded > 0) & (scoring.numbers == 5)) = 51;        
    end
    
    keep_EpochsBeginsSamples = hypnEpochsBeginsSamples((scoring.numbers == 5));
    keep_EpochsEndsSamples = hypnEpochsEndsSamples((scoring.numbers == 5));
    
    keep_samples = zeros(1,size(data.trial{1},2));
    for iKeep = 1:numel(keep_EpochsBeginsSamples)
        keep_samples(:,keep_EpochsBeginsSamples(iKeep):keep_EpochsEndsSamples(iKeep)) = 1;
    end
    
     
    if istrue(cfg.RemoveSignalThatIsNotREM)
        cut_EpochsBeginsSamples = hypnEpochsBeginsSamples(~(scoring.numbers == 5));
        cut_EpochsEndsSamples = hypnEpochsEndsSamples(~(scoring.numbers == 5));
        for iCut = 1:numel(cut_EpochsBeginsSamples)
            data.trial{1}(:,cut_EpochsBeginsSamples(iCut):cut_EpochsEndsSamples(iCut)) = 0;
        end
        data.trial{1}(:,hypnEpochsEndsSamples(end):end) = 0;
    end
    
    
    hasREM = true;
    if isempty(find(scoring.numbers == 5, 1))
        
        ft_warning(['No usable REM sleep'])
    
%         remDens = zeros(1,length(scoring.numbers));
%         meanREMSpeed = zeros(1,length(scoring.numbers));
%         meanREMWidth = zeros(1,length(scoring.numbers));
%         exactlyMarkedREMsVctSampling = usedFixedSamplingRate;
%         exactlyMarkedREMs = zeros(length(scoring.numbers)*epochLength*exactlyMarkedREMsVctSampling, 1);
        hasREM = false;   
    else
        
                
        % %           [hip, signal, REMsAddress, REMsData] = takeREMDataNotTwins(edfsPath, currFileName, false, false, 'nieee', sigInEDF);
        %
        %            %[hip, signal] = takeREMDataNotTwins(edfsPath, currFileName, false, false, 'nieee', sigInEDF);
        %            signal = takeEOGsig(edfFileAddress, sigInEDF);
        
        cfg_sd = [];
        cfg_sd.feedback = core_cfg.feedback;
        cfg_sd.channel = channel_label_EOG_left;
        data_EOG_left = ft_selectdata(cfg_sd,data);
        cfg_sd.channel = channel_label_EOG_right;
        data_EOG_right = ft_selectdata(cfg_sd,data);
        
        data = [];
        
        if hasOnlyOneChannel
            data_EOG_right.trial{1} = -1*data_EOG_right.trial{1};
        end
        

        
        
        fprintf('excluding artifacts\n');
        eogClass =    eogClinicClassifyFlexThr(data_EOG_left.trial{1}', data_EOG_right.trial{1}', fsample, scoring.numbers, false);%looks for artifacts!
        fprintf('detect REMs\n');
        [remDens, exactlyMarkedREMs, meanREMSpeed, meanREMWidth, EOG1, EOG2, EOG1_lp, EOG2_lp, EOG1_lphp, EOG2_lphp, EOG1_lp_dev, EOG2_lp_dev, EOG1_lphp_dev, EOG2_lphp_dev, eogNoise, thrBS, thrRS] ...
            = paperSciRemDetectorAdjustableThreshold(data_EOG_left.trial{1}', data_EOG_right.trial{1}', fsample, scoring.epochlength, eogClass, false, 100000000000, cfg.RelaxedThresholdFactorOfBasicThreshold, cfg.BasicThresholdFactorOfEOGnoise, cfg.amplitudemax);
        
        EOG1_lp_dev = [0; EOG1_lp_dev];
        EOG2_lp_dev = [0; EOG2_lp_dev];
        EOG1_lphp_dev = [0; EOG1_lphp_dev];
        EOG2_lphp_dev = [0; EOG2_lphp_dev];
        eogClass = [];
        data_EOG_right = [];
        data_EOG_left = [];
        
        EOG1 = [];EOG2 = []; EOG1_lp = []; EOG2_lp = []; EOG1_lphp = []; EOG2_lphp = [];

        
        
        [begins, ends] = getBeginsAndCorrespondingEndsIndicesAboveThreshold(exactlyMarkedREMs*2,1);
        
        
         
         begins = begins(1:numel(ends));
         
         if isempty(begins)
             indicesValid = [];
         else
             indicesValid = find(keep_samples(begins) == 1 | (keep_samples(ends) == 1)); %consider border effects filter
         end
         begins = begins(indicesValid);
         ends = ends(indicesValid);
        
         nDetected = length(ends);
        %         indicesValidSamples = find((begins >= firstSmplsMinFreqPostFreqBorderBufferLength) & (ends <= lastSmplsMinFreqPostFreqBorderBufferLength)); %consider border effects filter
        %
        %         begins = begins(indicesValidSamples);
        %         ends = ends(indicesValidSamples);
        %
        
        mean_detected_seconds_duration = NaN;
        mean_detected_REMSpeed_uV_per_second = NaN;
        mean_rightness_EOGRdeviation_minus_EOGLdeviation = NaN;
       
        detected = [];
        if (nDetected > 0)
            detectedLengthSamples = ends - begins + 1;
            
            %             indicesCandiates = find((tempCandidatesLengths >= smplsMinDetectionLength) & (tempCandidatesLengths <= smplsMaxDetectionLength));
            %
            %             nDetected = length(indicesCandiates);
            

            detected.id = (1:nDetected)';
            detected.seconds_begin = begins(1:nDetected)/usedFixedSamplingRate;
            detected.seconds_center = (begins(1:nDetected)+floor(detectedLengthSamples/2))/usedFixedSamplingRate;
            detected.seconds_end = ends(1:nDetected)/usedFixedSamplingRate;
            detected.seconds_duration = detectedLengthSamples(1:nDetected)/usedFixedSamplingRate;
            detected.REMSpeed_uV_per_second = zeros(nDetected,1);
            epochs = {};
            for iDetected = 1:nDetected
                detected.REMRightness(iDetected) = mean(EOG2_lp_dev(begins(iDetected):ends(iDetected)))-mean(EOG1_lp_dev(begins(iDetected):ends(iDetected)));                
                detected.REMSpeed_uV_per_second(iDetected) = mean([abs(mean(EOG1_lp_dev(begins(iDetected):ends(iDetected)))) abs(mean(EOG2_lp_dev(begins(iDetected):ends(iDetected))))])*fsample;                
                tempSample = begins(iDetected);
                tempInd = ((hypnEpochsBeginsSamples <= tempSample) & (tempSample <= hypnEpochsEndsSamples));
                if ~any(tempInd)
                    tempInd = ((hypnEpochsBeginsSamples <= tempSample+1) & (tempSample <= hypnEpochsEndsSamples));
                end
                if ~any(tempInd)
                    tempInd = ((hypnEpochsBeginsSamples <= tempSample) & (tempSample-1 <= hypnEpochsEndsSamples));
                end
                epochs(iDetected,:) = ...
                   [hypnStages(tempInd,1) ...
                    hypnStages(tempInd,2) ...
                    hypnStages(tempInd,3)];
            end
            
            detected.stage = epochs(:,1);
            detected.stage_alt = epochs(:,2);
            detected.stage_alt2 = epochs(:,3);
            epochs = [];
            
            detected.channel_label_EOG_left = cellstr(repmat(channel_label_EOG_left, nDetected, 1));
            detected.channel_label_EOG_right = cellstr(repmat(channel_label_EOG_right, nDetected, 1));
            detected.channel = cellstr(repmat([channel_label_EOG_left '+' channel_label_EOG_right], nDetected, 1));
            
            mean_detected_seconds_duration = mean(detected.seconds_duration);
            mean_detected_REMSpeed_uV_per_second = mean(detected.REMSpeed_uV_per_second);
            mean_rightness_EOGRdeviation_minus_EOGLdeviation = mean(detected.REMRightness);
            

        %else
        end
        
    end
    
    
    nEpochs_REM = sum(scoring.numbers == 5)+sum(scoring.numbers == 51);
    nEpochs_REM_withoutMA = sum(scoring.numbers == 5)-sum((scoring.excluded == 1) & (scoring.numbers == 5));
    
    
res_channel = [];
res_channel.ori = functionname;
res_channel.type = 'rems_channel';
res_channel.cfg = cfg;

tempvarnames = {'channel',...
            'count','density_per_REM_epoch','density_per_REM_without_MA_epoch',...
            'mean_duration_seconds','mean_speed_uV_per_second', 'mean_rightness_EOGRdeviation_minus_EOGLdeviation',...
            'used_basic_threshold_factor_of_EOG_noise','used_relaxed_threshold_factor_of_basic_threshold','observed_EOG_noise','observed_basic_threshold_of_EOG_noise','observed_relaxed_threshold_of_basic_threshold'};
if ~hasREM
    res_channel.table = cell2table(cell(0,numel(tempvarnames)), 'VariableNames', tempvarnames);
%     res_channel.table = table(...
%         [channel_label_EOG_left '+' channel_label_EOG_right],...
%         NaN,NaN,NaN,...
%         NaN,NaN,...
%         cfg.BasicThresholdFactorOfEOGnoise,cfg.RelaxedThresholdFactorOfBasicThreshold,...
%         'VariableNames',tempvarnames);
else
    res_channel.table = table(...
        {[channel_label_EOG_left '+' channel_label_EOG_right]},...
        nDetected,nDetected/nEpochs_REM,nDetected/nEpochs_REM_withoutMA,...
        mean_detected_seconds_duration,mean_detected_REMSpeed_uV_per_second,mean_rightness_EOGRdeviation_minus_EOGLdeviation,...
        cfg.BasicThresholdFactorOfEOGnoise,cfg.RelaxedThresholdFactorOfBasicThreshold, eogNoise, thrBS, thrRS,...
        'VariableNames',tempvarnames);
end
    

res_event = [];
res_event.ori = functionname;
res_event.type = 'rems_event';
res_event.cfg = cfg;
tempvarnames = {'channel',...
            'duration_seconds',...
            'speed_uV_per_second_seconds',...
            'rightness_EOGRdeviation_minus_EOGLdeviation',...
            'used_stages_for_detection',...
            'seconds_begin','seconds_end','seconds_center',...
            'id_within_channel',...
            'stage','stage_alt','stage_alt2'};
        
if ~hasREM
    res_event.table = cell2table(cell(0,numel(tempvarnames)), 'VariableNames', tempvarnames);
else
    if isempty(detected)
        res_event.table = cell2table(cell(0,numel(tempvarnames)), 'VariableNames', tempvarnames);
    else  
        res_event.table = table(cellstr(repmat([channel_label_EOG_left '+' channel_label_EOG_right],numel(detected.seconds_duration),1)),...
            detected.seconds_duration(:), ...
            detected.REMSpeed_uV_per_second(:),...
            detected.REMRightness(:),...
            cellstr(repmat('REM',numel(detected.seconds_duration),1)),...
            detected.seconds_begin(:)+scoring.dataoffset,detected.seconds_end(:)+scoring.dataoffset,detected.seconds_center(:)+scoring.dataoffset,...
            detected.id(:),...
            detected.stage(:),detected.stage_alt(:),detected.stage_alt2(:),...
            'VariableNames',tempvarnames);
    end
end


    
    
    % change computed REM density in matlab format into easy readable csv
    % format
        
    meanREMSpeed_as_microVolt_per_second_per_epoch = meanREMSpeed./(meanREMWidth./usedFixedSamplingRate);
    meanREMWidth_as_duration_seconds_per_epoch = meanREMWidth./usedFixedSamplingRate;
    
    %path_file_exactlyMarkedREMs = [pathOutputFolder filesep ouputFilesPrefixString 'REMsMaAd_exactlyMarkedREMs_' 'datanum_' num2str(iData) '.mat'];
    %path_file_standardREMdens = [pathOutputFolder filesep ouputFilesPrefixString 'REMsMaAd_standardREMdens_' 'datanum_' num2str(iData) '.mat'];
    
   % save(path_file_exactlyMarkedREMs,'exactlyMarkedREMs');
   % save(path_file_standardREMdens,'remDens');
    
    MarkedREMsSampling = usedFixedSamplingRate;
    
    [standardREMdensity1stCycle funName REMepochsNb1stCycle] = ...
        firstCycleREM_densFromExactVct(scoring.numbers, scoring.epochlength, exactlyMarkedREMs, MarkedREMsSampling);

    REM3Count1stCycle   = firstCycleMiniEp3WithREMsFromExactVct(scoring.numbers, scoring.epochlength, exactlyMarkedREMs, MarkedREMsSampling);
    REM1Count1stCycle =   firstCycleMiniEp1WithREMsFromExactVct(scoring.numbers, scoring.epochlength, exactlyMarkedREMs, MarkedREMsSampling);
    AllREMsCount1stCycle =  firstCycleAllREMsFromExactVct(scoring.numbers, scoring.epochlength, exactlyMarkedREMs, MarkedREMsSampling);
    
    allREMepochsNb =         allNightREMepochs(scoring.numbers, scoring.epochlength, exactlyMarkedREMs, MarkedREMsSampling);
    REMlatency =                giveREMlatency(scoring.numbers, scoring.epochlength, exactlyMarkedREMs, MarkedREMsSampling);
    
    standardREMdensityAllNight = allNightREM_densFromExactVct(scoring.numbers, scoring.epochlength, exactlyMarkedREMs, MarkedREMsSampling);
    REM3Count           = allNightMiniEp3WithREMsFromExactVct(scoring.numbers, scoring.epochlength, exactlyMarkedREMs, MarkedREMsSampling);
    REM1Count   =         allNightMiniEp1WithREMsFromExactVct(scoring.numbers, scoring.epochlength, exactlyMarkedREMs, MarkedREMsSampling);
    AllREMsCount =                allNightAllREMsFromExactVct(scoring.numbers, scoring.epochlength, exactlyMarkedREMs, MarkedREMsSampling);
    
    REMsCountInBurst =      allNightREMsInMidBurstsFromExactVct(scoring.numbers, scoring.epochlength, exactlyMarkedREMs, MarkedREMsSampling);
    burstsNumber =     allNightAllREMsNbOfMidBurstsFromExactVct(scoring.numbers, scoring.epochlength, exactlyMarkedREMs, MarkedREMsSampling);
    REMsInBurstsPrc =    allNightREMsPrcInMidBurstsFromExactVct(scoring.numbers, scoring.epochlength, exactlyMarkedREMs, MarkedREMsSampling);
    
    REM1CountInBurst =       allNightREM1InMidBurstsFromExactVct(scoring.numbers, scoring.epochlength, exactlyMarkedREMs, MarkedREMsSampling);
    REM1burstsNumber =     allNightREM1NbOfMidBurstsFromExactVct(scoring.numbers, scoring.epochlength, exactlyMarkedREMs, MarkedREMsSampling);
    REM1InBurstsPrc =     allNightREM1PrcInMidBurstsFromExactVct(scoring.numbers, scoring.epochlength, exactlyMarkedREMs, MarkedREMsSampling);
    
    
%     fid = fopen([pathOutputFolder filesep ouputFilesPrefixString 'REMsMaAd_' 'datanum_' num2str(iData) '.txt'],'wt');
%     
%     fprintf(fid, '---------- First cycle data ---------\n\n');
%     fprintf(fid, 'REM epochs number: %d\n', REMepochsNb1stCycle);
%     fprintf(fid, '\n++++++ REM density +++++++\n\n');
%     
%     fprintf(fid, '%2.1f sec miniep: %2.2f\n', scoring.epochlength/10, standardREMdensity1stCycle);
%     fprintf(fid, '1 sec miniep: %2.2f\n', REM1Count1stCycle/REMepochsNb1stCycle);
%     fprintf(fid, 'all REMs    : %2.2f\n', AllREMsCount1stCycle/REMepochsNb1stCycle);
%     fprintf(fid, '\n++++++ REM activity ++++++\n\n');
%     fprintf(fid, '%2.1f sec miniep: %d\n', scoring.epochlength/10, REM3Count1stCycle);
%     fprintf(fid, '1 sec miniep: %d\n', REM1Count1stCycle);
%     fprintf(fid, 'all REMs    : %d\n', AllREMsCount1stCycle);
%     
%     fprintf(fid, '\n---------- All night data ---------\n\n');
%     fprintf(fid, 'REM latency epochs number: %d\n', REMlatency);
%     fprintf(fid, 'REM epochs number: %d\n', allREMepochsNb);
%     fprintf(fid, '\n++++++ REM density +++++++\n\n');
%     fprintf(fid, '%2.1f sec miniep: %2.2f\n', scoring.epochlength/10, standardREMdensityAllNight);
%     fprintf(fid, '1 sec miniep: %2.2f\n', REM1Count/allREMepochsNb);
%     fprintf(fid, 'all REMs    : %2.2f\n', AllREMsCount/allREMepochsNb);
%     fprintf(fid, '\n++++++ REM activity ++++++\n\n');
%     fprintf(fid, '%2.1f sec miniep: %d\n', scoring.epochlength/10, REM3Count);
%     fprintf(fid, '1 sec miniep: %d\n', REM1Count);
%     fprintf(fid, 'all REMs    : %d\n', AllREMsCount);
%     fprintf(fid, '\n++++++  REM bursts  ++++++\n\n');
%     fprintf(fid, '>>> all REMs <<<\n');
%     fprintf(fid, 'activity: %d\n', REMsCountInBurst);
%     fprintf(fid, 'density : %2.2f\n', REMsCountInBurst/allREMepochsNb);
%     fprintf(fid, 'REM bursts number: %d\n', burstsNumber);
%     fprintf(fid, 'Percentage of REMs in bursts: %2.2f%%\n', REMsInBurstsPrc);
%     fprintf(fid, 'Nb of REMs in average burst : %2.2f\n', REMsCountInBurst/burstsNumber);
%     
%     fprintf(fid, '\n>>> 1 sec miniep <<<\n');
%     fprintf(fid, 'activity: %d\n', REM1CountInBurst);
%     fprintf(fid, 'density : %2.2f\n', REM1CountInBurst/allREMepochsNb);
%     fprintf(fid, 'REM bursts number: %d\n', REM1burstsNumber);
%     fprintf(fid, 'Percentage of REMs in bursts: %2.2f%%\n', REM1InBurstsPrc);
%     fprintf(fid, 'Nb of REMs in average burst : %2.2f\n', REM1CountInBurst/REM1burstsNumber);
%     
%     fclose(fid);
    
res_summary = [];
res_summary.ori = functionname;
res_summary.type = 'rems_summary';
res_summary.cfg = cfg;
tempvarnames = {...
        'AllREMsCount',...
        'REMsInBurstsPrc','avgREMsInBurst',...
        'REMactivityInBurstAllNight','REMdensityInBurstAllNight','burstsNumberAllNight', ...
        'REMepochsNb1stCycle','standardREMdensity1stCycle','REMdensity1stCycle1sMiniEpisode','allREMdensity1stCycle', ...
        'standardREMactivity1stCycle','REMactivity1stCycle1sMiniEpisode','allREMactivity1stCycle', ...
        'REMlatencyEpNb','REMepochsNb','standardREMdensityAllNight','REMdensityOf1sMiniEpisodeFromAllNight','allREMdensityAllNight', ...
        'allREMactivity3sMiniEpisode','allREMactivity1sMiniEpisode', 'REMactivityInBurst1sMiniEpisode','REMdensityInBurst1sMiniEpisode','burstsNumber1sMiniEpisode','burstsPrc1sMiniEpisode','NbOfREMsinAVGburst1sMiniEpisode'};
        
if ~hasREM
    res_summary.table = cell2table(cell(0,numel(tempvarnames)), 'VariableNames', tempvarnames);
else
    if isempty(detected)
        res_summary.table = cell2table(cell(0,numel(tempvarnames)), 'VariableNames', tempvarnames);
    else  
        res_summary.table = table(...
            AllREMsCount,...
            REMsInBurstsPrc, REMsCountInBurst/burstsNumber,...
            REMsCountInBurst, REMsCountInBurst/allREMepochsNb, burstsNumber, ...
            REMepochsNb1stCycle, standardREMdensity1stCycle, REM1Count1stCycle/REMepochsNb1stCycle,AllREMsCount1stCycle/REMepochsNb1stCycle,...
            REM3Count1stCycle, REM1Count1stCycle, AllREMsCount1stCycle, ...
            REMlatency, allREMepochsNb, standardREMdensityAllNight, REM1Count/allREMepochsNb, AllREMsCount/allREMepochsNb, ...
            REM3Count,REM1Count,REM1CountInBurst,REM1CountInBurst/allREMepochsNb,REM1burstsNumber,REM1InBurstsPrc,REM1CountInBurst/REM1burstsNumber,...
            'VariableNames',tempvarnames);
    end
end
    
    
    
    
    
    %%% END Alorythm by Marek_Adamczyk adapted by Frederik D. Weber %%%


fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end