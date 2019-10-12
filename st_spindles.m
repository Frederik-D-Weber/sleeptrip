function [res_channel, res_event, res_filter] = st_spindles(cfg, data)

% ST_SPINDLES detect sleep spindles and their properties like their count, density, amplitude, duration, frequency etc.
% results are stored it in a result structure
%
% Use as
%   [res_channel, res_event, res_filter] = st_spindles(cfg)
%   [res_channel, res_event, res_filter] = st_spindles(cfg, data)
%
% Required configuration parameters are:
%   cfg.scoring         = structure provided by ST_READ_SCORING
%   cfg.centerfrequency = number of the center frequency of the spindles in Hz,
%                         see cfg.leftofcenterfreq and cfg.rightofcenterfreq
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
%  cfg.channel            = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION
%  cfg.stages             = either a string or a Nx1 cell-array with strings that indicate sleep stages
%                            if no data structure is provided as input next to configuration also
%                            provide this parameter, possible stages are a subset of
%                            {'W', 'N1', 'N2', 'N3', 'S4', 'R'}
%                            (default = {'N2', 'N3', 'S4'} as the stages in which spindles typically appear)
%  cfg.thresholdstages    = either a string or a Nx1 cell-array with
%                            strings that indicate sleep stages to detemine the
%                            thresholds on. This must be a subset of cfg.stages 
%                            (default = cfg.stages)
%  cfg.leftofcenterfreq   = scalar, left of of the center frequency of the spindles in Hz,
%                            see cfg.rightofcenterfreq (default = 1)
%  cfg.rightofcenterfreq  = scalar, right of of the center frequency of the spindles in Hz,
%                            see cfg.leftofcenterfreq (default = 1)
%  cfg.thresholdaggmeth   = Method for aggregating the threshold of the channels
%                            e.g (where to take the standard deviation or mean from (see parameter cfg.thresholdformbase)
%                            either 'meanoverchan' or 'valuesoverchan' or 'respectivechan'
%                            i.e.
%                                 respectivechan uses SD of each channel for detection in that channel
%                                 meanoverchan uses mean of SD over all channels for each channel
%                                 valuesoverchan uses SD over all values in all channels for each channel
%                             (default = 'respectivechan')
%  cfg.envelopemeth         = the method for creating the envelope either
%                                'hilbertEnv' for the absolute value of the hilbert tranform of the signal or
%                                'smoothedRMSwd' for a smoothed moving symetric root mean squared of window of length cfg.rmstimewndw further symetrically smoothed with window length cfg.movavgtimewndw.
%                             (default = 'smoothedRMSwd')
%  cfg.factorthreshbeginend = factor to multiply the threshold basis gained by paramter decribed in cfg.thresholdformbase
%                              (e.g. standard deviations or means of signal or envelope)
%                              for threshold of the signal for envelope for begin and end of event
%                              (default = 1.5)
%  cfg.factorthresholdcriterion = factor to multiply the threshold basis gained by paramter decribed in cfg.thresholdformbase
%                              (e.g. standard deviations or means of signal or envelope)
%                              for threshold of the signal for envelope for for minimum criterion to reach to count as an event.
%                              must be greater or equal than cfg.factorthreshbeginend
%                              (default =
%                                1.5  for cfg.envelopemeth = 'smoothedRMSwd'
%                                2.25 for cfg.envelopemeth = 'hilbertEnv')
%  cfg.minduration          = Minimal duration within detection criteria in seconds i.e. RMS above threshold (default = 0.5)
%  cfg.maxduration          = Maximal duration within detection criteria in seconds (default = 2)
%  cfg.mergewithin          = proximity of boundaries of events in seconds in order to be merged within margins of cfg.minduration and cfg.maxduration.
%                             with 0 meaning no merging and all further selection filters like cfg.minamplitude and cfg.maxamplitude applied after merging
%                             (default = 0)
%  cfg.minamplitude         = minimum absolute potential difference (i.e.
%                              amplitude) to select as valid event (i.e. peak to peak or peak to trough)
%                              (default 0)
%  cfg.maxamplitude         = maximum absolute potential difference (i.e.
%                               amplitude) to select as valid event (i.e.peak to peak or peak to trough)
%                              (default 120)
%  cfg.filterSDamp          = number of standard deviations to fiter for
%                             amplitude within the range of +- those standard deviations around the
%                             mean e.g. +-5 SD. This is applied after cfg.minamplitude and cfg.maxamplitude criterions applied already (default 5)
%  cfg.filterSDdur          = number of standard deviations to fiter for
%                             duration within the range of +- those standard deviations around the
%                             mean e.g. +-5 SD. This is applied after cfg.minamplitude and cfg.maxamplitude criterions applied already (default 5)
%  cfg.filterSDfreq         = number of standard deviations to fiter for
%                             core frequency within the range of +- those standard deviations around the
%                             mean e.g. +-5 SD. This is applied after cfg.minamplitude and cfg.maxamplitude criterions applied already (default 5)
%  cfg.thresholdformbase    = gives the method in which the threshold is formed either 'mean' or 'std' (default = 'std')
%  cfg.thresholdsignal      = what to use for threshold generation, either
%                             'filtered_signal' or 'envelope' (default = 'filtered_signal')
%  cfg.rmstimewndw          = time window for RMS in seconds (default = 0.2)
%  cfg.movavgtimewndw       = time window for moving average in smoothing RMS in seconds (default = 0.2)
%  cfg.envelopethresholds   = overwrite the envelope begin-end and min
%                             threshold to those, e.g.
%                             cfg.envelopethresholds = [4, 6];
%                             then cfg.factorthresholdbeginend and
%                             cfg.factorthresholdcriterion do not apply to those
%                             thresholds and are ignored.
%  cfg.downsamplefs       = downsample the data to this frequency in Hz before doing the anlysis (default = 100/128)
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

tic
memtic
st = dbstack;
functionname = st.name;

fprintf([functionname ' function started\n']);

if ~isfield(cfg, 'scoring')
    cfg.scoring = st_read_scoring(cfg);
end

downsamplefsNotSet = false;
if ~isfield(cfg, 'downsamplefs')
    downsamplefsNotSet = true;
end

% set defaults
cfg.channel  = ft_getopt(cfg, 'channel', 'all', 1);
cfg.stages  = ft_getopt(cfg, 'stages', {'N2', 'N3', 'S4'});
cfg.thresholdstages  = ft_getopt(cfg, 'thresholdstages', cfg.stages);
cfg.leftofcenterfreq  = ft_getopt(cfg, 'leftofcenterfreq', 1);
cfg.rightofcenterfreq  = ft_getopt(cfg, 'rightofcenterfreq', 1);
cfg.thresholdaggmeth  = ft_getopt(cfg, 'thresholdaggmeth', 'respectivechan');
cfg.envelopemeth  = ft_getopt(cfg, 'envelopemeth','smoothedRMSwd');
cfg.factorthreshbeginend  = ft_getopt(cfg, 'factorthreshbeginend', 1.5);
cfg.minduration  = ft_getopt(cfg, 'minduration', 0.5);
cfg.maxduration  = ft_getopt(cfg, 'maxduration', 2.0);
cfg.mergewithin  = ft_getopt(cfg, 'mergewithin', 0);
cfg.minamplitude  = ft_getopt(cfg, 'minamplitude', 0);
cfg.maxamplitude  = ft_getopt(cfg, 'maxamplitude', 120);
cfg.filterSDamp  = ft_getopt(cfg, 'filterSDamp', 5);
cfg.filterSDdur  = ft_getopt(cfg, 'filterSDdur', 5);
cfg.filterSDfreq  = ft_getopt(cfg, 'filterSDfreq', 5);
cfg.downsamplefs     = ft_getopt(cfg, 'downsamplefs', 100);
cfg.thresholdformbase  = ft_getopt(cfg, 'thresholdformbase', 'std');
cfg.thresholdsignal  = ft_getopt(cfg, 'thresholdsignal', 'filtered_signal');
cfg.rmstimewndw  = ft_getopt(cfg, 'rmstimewndw', 0.2);
cfg.movavgtimewndw  = ft_getopt(cfg, 'movavgtimewndw', 0.2);

if ~isfield(cfg, 'factorthresholdcriterion')
    switch cfg.envelopemeth
        case 'smoothedRMSwd'
            cfg.factorthresholdcriterion = 1.5;
        case 'hilbertEnv'
            cfg.factorthresholdcriterion = 2.25;
        otherwise
            cfg.factorthresholdcriterion = cfg.factorthreshbeginend;
    end
end

if cfg.factorthreshbeginend > cfg.factorthresholdcriterion
    ft_error('cfg.factorthreshbeginend cannot be bigger than cfg.factorthresholdcriterion')
end

UseAbsoluteEnvelopeThreshold = 'no';
if isfield(cfg, 'envelopethresholds')
    UseAbsoluteEnvelopeThreshold = 'yes'; % If abosolute positive envelope threshold for all channels should be used. In this case the factorSDbeginEnd and factorSDcriterion is not applied. either yes or no default no
    AbsoluteEnvelopeThresholdBeginEnd = cfg.envelopethresholds(1); % Abosolute positive envelope potential threshold for all channels for begin and end of event default 4 (for microVolts potential)
    AbsoluteEnvelopeThresholdCriterion = cfg.envelopethresholds(2); % Abosolute positive envelope potential threshold for all channels for minimum criterion to reach to count as an event. must be greater or equal than AbsoluteEnvelopeThresholdBeginEnd default 4 (for microVolts potential)
    if AbsoluteEnvelopeThresholdBeginEnd > AbsoluteEnvelopeThresholdCriterion
        error(['cfg.envelopethresholds must be a vector with two numbers, and the first cannot be greater than the second'])
    end
end


if ~all(ismember(cfg.thresholdstages, cfg.stages))
	ft_error(['cfg.thresholdstages must be a subset of cfg.stages'])
end
thresholdingInDifferentStages = ~all(strcmp(cfg.stages,cfg.thresholdstages));

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

preCenterFreqFilterTo_FpassLeft = cfg.leftofcenterfreq;
postCenterFreqFilterTo_FpassRight = cfg.rightofcenterfreq;  % in Hz



if ~(strcmp(cfg.envelopemeth,'hilbertEnv') || strcmp(cfg.envelopemeth,'smoothedRMSwd'))
    error(['cfg.envelopemeth ' cfg.envelopemeth ' in parameters is unknown, use either hilbertEnv or smoothedRMSwd'])
end

if ~(strcmp(cfg.thresholdformbase,'mean') || strcmp(cfg.thresholdformbase,'std'))
    error(['cfg.thresholdformbase ' cfg.thresholdformbase ' in parameters is unknown, use either mean or std, default std'])
end

if ~(strcmp(cfg.thresholdsignal,'filtered_signal') || strcmp(cfg.thresholdsignal,'envelope'))
    error(['cfg.thresholdsignal ' cfg.thresholdsignal ' in parameters is unknown, use either filtered_signal or envelope, default filtered_signal'])
end

factorThresholdCriterion = cfg.factorthresholdcriterion;

if cfg.scoring.epochlength < (1/(cfg.centerfrequency - cfg.leftofcenterfreq))
    error(['the epoch length ' num2str(cfg.scoring.epochlength) 's must not be greater in order to support the requrested maximum of ' num2str((cfg.centerfrequency - cfg.leftofcenterfreq)) ' Hz!'])
end


if strcmp(UseAbsoluteEnvelopeThreshold,'yes')
    cfg.factorthreshbeginend  = 1;
    AbsoluteEnvelopeThresholdCriterionRatio = AbsoluteEnvelopeThresholdCriterion/AbsoluteEnvelopeThresholdBeginEnd;
    factorThresholdCriterion = AbsoluteEnvelopeThresholdCriterionRatio;
end




% set core parameters
load_core_cfg
% core_cfg

if (strcmp(core_cfg.bpfilttype,'FIRdesigned') || strcmp(core_cfg.bpfilttype,'IIRdesigned')) && (~strcmp(core_cfg.bpinstabilityfix,'no'))
    error(['band pass filter instability fix not supported for FIRdesigned or IIRdesigned'])
end

if (strcmp(core_cfg.hpfilttype,'FIRdesigned') || strcmp(core_cfg.hpfilttype,'IIRdesigned')) && (~strcmp(core_cfg.hpinstabilityfix,'no'))
    error(['high pass filter instability fix not supported for FIRdesigned or IIRdesigned'])
end

Apass_bp = Apass;
AstopLeft_bp = AstopLeft;
AstopRight_bp = AstopRight;

Apass_hp = Apass;
AstopLeft_hp = AstopLeft;

PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff = 0.5;



%filtfilt -- The length of the input x must be more than three times the filter order (N) defined as max(length(b)-1,length(a)-1).
%For best results, make sure the sequence you are filtering has length at least three times the filter order and tapers to zero on both edges.i.e.:
minSignalLengthSamples = cfg.scoring.epochlength*cfg.downsamplefs;
maxFilterOrder = floor((minSignalLengthSamples) / 3) - 1;

if strcmp(UseFixedFilterOrder_bp,'yes') && strcmp(useTwoPassFiltering_bp,'yes') && (FilterOrder_bp > maxFilterOrder)
    error(['filter order for band pass not optimal, use a maximal filter order three times of samples in one epoch, i.e. filter order of ' num2str(maxFilterOrder)])
elseif (FilterOrder_bp < maxFilterOrder) && strcmp(MaximizeFilterOrderIfFixedFilterOrderIsUsed,'yes')
    FilterOrder_bp = maxFilterOrder;
end

if strcmp(UseFixedFilterOrder_hp,'yes') && strcmp(useTwoPassFiltering_hp,'yes') && (FilterOrder_hp > maxFilterOrder)
    error(['filter order for high pass not optimal, use a maximal filter order three times of samples in one epoch, i.e. filter order of ' num2str(maxFilterOrder)])
elseif (FilterOrder_hp < maxFilterOrder)  && strcmp(MaximizeFilterOrderIfFixedFilterOrderIsUsed,'yes')
    FilterOrder_hp = maxFilterOrder;
end

if strcmp(UseFixedFilterOrder_bp,'yes') && logical(mod(FilterOrder_bp,2))
    error('band pass filter order must be an even number')
end

if strcmp(UseFixedFilterOrder_hp,'yes') && logical(mod(FilterOrder_hp,2))
    error('high pass order must be an even number')
end


% Note that a one- or two-pass filter has consequences for the
% strength of the filter, i.e. a two-pass filter with the same filter
% order will attenuate the signal twice as strong.
if strcmp(useTwoPassFiltering_bp,'yes') && strcmp(UseTwoPassAttenuationCorrection_bp,'yes')
    Apass_bp = Apass_bp/2;
    AstopLeft_bp = AstopLeft_bp/2;
    AstopRight_bp = AstopRight_bp/2;
end

if strcmp(useTwoPassFiltering_hp,'yes') && strcmp(UseTwoPassAttenuationCorrection_hp,'yes')
    Apass_hp = Apass_hp/2;
    AstopLeft_hp = AstopLeft_hp/2;
end

if strcmp(UseFixedFilterOrder_hp,'yes') && strcmp(core_cfg.hpfilttype,'FIRdesigned')
    error('UseFixedFilterOrder_hp not allowed for high pass filters of type FIRdesigned')
end



minFreq = cfg.centerfrequency - preCenterFreqFilterTo_FpassLeft;
maxFreq = cfg.centerfrequency + postCenterFreqFilterTo_FpassRight;

if ~(maxFreq*3 < fsample)
    error(['sample frequency of ' num2str(fsample) ' Hz must be MORE THAN three-fold (i.e. 3-fold) the maximal band frequency of ' num2str(maxFreq) ' Hz,\n even lower than requested by Nyquist-Shannon sample theorem.\n consider excludding higher frequency bands!'])
end

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
if isfield(cfg, 'montage'),     cfg_int.reref = cfg.montage;           end

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


cfg_int.stages   = cfg.stages;
cfg_int.scoring  = cfg.scoring;

hasROIs = true;
hasROIs_threshold = true;


if hasdata
    data_t = st_preprocessing(cfg_int, data);
    data_t = st_select_scoring(cfg_int, data_t);
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
    [cfg_int] = st_select_scoring(cfg_int);
    cfg_int.continuous   = 'yes';
    if isempty(cfg_int.trl)
        hasROIs = false;
        hasROIs_threshold = false;

        % read in dummy data
        cfg_int.trl = [1 round(fsample*60) 0];
        data = st_preprocessing(cfg_int);
        %         data.time = {};
        %         data.trial = {};
        data.sampleinfo = [0 -1];
    else
        data = st_preprocessing(cfg_int);
    end
    
end

if ~hasROIs
    for iT = 1:numel(data.trial)
        data.trial{iT}(:) = NaN;
    end
end

if (data.fsample == 128) && downsamplefsNotSet
    ft_warning('leaving 128 Hz sampling rate as default sampling rate')
    cfg.downsamplefs = 128;
end

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


use_hp = false;
use_lp = false;
use_bp = true;
%StopToPassTransitionWidth_hp_temp = StopToPassTransitionWidth_hp_predownsample;
adapt_filter_settings_to_toolbox

cfg_int = core_cfg;

cfg_int = [];
cfg_int = core_cfg;
cfg_int.bpfilter = 'yes';
FpassLeft = minFreq; %left pass frequency in Hz
FpassRight = maxFreq; %right pass frequency in Hz
FstopLeft = FpassLeft - StopToPassTransitionWidth_bp; %left stop frequency in Hz
FstopRight = FpassRight + PassToStopTransitionWidth_bp; %left stop frequency in Hz

usedFilterOrder_bp = NaN;
bp_hdm = [];
if strcmp(core_cfg.bpfilttype,'IIRdesigned') || strcmp(core_cfg.bpfilttype,'FIRdesigned')
    bp_d = [];
    bp_hd = [];
    fprintf('designing band pass filter\n');
    if strcmp(UseFixedFilterOrder_bp,'yes')
        bp_d = fdesign.bandpass('N,Fst1,Fp1,Fp2,Fst2',FilterOrder_bp,FstopLeft,FpassLeft,FpassRight,FstopRight,fsample);
        bp_hd = design(bp_d,'equiripple');
    else
        bp_d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',FstopLeft,FpassLeft,FpassRight,FstopRight,AstopLeft_bp,Apass_bp,AstopRight_bp,fsample);
        %bp_hd = design(bp_d,'equiripple','MinOrder', 'even');
        bp_hd = design(bp_d,'equiripple'); % order is even with 442 for 100 Hz 884 for 200 Hz
    end
    usedFilterOrder_bp = bp_hd.order;
    cfg_int.bpfilterdesign = bp_hd;
    bp_hdm = measure(bp_hd);
end
if strcmp(UseFixedFilterOrder_bp,'yes')
    cfg_int.bpfiltord     = FilterOrder_bp;
end
cfg_int.bpfreq        = [FpassLeft FpassRight];%dummy values are overwritten by low level function
cfg_int.feedback = core_cfg.feedback;
fprintf('reprocess and apply band filter\n');
data = st_preprocessing(cfg_int,data);










%     %ROI2
%     epochLengthSamples = cfg.scoring.epochlength * fsample;
%     [roiBegins, roiEnds] = getROIsByHypnogram(hypnogramPath,epochLengthSamples,sleepStagesOfInterest);
%
%     if (signalOffsetSamples ~= 0)
%         signalOffsetSamples_downsampled = floor(signalOffsetSeconds*fsample);
%         roiBegins = roiBegins + signalOffsetSamples_downsampled;
%         roiEnds = roiEnds + signalOffsetSamples_downsampled;
%     end
%
%     if strcmp(IgnoreDataSetHeader,'no')
%         roiBegins = roiBegins(1:indexLastIncludedROIinData);
%         roiEnds = roiEnds(1:indexLastIncludedROIinData);
%         nSampleLength = floor(nSampleLength*(FrqOfSmplWishedPar/preDownsampleFreq));
%         if (roiEnds(end) > nSampleLength)
%             roiEnds(indexLastIncludedROIinData) = nSampleLength;
%         end;
%     end
%
%
%     trlSampleBeginsAndEnds = [roiBegins roiEnds];

sampleLengthsAcrossROIs = sum(trlSampleLengths);
lengthsAcrossROIsSeconds = sampleLengthsAcrossROIs/fsample; % in seconds


smplsMinDetectionLength  = round(cfg.minduration*fsample);
smplsMaxDetectionLength  = round(cfg.maxduration*fsample);

minFreqPostFreqBorderBufferLength = round(fsample/minFreq);

%smplsRMSTimeWndw = round(RMSTimeWndw*fsample);

smplsRMSTimeWndw = 0;
%assure odd samplesize for moving average time window
if (mod(floor(cfg.rmstimewndw*fsample),2) == 0)
    smplsRMSTimeWndw = floor(cfg.rmstimewndw*fsample)+1;
else
    smplsRMSTimeWndw = floor(cfg.rmstimewndw*fsample);
end


smplsMovAvgTimeWndw = 0;
%assure odd samplesize for moving average time window
if (mod(floor(cfg.movavgtimewndw*fsample),2) == 0)
    smplsMovAvgTimeWndw = floor(cfg.movavgtimewndw*fsample)+1;
else
    smplsMovAvgTimeWndw = floor(cfg.movavgtimewndw*fsample);
end


nChannels = length(data.label);

SDfrqBndPssSignal = [];
if ~strcmp(UseAbsoluteEnvelopeThreshold,'yes')
    
    tempAllDataValues = [];
    
    for iTr = 1:size(data.trial,2)
        if iTr == 1
            tempAllDataValues = data.trial{iTr};
        else
            tempAllDataValues = cat(2,tempAllDataValues,data.trial{iTr});
        end
    end;
    
    if thresholdingInDifferentStages && hasROIs_threshold
        cfg_int = [];
        cfg_int.stages   = cfg.thresholdstages;
        cfg_int.scoring  = cfg.scoring;
        [begins, ends]   = st_select_scoring(cfg_int, data);
        if isempty(begins)
            hasROIs_threshold = false;
            ft_warning('cfg.thresholdstages not present in the data! This will result in no detection')
        end
        
        if ~hasROIs_threshold
            for iT = 1:numel(data.trial)
                data.trial{iT}(:) = NaN;
            end
        end
        
        tempAllDataValuesTimes = [];
        for iTr = 1:size(data.trial,2)
                tempAllDataValuesTimes = cat(2,tempAllDataValuesTimes,data.time{iTr});
        end
        
        indInStage = repmat(false,1,size(tempAllDataValues,2));
        for iBegins = 1:numel(begins)
            indInStage = indInStage | ((begins(iBegins) >= tempAllDataValuesTimes) & (tempAllDataValuesTimes <= ends(iBegins)));
        end
        
        tempAllDataValues = tempAllDataValues(:,indInStage);
        indInStage = [];
        tempAllDataValuesTimes = [];
        
    end
    
    if strcmp(cfg.thresholdsignal,'filtered_signal')
    elseif strcmp(cfg.thresholdsignal,'envelope')
        for iChan_temp = 1:size(tempAllDataValues,1)
            iNaN = isnan(tempAllDataValues(iChan_temp,:));
            tempAllDataValues(iChan_temp,iNaN) = 0;
            tempAllDataValues(iChan_temp,:) = abs(hilbert(tempAllDataValues(iChan_temp,:)));
            tempAllDataValues(iChan_temp,iNaN) = NaN;
        end
    end
    
    
    if strcmp('meanoverchan',cfg.thresholdaggmeth)
        if strcmp(cfg.thresholdformbase,'mean')
            SDfrqBndPssSignal = nanmean(nanmean(tempAllDataValues,2));
        elseif strcmp(cfg.thresholdformbase,'std')
            SDfrqBndPssSignal = nanmean(nanstd(tempAllDataValues,0,2));
        end
    elseif strcmp('valuesoverchan',cfg.thresholdaggmeth)
        if strcmp(cfg.thresholdformbase,'mean')
            SDfrqBndPssSignal = nanmean(tempAllDataValues(:));
        elseif strcmp(cfg.thresholdformbase,'std')
            SDfrqBndPssSignal = nanstd(tempAllDataValues(:));
        end
    elseif strcmp('respectivechan',cfg.thresholdaggmeth)
        if strcmp(cfg.thresholdformbase,'mean')
            SDfrqBndPssSignal = nanmean(tempAllDataValues,2);
        elseif strcmp(cfg.thresholdformbase,'std')
            SDfrqBndPssSignal = nanstd(tempAllDataValues,0,2);
        end
        if (1 == numel(SDfrqBndPssSignal)) && any(isnan(SDfrqBndPssSignal))
           SDfrqBndPssSignal = repmat(NaN,nChannels,1);   
        end
    else
        error('cfg.thresholdaggmeth for calculating mean or standard deviation over channels is unknown');
    end
    tempAllDataValues = [];%clear
end


dataset_factorThresholdBeginEnd = cfg.factorthreshbeginend ;
dataset_factorThresholdCriterion = factorThresholdCriterion;




maxPeaksOrTroughsPerSpindle = ceil(maxFreq*cfg.maxduration+1);


%cell(1,nChannels)
ch_detectedLengthSamples = [];
ch_detectedBeginSample = [];
ch_detectedEndSample = [];
ch_detectedPeak2Peaks = [];
ch_detectedPeaksSamples  = [];
ch_detectedTroughsSamples  = [];
ch_detectedSignalTroughsSamples = [];
ch_detectedSignalPeaksSamples = [];
ch_detectedSDofFilteredSignal = [];
ch_detectedMergeCount = [];
ch_nDetected = [];
ch_nMerged = [];
ch_densityPerEpoch = [];
ch_SDfrqBndPssSignal = [];
ch_detectedEnvelopeMaxs = [];
ch_detectedEnvelopeMaxSamples = [];

ch_detected_linear_regression_freq_slope = [];
ch_detected_linear_regression_freq_offset = [];
ch_detected_linear_regression_freq_R_squared = [];

ch_detected_inst_freq_troughs = [];
ch_detected_inst_freq_peaks = [];
ch_detectedTroughsPotential = [];
ch_detectedPeaksPotential = [];

ch_contigSegment = [];

%if hasROIs && hasROIs_threshold
for iChan = 1:nChannels
    %iChan = 1;
    
    fprintf('process channel %s ...\n',data.label{iChan});
    if strcmp(UseAbsoluteEnvelopeThreshold,'yes')
        ch_SDfrqBndPssSignal{iChan} = AbsoluteEnvelopeThresholdBeginEnd;
    else
        if strcmp('respectivechan',cfg.thresholdaggmeth)
            ch_SDfrqBndPssSignal{iChan} = SDfrqBndPssSignal(iChan);
        else
            ch_SDfrqBndPssSignal{iChan} = SDfrqBndPssSignal;
        end
    end
    
    %iChan = 1;
    cfg_int = [];
    cfg_int.feedback = core_cfg.feedback;
    cfg_int.channel = ft_channelselection(data.label{iChan}, data.label);
    chData = ft_selectdata(cfg_int,data);
    
    
    trl_detectedLengthSamples = [];
    trl_detectedBeginSample = [];
    trl_detectedEndSample = [];
    trl_detectedPeak2Peaks = [];
    trl_detectedPeaksSamples  = [];
    trl_detectedTroughsSamples  = [];
    trl_detectedSignalTroughsSamples = [];
    trl_detectedSignalPeaksSamples = [];
    trl_detectedSDofFilteredSignal = [];
    trl_detectedMergeCount = [];
    trl_nDetected = 0;
    trl_nMerged = 0;
    trl_detectedEnvelopeMaxs = [];
    trl_detectedEnvelopeMaxSamples = [];
    
    trl_detected_linear_regression_freq_slope = [];
    trl_detected_linear_regression_freq_offset = [];
    trl_detected_linear_regression_freq_R_squared = [];
    
    trl_detected_inst_freq_troughs = [];
    trl_detected_inst_freq_peaks = [];
    trl_detectedTroughsPotential = [];
    trl_detectedPeaksPotential = [];
    
    trl_contigSegment = [];
    
    for iTr = 1:size(chData.trial,2)
        fprintf('channel %s, subpart %i, preselect events in envelope\n',data.label{iChan},iTr);
        
        frqBndPssSignal = chData.trial{iTr};
        
        frqBndPssSignal_hilbert = frqBndPssSignal;
        iNaN = isnan(frqBndPssSignal_hilbert);
        frqBndPssSignal_hilbert(iNaN) = 0;
        frqBndPssSignal_hilbert = hilbert(frqBndPssSignal_hilbert);
        frqBndPssSignal_hilbert(iNaN) = NaN;
        
        thresholdForDetectionBeginEnd = ch_SDfrqBndPssSignal{iChan}*dataset_factorThresholdBeginEnd;
        thresholdForDetectionCriterion = ch_SDfrqBndPssSignal{iChan}*dataset_factorThresholdCriterion;
        
        lengthSignal = length(frqBndPssSignal);
        
        envelope = [];
        if strcmp(cfg.envelopemeth,'hilbertEnv')
            envelope = abs(frqBndPssSignal_hilbert)';
        elseif strcmp(cfg.envelopemeth,'smoothedRMSwd')
            envelope = smoothRMSwd(frqBndPssSignal,smplsRMSTimeWndw);
            if exist('smooth','file') == 2
                envelope = smooth(envelope,smplsMovAvgTimeWndw,'moving');
            else
                envelope = smoothwd(envelope,smplsMovAvgTimeWndw)';
            end
        end
        
        
        [begins, ends] = getBeginsAndCorrespondingEndsIndicesAboveThreshold(envelope,thresholdForDetectionBeginEnd);
        
        
        %indicesValidSamples = find((begins >= firstSmplsRMStimeWndw) & (ends <= lastSmplsRMStimeWndw)); %consider border effects of RMS
        firstSmplsMinFreqPostFreqBorderBufferLength = minFreqPostFreqBorderBufferLength;
        lastSmplsMinFreqPostFreqBorderBufferLength = lengthSignal - minFreqPostFreqBorderBufferLength;
        
        indicesValidSamples = [];
        if strcmp(cfg.envelopemeth,'hilbertEnv')
            indicesValidSamples = find((begins >= firstSmplsMinFreqPostFreqBorderBufferLength) & (ends <= lastSmplsMinFreqPostFreqBorderBufferLength)); %consider border effects of filter
        elseif strcmp(cfg.envelopemeth,'smoothedRMSwd')
            firstSmplsRMStimeWndw = (smplsRMSTimeWndw+1)/2 ;
            lastSmplsRMStimeWndw = lengthSignal - ((smplsRMSTimeWndw-1)/2);
            indicesValidSamples = find((begins >= firstSmplsRMStimeWndw) & (ends <= lastSmplsRMStimeWndw) & (begins >= firstSmplsMinFreqPostFreqBorderBufferLength) & (ends <= lastSmplsMinFreqPostFreqBorderBufferLength)); %consider border effects of RMS and filter
        end
        
        begins = begins(indicesValidSamples);
        ends = ends(indicesValidSamples);
        
        tempCandidatesLengths = ends - begins + 1;
        
        indicesCandiates = find((tempCandidatesLengths >= smplsMinDetectionLength) & (tempCandidatesLengths <= smplsMaxDetectionLength));
        
        smplsMergeEventsInProximityWithinDetectionMargins = round(cfg.mergewithin*fsample);
        
        nDetected = length(indicesCandiates);
        tempNmergedOuter = 0;
        detectedMergeCount = zeros(1,nDetected);
        
        if (nDetected > 0)
            %detectedLengthSamples = tempCandidatesLengths(indicesCandiates);
            detectedBeginSample = begins(indicesCandiates);
            detectedEndSample = ends(indicesCandiates);
            
            if smplsMergeEventsInProximityWithinDetectionMargins > 0
                fprintf('channel %s, subpart %i, merge events\n',data.label{iChan},iTr);
                %%%%%%%%
                tempNmergedInner = 1;
                
                while (nDetected > 1) && (tempNmergedInner > 0)
                    tempNmergedInner = 0;
                    
                    candidateDistanceSamples = detectedBeginSample(2:end) - detectedEndSample(1:(end-1));
                    
                    [sortCandidateDistanceSamples,preMergeDetectionIndices] = sort(candidateDistanceSamples);
                    
                    %preMergeDetectionIndices = uint64(preMergeDetectionIndices);
                    %preMergeDetectionIndices = 1:nDetected;
                    
                    tempMergeDetectionIndices = [];
                    
                    iIterCand = 1;
                    
                    
                    while (iIterCand <= length(preMergeDetectionIndices))
                        %iIterCand = 1
                        %if iIterCand < length(preMergeDetectionIndices)
                        if (  ((detectedEndSample(preMergeDetectionIndices(iIterCand)) + smplsMergeEventsInProximityWithinDetectionMargins) >= detectedBeginSample(preMergeDetectionIndices(iIterCand)+1)) && ...
                                (((detectedEndSample(preMergeDetectionIndices(iIterCand)+1)-detectedBeginSample(preMergeDetectionIndices(iIterCand)) + 1)) <= smplsMaxDetectionLength) )
                            tempNmergedInner = tempNmergedInner + 1;
                            tempNmergedOuter = tempNmergedOuter + 1;
                            tempMergeDetectionIndices(tempNmergedInner) = preMergeDetectionIndices(iIterCand);
                            preMergeDetectionIndices(preMergeDetectionIndices == (preMergeDetectionIndices(iIterCand)+1)) = [];
                            preMergeDetectionIndices(preMergeDetectionIndices == (preMergeDetectionIndices(iIterCand)-1)) = [];
                        end
                        iIterCand = iIterCand + 1;
                    end
                    
                    detectedEndSample(tempMergeDetectionIndices) = detectedEndSample(tempMergeDetectionIndices+1);
                    detectedMergeCount(tempMergeDetectionIndices) = detectedMergeCount(tempMergeDetectionIndices) + detectedMergeCount(tempMergeDetectionIndices+1) + 1;
                    detectedMergeCount(tempMergeDetectionIndices+1) = [];
                    detectedBeginSample(tempMergeDetectionIndices+1) = [];
                    detectedEndSample(tempMergeDetectionIndices+1) = [];
                    
                    nDetected = nDetected - tempNmergedInner;
                    
                end
            end
            detectedLengthSamples = detectedEndSample - detectedBeginSample + 1 ;
            
            
            %%%%%%%%
            
            
            % %                 tempPostMergeDetectedLengthSamples = zeros(1,nDetected);
            % %                 tempPostMergeDetectedBeginSample = zeros(1,nDetected);
            % %                 tempPostMergeDetectedEndSample = zeros(1,nDetected);
            
            minPeakDistanceSamples = ceil(((1/(maxFreq)) * fsample)/2); % half of max freq of interest in samples
            
            detectedPeak2Peaks = zeros(1,nDetected);
            detectedPeaksSamples  = zeros(1,nDetected);
            detectedTroughsSamples  = zeros(1,nDetected);
            detectedSignalTroughsSamples = zeros(nDetected,maxPeaksOrTroughsPerSpindle);
            detectedSignalPeaksSamples = zeros(nDetected,maxPeaksOrTroughsPerSpindle);
            detectedSignalTroughsSamples(:,:) = NaN;
            detectedSignalPeaksSamples(:,:) = NaN;
            detectedSDofFilteredSignal = zeros(1,nDetected);
            
            detectedEnvelopeMaxs = zeros(1,nDetected);
            detectedEnvelopeMaxSamples = zeros(1,nDetected);
            
            detected_linear_regression_freq_slope = zeros(1,nDetected);
            detected_linear_regression_freq_offset = zeros(1,nDetected);
            detected_linear_regression_freq_R_squared = zeros(1,nDetected);
            
            detected_inst_freq_troughs = zeros(nDetected,maxPeaksOrTroughsPerSpindle);
            detected_inst_freq_peaks = zeros(nDetected,maxPeaksOrTroughsPerSpindle);
            detected_inst_freq_troughs(:,:) = NaN;
            detected_inst_freq_peaks(:,:) = NaN;
            
            detectedTroughsPotential = zeros(nDetected,maxPeaksOrTroughsPerSpindle);
            detectedPeaksPotential = zeros(nDetected,maxPeaksOrTroughsPerSpindle);
            detectedTroughsPotential(:,:) = NaN;
            detectedPeaksPotential(:,:) = NaN;
            
            % %                 detectedMergeCount = zeros(1,nDetected);
            
            fprintf('channel %s, subpart %i, annotate events\n',data.label{iChan},iTr);
            for iIterCand = 1:nDetected
                % %                 preMergeDetectionIndices = 1:nDetected;
                % %                 postMergeDetectionIndicesBegins = [];
                % %                 postMergeDetectionIndicesEnds = [];
                % %                 iIterCand = 1;
                % %                 tempNmerged = 0;
                % %                 tempMergeCount = 0;
                % %                 %tempMergeDistances = [];%xxxx
                % %                 while (iIterCand <= length(preMergeDetectionIndices))
                %iIterCand = 1
                
                
                % %                     if iIterCand < length(preMergeDetectionIndices)
                % %                        if (  ((detectedEndSample(preMergeDetectionIndices(iIterCand)) + smplsMergeEventsInProximityWithinDetectionMargins) >= detectedBeginSample(preMergeDetectionIndices(iIterCand+1))) && ...
                % %                                (((detectedEndSample(preMergeDetectionIndices(iIterCand+1))-detectedBeginSample(preMergeDetectionIndices(iIterCand)) + 1)) <= smplsMaxDetectionLength) )
                % %                            postMergeDetectionIndicesBegins(iIterCand) = preMergeDetectionIndices(iIterCand);
                % %                            postMergeDetectionIndicesEnds(iIterCand) = preMergeDetectionIndices(iIterCand+1);
                % %                            tempNmerged = tempNmerged + 1;
                % %                            tempMergeCount = tempMergeCount + 1;
                % %                            %tempMergeDistances(tempMergeCount) = detectedBeginSample(preMergeDetectionIndices(iIterCand+1)) - detectedEndSample(preMergeDetectionIndices(iIterCand));%xxxx
                % %                            preMergeDetectionIndices(iIterCand+1) = [];
                % %                            continue;
                % %                        end
                % %                     end
                % %
                % %                     if tempMergeCount < 1
                % %                         postMergeDetectionIndicesBegins(iIterCand) = preMergeDetectionIndices(iIterCand);
                % %                         postMergeDetectionIndicesEnds(iIterCand) = preMergeDetectionIndices(iIterCand);
                % %                     end
                %postMergeDetectionIndices(iIterCand) = preMergeDetectionIndices(iIterCand);
                
                %currentRawDataSampleOffset = rawDataSampleOffset + detectedBeginSample(iIterCand) - 1;
                currentRawDataSampleOffset = detectedBeginSample(iIterCand) - 1;

                candSignal = frqBndPssSignal(detectedBeginSample(iIterCand):detectedEndSample(iIterCand));
                candSignal_hilbert = frqBndPssSignal_hilbert(detectedBeginSample(iIterCand):detectedEndSample(iIterCand));
                
                
                tempCandSignalminSample = find(min(candSignal) == candSignal);
                tempCandSignalmaxSample = find(max(candSignal) == candSignal);
                
                
                candSignalminSample = currentRawDataSampleOffset + tempCandSignalminSample;
                candSignalmaxSample = currentRawDataSampleOffset + tempCandSignalmaxSample;
                
                candSignalmin = candSignal(tempCandSignalminSample);
                candSignalmax = candSignal(tempCandSignalmaxSample);
                candPeak2Peak = candSignalmax - candSignalmin;
                
                %candSignal = candSignal + 5e-6;
                
                [tempCandSignalPeaks, tempCandSignalPeaksSamples] = findpeaks(candSignal,'MINPEAKHEIGHT',0,'MINPEAKDISTANCE',minPeakDistanceSamples);
                candSignalPeaks = candSignal(tempCandSignalPeaksSamples);
                candSignalPeaksSamples = currentRawDataSampleOffset + tempCandSignalPeaksSamples;
                
                
                tempAlignFactor = 1.0001*max(candSignal);
                candSignalInv = tempAlignFactor - candSignal;
                [tempCandSignalTroughs, tempCandSignalTroughsSamples] = findpeaks(candSignalInv,'MINPEAKHEIGHT',tempAlignFactor,'MINPEAKDISTANCE',minPeakDistanceSamples);
                candSignalTroughs = candSignal(tempCandSignalTroughsSamples);
                %validBelowZero = find(candSignalTroughs < 0);
                candSignalTroughsSamples = currentRawDataSampleOffset + tempCandSignalTroughsSamples;
                
                
                nCandSignalTroughs = length(candSignalTroughsSamples);
                nCandSignalPeaks = length(candSignalPeaksSamples);
                
                
                candSignalEnvelope = envelope(detectedBeginSample(iIterCand):detectedEndSample(iIterCand));
                tempCandSignalEnvelopemaxSample = find(max(candSignalEnvelope) == candSignalEnvelope);
                candSignalEnvelopemax = candSignalEnvelope(tempCandSignalEnvelopemaxSample);
                candSignalEnvelopemaxSample = currentRawDataSampleOffset + tempCandSignalEnvelopemaxSample;
                
                
                %chirp
                %                     tempConsecTroughsPeaks = sort([tempCandSignalPeaksSamples tempCandSignalTroughsSamples]/fsample);
                %                     temp_inst_freq2 = 1./(diff(tempConsecTroughsPeaks)*2);
                %
                %                     tempConsecTroughs = sort([tempCandSignalTroughsSamples]/fsample);
                %                     temp_inst_freq2_troughs = 1./(diff(tempConsecTroughs));
                %
                %                     tempConsecPeaks = sort([tempCandSignalPeaksSamples]/fsample);
                %                     temp_inst_freq2_peaks = 1./(diff(tempConsecPeaks));
                
                
                %inst_freq = (diff(unwrap(angle(hilbert(candSignal))))/(2*pi))*fsample;
                temp_inst_freq = (diff(unwrap(angle(candSignal_hilbert)))/(2*pi))*fsample;
                
                temp_time = (1:(length(candSignal)))./fsample;
                
                %[r_time_freq,p_time_freq] = corrcoef([temp_time(2:end)' temp_inst_freq']);
                
                
                temp_time_freq_regression = fitlm(temp_time(2:end),temp_inst_freq);
                %                     temp_time_freq_regression.Coefficients.Estimate(1)
                %                     temp_time_freq_regression.Coefficients.Estimate(2)
                %                     temp_time_freq_regression.Rsquared.Ordinary
                
                
                temp_linear_regression_freq_slope = temp_time_freq_regression.Coefficients.Estimate(2);
                temp_linear_regression_freq_offset = temp_time_freq_regression.Coefficients.Estimate(1);
                temp_linear_regression_freq_R_squared = temp_time_freq_regression.Rsquared.Ordinary;
                
                temp_inst_freq_troughs = temp_inst_freq(tempCandSignalTroughsSamples);
                temp_inst_freq_peaks = temp_inst_freq(tempCandSignalPeaksSamples);
                
                
                %                     figure
                %                     subplot(2,1,1);
                %                     plot(temp_time(2:end),temp_inst_freq); hold on;
                %                     plot(temp_time(2:end),temp_time_freq_regression.Coefficients.Estimate(2) * temp_time(2:end) + temp_time_freq_regression.Coefficients.Estimate(1))
                %                     hold off;
                %
                %                     subplot(2,1,2);
                %                     plot(temp_time(2:end),candSignal(2:end))
                
                
                
                
                
                %                     figure
                %                     plot((currentRawDataSampleOffset+1):(currentRawDataSampleOffset+(detectedEndSample(iIterCand)-detectedBeginSample(iIterCand))+1),candSignal)
                %                     hold on
                %                     %plot(candSignalPeaksSamples,candSignalPeaks,'g*')
                %                     %plot(begins(iCand):minPeakDistanceSamples:ends(iCand),0,'r*')
                %                     plot(candSignalminSample,candSignalmin,'b*',candSignalmaxSample,candSignalmax,'b*')
                %                     plot(candSignalTroughsSamples,candSignalTroughs,'r*',candSignalPeaksSamples,candSignalPeaks,'g*')
                
                detectedPeak2Peaks(iIterCand) = candPeak2Peak;
                detectedPeaksSamples(iIterCand) = candSignalmaxSample;
                detectedTroughsSamples(iIterCand) = candSignalminSample;
                detectedSignalTroughsSamples(iIterCand,1:nCandSignalTroughs) = candSignalTroughsSamples;
                detectedSignalPeaksSamples(iIterCand,1:nCandSignalPeaks) = candSignalPeaksSamples;
                detectedSDofFilteredSignal(iIterCand) = nanstd(candSignal);
                
                detectedEnvelopeMaxs(iIterCand) = candSignalEnvelopemax;
                detectedEnvelopeMaxSamples(iIterCand) = candSignalEnvelopemaxSample;
                
                
                
                detected_linear_regression_freq_slope(iIterCand) = temp_linear_regression_freq_slope;
                detected_linear_regression_freq_offset(iIterCand) = temp_linear_regression_freq_offset;
                detected_linear_regression_freq_R_squared(iIterCand) = temp_linear_regression_freq_R_squared;
                
                detected_inst_freq_troughs(iIterCand,1:nCandSignalTroughs) = temp_inst_freq_troughs;
                detected_inst_freq_peaks(iIterCand,1:nCandSignalPeaks) = temp_inst_freq_peaks;
                
                detectedTroughsPotential(iIterCand,1:nCandSignalTroughs) = candSignalTroughs;
                detectedPeaksPotential(iIterCand,1:nCandSignalPeaks) = candSignalPeaks;
                
                
                
                % %                     detectedMergeCount(iIterCand) = tempMergeCount;
                
                % %                     iIterCand = iIterCand + 1;
                % %                     tempMergeCount = 0;
                %tempMergeDistances = [];%xxxx
            end
            
            
            
            
            trl_detectedBeginSample = cat(2,trl_detectedBeginSample,data.time{iTr}(detectedBeginSample));
            trl_detectedEndSample = cat(2,trl_detectedEndSample,data.time{iTr}(detectedEndSample));
            trl_detectedLengthSamples = cat(2,trl_detectedLengthSamples,detectedLengthSamples/fsample);
            
            
            
            trl_detectedPeak2Peaks = cat(2,trl_detectedPeak2Peaks,detectedPeak2Peaks);
            trl_detectedPeaksSamples = cat(2,trl_detectedPeaksSamples,data.time{iTr}(detectedPeaksSamples));
            trl_detectedTroughsSamples  = cat(2,trl_detectedTroughsSamples ,data.time{iTr}(detectedTroughsSamples));
            tempsamples = detectedSignalTroughsSamples;
            tempsamples(~isnan(tempsamples)) = data.time{iTr}(tempsamples(~isnan(tempsamples)));
            trl_detectedSignalTroughsSamples = cat(1,trl_detectedSignalTroughsSamples,tempsamples);
            tempsamples = detectedSignalPeaksSamples;
            tempsamples(~isnan(tempsamples)) = data.time{iTr}(tempsamples(~isnan(tempsamples)));
            trl_detectedSignalPeaksSamples = cat(1,trl_detectedSignalPeaksSamples,tempsamples);
            trl_detectedSDofFilteredSignal = cat(2,trl_detectedSDofFilteredSignal,detectedSDofFilteredSignal);
            trl_detectedMergeCount = cat(2,trl_detectedMergeCount,detectedMergeCount);
            trl_nDetected = trl_nDetected + nDetected ;
            trl_nMerged = trl_nMerged + tempNmergedOuter; % unused
            trl_detectedEnvelopeMaxs = cat(2,trl_detectedEnvelopeMaxs,detectedEnvelopeMaxs);
            trl_detectedEnvelopeMaxSamples = cat(2,trl_detectedEnvelopeMaxSamples,data.time{iTr}(detectedEnvelopeMaxSamples));
            
            trl_detected_linear_regression_freq_slope = cat(2,trl_detected_linear_regression_freq_slope,detected_linear_regression_freq_slope);
            trl_detected_linear_regression_freq_offset = cat(2,trl_detected_linear_regression_freq_offset,detected_linear_regression_freq_offset);
            trl_detected_linear_regression_freq_R_squared = cat(2,trl_detected_linear_regression_freq_R_squared,detected_linear_regression_freq_R_squared);
            
            trl_detected_inst_freq_troughs = cat(1,trl_detected_inst_freq_troughs,detected_inst_freq_troughs);
            trl_detected_inst_freq_peaks = cat(1,trl_detected_inst_freq_peaks,detected_inst_freq_peaks);
            
            trl_detectedTroughsPotential = cat(1,trl_detectedTroughsPotential,detectedTroughsPotential);
            trl_detectedPeaksPotential = cat(1,trl_detectedPeaksPotential,detectedPeaksPotential);
            
            trl_contigSegment = cat(2,trl_contigSegment,repmat(iTr,1,nDetected));
            
            
        end
    end
    
    fprintf('channel %s, select events\n',data.label{iChan});
    tempIndexWithinThresholds = (trl_detectedPeak2Peaks >= cfg.minamplitude) & (trl_detectedPeak2Peaks <= cfg.maxamplitude) & (trl_detectedEnvelopeMaxs >= thresholdForDetectionCriterion);
    
    %filter events by SD of amplitude
    std_ampl_filter = nanstd(trl_detectedPeak2Peaks(tempIndexWithinThresholds));
    mean_ampl_filter = nanmean(trl_detectedPeak2Peaks);
    tempIndexWithinThresholds_filter_amp = (trl_detectedPeak2Peaks > (mean_ampl_filter + cfg.filterSDamp*std_ampl_filter)) | (trl_detectedPeak2Peaks < (mean_ampl_filter - cfg.filterSDamp*std_ampl_filter));

    std_dur_filter = nanstd(trl_detectedLengthSamples(tempIndexWithinThresholds));
    mean_dur_filter = nanmean(trl_detectedLengthSamples);
    tempIndexWithinThresholds_filter_dur = (trl_detectedLengthSamples > (mean_dur_filter + cfg.filterSDdur*std_dur_filter)) | (trl_detectedLengthSamples < (mean_dur_filter - cfg.filterSDdur*std_dur_filter));

    std_freq_filter = nanstd(trl_detected_linear_regression_freq_slope(tempIndexWithinThresholds));
    mean_freq_filter = nanmean(trl_detected_linear_regression_freq_slope);
    tempIndexWithinThresholds_filter_freq = (trl_detected_linear_regression_freq_slope > (mean_freq_filter + cfg.filterSDfreq*std_freq_filter)) | (trl_detected_linear_regression_freq_slope < (mean_freq_filter - cfg.filterSDfreq*std_freq_filter));

    tempIndexWithinThresholds = tempIndexWithinThresholds & ~tempIndexWithinThresholds_filter_amp & ~tempIndexWithinThresholds_filter_dur & ~tempIndexWithinThresholds_filter_freq;
    

    tempIndexWithinThresholds = find(tempIndexWithinThresholds);
    
    ch_detectedLengthSamples{iChan} = trl_detectedLengthSamples(tempIndexWithinThresholds);
    ch_detectedBeginSample{iChan} = trl_detectedBeginSample(tempIndexWithinThresholds);
    ch_detectedEndSample{iChan} = trl_detectedEndSample(tempIndexWithinThresholds);
    ch_detectedPeak2Peaks{iChan} = trl_detectedPeak2Peaks(tempIndexWithinThresholds);
    ch_detectedPeaksSamples{iChan} = trl_detectedPeaksSamples(tempIndexWithinThresholds);
    ch_detectedTroughsSamples{iChan}  = trl_detectedTroughsSamples(tempIndexWithinThresholds);
    ch_detectedSignalTroughsSamples{iChan} = trl_detectedSignalTroughsSamples(tempIndexWithinThresholds,:);
    ch_detectedSignalPeaksSamples{iChan} = trl_detectedSignalPeaksSamples(tempIndexWithinThresholds,:);
    ch_detectedSDofFilteredSignal{iChan} = trl_detectedSDofFilteredSignal(tempIndexWithinThresholds);
    ch_detectedMergeCount{iChan} = trl_detectedMergeCount(tempIndexWithinThresholds);
    ch_nDetected{iChan} = length(tempIndexWithinThresholds);
    ch_nMerged{iChan} = sum(trl_detectedMergeCount(tempIndexWithinThresholds));
    ch_densityPerEpoch{iChan} = length(tempIndexWithinThresholds)/(lengthsAcrossROIsSeconds/cfg.scoring.epochlength);
    ch_detectedEnvelopeMaxs{iChan} = trl_detectedEnvelopeMaxs(tempIndexWithinThresholds);
    ch_detectedEnvelopeMaxSamples{iChan} = trl_detectedEnvelopeMaxSamples(tempIndexWithinThresholds);
    
    ch_detected_linear_regression_freq_slope{iChan} = trl_detected_linear_regression_freq_slope(tempIndexWithinThresholds);
    ch_detected_linear_regression_freq_offset{iChan} = trl_detected_linear_regression_freq_offset(tempIndexWithinThresholds);
    ch_detected_linear_regression_freq_R_squared{iChan} = trl_detected_linear_regression_freq_R_squared(tempIndexWithinThresholds);
    
    ch_detected_inst_freq_troughs{iChan} = trl_detected_inst_freq_troughs(tempIndexWithinThresholds,:);
    ch_detected_inst_freq_peaks{iChan} = trl_detected_inst_freq_peaks(tempIndexWithinThresholds,:);
    
    ch_detectedTroughsPotential{iChan} = trl_detectedTroughsPotential(tempIndexWithinThresholds,:);
    ch_detectedPeaksPotential{iChan} = trl_detectedPeaksPotential(tempIndexWithinThresholds,:);
    
    ch_contigSegment{iChan} = trl_contigSegment(tempIndexWithinThresholds);
end
%end
fprintf('write results\n');

if PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff == 0
    usedFilterOrder_hp_preDS = 0;
    hp_preDS_hdm.Fs = NaN;
    hp_preDS_hdm.Astop = NaN;
    hp_preDS_hdm.Fstop = NaN;
    hp_preDS_hdm.F6dB = NaN;
    hp_preDS_hdm.F3dB = NaN;
    hp_preDS_hdm.TransitionWidth = NaN;
    hp_preDS_hdm.Fpass = NaN;
    hp_preDS_hdm.Apass = NaN;
end


if ~(strcmp(core_cfg.hpfilttype,'FIRdesigned') || strcmp(core_cfg.hpfilttype,'IIRdesigned'))
    
    usedFilterOrder_hp_preDS = NaN;
    hp_preDS_hdm.Fs = preDownsampleFreq;
    hp_preDS_hdm.Astop = NaN;
    hp_preDS_hdm.Fstop = NaN;
    hp_preDS_hdm.F6dB = NaN;
    hp_preDS_hdm.F3dB = PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff;
    hp_preDS_hdm.TransitionWidth = NaN;
    hp_preDS_hdm.Fpass = NaN;
    hp_preDS_hdm.Apass = NaN;
    
    if strcmp(core_cfg.hpfilttype,'but')
        if strcmp(UseFixedFilterOrder_hp,'yes')
            usedFilterOrder_hp = FilterOrder_hp;
        else
            usedFilterOrder_hp = 6;
        end
    end
end

if ~(strcmp(core_cfg.bpfilttype,'FIRdesigned') || strcmp(core_cfg.bpfilttype,'IIRdesigned'))
    
    usedFilterOrder_bp = NaN;
    bp_hdm.Fs = fsample;
    bp_hdm.Astop1 = NaN;
    bp_hdm.Fstop1 = NaN;
    bp_hdm.TransitionWidth1 = NaN;
    bp_hdm.F3dB1 = minFreq;
    bp_hdm.F6dB1 = NaN;
    bp_hdm.Fpass1 = NaN;
    bp_hdm.Apass = NaN;
    bp_hdm.Fpass2 = minFreq;
    bp_hdm.F3dB2 = NaN;
    bp_hdm.F6dB2 = NaN;
    bp_hdm.TransitionWidth2 = NaN;
    bp_hdm.Astop2 = NaN;
    bp_hdm.Fstop2 = NaN;
    if strcmp(core_cfg.bpfilttype,'but')
        if strcmp(UseFixedFilterOrder_bp,'yes')
            usedFilterOrder_bp = FilterOrder_bp;
        else
            usedFilterOrder_bp = 4;
        end
    end
end


bp_f_type_detail = '';
switch core_cfg.bpfilttype
    case 'but'
        bp_f_type_detail = 'IIR_Butterworth_ml_butter';
    case 'fir'
        bp_f_type_detail = 'FIR_window_Hamming_ml_fir1';
    case 'FIRdesigned'
        bp_f_type_detail = 'FIR_equiripple_signal_toolbox';
    case 'IIRdesigned'
        bp_f_type_detail = 'IIR_Butterworth_signal_toolbox';
end
hp_f_type_detail = '';
switch core_cfg.hpfilttype
    case 'but'
        hp_f_type_detail = 'IIR_Butterworth_ml_butter';
    case 'fir'
        hp_f_type_detail = 'FIR_window_Hamming_ml_fir1';
    case 'FIRdesigned'
        hp_f_type_detail = 'FIR_equiripple_signal_toolbox';
    case 'IIRdesigned'
        hp_f_type_detail = 'IIR_Butterworth_signal_toolbox';
end





res_filter = [];
res_filter.ori = functionname;
res_filter.type = 'spindles_filter';
res_filter.cfg = cfg;
res_filter.table = table(...
    {replaceIfEmpty(core_cfg.hpfilttype)},{replaceIfEmpty(hp_f_type_detail)},{replaceIfEmpty(core_cfg.hpfiltdir)},replaceIfEmpty(usedFilterOrder_hp_preDS),replaceIfEmpty(hp_preDS_hdm.Fs),replaceIfEmpty(hp_preDS_hdm.Astop),replaceIfEmpty(hp_preDS_hdm.Fstop),replaceIfEmpty(hp_preDS_hdm.F6dB),replaceIfEmpty(hp_preDS_hdm.F3dB),replaceIfEmpty(hp_preDS_hdm.TransitionWidth),replaceIfEmpty(hp_preDS_hdm.Fpass),replaceIfEmpty(hp_preDS_hdm.Apass),...
    {replaceIfEmpty(core_cfg.bpfilttype)},{replaceIfEmpty(bp_f_type_detail)},{replaceIfEmpty(core_cfg.bpfiltdir)},replaceIfEmpty(usedFilterOrder_bp),replaceIfEmpty(bp_hdm.Fs),replaceIfEmpty(bp_hdm.Astop1),replaceIfEmpty(bp_hdm.Fstop1),replaceIfEmpty(bp_hdm.TransitionWidth1),replaceIfEmpty(bp_hdm.F3dB1),replaceIfEmpty(bp_hdm.F6dB1),replaceIfEmpty(bp_hdm.Fpass1),replaceIfEmpty(bp_hdm.Apass),replaceIfEmpty(bp_hdm.Fpass2),replaceIfEmpty(bp_hdm.F3dB2),replaceIfEmpty(bp_hdm.F6dB2),replaceIfEmpty(bp_hdm.TransitionWidth2),replaceIfEmpty(bp_hdm.Fstop2),replaceIfEmpty(bp_hdm.Astop2),...
    'VariableNames',{...
    'hp_preDS_filter','hp_preDS_filter_type','hp_dir_and_passing','usedFilterOrder_hp_preDS','hp_preDS_Fs_Hz','hp_preDS_Astop_dB','hp_preDS_Fstop_Hz','hp_preDS_F6dB_Hz','hp_preDS_F3dB_Hz','hp_preDS_TransitionWidth_Hz','hp_preDS_Fpass_Hz','hp_preDS_Apass_dB',...
    'bp_filter','bp_filter_type','bp_dir_and_passing','usedFilterOrder_bp','bp_Fs_Hz','bp_Astop1_dB','bp_Fstop1_Hz','bp_TransitionWidth1_Hz','bp_F3dB1_Hz','bp_F6dB1_Hz','bp_Fpass1_Hz','bp_Apass_dB','bp_Fpass2_Hz','bp_F3dB2_Hz','bp_F6dB2_Hz','bp_TransitionWidth2_Hz','bp_Fstop2_Hz','bp_Astop2_dB'}...
    );


epochs = cellfun(@sleepStage2str,cfg.scoring.epochs','UniformOutput',0);
hypnStages = [cellfun(@sleepStage2str,epochs,'UniformOutput',0) ...
    cellfun(@sleepStage2str_alt,epochs,'UniformOutput',0) ...
    cellfun(@sleepStage2str_alt2,epochs,'UniformOutput',0)];

nRowsCh = numel(ch_nDetected);
nRowsEv = 0;
indEv = {};
for iChan = 1:nChannels
    iStart = nRowsEv + 1;
    nRowsEv = nRowsEv + ch_nDetected{iChan};
    iEnd = iStart + ch_nDetected{iChan} - 1;
    indEv{iChan} = iStart:iEnd;
end

output = cell(nRowsEv,28);

chs = cell(nRowsCh,1);
% ch_nDetected
% ch_densityPerEpoch
ch_detectedLengthSampless = zeros(nRowsCh,1);
ch_detectedPeak2Peakss = zeros(nRowsCh,1);
ch_freqbypeaks = zeros(nRowsCh,1);
epochlengths = repmat(cfg.scoring.epochlength,nRowsCh,1);

% ch_nMerged
lengthsAcrossROIsSecondss = repmat(lengthsAcrossROIsSeconds,nRowsCh,1);
% ch_SDfrqBndPssSignal


dataset_factorThresholdBeginEnds = repmat(dataset_factorThresholdBeginEnd,nRowsCh,1);
dataset_factorThresholdCriterions = repmat(dataset_factorThresholdCriterion,nRowsCh,1);
minFreqs = repmat(minFreq,nRowsCh,1);
maxFreqs = repmat(maxFreq,nRowsCh,1);
centerFreqFilters = repmat(cfg.centerfrequency,nRowsCh,1);
ch_detectedSDofFilteredSignals = zeros(nRowsCh,1);
    
tempNtroughsMeans = zeros(nRowsCh,1);
tempNpeaksMeans = zeros(nRowsCh,1);
ch_detected_linear_regression_freq_slopes = zeros(nRowsCh,1);
ch_detected_linear_regression_freq_offsets = zeros(nRowsCh,1);
ch_detected_linear_regression_freq_R_squareds = zeros(nRowsCh,1); 

hypnEpochsEndsSamples = cumsum(repmat(cfg.scoring.epochlength,numel(cfg.scoring.epochs),1));
hypnEpochsBeginsSamples = hypnEpochsEndsSamples - cfg.scoring.epochlength;

hypnEpochsEndsSamples = hypnEpochsEndsSamples + cfg.scoring.dataoffset;
hypnEpochsBeginsSamples = hypnEpochsBeginsSamples + cfg.scoring.dataoffset;



% for iChan = 1:nChannels
%     
%     ch = data.label{iChan};
%     
% end


for iChan = 1:nChannels
    
    ch = data.label{iChan};
    
    epochs = cell(length(ch_detectedTroughsSamples{iChan}),3);
    for iDet = 1:length(ch_detectedTroughsSamples{iChan})
        tempSample = ch_detectedTroughsSamples{iChan}(iDet);
        tempInd = ((hypnEpochsBeginsSamples <= tempSample) & (tempSample < hypnEpochsEndsSamples));
        if ~any(tempInd)
            tempInd = ((hypnEpochsBeginsSamples <= tempSample+1/fsample) & (tempSample < hypnEpochsEndsSamples));
        end
        if ~any(tempInd)
            tempInd = ((hypnEpochsBeginsSamples <= tempSample) & (tempSample-1/fsample < hypnEpochsEndsSamples));
        end
         if any(tempInd)
            epochs{iDet,1} = hypnStages(tempInd,1);
            epochs{iDet,2} = hypnStages(tempInd,2);
            epochs{iDet,3} = hypnStages(tempInd,3);
        else
            epochs{iDet,1} = '?';
            epochs{iDet,2} = '?';
            epochs{iDet,3} = '?';  
         end
    end;
    
    chs{iChan} = ch;
    tempLengthMeanSeconds = nanmean(ch_detectedLengthSamples{iChan}');
    ch_detectedLengthSampless(iChan) = nanmean(ch_detectedLengthSamples{iChan}');
    
    ch_detectedPeak2Peakss(iChan) = nanmean(ch_detectedPeak2Peaks{iChan});
    tempNtroughsMean = (length(find(~isnan(ch_detectedSignalTroughsSamples{iChan}))))/size(ch_detectedSignalTroughsSamples{iChan},1);
    tempNpeaksMean = (length(find(~isnan(ch_detectedSignalPeaksSamples{iChan}))))/size(ch_detectedSignalPeaksSamples{iChan},1);
    ch_freqbypeaks(iChan) = ((tempNtroughsMean + tempNpeaksMean)/2)/tempLengthMeanSeconds;
    
    ch_detectedSDofFilteredSignals(iChan) = nanmean(ch_detectedSDofFilteredSignal{iChan});
    
    tempNtroughsMeans(iChan) = tempNtroughsMean;
    tempNpeaksMeans(iChan) = tempNpeaksMean;

    
    ch_detected_linear_regression_freq_slopes(iChan) = nanmean(ch_detected_linear_regression_freq_slope{iChan});
    ch_detected_linear_regression_freq_offsets(iChan) = nanmean(ch_detected_linear_regression_freq_offset{iChan});
    ch_detected_linear_regression_freq_R_squareds(iChan) = nanmean(ch_detected_linear_regression_freq_R_squared{iChan});

    
    
    
    
    
    if ch_nDetected{iChan} > 0
        
%         output = cell(ch_nDetected{iChan},11);
        
        output(indEv{iChan},1) = cellstr(repmat(ch, ch_nDetected{iChan}, 1));
        
        output(indEv{iChan},2) = num2cell(ch_detectedLengthSamples{iChan}');
        output(indEv{iChan},3) = num2cell(ch_detectedBeginSample{iChan}');
        output(indEv{iChan},4) = num2cell(ch_detectedEndSample{iChan}');
        output(indEv{iChan},5) = num2cell(ch_detectedPeaksSamples{iChan}');
        output(indEv{iChan},6) = num2cell(ch_detectedTroughsSamples{iChan}');
        
        output(indEv{iChan},7) = num2cell(ch_detectedPeak2Peaks{iChan});
        output(indEv{iChan},8) = num2cell(ch_detectedEnvelopeMaxs{iChan});
        
%         output(indEv{iChan},9) = cellstr(repmat(datasetsPath, ch_nDetected{iChan}, 1));
%         output(indEv{iChan},10) = cellstr(repmat(hypnogramPath, ch_nDetected{iChan}, 1));
        output(indEv{iChan},9) = cellstr(repmat(strjoin(cfg.stages,' '), ch_nDetected{iChan}, 1));
        output(indEv{iChan},10) = num2cell(repmat(ch_SDfrqBndPssSignal{iChan}, ch_nDetected{iChan}, 1));
        
        output(indEv{iChan},11) = num2cell(repmat(minFreq, ch_nDetected{iChan}, 1));
        output(indEv{iChan},12) = num2cell(repmat(maxFreq, ch_nDetected{iChan}, 1));
        output(indEv{iChan},13) = num2cell(repmat(cfg.centerfrequency, ch_nDetected{iChan}, 1));
        
        %             output(:,16) = num2cell(ch_detectedLengthSamples{iChan}'/fsample);
        %             output(:,17) = num2cell(ch_detectedBeginSample{iChan}'/fsample);
        %             output(:,18) = num2cell(ch_detectedEndSample{iChan}'/fsample);
        %             output(:,19) = num2cell(ch_detectedPeaksSamples{iChan}'/fsample);
        %             output(:,20) = num2cell(ch_detectedTroughsSamples{iChan}'/fsample);
        %             output(:,21) = num2cell(ch_detectedEnvelopeMaxSamples{iChan}'/fsample);
        
        output(indEv{iChan},14) = num2cell((1:ch_nDetected{iChan})');
        output(indEv{iChan},15) = [epochs{:,1}];
        output(indEv{iChan},16) = [epochs{:,2}];
        output(indEv{iChan},17) = [epochs{:,3}];%cellstr(epochs(:,3));
        
        output(indEv{iChan},18) = num2cell(ch_detectedSDofFilteredSignal{iChan}');
        output(indEv{iChan},19) = num2cell(ch_detectedMergeCount{iChan}');
        
        output(indEv{iChan},20) = num2cell(ch_detected_linear_regression_freq_slope{iChan});
        output(indEv{iChan},21) = num2cell(ch_detected_linear_regression_freq_offset{iChan});
        output(indEv{iChan},22) = num2cell(ch_detected_linear_regression_freq_R_squared{iChan});
        
        output(indEv{iChan},23) = num2cell(ch_contigSegment{iChan});
        
                    tempLengthSeconds = output(indEv{iChan},2);

        evfreqs = zeros(numel(indEv{iChan}),1);
        tempNTroughs = zeros(numel(indEv{iChan}),1);
        tempNPeaks = zeros(numel(indEv{iChan}),1);

                for iLine=1:(size(output(indEv{iChan},:),1))
                    tempNTroughs(iLine) = length(find(~isnan(ch_detectedSignalTroughsSamples{iChan}(iLine,:))));
                    tempNPeaks(iLine) = length(find(~isnan(ch_detectedSignalPeaksSamples{iChan}(iLine,:))));
                    evfreqs(iLine) = ((tempNTroughs(iLine) + tempNPeaks(iLine))/2)/tempLengthSeconds{iLine};
                end
        output(indEv{iChan},24) = num2cell(evfreqs);
        output(indEv{iChan},25) = num2cell(tempNTroughs);
        output(indEv{iChan},26) = num2cell(tempNPeaks);
        output(indEv{iChan},27) = num2cell(ch_detectedEnvelopeMaxSamples{iChan});
        output(indEv{iChan},28) = cellstr(repmat(strjoin(cfg.thresholdstages,' '), ch_nDetected{iChan}, 1));




        
    
%         for iLine=1:(size(output,1))
%             fprintf(fide,'%s,',output{iLine,1});
%             tempLengthSeconds = output{iLine,16};
%             fprintf(fide,'%f,',tempLengthSeconds); %
%             fprintf(fide,'%e,',output{iLine,7});
%             tempNtroughs = ;
%             tempNpeaks = ;
%             fprintf(fide,'%f,',;
%             
%             fprintf(fide,'%i,',output{iLine,2});
%             fprintf(fide,'%i,',output{iLine,3});
%             fprintf(fide,'%i,',output{iLine,4});
%             fprintf(fide,'%i,',output{iLine,5});
%             fprintf(fide,'%i,',output{iLine,6});

%             fprintf(fide,'%e,',output{iLine,8});
%             
%             fprintf(fide,'%s,',output{iLine,9});
%             fprintf(fide,'%s,',output{iLine,10});
%             fprintf(fide,'%s,',output{iLine,11});
%             fprintf(fide,'%e,',output{iLine,12});
%             fprintf(fide,'%f,',output{iLine,13});
%             fprintf(fide,'%f,',output{iLine,14});
%             fprintf(fide,'%f,',output{iLine,15});
%             
%             
% %             fprintf(fide,'%f,',output{iLine,17});
% %             fprintf(fide,'%f,',output{iLine,18});
% %             fprintf(fide,'%f,',output{iLine,19});
% %             fprintf(fide,'%f,',output{iLine,20});
% %             fprintf(fide,'%f,',output{iLine,21});
%             fprintf(fide,'%i,',output{iLine,16});
%             fprintf(fide,'%s,',output{iLine,17});
%             fprintf(fide,'%s,',output{iLine,18});
%             fprintf(fide,'%s,',output{iLine,19});
%             
%             fprintf(fide,'%i,',output{iLine,26});
%             
%             fprintf(fide,'%e,',output{iLine,20});
%             fprintf(fide,'%i,',output{iLine,21});
%             fprintf(fide,'%f,',output{iLine,22});
%             fprintf(fide,'%f,',output{iLine,23});
%             fprintf(fide,'%f,',output{iLine,25});
%             
%             fprintf(fide,'%i,',tempNtroughs);
%             fprintf(fide,'%i\n',tempNpeaks);
%         end
        
        
%         for iLine=1:(size(output,1))
%             fprintf(fidp,'%s,',ch);
%             fprintf(fidp,'%i,',output{iLine,20});
%             fprintf(fidp,'%f,',ch_detectedSignalPeaksSamples{iChan}(iLine,1:end-1)/fsample );
%             fprintf(fidp,'%f\n',ch_detectedSignalPeaksSamples{iChan}(iLine,end)/fsample );
%         end
%         
%         for iLine=1:(size(output,1))
%             fprintf(fidt,'%s,',ch);
%             fprintf(fidt,'%i,',output{iLine,20});
%             fprintf(fidt,'%f,',ch_detectedSignalTroughsSamples{iChan}(iLine,1:end-1)/fsample );
%             fprintf(fidt,'%f\n',ch_detectedSignalTroughsSamples{iChan}(iLine,end)/fsample );
%         end
        


        

        

        
    end
    
    
end

res_channel = [];
res_channel.ori = functionname;
res_channel.type = 'spindles_channel';
res_channel.cfg = cfg;

tempvarnames = {...
    'channel',...
    'count','density_per_epoch','mean_duration_seconds','mean_amplitude_trough2peak_potential',...
    'mean_frequency_by_mean_pk_trgh_cnt_per_dur','epoch_length_seconds','merged_count','lengths_ROI_seconds',...
    'used_stages_for_detection','used_stages_for_thresholds','used_threshold_basis',...
    'used_factor_for_threshold_basis_begin_end','used_factor_for_threshold_basis_criterion',...
    'used_min_detection_pass_or_cutoff_freq','used_max_detection_pass_or_cutoff_freq',...
    'used_center_freq','mean_SD_of_filtered_signal',...
    'mean_troughs_per_event','mean_peaks_per_event',...
    'mean_linear_regression_freq_slope','mean_linear_regression_freq_offset','mean_linear_regression_freq_R_squared'};

% if isempty(output)
%     res_channel.table = cell2table(cell(0,numel(tempvarnames)), 'VariableNames', tempvarnames);
% else
    res_channel.table = table(...
        chs,...
        [ch_nDetected{:}]', [ch_densityPerEpoch{:}]', ch_detectedLengthSampless,ch_detectedPeak2Peakss,...
        ch_freqbypeaks,epochlengths,[ch_nMerged{:}]',lengthsAcrossROIsSecondss,...
        cellstr(repmat(strjoin(cfg.stages,' '), nRowsCh, 1)),cellstr(repmat(strjoin(cfg.thresholdstages,' '), nRowsCh, 1)),[ch_SDfrqBndPssSignal{:}]',...
        dataset_factorThresholdBeginEnds,dataset_factorThresholdCriterions,...
        minFreqs,maxFreqs,...
        centerFreqFilters,ch_detectedSDofFilteredSignals,...
        tempNtroughsMeans,tempNpeaksMeans,...
        ch_detected_linear_regression_freq_slopes,ch_detected_linear_regression_freq_offsets,ch_detected_linear_regression_freq_R_squareds,...
        'VariableNames',tempvarnames);
% end

res_event = [];
res_event.ori = functionname;
res_event.type = 'spindles_event';
res_event.cfg = cfg;
tempvarnames = {...
    'channel','duration_seconds','amplitude_peak2trough_max','envelope_max','frequency_by_mean_pk_trgh_cnt_per_dur',...
    'seconds_begin','seconds_end','seconds_peak_max','seconds_trough_max','seconds_envelope_max',...  %'duration_samples','sample_begin','sample_end','sample_peak_max','sample_trough_max','envelope_max',...
    'used_stages_for_detection','used_stages_for_thresholds','used_threshold_basis','used_min_pass_or_cutoff_detection_freq','used_max_detection_pass_or_cutoff_freq','used_center_freq',...
    'id_within_channel',...
    'stage','stage_alt','stage_alt2',...
    'contig_segment',...
    'SD_of_filtered_signal','merged_count',...
    'linear_regression_freq_slope','linear_regression_freq_offset','linear_regression_freq_R_squared',...
    'number_troughs','number_peaks'};

if isempty(output)
    res_event.table = cell2table(cell(0,numel(tempvarnames)), 'VariableNames', tempvarnames);
else
    res_event.table = table(...
        output(:,1),[output{:,2}]',[output{:,7}]',[output{:,8}]',[output{:,24}]',...
        [output{:,3}]',[output{:,4}]',[output{:,5}]',[output{:,6}]',[output{:,27}]',...
        output(:,9),output(:,28),[output{:,10}]',[output{:,11}]',[output{:,12}]',[output{:,13}]',...
        [output{:,14}]',...
        output(:,15),output(:,16),output(:,17),...
        [output{:,23}]',...
        [output{:,18}]',[output{:,19}]',...
        [output{:,20}]',[output{:,21}]',[output{:,22}]',...
        [output{:,25}]',[output{:,26}]',...
        'VariableNames',tempvarnames);
end
data = [];%clear
chData = [];%clear

fprintf([functionname ' function finished\n']);
toc
memtoc
end

function v = replaceIfEmpty(v)
if isempty(v)
    if isnumeric(v)
        v = NaN;
    else
        v = '';
    end
end
end
