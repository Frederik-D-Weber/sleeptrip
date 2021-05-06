function [res_power_bin, res_power_band] = st_power(cfg, data)

% ST_POWER calculates power and power density as well as engergy and energy
% density based on Welch's method with FFTs to overlapping data segments of a length that   
% are windowed by a Hanning window and overlap
% results are stored it in a result structure
%
% Use as
%   [res_power_bin res_power_band] = st_power(cfg, data)
%   [res_power_bin res_power_band] = st_power(cfg)
%
% Required configuration parameters are:
%   cfg.scoring  = structure provided by ST_READ_SCORING
%   cfg.stages   = either a string or a Nx1 cell-array with strings that indicate sleep stages
%                  if no data structure is provided as input next to configuration also
%                  provide this parameter, possible stages are a subset of
%                  {'W', 'N1', 'N2', 'N3', 'S4', 'R'}
%
%  if no data structure is provided as parameter then you need to define
%   cfg.dataset  = string with the filename
%
% Optional configuration parameters are:
%   cfg.quick        = string, either 'yes' or 'no' (default = 'no'). if yes this
%                       will speed up calcualtion considerably by using matlab spectrogram function and 
%                       force cfg.windowproportion = 1/scoring.epochlength
%                       (thus 0.5 seconds of power per sleep stage are lost
%                       at the tails of each sleep stage) and force
%                       as well cfg.segmentlength = scoring.epochlength 
%                       and no overlap of segments cfg.segmentoverlap = 0
%
%   cfg.segmentlength  = scalar, segment length in seconds,
%                        should be lower or equal to the epoch length (default = 10, or the epoch length)
%   cfg.segmentoverlap = fraction of the overlap of segments during
%                        segmentation in the interval [0,1) (default = 0.1)
% 
%   cfg.channel  = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION
%   cfg.foilim   = [begin end], frequency band of interest, default = [0.5 30]
%   cfg.bands    = structure of frequency band definitions in the form of
%                  the default one, e.g. 
%                  cfg.bands = ...
%                             {{'SO',0.5,1},...
%                             {'SWA',0.5,4},...
%                             {'delta',1,4},...
%                             {'theta',4,8},...
%                             {'alpha',8,11},...
%                             {'slow_spindle',9,12},...
%                             {'fast_spindle',12,15},...
%                             {'spindle',9,15},...
%                             {'beta',15,30}}
%   cfg.windowproportion = the fraction of hanning window (only) that is applied. 
%                          e.g. 1 means 100% of hanning window applied and 
%                          0.5 means 50% of hanning window with symmetrically 
%                          25% of each segment tail (left and right) given a hanning shape. 
%                          cfg.segmentoverlap should be adpated to at least
%                          half this values (which is the optimal for speed or calulation, higher overlaps will increase computation time)
%                          (default = 0.2)
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

ttic = tic;
mtic = memtic;
functionname = getfunctionname();
fprintf([functionname ' function started\n']);


if ~isfield(cfg, 'scoring')
    cfg.scoring = st_read_scoring(cfg);
end

downsamplefsNotSet = false;
if ~isfield(cfg, 'downsamplefs')
    downsamplefsNotSet = true;
end

if isfield(cfg, 'montage') && isfield(cfg, 'channel')
    if ~any(ismember(cfg.channel,cfg.montage.labelnew))
        ft_error('The selected channels with cfg.channel does not match any of the ones defined in cfg.montage.\nPlease make them match or select all channels in the cfg.montage and NOT setting the cfg.channel')
    end
end

% set defaults
cfg.channel  = ft_getopt(cfg, 'channel', 'all', 1);
cfg.foilim   = ft_getopt(cfg, 'foilim', [0.5 30]);
cfg.bands    = ft_getopt(cfg, 'bands', ...
    {{'SO',0.5,1},...
    {'SWA',0.5,4},...
    {'delta',1,4},...
    {'theta',4,8},...
    {'alpha',8,11},...
    {'slow_spindle',9,12},...
    {'fast_spindle',12,15},...
    {'spindle',9,15},...
    {'beta',15,30}});
cfg.windowproportion = ft_getopt(cfg, 'windowproportion', 0.2);
cfg.segmentlength    = ft_getopt(cfg, 'segmentlength', 10);
cfg.segmentoverlap   = ft_getopt(cfg, 'segmentoverlap', 0.1);    
cfg.downsamplefs     = ft_getopt(cfg, 'downsamplefs', 100); 
cfg.quick            = ft_getopt(cfg, 'quick', 'no');                                      



if (cfg.windowproportion/2) ~= cfg.segmentoverlap
    ft_warning('the cfg.windowproportion should be exactly double of cfg.segmentoverlap to have good coverage and optimal computation time.')
end

%check if the frequency bands are in the requested frequency range                      
for band = cfg.bands
    b = band{:};
    if (b{2} < cfg.foilim(1)) || (b{3} > cfg.foilim(2))
        ft_error('frequency band %s %f to %f not in range of cfg.foilim = [%f %f].\nPlease change either the range or the frequency bands.',b{1},b{2},b{3},cfg.foilim(1),cfg.foilim(2))
    end
end
                                         

minBandFreq = cfg.foilim(1);
maxBandFreq = cfg.foilim(2);

PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff = minBandFreq/2;

if (cfg.segmentlength > cfg.scoring.epochlength) && (~istrue(cfg.quick))
    ft_warning(['Parameter cfg.segmentlength ' num2str(cfg.segmentlength) ' s should not be greater than epoch length ' num2str(cfg.scoring.epochlength) ' s!'])
    ft_warning(['Parameter cfg.segmentlength ' num2str(cfg.segmentlength) ' was set to epoch length ' num2str(cfg.scoring.epochlength) ' s!'])
    cfg.segmentlength = cfg.scoring.epochlength;
end

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

if ~(maxBandFreq*3 < fsample)
    error(['sample frequency of ' num2str(fsample) ' Hz must be MORE THAN three-fold (i.e. 3-fold) the maximal band frequency of ' num2str(maxBandFreq) ' Hz,\n even lower than requested by Nyquist-Shannon sample theorem.\n consider excludding higher frequency bands!'])
end

if ~(maxBandFreq*3 < cfg.downsamplefs)
    error(['cfg.downsamplefs resample frequency of ' num2str(cfg.downsamplefs) ' Hz must be MORE THAN three-fold (i.e. 3-fold) the maximal band frequency of ' num2str(maxBandFreq) ' Hz,\n even lower than requested by Nyquist-Shannon sample theorem.\n consider excludding higher frequency bands!'])
end

if cfg.segmentlength < (1/minBandFreq)
    error(['Parameter cfg.segmentlength of ' num2str(cfg.segmentlength) ' s is to short\n to support minimum band frequency of ' num2str(minBandFreq) ' Hz,\n choose either mimimum band ferquency of ' num2str((1/cfg.segmentlength)) ' Hz or adapt sequence cfg.segmentlength to ' num2str((1/minBandFreq)) ' s!'])
end

% set core parameters
load_core_cfg

Apass_hp = Apass;
AstopLeft_hp = AstopLeft;

useTwoPassFiltering_hp = 'no';

if ~isempty(strfind(core_cfg.hpfiltdir,'two'))
    useTwoPassFiltering_hp = 'yes';
end

%filtfilt --  The length of the input x must be more than three times the filter order (N) defined as max(length(b)-1,length(a)-1).
%For best results, make sure the sequence you are filtering has length at least three times the filter order and tapers to zero on both edges.i.e.:
minSignalLengthSamples = cfg.scoring.epochlength*fsample;
maxFilterOrder = floor((minSignalLengthSamples) / 3) - 1;


if strcmp(UseFixedFilterOrder_hp,'yes') && strcmp(useTwoPassFiltering_hp,'yes') && (FilterOrder_hp > maxFilterOrder)
    error(['filter order for high pass not optimal, use a maximal filter order three times of samples in one epoch, i.e. filter order of ' num2str(maxFilterOrder)])
elseif (FilterOrder_hp < maxFilterOrder) && strcmp(MaximizeFilterOrderIfFixedFilterOrderIsUsed,'yes')
    FilterOrder_hp = maxFilterOrder;
end



if strcmp(UseFixedFilterOrder_hp,'yes') && logical(mod(FilterOrder_hp,2))
    error('high pass order must be an even number')
end


% Note that a one- or two-pass filter has consequences for the
% strength of the filter, i.e. a two-pass filter with the same filter
% order will attenuate the signal twice as strong.

if strcmp(useTwoPassFiltering_hp,'yes') && strcmp(UseTwoPassAttenuationCorrection_hp,'yes')
    Apass_hp = Apass_hp/2;
    AstopLeft_hp = AstopLeft_hp/2;
end

if strcmp(UseFixedFilterOrder_hp,'yes') && strcmp(core_cfg.hpfilttype,'FIRdesigned')
    error('UseFixedFilterOrder_hp not allowed for high pass filters of type FIRdesigned')
end


ft_power_cfg_taper = 'hanning'; %choose either dpss or hanning as parameter for tapering method during power sprectrum analysis resulting in the use of either multi-taper method (MTM) based on discrete prolate spheroidal sequences (DPSS or Slepian sequences) as tapers or a single-taper based on Hanning or on Hamming windows. in case of hanning the ft_power_cfg_tapsmofrq parameter is ignored default hanning
ft_power_cfg_tapsmofrq = 0.1; % taper smoothing frequency for spectral smoothing trough multi-tapering (only). Note that 1 Hz smoothing meansplus-minus 1 Hz (2 Hz smoothing box) and lower values speed up calculation default 0.1

if ~(strcmp(ft_power_cfg_taper,'hanning') || strcmp(ft_power_cfg_taper,'dpss'))
    error(['Parameter ft_power_cfg_taper = ' ft_power_cfg_taper ' is not supported, for now only ft_power_cfg_taper = hanning or ft_power_cfg_taper = dpss are supported !'])
end

NENBW = [];
%NENBW.hanning = 1.5;
%NENBW.hamming = 1.3628;


fprintf([functionname ' function initialized\n']);

% signalMultiplicator = 1;
% signalOffsetSamples = cfg.scoring.dataoffset*fsample;

% %ROI
% epochLengthSamples = cfg.scoring.epochlength * fsample;
% [roiBegins, roiEnds] = getROIsByHypnogram(hypnogramPath,epochLengthSamples,cfg.stages);
% 
% if length(roiEnds) < 1
%     error(['no ROI in data left for analysis']);
% end
% 
% if (signalOffsetSamples ~= 0)
%     signalOffsetSeconds = signalOffsetSamples/fsample;
%     roiBegins = roiBegins + signalOffsetSamples;
%     roiEnds = roiEnds + signalOffsetSamples;
% end
% 
% indexLastIncludedROIinData = length(roiBegins);
% nSampleLength = -1;
% if strcmp(IgnoreDataSetHeader,'no')
%     nSampleLength = hdr.nSamples*hdr.nTrials + hdr.nSamplesPre;
%     if (roiEnds(end) > nSampleLength)
%         nMissingSamples  = roiEnds(end) - nSampleLength;
%         warning ([ num2str(nMissingSamples) ' scored (hypnogram) samples (' num2str(nMissingSamples/fsample) ' seconds) missing in the dataset ' num2str(iData)]);
%         
%         indexLastIncludedROIinData = find(roiBegins < nSampleLength,1,'last');
%         roiBegins = roiBegins(1:indexLastIncludedROIinData);
%         roiEnds = roiEnds(1:indexLastIncludedROIinData);
%         
%         if (roiEnds(end) > nSampleLength)
%             roiEnds(indexLastIncludedROIinData) = nSampleLength;
%         end;
%         
%     end
% end
% 
% if length(roiEnds) < 1
%     error(['no ROI in data left for analysis']);
% end


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
        fprintf('power: designing high pass filter for pre downsampling filtering \n');
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
fprintf('power: preprocess and pre filter data\n');


cfg_int.stages   = cfg.stages;
cfg_int.scoring  = cfg.scoring;








hasROIs = true;

if hasdata
    data_t = st_preprocessing(cfg_int, data);
    if ~istrue(cfg.quick)
    data_t = st_select_scoring(cfg_int, data_t);
    end
    if isempty(data_t.trial)
        hasROIs = false;
        % read in dummy data
        cfg_int.trl = [1 round(cfg.segmentlength*fsample)+1 0];
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
    if ~istrue(cfg.quick)
    [cfg_int] = st_select_scoring(cfg_int);
    end
    cfg_int.continuous   = 'yes';
    if isempty(cfg_int.trl) 
        hasROIs = false;
        % read in dummy data
        cfg_int.trl = [1 round(cfg.segmentlength*fsample)+1 0];
        data = st_preprocessing(cfg_int); 
%         data.time = {};
%         data.trial = {};
        data.sampleinfo = [0 -1];
    else
        data = st_preprocessing(cfg_int); 
    end

end

cfg_chan = [];
cfg_chan.channel = ft_channelselection(cfg.channel, data.label);
data = ft_selectdata(cfg_chan, data);

if (data.fsample == 128) && downsamplefsNotSet
    ft_warning('leaving 128 Hz sampling rate as default sampling rate')
    cfg.downsamplefs = 128;
end

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

scoring = cfg.scoring;

if istrue(cfg.quick)
    if scoring.dataoffset ~= 0
        maxtimedata = max(data.time{1});
        
        if scoring.dataoffset < 0
            iCut = ceil(abs(scoring.dataoffset)/scoring.epochlength);
            scoring.epochs = scoring.epochs{(iCut+1):end};
            scoring.excluded = scoring.excluded((iCut+1):end);
            scoring.dataoffset = scoring.dataoffset + iCut*scoring.epochlength;
        end
        
        if scoring.dataoffset > (maxtimedata - scoring.epochlength)
            cfg_int = [];
            hasROIs = false;
            % make in dummy data
            cfg_sd = [];
            cfg_sd.latency = [0 round(cfg.segmentlength)+1/fsample];
            data = ft_selectdata(cfg_sd,data);
        data.sampleinfo = [0 -1];
        else
        cfg_sd = [];
        cfg_sd.latency = [scoring.dataoffset maxtimedata];
        data = ft_selectdata(cfg_sd,data);
        end
    end
    
    
cfg.windowproportion = 1/scoring.epochlength;
cfg.segmentlength = scoring.epochlength;
cfg.segmentoverlap = 0;

cfg_tfr = [];
cfg_tfr.approach = 'spectrogram';
cfg_tfr.length  = cfg.segmentlength;
cfg_tfr.overlap  = cfg.segmentoverlap;
cfg_tfr.transform  = 'none';
%cfg_tfr.taper  = 'dpss';
cfg_tfr.taper  = 'hanning_proportion';
cfg_tfr.windowproportion = cfg.windowproportion;
cfg_tfr.foi = min(cfg.foilim):(1/scoring.epochlength):max(cfg.foilim);
%cfg_tfr.channel = cfg.channel_eog;
tfa = st_tfr_continuous(cfg_tfr, data);

%chan_freq_time

pFreq = tfa.freq;
pPower = tfa.powspctrm;

nEpochsData = size(pPower,3);

if numel(scoring.epochs) > nEpochsData
    scoring.epochs = scoring.epochs(1:nEpochsData);
    scoring.excluded = scoring.excluded(1:nEpochsData);
end
if numel(scoring.epochs) < nEpochsData
pPower(:,:,(numel(scoring.epochs)+1):end) = [];
end
idxValidStages = ~scoring.excluded' & ismember(scoring.epochs,cfg.stages)';
pPower(:,:,find(~idxValidStages)) = [];

%chan_freq_time  to %time_chan_freq
pPower = permute(pPower,[3 1 2]);

NSegments = sum(idxValidStages);

if hasROIs
    NSamplesPerSegment = round(cfg.segmentlength*fsample);
    NconsecutiveROIs = 0;
else
    NSegments = 0;
    NSamplesPerSegment = round(cfg.segmentlength*data.fsample);
    NconsecutiveROIs = 0;
    trlSampleLengths = 0;
end
sampleLengthsAcrossROIs = sum(trlSampleLengths);
lengthsAcrossROIsSeconds = sampleLengthsAcrossROIs/fsample; % in seconds
guaranteedROIsegmentCoverage = lengthsAcrossROIsSeconds - ((cfg.segmentlength * (1 - cfg.segmentoverlap)) * NconsecutiveROIs);
freqResolutionCalculation = (fsample/NSamplesPerSegment);
ENBW = NaN;
if strcmp(ft_power_cfg_taper,'hanning')
    %ENBW = NENBW.hanning * freqResolutionCalculation;
    %     elseif strcmp(ft_power_cfg_taper,'hamming')
    %        ENBW = NENBW.hamming * freqResolutionCalculation;
    %ft_power_cfg_taper = 'hamming'
    
    windowFunction = ft_power_cfg_taper;
    %windowFunctionValues =  window(windowFunction, NSamplesPerSegment);
    %plot(windowFunctionValues);
    %windowFunctionValues = windowFunctionValues ./ norm(windowFunctionValues);
    
    
    temp_windowFunctionValues = window(windowFunction, floor(cfg.windowproportion*NSamplesPerSegment));
    %plot(temp_windowFunctionValues);
    
    windowFunctionValuesLeft = temp_windowFunctionValues(1:floor(end/2));
    windowFunctionValuesRight = temp_windowFunctionValues(floor(end/2)+1:end);
    windowFunctionValues = ones(NSamplesPerSegment,1);
    windowFunctionValues(1:length(windowFunctionValuesLeft)) = windowFunctionValuesLeft;
    windowFunctionValues(end-length(windowFunctionValuesRight)+1:end) = windowFunctionValuesRight;
    
    %plot(windowFunctionValues);
    S1 = sum(windowFunctionValues);
    S2 = sum(windowFunctionValues.^2);
    
    NENBW = NSamplesPerSegment*(S2/(S1^2));
    ENBW = NENBW * freqResolutionCalculation;
    
end
    
else



%if hasROIs
    cfg_int = [];
    cfg_int.length    = cfg.segmentlength;%single number (in unit of time, typically seconds) of the required snippets
    cfg_int.overlap   = cfg.segmentoverlap;%single number (between 0 and 1 (exclusive)) specifying the fraction of overlap between snippets (0 = no overlap)
    cfg_int.feedback = core_cfg.feedback;
    data = ft_redefinetrial(cfg_int,data);
%end

NSegments = length(data.time);

if hasROIs
    NSamplesPerSegment = length(data.time{1,1});
    NconsecutiveROIs = size(data.trial,2);
else
    NSegments = 0;
    NSamplesPerSegment = round(cfg.segmentlength*data.fsample);
    NconsecutiveROIs = 0;
    trlSampleLengths = 0;
end
sampleLengthsAcrossROIs = sum(trlSampleLengths);
lengthsAcrossROIsSeconds = sampleLengthsAcrossROIs/fsample; % in seconds
guaranteedROIsegmentCoverage = lengthsAcrossROIsSeconds - ((cfg.segmentlength * (1 - cfg.segmentoverlap)) * NconsecutiveROIs);
freqResolutionCalculation = (fsample/NSamplesPerSegment);
ENBW = NaN;
if strcmp(ft_power_cfg_taper,'hanning')
    %ENBW = NENBW.hanning * freqResolutionCalculation;
    %     elseif strcmp(ft_power_cfg_taper,'hamming')
    %        ENBW = NENBW.hamming * freqResolutionCalculation;
    %ft_power_cfg_taper = 'hamming'
    
    windowFunction = ft_power_cfg_taper;
    %windowFunctionValues =  window(windowFunction, NSamplesPerSegment);
    %plot(windowFunctionValues);
    %windowFunctionValues = windowFunctionValues ./ norm(windowFunctionValues);
    
    
    temp_windowFunctionValues = window(windowFunction, floor(cfg.windowproportion*NSamplesPerSegment));
    %plot(temp_windowFunctionValues);
    
    windowFunctionValuesLeft = temp_windowFunctionValues(1:floor(end/2));
    windowFunctionValuesRight = temp_windowFunctionValues(floor(end/2)+1:end);
    windowFunctionValues = ones(NSamplesPerSegment,1);
    windowFunctionValues(1:length(windowFunctionValuesLeft)) = windowFunctionValuesLeft;
    windowFunctionValues(end-length(windowFunctionValuesRight)+1:end) = windowFunctionValuesRight;
    
    %plot(windowFunctionValues);
    S1 = sum(windowFunctionValues);
    S2 = sum(windowFunctionValues.^2);
    
    NENBW = NSamplesPerSegment*(S2/(S1^2));
    ENBW = NENBW * freqResolutionCalculation;
end



cfg_int = [];
cfg_int.method = 'mtmfft';
cfg_int.output = 'pow';
%cfg.pad = 'maxperlen';
cfg_int.foilim  = [minBandFreq maxBandFreq];
%cfg.foi  = [minFreq:FreqStepSize:maxFreq];
cfg_int.taper = ft_power_cfg_taper;%'hanning'
if strcmp(cfg_int.taper,'hanning') %%&& (cfg.windowproportion ~= 1)
    cfg_int.taper = 'hanning_proportion';
    cfg_int.tapervalues = windowFunctionValues;
end
%if (~(strcmp(ft_power_cfg_taper,'hanning') || strcmp(ft_power_cfg_taper,'hamming')) && strcmp(ft_power_cfg_taper,'dpss'))
if (~strcmp(ft_power_cfg_taper,'hanning')) && strcmp(ft_power_cfg_taper,'dpss')
    cfg_int.tapsmofrq = ft_power_cfg_tapsmofrq;%0.1;
end

cfg_int.channel = ft_channelselection(cfg.channel, data.label);

cfg_int.keeptrials = 'yes';
cfg_int.feedback = core_cfg.feedback;
fprintf('power: filter %i to %i Hz\n',minBandFreq,maxBandFreq);

tfa = ft_freqanalysis(cfg_int,data);
pFreq = tfa.freq;
pPower = tfa.powspctrm;
end


if ~hasROIs
    pPower(:) = 0;
end

newChannelLabels = tfa.label;
nChannels = numel(newChannelLabels);

tfa = [];%clear


%W = (trlSampleLengths./sampleLengthsAcrossROIs); %Nx1

band_ch_meanPowerSumOverSegments = [];
band_ch_meanPowerMeanOverSegments = [];


for iBand = 1:(numel(cfg.bands))
    %iBand = 1;
    band = cfg.bands{iBand};
    temp_freq = [band{2} band{3}];
    bandMinFreq = min(temp_freq);
    bandMaxFreq = max(temp_freq);
    bandfoiIndex = (pFreq >= bandMinFreq) & (pFreq <= bandMaxFreq);
    
    
    for iChan = 1:nChannels
        %iChan = 1;
        fprintf('power: process band %i to %i Hz in channel %s\n',bandMinFreq,bandMaxFreq,newChannelLabels{iChan});
        
        trl_meanPower = [];
        for iTr = 1:size(pPower,1)
            trl_meanPower(iTr,:) = pPower(iTr,iChan,bandfoiIndex);
        end
        %band_ch_meanPower{iBand,iChan} = sum(W .* mean(trl_meanPower,2))/sum(W);%weighted mean
        band_ch_meanPowerSumOverSegments{iBand,iChan} = sum(mean(trl_meanPower,2));%sum over meaned power in band of channel
        band_ch_meanPowerMeanOverSegments{iBand,iChan} = mean(mean(trl_meanPower,2));%mean over meaned power in band of channel
    end
    
end

fprintf('power: write results\n');

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

if ~(strcmp(core_cfg.hpfilttype,'FIRdesigned') || strcmp(core_cfg.hpfilttype,'IIRdesigned')) || isempty(hp_preDS_hdm.Fs);
    
    usedFilterOrder_hp_preDS = NaN;
    hp_preDS_hdm.Fs = fsample;
    hp_preDS_hdm.Astop = NaN;
    hp_preDS_hdm.Fstop = NaN;
    hp_preDS_hdm.F6dB = NaN;
    hp_preDS_hdm.F3dB = PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff;
    hp_preDS_hdm.TransitionWidth = NaN;
    hp_preDS_hdm.Fpass = NaN;
    hp_preDS_hdm.Apass = NaN;
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





% res_filters = table(...
%     {core_cfg.hpfilttype},{hp_f_type_detail},{core_cfg.hpfiltdir},usedFilterOrder_hp_preDS,hp_preDS_hdm.Fs,hp_preDS_hdm.Astop,hp_preDS_hdm.Fstop,hp_preDS_hdm.F6dB,hp_preDS_hdm.F3dB,hp_preDS_hdm.TransitionWidth,hp_preDS_hdm.Fpass,hp_preDS_hdm.Apass, ...
%     'VariableNames',{...
%     'hp_preDS_filter','hp_preDS_filter_type','hp_dir_and_passing','usedFilterOrder_hp_preDS','hp_preDS_Fs_Hz','hp_preDS_Astop_dB','hp_preDS_Fstop_Hz','hp_preDS_F6dB_Hz','hp_preDS_F3dB_Hz','hp_preDS_TransitionWidth_Hz','hp_preDS_Fpass_Hz','hp_preDS_Apass_dB'}...
%     );



nRows = numel(cfg.bands)*nChannels;

bandNames = cell(nRows,1);
chs = cell(nRows,1);

band_ch_meanPowerMeanOverSegmentss = zeros(nRows,1);
band_ch_meanPowerMeanOverSegments_density = zeros(nRows,1);
band_ch_meanPowerSumOverSegmentss = zeros(nRows,1);
band_ch_meanPowerSumOverSegmentss_density = zeros(nRows,1);

bandMinFreqs = zeros(nRows,1);
bandMaxFreqs = zeros(nRows,1);
stagess = repmat(strjoin(cfg.stages),nRows,1);
epochlengths = repmat(cfg.scoring.epochlength,nRows,1);
segmentlengths = repmat(cfg.segmentlength,nRows,1);
segmentoverlaps = repmat(cfg.segmentoverlap,nRows,1);
NSegmentss = repmat(NSegments,nRows,1);
NconsecutiveROIss = repmat(NconsecutiveROIs,nRows,1);
guaranteedROIsegmentCoverages = repmat(guaranteedROIsegmentCoverage,nRows,1);
freqResolutionCalculations = repmat(freqResolutionCalculation,nRows,1);
lengthsAcrossROIsSecondss = repmat(lengthsAcrossROIsSeconds,nRows,1);





for iBand = 1:(numel(cfg.bands))
    %iBand = 1;
    band = cfg.bands{iBand};
    
    bandName = band{1};
    temp_freq = [band{2} band{3}];
    bandMinFreq = min(temp_freq);
    bandMaxFreq = max(temp_freq);
    
    for iChan = 1:nChannels
        ch = data.label{iChan};
        iRow = (iBand - 1)*nChannels + iChan;
        
        bandNames{iRow} = bandName;
        chs{iRow} = ch;
        
        band_ch_meanPowerMeanOverSegmentss(iRow) = band_ch_meanPowerMeanOverSegments{iBand,iChan};
        band_ch_meanPowerMeanOverSegments_density(iRow) = band_ch_meanPowerMeanOverSegments{iBand,iChan}/ENBW;
        band_ch_meanPowerSumOverSegmentss(iRow) = band_ch_meanPowerSumOverSegments{iBand,iChan};
        band_ch_meanPowerSumOverSegmentss_density(iRow) = band_ch_meanPowerSumOverSegments{iBand,iChan}/ENBW;
        
        bandMinFreqs(iRow) = bandMinFreq;
        bandMaxFreqs(iRow) = bandMaxFreq;
    end
    
end

res_power_band = [];
res_power_band.ori = functionname;
res_power_band.type = 'power_band';
res_power_band.cfg = cfg;
res_power_band.table = table(...
    bandNames,chs,...
    band_ch_meanPowerMeanOverSegmentss,band_ch_meanPowerMeanOverSegments_density,...
    band_ch_meanPowerSumOverSegmentss,band_ch_meanPowerSumOverSegmentss_density,...
    bandMinFreqs,bandMaxFreqs,...
    stagess,epochlengths,segmentlengths,segmentoverlaps,...
    NSegmentss,NconsecutiveROIss, guaranteedROIsegmentCoverages, freqResolutionCalculations,...
    lengthsAcrossROIsSecondss,...
    'VariableNames',{...
    'band','channel',...
    'mean_band_of_mean_power_over_segments','mean_band_of_mean_powerDensity_over_segments',...
    'mean_band_of_arb_energy_over_segments','mean_band_of_arb_energyDensity_over_segments',...
    'min_freq_Hz','max_freq_Hz',....
    'sleep_stages','epoch_length_seconds','segment_length_seconds','segments_overlap_proportion',...
    'segment_count','consecutive_ROI_count','guaranteed_ROI_segment_coverage_seconds','frequency_true_resolution_calculation',...%'frequency_stepsize_output',...
    'lengths_ROI_seconds'}...
    );



nFreqs = length(pFreq);

nRows = nChannels*nFreqs;

chs = cell(nRows,1);
pFreqs = zeros(nRows,1);

band_ch_meanPowerMeanOverSegmentss = zeros(nRows,1);
band_ch_meanPowerMeanOverSegments_density = zeros(nRows,1);
band_ch_meanPowerSumOverSegmentss = zeros(nRows,1);
band_ch_meanPowerSumOverSegmentss_density = zeros(nRows,1);

stagess = repmat(strjoin(cfg.stages),nRows,1);
epochlengths = repmat(cfg.scoring.epochlength,nRows,1);
segmentlengths = repmat(cfg.segmentlength,nRows,1);
segmentoverlaps = repmat(cfg.segmentoverlap,nRows,1);
NSegmentss = repmat(NSegments,nRows,1);
NconsecutiveROIss = repmat(NconsecutiveROIs,nRows,1);
guaranteedROIsegmentCoverages = repmat(guaranteedROIsegmentCoverage,nRows,1);
freqResolutionCalculations = repmat(freqResolutionCalculation,nRows,1);
lengthsAcrossROIsSecondss = repmat(lengthsAcrossROIsSeconds,nRows,1);


for iChan = 1:nChannels
    ch = data.label{iChan};
    tPowSumed = squeeze(sum(pPower(:,iChan,:),1));
    tPowMeaned = squeeze(mean(pPower(:,iChan,:),1));
    for iFreq = 1:nFreqs
        iRow = (iChan - 1)*nFreqs + iFreq;
        
        chs{iRow} = ch;
        pFreqs(iRow) = pFreq(iFreq);
        
        band_ch_meanPowerMeanOverSegmentss(iRow) = tPowMeaned(iFreq);
        band_ch_meanPowerMeanOverSegments_density(iRow) = tPowMeaned(iFreq)/ENBW;
        band_ch_meanPowerSumOverSegmentss(iRow) = tPowSumed(iFreq);
        band_ch_meanPowerSumOverSegmentss_density(iRow) = tPowSumed(iFreq)/ENBW;
    end
end

res_power_bin = [];
res_power_bin.ori = functionname;
res_power_bin.type = 'power_bin';
res_power_bin.cfg = cfg;
res_power_bin.table = table(...
    chs,pFreqs,...
    band_ch_meanPowerMeanOverSegmentss,band_ch_meanPowerMeanOverSegments_density,...
    band_ch_meanPowerSumOverSegmentss,band_ch_meanPowerSumOverSegmentss_density,...
    stagess,epochlengths,segmentlengths,segmentoverlaps,...
    NSegmentss,NconsecutiveROIss, guaranteedROIsegmentCoverages, freqResolutionCalculations,...
    lengthsAcrossROIsSecondss,...
    'VariableNames',{...
    'channel','freq',...
    'mean_power_over_segments','mean_powerDensity_over_segments',...
    'arb_energy_over_segments','arb_energyDensity_over_segments',...
    'sleep_stages','epoch_length_seconds','segment_length_seconds','segments_overlap_proportion',...
    'segment_count','consecutive_ROI_count','guaranteed_ROI_segment_coverage_seconds','frequency_true_resolution_calculation',...%'frequency_stepsize_output',
    'lengths_ROI_Seconds'}...
    );



data = [];%clear
pPower = [];%clear
pFreq = [];%clear



fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end

