function [res_power_fluct_bin, res_power_fluct_signal] = st_power_fluct(cfg, data)

% ST_POWER_FLUCT calculates power fluctuations
% results are stored it in a result structure
%
% Use as
%   [res_power_fluct_bin res_power_fluct_signal] = st_power_fluct(cfg, data)
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
%
% %   cfg.segmentlength  = scalar, segment length in seconds,
% %                        should be lower than the epoch length (default = 10)
% %   cfg.segmentoverlap = fraction of the overlap of segments during
% %                        segmentation in the interval [0,1) (default = 0.1)
%
%   cfg.channel  = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION
%                           Note that the channels are selected 
%                           AFTER cfg.montage is applied (if applied)
% %   cfg.foilim   = [begin end], frequency band of interest, default = [0.5 30]
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
% %   cfg.windowproportion = the fraction of hanning window (only) that is applied.
% %                          e.g. 1 means 100% of hanning window applied and
% %                          0.5 means 50% of hanning window with symmetrically
% %                          25% of each segment tail (left and right) given a hanning shape.
% %                          cfg.segmentoverlap should be adpated to at least
% %                          half this values (which is the optimal for speed or calulation, higher overlaps will increase computation time)
% %                          (default = 0.2)
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


% set defaults
cfg.channel  = ft_getopt(cfg, 'channel', 'all', 1);
cfg.foilim   = ft_getopt(cfg, 'foilim', [10 15]);
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
%cfg.windowproportion = ft_getopt(cfg, 'windowproportion', 0.2);
%cfg.segmentlength    = ft_getopt(cfg, 'segmentlength', 10);
%cfg.segmentoverlap   = ft_getopt(cfg, 'segmentoverlap', 0.1);
cfg.downsamplefs     = ft_getopt(cfg, 'downsamplefs', 100);
cfg.FreqSteps = ft_getopt(cfg, 'FreqSteps', 0.2);
cfg.TimeSteps = ft_getopt(cfg, 'TimeSteps', 0.1);

cfg.fflucmin = ft_getopt(cfg, 'fflucmin', 0.001);
cfg.fflucmax = ft_getopt(cfg, 'fflucmax', 0.12);
cfg.fflucstep = ft_getopt(cfg, 'fflucstep', 0.001);
cfg.flucStepSize = ft_getopt(cfg, 'flucStepSize', 0.5);
cfg.cycles = ft_getopt(cfg, 'cycles', 4);
cfg.flucCycles = ft_getopt(cfg, 'flucCycles', 4);
cfg.minSecBouts = ft_getopt(cfg, 'minSecBouts', 120);



%check if the frequency bands are in the requested frequency range
ftol = 0.001;
minBandFreq = cfg.foilim(1);
maxBandFreq = cfg.foilim(2);
for band = cfg.bands
    b = band{:};
    if (b{2} < minBandFreq) || (b{3} > maxBandFreq)
        ft_error('frequency band %s %f to %f not in range of freq limits = [%f %f].\nPlease change either the range or the freq input.',b{1},b{2},b{3},minBandFreq,maxBandFreq)
    end
end

PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff = minBandFreq/2;

% if cfg.segmentlength > cfg.scoring.epochlength
%     error(['Parameter cfg.segmentlength ' num2str(cfg.segmentlength) ' s must not be greater than epoch length ' num2str(cfg.scoring.epochlength) ' s!'])
% end

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

% if cfg.segmentlength < (1/minBandFreq)
%     error(['Parameter cfg.segmentlength of ' num2str(cfg.segmentlength) ' s is to short\n to support minimum band frequency of ' num2str(minBandFreq) ' Hz,\n choose either mimimum band ferquency of ' num2str((1/cfg.segmentlength)) ' Hz or adapt sequence cfg.segmentlength to ' num2str((1/minBandFreq)) ' s!'])
% end

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



fprintf([functionname ' function initialized\n']);



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
    data_t = st_select_scoring(cfg_int, data_t);
    if isempty(data_t.trial)
        hasROIs = false;
        % read in dummy data
        cfg_int.trl = [1 round(cfg.segmentlength*fsample)+1 0];
        data_t = st_preprocessing(cfg_int, data); 
%         data.time = {};
%         data.trial = {};
        data_t.sampleinfo = [0 -1];
    else
        %data = data_t;
        %data_t = [];
    end
else
    cfg_int.dataset  = cfg.dataset;
    [cfg_int] = st_select_scoring(cfg_int);
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







cfg_int = [];
cfg_int.method    = 'wavelet';%single number (in unit of time, typically seconds) of the required snippets
cfg_int.output   = 'pow';%single number (between 0 and 1 (exclusive)) specifying the fraction of overlap between snippets (0 = no overlap)
cfg_int.foi = [cfg.foilim(1):cfg.FreqSteps:cfg.foilim(2)];%
%cfg.foi = [10:0.25:15];

cfg_int.width = cfg.cycles;%7
cfg_int.pad = 'maxperlen';
cfg_int.feedback = 'no';
cfg_int.keeptrials = 'yes';
cfg_int.toi = [min(cellfun(@min,data_t.time)):cfg.TimeSteps:max(cellfun(@max,data_t.time))];
%cfg_int.channel = {'C4:A1-A2','P4:A1-A2'};
freq_t = ft_freqanalysis(cfg_int,data_t);
cfg_int.toi = [0:cfg.TimeSteps:max(cellfun(@max,data.time))];
freq_complete = ft_freqanalysis(cfg_int,data);



flucs_powers = {};
flucs_powers_freqs = {};
flucs_powers_norm = {};

%     datas = {};
%     datas_complete = {};


number_of_bouts = {};
mean_used_of_bout_length = {};

norm_mean = {};

power_signal_freq2 = {};
power_signal_freq2_complete = {};

power_signal_freq2_fig = {};
power_signal_freq2_complete_fig = {};


% 27 is 2Hz band 3 Hz below fast spindle peak
% 28 is 2Hz band 3 Hz above fast spindle peak
% 29 is fast spindle peak +- 1Hz (2Hz band)
% 30 is SWA 0.5 to 4 Hz
% 31 is sigma 10 to 15 Hz
% 32 is theta 6 to 10 Hz
% 33 is above sigma 16 to 20 Hz
% 34 is beta 20 to 24 Hz
% 35 is beta 4 to 8 Hz
for iFreqBand = 1:numel(cfg.bands)
    %for iFreqMin = [1:1:24 27 28 29 30 31 32 33 34]
    band = cfg.bands{iFreqBand};
    minFreq_temp = band{2};
    maxFreq_temp = band{3};
    
    
    cfg_temp = [];
    cfg_temp.frequency  = [minFreq_temp maxFreq_temp];
    freq_t_temp = ft_selectdata(cfg_temp,freq_t);
    %size(freq1_temp.powspctrm)
    %size(freq1.powspctrm)
    
    
    freq2 = freq_t_temp;
    %freq2.powspctrm = squeeze(nanmean(freq_t_temp.powspctrm,3));
    freq2.powspctrm = nanmean(freq_t_temp.powspctrm,3);
    sz = size(freq2.powspctrm);
    freq2.powspctrm = reshape(freq2.powspctrm(:),sz([1 2 4]));


    freq2.trial = {};
    freq2.time2 = {};
    freq2_fig = freq2;
    
    freq_complete_temp = ft_selectdata(cfg_temp,freq_complete);
    freq_complete = [];
    
    freq2_complete = freq_complete_temp;
    freq2_complete.powspctrm = nanmean(freq_complete_temp.powspctrm,3);
    freq_complete_temp = [];
    sz = size(freq2_complete.powspctrm);
    freq2_complete.powspctrm = reshape(freq2_complete.powspctrm(:),sz([2 4]));
    freq2_complete.trial = {};
    freq2_complete.time2 = {};
    freq2_complete_fig = freq2_complete;
    
    for iChn = 1:size(freq2.powspctrm,2)
        %iChn = 1
        %tempValuesToMean = freq2.powspctrm(:,iChn,:);%all parts
        %mice: normalize pwr to NonREM of the first 100 min after lights off
        %mice: normalize pwr to NonREM of the first 40 min after lights off
        tempValuesToMean = freq2.powspctrm(:,iChn,:);%all parts
        tempValuesToMean = tempValuesToMean(:);
        tempValuesToMean = tempValuesToMean(1:min(floor((1/cfg.TimeSteps)*60*100),length(tempValuesToMean(:))));
        temp_mean = nanmean(tempValuesToMean(:));
        norm_mean{iChn,iFreqBand} = temp_mean;
        for iRpt = 1:size(freq2.powspctrm,1)
            %iRpt = 1
            freq2_fig.powspctrm(iRpt,iChn,:) = freq2_fig.powspctrm(iRpt,iChn,:)./temp_mean;
            
            %freq2.powspctrm(iRpt,iChn,:) = smooth(freq2.powspctrm(iRpt,iChn,:),4/cfg.TimeSteps,'moving');%smooth the signal to 4 seconds window
            freq2.powspctrm(iRpt,iChn,:) = ( freq2.powspctrm(iRpt,iChn,:) -  temp_mean) ./temp_mean;
        end
        
        freq2_complete_fig.powspctrm(iChn,:) = freq2_complete_fig.powspctrm(iChn,:)./temp_mean;
        
        %freq2_complete.powspctrm(iChn,:) = smooth(freq2_complete.powspctrm(iChn,:),4/cfg.TimeSteps,'moving');%smooth the signal to 4 seconds window
        freq2_complete.powspctrm(iChn,:) = ( freq2_complete.powspctrm(iChn,:) -  temp_mean) ./temp_mean;
        
    end
    for iRpt = 1:size(freq2.powspctrm,1)
        %iRpt = 1;
        
        %size(freq2.powspctrm)
        freq2.trial{iRpt} = freq2.powspctrm(iRpt,:,:);
        sz = size(freq2.trial{iRpt});
        freq2.trial{iRpt} = reshape(freq2.trial{iRpt}(:),sz([2 3]));
        tempIndNan = isnan(freq2.trial{iRpt}(1,:));
        freq2.trial{iRpt}(:,tempIndNan) = [];
        freq2.time2{iRpt} = freq2.time(~tempIndNan);
        
        freq2_fig.trial{iRpt} = freq2_fig.powspctrm(iRpt,:,:);
        sz = size(freq2_fig.trial{iRpt});
        freq2_fig.trial{iRpt} = reshape(freq2_fig.trial{iRpt}(:),sz([2 3]));
        tempIndNan = isnan(freq2_fig.trial{iRpt}(1,:));
        freq2_fig.trial{iRpt}(:,tempIndNan) = [];
        freq2_fig.time2{iRpt} = freq2_fig.time(~tempIndNan);
        
    end
    iRpt = 1;
    freq2_complete.trial{iRpt} = squeeze(freq2_complete.powspctrm(:,:));
    tempIndNan = isnan(freq2_complete.trial{iRpt}(1,:));
    freq2_complete.trial{iRpt}(:,tempIndNan) = [];
    freq2_complete.time2{iRpt} = freq2_complete.time(~tempIndNan);
    
    freq2_complete_fig.trial{iRpt} = squeeze(freq2_complete_fig.powspctrm(:,:));
    tempIndNan = isnan(freq2_complete_fig.trial{iRpt}(1,:));
    freq2_complete_fig.trial{iRpt}(:,tempIndNan) = [];
    freq2_complete_fig.time2{iRpt} = freq2_complete_fig.time(~tempIndNan);
    
    power_signal_freq2{iFreqBand} = freq2;
    power_signal_freq2_complete{iFreqBand} = freq2_complete;
    freq2_complete = [];
    
    power_signal_freq2_fig{iFreqBand} = freq2_fig;
    power_signal_freq2_complete_fig{iFreqBand} = freq2_complete_fig;
    
    freq2_complete_fig = []
    
    
    freq2.time = freq2.time2;
    
    freq2 = rmfield(freq2,'dimord');
    freq2 = rmfield(freq2,'cumtapcnt');
    freq2 = rmfield(freq2,'freq');
    freq2 = rmfield(freq2,'powspctrm');
    freq2 = rmfield(freq2,'time2');
    
    freq2.fsample = 1/cfg.FreqSteps;
    
    freq2_fig.time = freq2_fig.time2;
    
    freq2_fig = rmfield(freq2_fig,'dimord');
    freq2_fig = rmfield(freq2_fig,'cumtapcnt');
    freq2_fig = rmfield(freq2_fig,'freq');
    freq2_fig = rmfield(freq2_fig,'powspctrm');
    freq2_fig = rmfield(freq2_fig,'time2');
    
    freq2_fig.fsample = 1/cfg.FreqSteps;
    
    freq2_fig = []   


    
    
    
    
    
    
    cfg_int = [];
    cfg_int.method    = 'wavelet';%single number (in unit of time, typically seconds) of the required snippets
    cfg_int.output   = 'pow';%single number (between 0 and 1 (exclusive)) specifying the fraction of overlap between snippets (0 = no overlap)
    cfg_int.foi = [cfg.fflucmin:cfg.fflucstep:cfg.fflucmax];
    cfg_int.width = cfg.flucCycles;
    cfg_int.pad = 'maxperlen';
    cfg_int.feedback = core_cfg.feedback;
    cfg_int.keeptrials = 'yes';
    cfg_int.toi = [0:cfg.flucStepSize:max(cellfun(@max,data_t.time))];
    freq3 = ft_freqanalysis(cfg_int,freq2);
    
    if length(freq3.label) > 1
        powspec = squeeze(nanmean(freq3.powspctrm,4));
    else
        powspec = (nanmean(freq3.powspctrm,4));
    end
    
    %     plot(freq3.freq,squeeze(powspec(1,1,:)))
    %     plot(freq3.freq,squeeze(powspec(2,1,:)))
    %     plot(freq3.freq,squeeze(powspec(3,1,:)))
    %     plot(freq3.freq,squeeze(powspec(4,1,:)))
    %     plot(freq3.freq,squeeze(powspec(5,1,:)))
    %     plot(freq3.freq,squeeze(powspec(6,1,:)))
    
    
    for iChanNum = 1:length(freq3.label)
        %temp_exampleChannel =         char(strrep(freq3.label(iChanNum),'_','_'));
        %titleName =  [ouputFilesPrefixString num2str(minFreq_temp) ' to ' num2str(maxFreq_temp) ' Hz Power fluctuations dataset ' num2str(iData) ' freq iter' num2str(iFreqBand) ' in ' temp_exampleChannel];
        
        
        %size(powspec)
        
        w_start = cellfun(@min,freq2.time);
        w_end = cellfun(@max,freq2.time);
        w = w_end - w_start;
        
        
        number_of_bouts{iFreqBand,iChanNum} = length(w2);
        mean_used_of_bout_length{iFreqBand,iChanNum} =  mean(w2);
        if length(w2) <2
            weighted_mean = powspec(w>=cfg.minSecBouts,iChanNum,:);
            sz = size(weighted_mean);
            weighted_mean = reshape(weighted_mean(:),sz([1 3]))';
        else
            W = w2./sum(w2);
            W = W';
            A = powspec(w>=cfg.minSecBouts,iChanNum,:);
            sz = size(A);
            A = reshape(A(:),sz([1 3]));

            %A = squeeze(powspec(1:2,iChanNum,:));
            
            
            weighted_mean = zeros(1,size(A,2));
            for iFreq = 1:size(A,2)
                weighted_mean(iFreq) = nansum(W.*A(:,iFreq))/sum(W(~isnan(A(:,iFreq))));
            end
        end
        
        %         fh = figure;
        %         plot(freq3.freq,weighted_mean);
        %
        %         title(titleName);
        %         xlabel('Frequncy (Hz)');
        %         ylabel('Power fluctuation relative to NonREM bouts');
        %
        %         figure_width = 6;     % Width in inches
        %         figure_height = 5;    % Height in inches
        %         pos = get(fh, 'Position');
        %         set(fh, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
        %         % Here we preserve the size of the image when we save it.
        %         set(fh,'InvertHardcopy','on');
        %         set(fh,'PaperUnits', 'inches');
        %         papersize = get(fh, 'PaperSize');
        %         left = (papersize(1)- figure_width)/2;
        %         bottom = (papersize(2)- figure_height)/2;
        %         myfiguresize = [left, bottom, figure_width, figure_height];
        %         set(fh,'PaperPosition', myfiguresize);
        %         set(fh,'PaperOrientation', 'portrait');
        %         %print([titleName '.png'],'-dpng','-r600');
        %         saveas(fh, [titleName '.fig']);
        %         saveas(fh, [titleName '.png']);
        %
        %         close(fh);
        
        
        
        flucs_powers{iFreqBand,iChanNum} = weighted_mean;
        flucs_powers_norm{iFreqBand,iChanNum} = weighted_mean./nanmean(weighted_mean);
        flucs_powers_freqs{iFreqBand,iChanNum} = freq3.freq;
        
        
        
        
        %plot(flucs_powers_freqs{iData,iFreqMin},flucs_powers_norm{iData,iFreqMin})
        
    end
    
end
res_power_fluct_bin.flucs_powers = flucs_powers;
res_power_fluct_bin.flucs_powers_norm = flucs_powers_norm;
res_power_fluct_bin.flucs_powers_freqs = flucs_powers_freqs;
res_power_fluct_bin.ori = functionname;
res_power_fluct_bin.type = 'power_fluct_bin';
res_power_fluct_bin.cfg = cfg;

res_power_fluct_signal.number_of_bouts = number_of_bouts;
res_power_fluct_signal.mean_used_of_bout_length = mean_used_of_bout_length;

res_power_fluct_signal.power_signal_freq2 = power_signal_freq2;
res_power_fluct_signal.power_signal_freq2_complete = power_signal_freq2_complete;

res_power_fluct_signal.power_signal_freq2_fig = power_signal_freq2_fig;
res_power_fluct_signal.power_signal_freq2_complete_fig = power_signal_freq2_complete_fig;
res_power_fluct_signal.norm_mean = norm_mean;
res_power_fluct_signal.ori = functionname;
res_power_fluct_signal.type = 'power_fluct_signal';
res_power_fluct_signal.cfg = cfg;

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end

