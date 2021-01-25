function [res_channel, res_event, res_filter] = st_slowwaves(cfg, data)

% ST_SLOWWAVES detect slow waves/K-complexes and their properties 
% like their count, density, amplitude, down-slope, up-slope, duration, frequency etc.
% results are stored it in a result structure
%
% Use as
%   [res_channel, res_event, res_filter] = st_slowwaves(cfg)
%   [res_channel, res_event, res_filter] = st_slowwaves(cfg, data)
%
% Required configuration parameters are:
%   cfg.scoring         = structure provided by ST_READ_SCORING
%
%  if no data structure is provided as parameter then you need to define
%   cfg.dataset  = string with the filename
%
% THE DEFAULT options are an adapted and improved version of the Mölle et
%             al. 2002-2013 papers
%
% Optional configuration parameters are:
%
%  cfg.channel            = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION
%                           Note that the channels are selected 
%                           AFTER cfg.montage is applied (if applied)
%  cfg.stages             = either a string or a Nx1 cell-array with strings that indicate sleep stages
%                             if no data structure is provided as input next to configuration also
%                             provide this parameter, possible stages are a subset of
%                             {'W', 'N1', 'N2', 'N3', 'S4', 'R'}
%                             (default = {'N2', 'N3', 'S4'} as the stages in which spindles typically appear)
%  cfg.thresholdstages    = either a string or a Nx1 cell-array with
%                            strings that indicate sleep stages to detemine the
%                            thresholds on. This must be a subset of cfg.stages 
%                            (default = cfg.stages)
%  cfg.minfreq            = Minimal frequency within detection criteria in Hz 
%                             i.e. select intervals of zero corssings with 
%                             minmal length of 1/f (default = 0.50)
%  cfg.maxfreq            = Maximal frequency within detection criteria in Hz 
%                             i.e. select intervals of zero corssings with 
%                             maximal length of 1/f (default = 1.11)
%  cfg.minfreqdetectfilt  = -3dB (single-)pass frequency in Hz for the high pass filter 
%                             after downsampling to look for zero crossings. 
%                             Should be greater than 0.2 (= core cfg: StopToPassTransitionWidth_hp)
%                             (default = 0.3)
%  cfg.maxfreqdetectfilt  = pass frequency in Hz for the low pass filter 
%                             after downsampling to look for zero crossings. 
%                             Should be greater than cfg.minfreqdetectfilt 
%                             (default = 3.5)
%  cfg.meanfactoramp      = Factor in means for threshold of amplitude (negative peak to positive peak) (default = 1)
%  cfg.meanfactordwnpk    = Factor in means for threshold of negative down peaks default (default = 1)
%  cfg.minamplitude       = Minimum absolute potential difference previous to take means to select as valid event (i.e. peak to peak or peak to trough)                         
%                             (default = 0)
%  cfg.maxamplitude       = Maximum absolute potential difference previous to take means to select as valid event (i.e. peak to peak or peak to trough)
%                             (default = 600)
%  cfg.minuppeak          = Minimum absolute potential previous to take means to select as valid event (e.g. 5E-6 for 5 micro Potential above zero)
%                             (default = 10)
%  cfg.maxdownpeak        = Maximum absolute potential previous to take means to select as valid event (e.g. -5E-6 for 5 micro Potential below zero)
%                             (default = -15)
%  cfg.centerfrequency    = number of the center frequency of the slow waves in Hz,
%                             see cfg.leftofcenterfreq and cfg.rightofcenterfreq
%  cfg.leftofcenterfreq   = scalar, left of of the center frequency of the spindles in Hz,
%                             see cfg.rightofcenterfreq (default = 1)
%  cfg.rightofcenterfreq  = scalar, right of of the center frequency of the spindles in Hz,
%                             see cfg.leftofcenterfreq (default = 1)
%  cfg.filterSDamp        = number of standard deviations to fiter for
%                             amplitude within the range of +- those standard deviations around the
%                             mean e.g. +-5 SD. This is applied after cfg.minamplitude and cfg.maxamplitude and cfg.minuppeak and cfg.mindownpeak criterions applied already (default 5)
%  cfg.filterSDdur        = number of standard deviations to fiter for
%                             duration within the range of +- those standard deviations around the
%                             mean e.g. +-5 SD. This is applied after cfg.minamplitude and cfg.maxamplitude and cfg.minuppeak and cfg.mindownpeak criterions applied already (default 5)
%  cfg.filterSDslopes     = number of standard deviations to fiter for
%                             upd and down slopes within the range of +- those standard deviations around the
%                             mean e.g. +-5 SD. This is applied after cfg.minamplitude and cfg.maxamplitude and cfg.minuppeak and cfg.mindownpeak criterions applied already (default 5)
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
cfg.channel     = ft_getopt(cfg, 'channel', 'all', 1);
cfg.stages      = ft_getopt(cfg, 'stages', {'N2', 'N3', 'S4'});
cfg.thresholdstages  = ft_getopt(cfg, 'thresholdstages', cfg.stages);
cfg.minfreq     = ft_getopt(cfg, 'minfreq', 0.5);
cfg.maxfreq     = ft_getopt(cfg, 'maxfreq', 1.11);
cfg.minfreqdetectfilt     = ft_getopt(cfg, 'minfreqdetectfilt', 0.3);
cfg.maxfreqdetectfilt     = ft_getopt(cfg, 'maxfreqdetectfilt', 3.5);
cfg.meanfactoramp    = ft_getopt(cfg, 'meanfactoramp', 1);
cfg.meanfactordwnpk  = ft_getopt(cfg, 'meanfactordwnpk', 1);
cfg.minamplitude     = ft_getopt(cfg, 'minamplitude', 0);
cfg.maxamplitude     = ft_getopt(cfg, 'maxamplitude', 600);
cfg.minuppeak        = ft_getopt(cfg, 'minuppeak', 10);
cfg.maxdownpeak      = ft_getopt(cfg, 'maxdownpeak', -15);
% cfg.leftofcenterfreq  = ft_getopt(cfg, 'leftofcenterfreq', 1);
% cfg.rightofcenterfreq  = ft_getopt(cfg, 'rightofcenterfreq', 1);
% cfg.thresholdaggmeth  = ft_getopt(cfg, 'thresholdaggmeth', 'respectivechan');
% cfg.envelopemeth  = ft_getopt(cfg, 'envelopemeth','smoothedRMSwd');
% cfg.factorthreshbeginend  = ft_getopt(cfg, 'factorthreshbeginend', 1.5);
% cfg.minduration  = ft_getopt(cfg, 'minduration', 0.5);
% cfg.maxduration  = ft_getopt(cfg, 'maxduration', 2.0);
cfg.filterSDamp     = ft_getopt(cfg, 'filterSDamp', 5);
cfg.filterSDdur     = ft_getopt(cfg, 'filterSDdur', 5);
cfg.filterSDslopes  = ft_getopt(cfg, 'filterSDslopes', 5);
cfg.downsamplefs    = ft_getopt(cfg, 'downsamplefs', 100);


if ~all(ismember(cfg.thresholdstages, cfg.stages))
	ft_error(['cfg.thresholdstages must be a subset of cfg.stages'])
end
thresholdingInDifferentStages = ~all(strcmp(cfg.stages,cfg.thresholdstages));



minFreq = cfg.minfreq;
maxFreq = cfg.maxfreq;

PreDetectionHighPassFilter_FpassLeft_or_F3dBcutoff = cfg.minfreqdetectfilt;
PreDetectionLowPassFilterFreq_FpassRight = cfg.maxfreqdetectfilt;


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

if isfield(cfg,'centerfrequency')
    if ~isfield(cfg,'leftofcenterfreq') || ~isfield(cfg,'rightofcenterfreq')
        error('if cfg.centerfrequency is defined then also cfg.leftofcenterfreq and cfg.rightofcenterfreq need to be defined') 
    end
    if cfg.scoring.epochlength < (1/(cfg.centerfrequency - cfg.leftofcenterfreq))
        error(['the epoch length ' num2str(cfg.scoring.epochlength) 's must not be greater in order to support the requrested maximum of ' num2str((cfg.centerfrequency - cfg.leftofcenterfreq)) ' Hz!'])
    end
    preCenterFreqFilterTo_FpassLeft = cfg.leftofcenterfreq;
    postCenterFreqFilterTo_FpassRight = cfg.rightofcenterfreq;  % in Hz
    minFreq = cfg.centerfrequency - preCenterFreqFilterTo_FpassLeft;
    maxFreq = cfg.centerfrequency + postCenterFreqFilterTo_FpassRight;
end


% set core parameters
load_core_cfg
% core_cfg

if (strcmp(core_cfg.lpfilttype,'FIRdesigned') || strcmp(core_cfg.lpfilttype,'IIRdesigned')) && (~strcmp(core_cfg.lpinstabilityfix,'no'))
    error(['band pass filter instability fix not supported for FIRdesigned or IIRdesigned'])
end

if (strcmp(core_cfg.hpfilttype,'FIRdesigned') || strcmp(core_cfg.hpfilttype,'IIRdesigned')) && (~strcmp(core_cfg.hpinstabilityfix,'no'))
    error(['high pass filter instability fix not supported for FIRdesigned or IIRdesigned'])
end

Apass_lp = Apass;
AstopRight_lp = AstopRight;

Apass_hp = Apass;
AstopLeft_hp = AstopLeft;

PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff = 0.2;



%filtfilt -- The length of the input x must be more than three times the filter order (N) defined as max(length(b)-1,length(a)-1).
%For best results, make sure the sequence you are filtering has length at least three times the filter order and tapers to zero on both edges.i.e.:
minSignalLengthSamples = cfg.scoring.epochlength*cfg.downsamplefs;
maxFilterOrder = floor((minSignalLengthSamples) / 3) - 1;

if strcmp(UseFixedFilterOrder_lp,'yes') && strcmp(useTwoPassFiltering_lp,'yes') && (FilterOrder_lp > maxFilterOrder)
    error(['filter order for low pass not optimal, use a maximal filter order three times of samples in one epoch, i.e. filter order of ' num2str(maxFilterOrder)])
elseif (FilterOrder_lp < maxFilterOrder) && strcmp(MaximizeFilterOrderIfFixedFilterOrderIsUsed,'yes')
    FilterOrder_lp = maxFilterOrder;
end

if strcmp(UseFixedFilterOrder_hp,'yes') && strcmp(useTwoPassFiltering_hp,'yes') && (FilterOrder_hp > maxFilterOrder)
    error(['filter order for high pass not optimal, use a maximal filter order three times of samples in one epoch, i.e. filter order of ' num2str(maxFilterOrder)])
elseif (FilterOrder_hp < maxFilterOrder)  && strcmp(MaximizeFilterOrderIfFixedFilterOrderIsUsed,'yes')
    FilterOrder_hp = maxFilterOrder;
end

if strcmp(UseFixedFilterOrder_lp,'yes') && logical(mod(FilterOrder_lp,2))
    error('low pass filter order must be an even number')
end

if strcmp(UseFixedFilterOrder_hp,'yes') && logical(mod(FilterOrder_hp,2))
    error('high pass filter order must be an even number')
end


% Note that a one- or two-pass filter has consequences for the
% strength of the filter, i.e. a two-pass filter with the same filter
% order will attenuate the signal twice as strong.
if strcmp(useTwoPassFiltering_lp,'yes') && strcmp(UseTwoPassAttenuationCorrection_lp,'yes')
    Apass_lp = Apass_lp/2;
    AstopRight_lp = AstopRight_lp/2;
end

if strcmp(useTwoPassFiltering_hp,'yes') && strcmp(UseTwoPassAttenuationCorrection_hp,'yes')
    Apass_hp = Apass_hp/2;
    AstopLeft_hp = AstopLeft_hp/2;
end

if strcmp(UseFixedFilterOrder_hp,'yes') && strcmp(core_cfg.hpfilttype,'FIRdesigned')
    error('UseFixedFilterOrder_hp not allowed for high pass filters of type FIRdesigned')
end



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
    if isempty(cfg_int.contendsample)
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


cfg_chan = [];
cfg_chan.channel = ft_channelselection(cfg.channel, data.label);
data = ft_selectdata(cfg_chan, data);


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






use_hp = true;
use_lp = true;
use_bp = false;
adapt_filter_settings_to_toolbox


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



    cfg_int = core_cfg;
    cfg_int.hpfilter = 'yes';
    FpassLeft = PreDetectionHighPassFilter_FpassLeft_or_F3dBcutoff; %left pass frequency in Hz
    FstopLeft = FpassLeft - StopToPassTransitionWidth_hp; %left stop frequency in Hz
    
    usedFilterOrder_hp = NaN;
    hp_hdm = [];
    if strcmp(core_cfg.hpfilttype,'IIRdesigned') || strcmp(core_cfg.hpfilttype,'FIRdesigned')
        
        hp_d = [];
        hp_hd = [];
        if strcmp(UseFixedFilterOrder_hp,'yes')
            hp_d = fdesign.highpass('N,F3db',FilterOrder_hp,FpassLeft,fsample);%designmethods(hp_preDS_d); help(hp_preDS_d,'equiripple'); help(hp_preDS_d,'butter')
        else
            hp_d = fdesign.highpass('Fst,Fp,Ast,Ap',FstopLeft,FpassLeft,AstopLeft_hp,Apass_hp,fsample);%designmethods(hp_preDS_d); help(hp_preDS_d,'equiripple'); help(hp_preDS_d,'butter')
        end
        fprintf('dataset %i: designing high pass filter for filtering SO band \n',iData);
        if strcmp(core_cfg.hpfilttype,'IIRdesigned')
            hp_hd = design(hp_d,'butter'); %isstable(hp_hd)
        elseif strcmp(core_cfg.hpfilttype,'FIRdesigned')
            hp_hd = design(hp_d,'equiripple','MinOrder', 'even');
        else
            error(['highpass filter type of ' core_cfg.hpfilttype ' unknown or not allowed'])
        end
        usedFilterOrder_hp = hp_hd.order;
        cfg_int.hpfilterdesign = hp_hd;
        hp_hdm = measure(hp_hd);
    end
    if strcmp(UseFixedFilterOrder_hp,'yes')
        cfg_int.hpfiltord     = FilterOrder_hp;
    end
    cfg_int.hpfreq        = [FpassLeft];%dummy values are overwritten by low level function
    cfg_int.feedback = core_cfg.feedback;
    fprintf('reprocess and apply high-pass filter\n');
    data = st_preprocessing(cfg_int,data);
    
    
    
    cfg_int = core_cfg;
    cfg_int.lpfilter = 'yes';
    FpassRight = PreDetectionLowPassFilterFreq_FpassRight; %right pass frequency in Hz
    FstopRight = FpassRight + PassToStopTransitionWidth_lp; %right stop frequency in Hz
    usedFilterOrder_lp = NaN;
    lp_hdm = [];
    if strcmp(core_cfg.lpfilttype,'IIRdesigned') || strcmp(core_cfg.lpfilttype,'FIRdesigned')
        lp_d = [];
        lp_hd = [];
        fprintf('dataset %i: designing low pass filter for filtering SO band \n',iData);
        if strcmp(UseFixedFilterOrder_lp,'yes')
            lp_d = fdesign.lowpass('N,Fp,Fst',FilterOrder_lp,FpassRight,FstopRight,fsample);
            lp_hd = design(lp_d,'equiripple');
        else
            lp_d = fdesign.lowpass('Fp,Fst,Ap,Ast',FpassRight,FstopRight,Apass_lp,AstopRight_lp,fsample);
            lp_hd = design(lp_d,'equiripple','MinOrder', 'even');
        end
        usedFilterOrder_lp = lp_hd.order;
        cfg_int.lpfilterdesign = lp_hd;
        lp_hdm = measure(lp_hd);
    end
    if strcmp(UseFixedFilterOrder_lp,'yes')
        cfg_int.lpfiltord     = FilterOrder_lp;
    end
    cfg_int.lpfreq        = [FpassRight];%dummy values are overwritten by low level function
    cfg_int.feedback = core_cfg.feedback;
    fprintf('reprocess and apply low-pass filter\n');
    data = st_preprocessing(cfg_int,data);
    





sampleLengthsAcrossROIs = sum(trlSampleLengths);
lengthsAcrossROIsSeconds = sampleLengthsAcrossROIs/fsample; % in seconds


smplsMinDetectionLength  = round(fsample/maxFreq);
smplsMaxDetectionLength  = round(fsample/minFreq);

minFreqPostFreqBorderBufferLength = round(fsample/minFreq);


nChannels = length(data.label);

ch_detectedLengthSamples = [];
ch_detectedBeginSample = [];
ch_detectedEndSample = [];
ch_detectedSignalMin = [];
ch_detectedSignalMax = [];
ch_detectedPeak2Peaks = [];
ch_detectedPeaksSamples  = [];
ch_detectedTroughsSamples  = [];
ch_detectedMaxSlopes  = [];
ch_detectedMaxSlopesSamples  = [];
ch_detectedMaxDownSlopes  = [];
ch_detectedMaxDownSlopesSamples  = [];
ch_detectedTroughToZeroCrossingSlopes = [];
ch_detectedTroughToZeroCrossingSlopesSamples = [];
ch_detectedZeroCrossingToTroughSlopes = [];
%ch_detectedZeroCrossingToTroughSlopesSamples = [];
ch_detectedSDofFilteredSignal = [];
ch_nDetected = [];
ch_densityPerEpoch = [];
ch_meanNegativePeak =  [];
ch_meanPeak2Peak = [];

ch_contigSegment = [];



for iChan = 1:nChannels
    
    %iChan = 1;
    
    fprintf('process channel %s ...\n',data.label{iChan});
    cfg_int = [];
    cfg_int.feedback = core_cfg.feedback;
    cfg_int.channel = ft_channelselection(data.label{iChan}, data.label);
    chData = ft_selectdata(cfg_int,data);
    
    
    
    trl_detectedLengthSamples = [];
    trl_detectedBeginSample = [];
    trl_detectedEndSample = [];
    trl_detectedSignalMin = [];
    trl_detectedSignalMax = [];
    trl_detectedPeak2Peaks = [];
    trl_detectedPeaksSamples  = [];
    trl_detectedTroughsSamples  = [];
    trl_detectedMaxSlopes  = [];
    trl_detectedMaxSlopesSamples  = [];
    trl_detectedMaxDownSlopes  = [];
    trl_detectedMaxDownSlopesSamples  = [];
    trl_detectedTroughToZeroCrossingSlopes = [];
    trl_detectedTroughToZeroCrossingSlopesSamples = [];
    trl_detectedZeroCrossingToTroughSlopes = [];
    %trl_detectedZeroCrossingToTroughSlopesSamples = [];
    trl_detectedSDofFilteredSignal = [];
    trl_nDetected = 0;
    
    trl_contigSegment = [];
    
    for iTr = 1:size(chData.trial,2)
        %fprintf('channel %s, subpart %i, preselect events of appropriate zero-crossings\n',data.label{iChan},iTr);
        fprintf('channel %s, subpart %i, detecting and annotating events\n',data.label{iChan},iTr);
                    
            frqBndPssSignal = chData.trial{iTr};
            %frqBndPssSignal_hilbert = hilbert(frqBndPssSignal);

            
            %thresholdForDetection = ch_SDfrqBndPssSignal*factorSD;
            
            lengthSignal = length(frqBndPssSignal);
            
            [begins, ends] = getBeginsAndCorrespondingEndsIndicesBelowThreshold(frqBndPssSignal,0);
            
            firstSmplsMinFreqPostFreqBorderBufferLength = minFreqPostFreqBorderBufferLength;
            lastSmplsMinFreqPostFreqBorderBufferLength = lengthSignal - minFreqPostFreqBorderBufferLength;
            indicesValidSamples = find((begins >= firstSmplsMinFreqPostFreqBorderBufferLength) & (ends <= lastSmplsMinFreqPostFreqBorderBufferLength)); %consider border effects filter
            
            begins = begins(indicesValidSamples);
            ends = ends(indicesValidSamples);
            
            nDetected = length(begins)-1;
            if (nDetected > 0)
                ends = begins(2:end)';
                begins = begins(1:end-1)';
                
                tempCandidatesLengths = ends - begins + 1;
                
                indicesCandiates = find((tempCandidatesLengths >= smplsMinDetectionLength) & (tempCandidatesLengths <= smplsMaxDetectionLength));
                
                nDetected = length(indicesCandiates);
            end
     
        
        if (nDetected > 0)
            %detectedLengthSamples = tempCandidatesLengths(indicesCandiates);
            detectedBeginSample = begins(indicesCandiates);
            detectedEndSample = ends(indicesCandiates);
            detectedLengthSamples = detectedEndSample - detectedBeginSample + 1 ;
            
            minPeakDistanceSamples = ceil(((1/(maxFreq)) * fsample)/2); % half of max freq of interest in samples
            
            detectedSignalMin = zeros(1,nDetected);
            detectedSignalMax = zeros(1,nDetected);
            detectedPeak2Peaks = zeros(1,nDetected);
            detectedPeaksSamples  = zeros(1,nDetected);
            detectedTroughsSamples  = zeros(1,nDetected);
            detectedMaxSlopes  = zeros(1,nDetected);
            detectedMaxSlopesSamples  = zeros(1,nDetected);
            detectedMaxDownSlopes  = zeros(1,nDetected);
            detectedMaxDownSlopesSamples  = zeros(1,nDetected);
            detectedTroughToZeroCrossingSlopes = zeros(1,nDetected);
            detectedTroughToZeroCrossingSlopesSamples = zeros(1,nDetected);
            detectedZeroCrossingToTroughSlopes  = zeros(1,nDetected);
            %detectedZeroCrossingToTroughSlopesSamples  = zeros(1,nDetected);
            
            detectedSDofFilteredSignal = zeros(1,nDetected);
            
            %fprintf('channel %s, subpart %i, annotate events\n',data.label{iChan},iTr);
            for iIterCand = 1:nDetected
                
                %currentRawDataSampleOffset = rawDataSampleOffset + detectedBeginSample(iIterCand) - 1;
                currentRawDataSampleOffset = detectedBeginSample(iIterCand) - 1;

                candSignal = frqBndPssSignal(detectedBeginSample(iIterCand):detectedEndSample(iIterCand));
                
                tempCandSignalmaxSample = find(max(candSignal) == candSignal,1,'first');
                tempCandSignalminSample = find(min(candSignal(1:tempCandSignalmaxSample)) == candSignal(1:tempCandSignalmaxSample),1,'first');
                              
                candSignalminSample = currentRawDataSampleOffset + tempCandSignalminSample;
                candSignalmaxSample = currentRawDataSampleOffset + tempCandSignalmaxSample;
                
                
                candSignalMaxSlope = max(diff(candSignal(tempCandSignalminSample:tempCandSignalmaxSample)))*fsample;% in potential/s
                 if isempty(candSignalMaxSlope)
                    candSignalMaxSlope = NaN;
                    candSignalMaxSlopeSample = tempCandSignalmaxSample;
                else
                    candSignalMaxSlopeSample = find(diff(candSignal)*fsample == candSignalMaxSlope);
                    candSignalMaxSlopeSample = intersect(candSignalMaxSlopeSample,tempCandSignalminSample:tempCandSignalmaxSample);
                    candSignalMaxSlopeSample = candSignalMaxSlopeSample(find(abs(candSignal(candSignalMaxSlopeSample)) == min(abs(0-(candSignal(candSignalMaxSlopeSample)))),1,'first'));
                    candSignalMaxSlopeSample =  currentRawDataSampleOffset + candSignalMaxSlopeSample;
                 end

                candSignalMaxDownSlope = min(diff(candSignal(1:tempCandSignalminSample)))*fsample;% in potential/s
                if isempty(candSignalMaxDownSlope)
                    candSignalMaxDownSlope = NaN;
                    candSignalMaxDownSlopeSample = tempCandSignalminSample;
                else
                    candSignalMaxDownSlopeSample = find(diff(candSignal)*fsample == candSignalMaxDownSlope);
                    candSignalMaxDownSlopeSample = intersect(candSignalMaxDownSlopeSample,1:tempCandSignalminSample);
                    candSignalMaxDownSlopeSample = candSignalMaxDownSlopeSample(find(abs(candSignal(candSignalMaxDownSlopeSample)) == min(abs(0-(candSignal(candSignalMaxDownSlopeSample)))),1,'first'));
                    candSignalMaxDownSlopeSample =  currentRawDataSampleOffset + candSignalMaxDownSlopeSample;
                end
                    
                    
                    
                candSignalmin = candSignal(tempCandSignalminSample);
                candSignalmax = candSignal(tempCandSignalmaxSample);
                candPeak2Peak = candSignalmax - candSignalmin;
                
                tempCandSignalMinMaxZeroCrossingSample = tempCandSignalminSample + find( abs(0-candSignal(tempCandSignalminSample:tempCandSignalmaxSample)) == min( abs(0-candSignal(tempCandSignalminSample:tempCandSignalmaxSample)) ),1,'first' ) - 1;
                candSignalMinMaxZeroCrossingSample = currentRawDataSampleOffset + tempCandSignalMinMaxZeroCrossingSample;
                candSignalTroughToZeroCrossingSlope = abs(candSignalmin)/((tempCandSignalMinMaxZeroCrossingSample - tempCandSignalminSample)/fsample); % in potential/s
                
                
                %tempCandSignalDownMinMaxZeroCrossingSample = 0 + find( abs(0-candSignal(1:tempCandSignalminSample)) == min( abs(0-candSignal(1:tempCandSignalminSample)) ),1,'first' ) - 1;
                %candSignalDownMinMaxZeroCrossingSample = currentRawDataSampleOffset + tempCandSignalDownMinMaxZeroCrossingSample;
                candSignalDownZeroCrossingToTroughSlope = -abs(candSignalmin)/((tempCandSignalminSample - 0)/fsample); % in potential/s
                
    
                detectedSignalMin(iIterCand) = candSignalmin;
                detectedSignalMax(iIterCand) = candSignalmax;
                detectedPeak2Peaks(iIterCand) = candPeak2Peak;
                detectedPeaksSamples(iIterCand) = candSignalmaxSample;
                detectedTroughsSamples(iIterCand) = candSignalminSample;
                detectedMaxSlopes(iIterCand) = candSignalMaxSlope;
                detectedMaxSlopesSamples(iIterCand) = candSignalMaxSlopeSample;
                detectedMaxDownSlopes(iIterCand) = candSignalMaxDownSlope;
                detectedMaxDownSlopesSamples(iIterCand) = candSignalMaxDownSlopeSample;
                detectedTroughToZeroCrossingSlopes(iIterCand) = candSignalTroughToZeroCrossingSlope;
                detectedTroughToZeroCrossingSlopesSamples(iIterCand) = candSignalMinMaxZeroCrossingSample;
                detectedZeroCrossingToTroughSlopes(iIterCand) = candSignalDownZeroCrossingToTroughSlope;
                %detectedZeroCrossingToTroughSlopesSamples(iIterCand) = candSignalDownMinMaxZeroCrossingSample;
                %detectedSignalTroughsSamples(iIterCand,1:nCandSignalTroughs) = candSignalTroughsSamples;
                %detectedSignalPeaksSamples(iIterCand,1:nCandSignalPeaks) = candSignalPeaksSamples;
                detectedSDofFilteredSignal(iIterCand) = nanstd(candSignal);
                

            end
            
        
            
            trl_detectedLengthSamples = cat(2,trl_detectedLengthSamples,detectedLengthSamples/fsample);
            trl_detectedBeginSample = cat(2,trl_detectedBeginSample,data.time{iTr}(detectedBeginSample));
            trl_detectedEndSample = cat(2,trl_detectedEndSample,data.time{iTr}(detectedEndSample));
            
            trl_detectedSignalMin = cat(2,trl_detectedSignalMin,detectedSignalMin);
            trl_detectedSignalMax = cat(2,trl_detectedSignalMax,detectedSignalMax);
            trl_detectedPeak2Peaks = cat(2,trl_detectedPeak2Peaks,detectedPeak2Peaks);
            trl_detectedPeaksSamples = cat(2,trl_detectedPeaksSamples,data.time{iTr}(detectedPeaksSamples));
            trl_detectedTroughsSamples  = cat(2,trl_detectedTroughsSamples,data.time{iTr}(detectedTroughsSamples));
            trl_detectedMaxSlopes = cat(2,trl_detectedMaxSlopes ,detectedMaxSlopes);
            trl_detectedMaxSlopesSamples = cat(2,trl_detectedMaxSlopesSamples,data.time{iTr}(detectedMaxSlopesSamples));
            trl_detectedMaxDownSlopes = cat(2,trl_detectedMaxDownSlopes ,detectedMaxDownSlopes);
            trl_detectedMaxDownSlopesSamples = cat(2,trl_detectedMaxDownSlopesSamples,data.time{iTr}(detectedMaxDownSlopesSamples));
            trl_detectedTroughToZeroCrossingSlopes = cat(2,trl_detectedTroughToZeroCrossingSlopes,detectedTroughToZeroCrossingSlopes);
            trl_detectedTroughToZeroCrossingSlopesSamples = cat(2,trl_detectedTroughToZeroCrossingSlopesSamples,data.time{iTr}(detectedTroughToZeroCrossingSlopesSamples));
            trl_detectedZeroCrossingToTroughSlopes = cat(2,trl_detectedZeroCrossingToTroughSlopes,detectedZeroCrossingToTroughSlopes);
            %trl_detectedZeroCrossingToTroughSlopesSamples = cat(2,trl_detectedZeroCrossingToTroughSlopesSamples,detectedZeroCrossingToTroughSlopesSamples);
            %trl_detectedSignalTroughsSamples = cat(1,trl_detectedSignalTroughsSamples,detectedSignalTroughsSamples);
            %trl_detectedSignalPeaksSamples = cat(1,trl_detectedSignalPeaksSamples,detectedSignalPeaksSamples);
            trl_detectedSDofFilteredSignal = cat(2,trl_detectedSDofFilteredSignal,detectedSDofFilteredSignal);
            
            trl_nDetected = trl_nDetected + nDetected;
            
            trl_contigSegment = cat(2,trl_contigSegment,repmat(iTr,1,nDetected));
            
        end
    end
    
    
     indInStage = [true];

    if thresholdingInDifferentStages && hasROIs_threshold
        indInStage = [false];
        cfg_int = [];
        cfg_int.stages   = cfg.thresholdstages;
        cfg_int.scoring  = cfg.scoring;
        
        [begins, ends] = st_select_scoring(cfg_int, data);
        for iBegins = 1:numel(begins)
            indInStage = indInStage | ((begins(iBegins) >= trl_detectedTroughsSamples) & (trl_detectedTroughsSamples <= ends(iBegins)));
        end
        
        
    end

    fprintf('channel %s, select events\n',data.label{iChan});
        tempIndexAboveMeanThreshold_premean = (trl_detectedPeak2Peaks >= cfg.minamplitude) ...
            & (trl_detectedPeak2Peaks <= cfg.maxamplitude) ...
            & (trl_detectedSignalMax >= cfg.minuppeak) ...
            & (trl_detectedSignalMin <= cfg.maxdownpeak ...
            & indInStage);
        
        ch_meanNegativePeak{iChan} =  nanmean(trl_detectedSignalMin(tempIndexAboveMeanThreshold_premean));
        ch_meanPeak2Peak{iChan} =  nanmean(trl_detectedPeak2Peaks(tempIndexAboveMeanThreshold_premean));
        
        tempIndexAboveMeanThreshold_temp = (trl_detectedTroughsSamples >= (ch_meanNegativePeak{iChan}*cfg.meanfactordwnpk)) ...
            & (trl_detectedPeak2Peaks >= (ch_meanPeak2Peak{iChan}*cfg.meanfactoramp));
        
        tempIndexWithinThresholds = (tempIndexAboveMeanThreshold_temp & tempIndexAboveMeanThreshold_premean);
        
        
    %filter events by SD of amplitude
    std_ampl_filter = nanstd(trl_detectedPeak2Peaks(tempIndexWithinThresholds));
    mean_ampl_filter = nanmean(trl_detectedPeak2Peaks);
    tempIndexWithinThresholds_filter_amp = (trl_detectedPeak2Peaks > (mean_ampl_filter + cfg.filterSDamp*std_ampl_filter)) | (trl_detectedPeak2Peaks < (mean_ampl_filter - cfg.filterSDamp*std_ampl_filter));

    std_dur_filter = nanstd(trl_detectedLengthSamples(tempIndexWithinThresholds));
    mean_dur_filter = nanmean(trl_detectedLengthSamples);
    tempIndexWithinThresholds_filter_dur = (trl_detectedLengthSamples > (mean_dur_filter + cfg.filterSDdur*std_dur_filter)) | (trl_detectedLengthSamples < (mean_dur_filter - cfg.filterSDdur*std_dur_filter));
        
    std_downslope_filter = nanstd(trl_detectedMaxDownSlopes(tempIndexWithinThresholds));
    mean_downslope_filter = nanmean(trl_detectedMaxDownSlopes);
    tempIndexWithinThresholds_filter_downslope = (trl_detectedMaxDownSlopes > (mean_downslope_filter + cfg.filterSDslopes*std_downslope_filter)) | (trl_detectedMaxDownSlopes < (mean_downslope_filter - cfg.filterSDslopes*std_downslope_filter));

    std_upslope_filter = nanstd(trl_detectedMaxSlopes(tempIndexWithinThresholds));
    mean_upslope_filter = nanmean(trl_detectedMaxSlopes);
    tempIndexWithinThresholds_filter_upslope = (trl_detectedMaxSlopes > (mean_upslope_filter + cfg.filterSDslopes*std_upslope_filter)) | (trl_detectedMaxSlopes < (mean_upslope_filter - cfg.filterSDslopes*std_upslope_filter));

    tempIndexWithinThresholds = tempIndexWithinThresholds & ~tempIndexWithinThresholds_filter_amp & ~tempIndexWithinThresholds_filter_dur & ~tempIndexWithinThresholds_filter_downslope & ~tempIndexWithinThresholds_filter_upslope;

    tempIndexWithinThresholds = find(tempIndexWithinThresholds);
    

    ch_detectedLengthSamples{iChan} = trl_detectedLengthSamples(tempIndexWithinThresholds);
    ch_detectedBeginSample{iChan} = trl_detectedBeginSample(tempIndexWithinThresholds);
    ch_detectedEndSample{iChan} = trl_detectedEndSample(tempIndexWithinThresholds);
    ch_detectedSignalMin{iChan} = trl_detectedSignalMin(tempIndexWithinThresholds);
    ch_detectedSignalMax{iChan} = trl_detectedSignalMax(tempIndexWithinThresholds);
    ch_detectedPeak2Peaks{iChan} = trl_detectedPeak2Peaks(tempIndexWithinThresholds);
    ch_detectedPeaksSamples{iChan} = trl_detectedPeaksSamples(tempIndexWithinThresholds);
    ch_detectedTroughsSamples{iChan}  = trl_detectedTroughsSamples(tempIndexWithinThresholds);
    ch_detectedMaxSlopes{iChan} = trl_detectedMaxSlopes(tempIndexWithinThresholds);
    ch_detectedMaxSlopesSamples{iChan} = trl_detectedMaxSlopesSamples(tempIndexWithinThresholds);
    ch_detectedMaxDownSlopes{iChan} = trl_detectedMaxDownSlopes(tempIndexWithinThresholds);
    ch_detectedMaxDownSlopesSamples{iChan} = trl_detectedMaxDownSlopesSamples(tempIndexWithinThresholds);
    
    ch_detectedTroughToZeroCrossingSlopes{iChan} = trl_detectedTroughToZeroCrossingSlopes(tempIndexWithinThresholds);
    ch_detectedTroughToZeroCrossingSlopesSamples{iChan} = trl_detectedTroughToZeroCrossingSlopesSamples(tempIndexWithinThresholds);
    ch_detectedZeroCrossingToTroughSlopes{iChan} = trl_detectedZeroCrossingToTroughSlopes(tempIndexWithinThresholds);
    %ch_detectedZeroCrossingToTroughSlopesSamples{iChan} = trl_detectedZeroCrossingToTroughSlopesSamples(tempIndexWithinThresholds);
    
    %ch_detectedSignalTroughsSamples{iChan} = trl_detectedSignalTroughsSamples;
    %ch_detectedSignalPeaksSamples{iChan} = trl_detectedSignalPeaksSamples;
    ch_detectedSDofFilteredSignal{iChan} = trl_detectedSDofFilteredSignal(tempIndexWithinThresholds);
    
    ch_nDetected{iChan} = length(tempIndexWithinThresholds);
    
    ch_densityPerEpoch{iChan} = length(tempIndexWithinThresholds)/(lengthsAcrossROIsSeconds/cfg.scoring.epochlength);
    ch_contigSegment{iChan} = trl_contigSegment(tempIndexWithinThresholds);
    
end


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
    
    usedFilterOrder_hp = NaN;
    hp_hdm.Fs = fsample;
    hp_hdm.Astop = NaN;
    hp_hdm.Fstop = NaN;
    hp_hdm.F6dB = NaN;
    hp_hdm.F3dB = FpassLeft;
    hp_hdm.TransitionWidth = NaN;
    hp_hdm.Fpass = NaN;
    hp_hdm.Apass = NaN;
    
    if strcmp(core_cfg.hpfilttype,'but')
        if strcmp(UseFixedFilterOrder_hp,'yes')
            usedFilterOrder_hp = FilterOrder_hp;
            usedFilterOrder_hp_preDS = FilterOrder_hp;
        else
            usedFilterOrder_hp = 6;
            usedFilterOrder_hp_preDS = 6;
        end
    end
end

if ~(strcmp(core_cfg.lpfilttype,'FIRdesigned') || strcmp(core_cfg.lpfilttype,'IIRdesigned'))
        
        usedFilterOrder_lp = NaN;
        lp_hdm.Fs = fsample;
        lp_hdm.Astop = NaN;
        lp_hdm.Fstop = NaN;
        lp_hdm.F6dB = NaN;
        lp_hdm.F3dB = FpassRight;
        lp_hdm.TransitionWidth = NaN;
        lp_hdm.Fpass = NaN;
        lp_hdm.Apass = NaN;
        
        if strcmp(core_cfg.lpfilttype,'but')
            if strcmp(UseFixedFilterOrder_lp,'yes')
                usedFilterOrder_lp = FilterOrder_lp;
            else
                usedFilterOrder_lp = 6;
            end
        end
        
end
    

lp_f_type_detail = '';
switch core_cfg.lpfilttype
    case 'but'
        lp_f_type_detail = 'IIR_Butterworth_ml_butter';
    case 'fir'
        lp_f_type_detail = 'FIR_window_Hamming_ml_fir1';
    case 'FIRdesigned'
        lp_f_type_detail = 'FIR_equiripple_signal_toolbox';
    case 'IIRdesigned'
        lp_f_type_detail = 'IIR_Butterworth_signal_toolbox';
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
res_filter.ori   = functionname;
res_filter.type  = 'slowwaves_filter';
res_filter.cfg   = cfg;
res_filter.table = table(...
    {replaceIfEmpty(core_cfg.hpfilttype)},{replaceIfEmpty(hp_f_type_detail)},{replaceIfEmpty(core_cfg.hpfiltdir)},replaceIfEmpty(usedFilterOrder_hp_preDS),replaceIfEmpty(hp_preDS_hdm.Fs),replaceIfEmpty(hp_preDS_hdm.Astop),replaceIfEmpty(hp_preDS_hdm.Fstop),replaceIfEmpty(hp_preDS_hdm.F6dB),replaceIfEmpty(hp_preDS_hdm.F3dB),replaceIfEmpty(hp_preDS_hdm.TransitionWidth),replaceIfEmpty(hp_preDS_hdm.Fpass),replaceIfEmpty(hp_preDS_hdm.Apass),...
    {replaceIfEmpty(core_cfg.hpfilttype)},{replaceIfEmpty(hp_f_type_detail)},{replaceIfEmpty(core_cfg.hpfiltdir)},replaceIfEmpty(usedFilterOrder_hp),replaceIfEmpty(hp_hdm.Fs),replaceIfEmpty(hp_hdm.Astop),replaceIfEmpty(hp_hdm.Fstop),replaceIfEmpty(hp_hdm.F6dB),replaceIfEmpty(hp_hdm.F3dB),replaceIfEmpty(hp_hdm.TransitionWidth),replaceIfEmpty(hp_hdm.Fpass),replaceIfEmpty(hp_hdm.Apass),...
    {replaceIfEmpty(core_cfg.lpfilttype)},{replaceIfEmpty(lp_f_type_detail)},{replaceIfEmpty(core_cfg.lpfiltdir)},replaceIfEmpty(usedFilterOrder_lp),replaceIfEmpty(lp_hdm.Fs),replaceIfEmpty(lp_hdm.Astop),replaceIfEmpty(lp_hdm.Fstop),replaceIfEmpty(lp_hdm.TransitionWidth),replaceIfEmpty(lp_hdm.F3dB),replaceIfEmpty(lp_hdm.F6dB),replaceIfEmpty(lp_hdm.Fpass),replaceIfEmpty(lp_hdm.Apass),...
    'VariableNames',{...
    'hp_preDS_filter','hp_preDS_filter_type','hp_preDS_dir_and_passing','usedFilterOrder_hp_preDS','hp_preDS_Fs_Hz','hp_preDS_Astop_dB','hp_preDS_Fstop_Hz','hp_preDS_F6dB_Hz','hp_preDS_F3dB_Hz','hp_preDS_TransitionWidth_Hz','hp_preDS_Fpass_Hz','hp_preDS_Apass_dB',...
    'hp_filter','hp_filter_type','hp_dir_and_passing','usedFilterOrder_hp','hp_Fs_Hz','hp_Astop_dB','hp_Fstop_Hz','hp_F6dB_Hz','hp_F3dB_Hz','hp_TransitionWidth_Hz','hp_Fpass_Hz','hp_Apass_dB',...
    'lp_filter','lp_filter_type','lp_dir_and_passing','usedFilterOrder_lp','lp_Fs_Hz','lp_Astop_dB','lp_Fstop_Hz','lp_F6dB_Hz','lp_F3dB_Hz','lp_TransitionWidth_Hz','lp_Fpass_Hz','lp_Apass_dB'}...
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

ch_detectedMaxDownSlopess = zeros(nRowsCh,1);
ch_detectedZeroCrossingToTroughSlopess = zeros(nRowsCh,1);

ch_detectedMaxSlopess = zeros(nRowsCh,1);
ch_detectedTroughToZeroCrossingSlopess = zeros(nRowsCh,1);

ch_frequency_by_durations = zeros(nRowsCh,1);
ch_frequency_by_trough_to_peak_latencys = zeros(nRowsCh,1);

ch_meanNegativePeaks = zeros(nRowsCh,1);
ch_meanPeak2Peaks = zeros(nRowsCh,1);

ch_detectedSignalMins = zeros(nRowsCh,1);
ch_detectedSignalMaxs = zeros(nRowsCh,1);


epochlengths = repmat(cfg.scoring.epochlength,nRowsCh,1);

% ch_nMerged
lengthsAcrossROIsSecondss = repmat(lengthsAcrossROIsSeconds,nRowsCh,1);
% ch_SDfrqBndPssSignal

minFreqs = repmat(minFreq,nRowsCh,1);
maxFreqs = repmat(maxFreq,nRowsCh,1);

ch_detectedSDofFilteredSignals = zeros(nRowsCh,1);
    

hypnEpochsEndsSamples = cumsum(repmat(cfg.scoring.epochlength,numel(cfg.scoring.epochs),1));
hypnEpochsBeginsSamples = hypnEpochsEndsSamples - cfg.scoring.epochlength;

hypnEpochsEndsSamples = hypnEpochsEndsSamples + cfg.scoring.dataoffset;
hypnEpochsBeginsSamples = hypnEpochsBeginsSamples + cfg.scoring.dataoffset;


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

    
    ch_detectedMaxDownSlopess(iChan) = nanmean(ch_detectedPeak2Peaks{iChan});
    
    ch_detectedMaxDownSlopess(iChan) = nanmean(ch_detectedMaxDownSlopes{iChan});
    ch_detectedZeroCrossingToTroughSlopess(iChan) = nanmean(ch_detectedZeroCrossingToTroughSlopes{iChan});
    
    ch_detectedMaxSlopess(iChan) = nanmean(ch_detectedMaxSlopes{iChan});
    ch_detectedTroughToZeroCrossingSlopess(iChan) = nanmean(ch_detectedTroughToZeroCrossingSlopes{iChan});
    
    
    ch_frequency_by_durations(iChan) = (1/tempLengthMeanSeconds);
    tempTroughPeakLengthMean = nanmean((ch_detectedPeaksSamples{iChan} - ch_detectedTroughsSamples{iChan})');
    ch_frequency_by_trough_to_peak_latencys(iChan) = (1/(tempTroughPeakLengthMean*2));
    
        
    
    ch_meanNegativePeaks(iChan) = nanmean(ch_meanNegativePeak{iChan});
    ch_meanPeak2Peaks(iChan) = nanmean(ch_meanPeak2Peak{iChan});
    
    ch_detectedSignalMins(iChan) = nanmean(ch_detectedSignalMin{iChan});
    ch_detectedSignalMaxs(iChan) = nanmean(ch_detectedSignalMax{iChan});

    ch_detectedSDofFilteredSignals(iChan) = nanmean(ch_detectedSDofFilteredSignal{iChan});
 
    
    if ch_nDetected{iChan} > 0
        
%         output = cell(ch_nDetected{iChan},11);



            output(indEv{iChan},1) = cellstr(repmat(ch, ch_nDetected{iChan}, 1));
%             output(indEv{iChan},2) = num2cell(ch_detectedLengthSamples{iChan}');
%             output(indEv{iChan},3) = num2cell(ch_detectedBeginSample{iChan}');
%             output(indEv{iChan},4) = num2cell(ch_detectedEndSample{iChan}');
%             output(indEv{iChan},5) = num2cell(ch_detectedPeaksSamples{iChan}');
%             output(indEv{iChan},6) = num2cell(ch_detectedTroughsSamples{iChan}');
            
            output(indEv{iChan},2) = num2cell(ch_detectedPeak2Peaks{iChan});
            
            tempTroughPeakLength = (ch_detectedPeaksSamples{iChan}' - ch_detectedTroughsSamples{iChan}');
            
            output(indEv{iChan},3) = num2cell(1./tempTroughPeakLength);
            output(indEv{iChan},4) = num2cell((1./(tempTroughPeakLength*2)));
            
            output(indEv{iChan},5) = cellstr(repmat(strjoin(cfg.stages,' '), ch_nDetected{iChan}, 1));
            output(indEv{iChan},6) = num2cell(repmat(ch_meanNegativePeak{iChan}, ch_nDetected{iChan}, 1));
            output(indEv{iChan},7) = num2cell(repmat(ch_meanPeak2Peak{iChan}, ch_nDetected{iChan}, 1));
            
            output(indEv{iChan},8) = num2cell(ch_detectedLengthSamples{iChan}');
            output(indEv{iChan},9) = num2cell(ch_detectedBeginSample{iChan}');
            output(indEv{iChan},10) = num2cell(ch_detectedEndSample{iChan}');
            output(indEv{iChan},11) = num2cell(ch_detectedPeaksSamples{iChan}');
            output(indEv{iChan},12) = num2cell(ch_detectedTroughsSamples{iChan}');
            
            output(indEv{iChan},13) = num2cell((1:ch_nDetected{iChan})');
            
            output(indEv{iChan},14) = num2cell(ch_detectedMaxSlopes{iChan}');
            output(indEv{iChan},15) = num2cell(ch_detectedTroughToZeroCrossingSlopes{iChan}');
            output(indEv{iChan},16) = num2cell(ch_detectedMaxSlopesSamples{iChan}');
            output(indEv{iChan},17) = num2cell(ch_detectedTroughToZeroCrossingSlopesSamples{iChan}');
            
            
            output(indEv{iChan},18) = num2cell(ch_detectedSignalMin{iChan}');
            output(indEv{iChan},19) = num2cell(ch_detectedSignalMax{iChan}');
            output(indEv{iChan},20) = [epochs{:,1}];
            output(indEv{iChan},21) = [epochs{:,2}];
            output(indEv{iChan},22) = [epochs{:,3}];
            output(indEv{iChan},23) = num2cell(ch_detectedSDofFilteredSignal{iChan});
            
            output(indEv{iChan},24) = num2cell(ch_detectedMaxDownSlopes{iChan}');
            output(indEv{iChan},25) = num2cell(ch_detectedZeroCrossingToTroughSlopes{iChan}');
            output(indEv{iChan},26) = num2cell(ch_detectedMaxDownSlopesSamples{iChan}');
            
            
            output(indEv{iChan},27) = num2cell(ch_contigSegment{iChan});
            output(indEv{iChan},28) = cellstr(repmat(strjoin(cfg.thresholdstages,' '), ch_nDetected{iChan}, 1));

            
%             for iLine=1:(size(output,1))
%                 fprintf(fide,'%i,',iData);
%                 fprintf(fide,'%s,',output{iLine,1});
%                 tempLengthSeconds = output{iLine,13};
%                 fprintf(fide,'%f,',tempLengthSeconds);
%                 fprintf(fide,'%e,',output{iLine,7});
%                 fprintf(fide,'%e,',output{iLine,31});
%                 fprintf(fide,'%e,',output{iLine,32});
%                 fprintf(fide,'%e,',output{iLine,19});
%                 fprintf(fide,'%e,',output{iLine,20});
%                 tempTroughPeakLength = (output{iLine,16} - output{iLine,17});
%                 fprintf(fide,'%f,', (1/tempLengthSeconds));
%                 fprintf(fide,'%f,', (1/(tempTroughPeakLength*2)));
%                 fprintf(fide,'%i,',output{iLine,2});
%                 fprintf(fide,'%i,',output{iLine,3});
%                 fprintf(fide,'%i,',output{iLine,4});
%                 fprintf(fide,'%i,',output{iLine,5});
%                 fprintf(fide,'%i,',output{iLine,6});
%                 fprintf(fide,'%s,',output{iLine,8});
%                 fprintf(fide,'%s,',output{iLine,9});
%                 fprintf(fide,'%s,',output{iLine,10});
%                 fprintf(fide,'%e,',output{iLine,11});
%                 fprintf(fide,'%e,',output{iLine,12});
%                 
%                 
%                 fprintf(fide,'%f,',output{iLine,14});
%                 fprintf(fide,'%f,',output{iLine,15});
%                 
%                 fprintf(fide,'%f,',output{iLine,16});
%                 fprintf(fide,'%f,',output{iLine,17});
%                 fprintf(fide,'%i,',output{iLine,18});
% 
%                 fprintf(fide,'%i,',output{iLine,33});                
%                 fprintf(fide,'%i,',output{iLine,21});
%                 fprintf(fide,'%i,',output{iLine,22});
%                 
%                 fprintf(fide,'%f,',output{iLine,34});
%                 fprintf(fide,'%f,',output{iLine,23});
%                 fprintf(fide,'%f,',output{iLine,24});
%                 
%                 fprintf(fide,'%e,',output{iLine,25});
%                 fprintf(fide,'%e,',output{iLine,26});
%                 
%                 fprintf(fide,'%s,',output{iLine,27});
%                 fprintf(fide,'%s,',output{iLine,28});
%                 fprintf(fide,'%s,',output{iLine,29});
%                 
%                 fprintf(fide,'%i,',output{iLine,34});
% 
%                 
%                 fprintf(fide,'%e\n',output{iLine,30});
%             end













    end
    
    
end

res_channel = [];
res_channel.ori = functionname;
res_channel.type = 'slowwaves_channel';
res_channel.cfg = cfg;

tempvarnames = {...
        'channel',...
        'count','density_per_epoch',...
        'mean_duration_seconds','mean_amplitude_peak2trough_potential',...
        'mean_slope_to_trough_min_potential_per_second','mean_slope_zerocrossing_to_trough_potential_per_second',...
        'mean_slope_trough_to_up_max_potential_per_second','mean_slope_trough_to_zerocrossing_potential_per_second',...
        'mean_frequency_by_duration','mean_frequency_by_trough_to_peak_latency',...
        'epoch_length_seconds','lengths_ROI_seconds',...
        'used_stages_for_detection','used_stages_for_thresholds','used_meanNegPeak_potential','used_meanP2T_potential',...
        'mean_SD_of_filtered_signal',...
        'mean_negative_peak_potential','mean_positive_peak_potential',...
        'used_min_freq','used_max_freq'};

% if isempty(output)
%     res_channel.table = cell2table(cell(0,numel(tempvarnames)), 'VariableNames', tempvarnames);
% else
    res_channel.table = table(...
        chs,...
        [ch_nDetected{:}]', [ch_densityPerEpoch{:}]', ...
        ch_detectedLengthSampless,ch_detectedPeak2Peakss,...
        ch_detectedMaxDownSlopess,ch_detectedZeroCrossingToTroughSlopess,...
        ch_detectedMaxSlopess,ch_detectedTroughToZeroCrossingSlopess,...
        ch_frequency_by_durations,ch_frequency_by_trough_to_peak_latencys,...
        epochlengths,lengthsAcrossROIsSecondss,...
        cellstr(repmat(strjoin(cfg.stages,' '), nRowsCh, 1)),cellstr(repmat(strjoin(cfg.thresholdstages,' '), nRowsCh, 1)),ch_meanNegativePeaks,ch_meanPeak2Peaks,...
        ch_detectedSDofFilteredSignals,...
        ch_detectedSignalMins,ch_detectedSignalMaxs,...
        minFreqs,maxFreqs,...
        'VariableNames',tempvarnames);
% end


res_event = [];
res_event.ori = functionname;
res_event.type = 'slowwaves_event';
res_event.cfg = cfg;
tempvarnames = {...
        'channel','duration_seconds','amplitude_peak2trough_max',...
        'slope_to_trough_min_potential_per_second','slope_zeroxing_to_trough_potential_per_second',...
        'slope_trough_to_up_max_potential_per_second','slope_trough_to_zeroxing_potential_per_second',...
        'frequency_by_duration','frequency_by_trough_to_peak_latency',...
        'used_stages_for_detection','used_stages_for_thresholds','used_meanNegPeak','used_meanP2T',...
        'seconds_begin','seconds_end','seconds_peak_max','seconds_trough_max',...
        'id_within_channel',...
        'seconds_slope_to_trough_min','seconds_slope_trough_to_up_max',...
        'seconds_zeroxing_of_trough_to_zeroxing_slope_trough_to_up_max',...
        'negative_peak_potential','positive_peak_potential',...
        'stage','stage_alt','stage_alt2',...
        'contig_segment',...
        'SD_of_filtered_signal'};

if isempty(output)
    res_event.table = cell2table(cell(0,numel(tempvarnames)), 'VariableNames', tempvarnames);
else
    res_event.table = table(...
        output(:,1),[output{:,8}]',[output{:,2}]',...
        [output{:,24}]',[output{:,25}]',...
        [output{:,14}]',[output{:,15}]',...
        [output{:,3}]',[output{:,4}]',...
        output(:,5),output(:,28),[output{:,6}]',[output{:,7}]',...
        [output{:,9}]',[output{:,10}]',[output{:,11}]',[output{:,12}]',...
        [output{:,13}]',...
        [output{:,26}]',[output{:,16}]',...
        [output{:,17}]',...
        [output{:,18}]',[output{:,19}]',...
        output(:,20),output(:,21),output(:,22),...
        [output{:,27}]',...
        [output{:,23}]',...
        'VariableNames',tempvarnames);
end
data = [];%clear
chData = [];%clear


fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
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
