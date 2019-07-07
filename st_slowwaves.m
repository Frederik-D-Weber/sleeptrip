function [res_channel res_events] = st_spindles(cfg, data)

% ST_SPINDLES detect sleep spindles and their properties like their count, density, amplitude, duration, frequency etc. 
% results are stored it in a result structure
%
% Use as
%   [res_channel res_events] = st_spindles(cfg, data)
%   [res_channel res_events] = st_spindles(cfg, data)
%
% Required configuration parameters are:
%   cfg.scoring  = structure provided by ST_READ_SCORING
%   cfg.stages   = either a string or a Nx1 cell-array with strings that indicate sleep stages
%                  if no data structure is provided as input next to configuration also
%                  provide this parameter
%
%  if no data structure is provided as parameter then you need to define
%   cfg.dataset  = string with the filename
%
% Optional configuration parameters are:
%
%   cfg.segmentlength  = scalar, segment length in seconds,
%                        should be lower than the epoch length (default = 10)
%   cfg.segmentoverlap = fraction of the overlap of segments during
%                        segmentation in the interval [0,1) (default = 0.9)
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
%                          (default = 0.2)
%
% Some additional parameters from FT_PREPROCESSING can also be included
% including thje reprocessing options that you can only use for EEG data are:
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

if ~isfield(cfg, 'scoring')
    cfg.scoring = st_read_scoring(cfg);
end

% set defaults
cfg.channel  = ft_getopt(cfg, 'channel', 'all', 1);                                  

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

% set core parameters
load_core_cfg



res_bands = [];
res_bands.ori = functionname;
res_bands.type = 'power_full';
res_bands.cfg = cfg;
res_bands.table = table(...
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

res = [];
res.ori = functionname;
res.type = 'power_full';
res.cfg = cfg;
res.table = table(...
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



fprintf('st_power function finished\n');
toc
memtoc
end

