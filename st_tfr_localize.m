function [time_freq_pow_centroids] = st_tfr_localize(cfg, freq)

% ST_TFR_LOCALIZE finds the power/value maximum poin in a tfr map (as from FT_FREQANALYSIS)
% or ST_CHANNEL_EVENT_TFR
%
% Use as
%   [time_freq_pow_centroids] = st_tfr_localize(cfg, freq)
%
% Optional configuration parameters are:
%
% cfg.timewin = a two value numeric vector for the time boundaries to
%               localize in e.g. e.g. [0 0.75]; (default = [min(freq.time) max(freq.time)])
% cfg.freqwin = a two value numeric vector for the frequency boundaries to
%               localize in  e.g.  [11 16] (default = [min(freq.freq) max(freq.freq)])
% 
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
st_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar freq
ft_preamble provenance freq
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% set defaults
cfg.timewin  = ft_getopt(cfg, 'timewin', [min(freq.time) max(freq.time)]);
cfg.freqwin  = ft_getopt(cfg, 'freqwin', [min(freq.freq) max(freq.freq)]);
                                    
fprintf([functionname ' function initialized\n']);

nChannel = numel(freq.label);
time_freq_pow_centroids = nan(nChannel,3);

for iCh = 1:nChannel
    ch = freq.label{iCh};
    cfg_sd = [];
    cfg_sd.latency    = cfg.timewin;
    cfg_sd.frequency  = cfg.freqwin;
    cfg_sd.channel  = ch;
    freq_ch_localize = ft_selectdata(cfg_sd,freq);
    powermax = max(freq_ch_localize.powspctrm(:));
    pos_cluster_mat = ([freq_ch_localize.powspctrm]>=powermax);
    pos_cluster_pow = freq_ch_localize.powspctrm(pos_cluster_mat);
    matclust = squeeze(pos_cluster_mat);
    matclust=matclust/sum(matclust(:));
    [m,n]=size(matclust);
    [I,J]=ndgrid(1:m,1:n);
    centroid=[dot(I(:),matclust(:)),  dot(J(:),matclust(:))];
    centroid = round(centroid);
    time_freq_pow_centroids(iCh,:) = [freq_ch_localize.time(centroid(2)), freq_ch_localize.freq(centroid(1)) powermax];
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
% ft_postamble previous data
% data = datanew
% ft_postamble provenance data
% ft_postamble history    data
% ft_postamble savevar    data

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end
