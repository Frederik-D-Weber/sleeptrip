function cfg_detector_set=st_get_default_detector_set(cfg)

% ST_GET_DEFAULT_DETECTOR_SET returns a default "detector set" (=collection of artifact detectors),
% which can be supplied to st_run_detector_set.
%
% Use as:
%     detector_set=st_get_default_detector_set(cfg)
%
% Required configuration parameters:
%     cfg.elec      = structure, containing channel information (FieldTrip format)
%     cfg.fsample   = numeric, sample rate of data (should match data.fsample)
%
% Optional configuration parameters:
%     cfg.include   = cell array of strings, with names of each detector to
%     include. (default: 'all' [string])
%
% The following detector names are recognized:
%   'absamp'
%   'lowfreq'
%   'highfreq'
%   'jump'
%   'flatline'
%   'deviant'
%   'similar'
%
% Output:
%     cfg_detector_set  = output cfg, containing details of individual
%     detectors
%
% See also ST_RUN_DETECTOR_SET ST_COMBINE_DETECTORS

% Copyright (C) 2022-, Roy Cox, Frederik D. Weber
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

%TO DO:
% - consider if filter orders can be based on fsample instead of defaults


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
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
    return
end

%----core function starts


%---input checks and defaults----
ft_checkconfig(cfg,'required',{'elec','fsample'});
cfg.include  = ft_getopt(cfg, 'include', 'all');

fprintf([functionname ' function initialized\n']);

fsample=cfg.fsample;

%initialize empty cell to populate individual detectors
all_detectors_cfg={};
%%
%--absolute amplitude detector
cfg_detector=[];
cfg_detector.label='highamp';

%st
cfg_detector.st.method='threshold';
cfg_detector.st.abs='yes';
cfg_detector.st.thresholddirection='above';
cfg_detector.st.thresholdvalue=300;
cfg_detector.st.paddingduration = 0.1; %extend/pad (in seconds). minor extension before/after event
cfg_detector.st.mergeduration=1; %merge artifacts if within this range (in seconds)

all_detectors_cfg{end+1}=cfg_detector;
%%
%--low frequency detector
freqRange=[0.3 15];

cfg_detector=[];
cfg_detector.label='lowfreq';

%ft
cfg_detector.ft.bpfilter     = 'yes';
cfg_detector.ft.bpfreq       = freqRange;
cfg_detector.ft.hilbert      = 'abs';
cfg_detector.ft.boxcar       = 0.2; %default taken from ft

%st
cfg_detector.st.method='threshold';
cfg_detector.st.zscore='yes';
cfg_detector.st.thresholddirection='above';
cfg_detector.st.thresholdvalue =8; %artifact threshold: > maxZ. Note: lower value risks detecting slow waves
cfg_detector.st.paddingduration = 3; %extend/pad (in seconds). extend a bit to capture slow activity before/after event
cfg_detector.st.mergeduration=1; %merge artifacts if within this range (in seconds)

all_detectors_cfg{end+1}=cfg_detector;

%%
%--high frequency detector
freqRange=[60 120];
minHigherBand=80;%lower bandpass edge is 60 Hz, so upper edge of at least 80 Hz for meaningful bandpass
freqNyq=fsample/2; %Nyquist freqency

if freqNyq>=minHigherBand
    cfg_detector=[];
    cfg_detector.label='highfreq';

    %ft
    cfg_detector.ft.bpfilter     = 'yes';
    cfg_detector.ft.bpfreq       = [freqRange(1) min([freqNyq-1 freqRange(2)])]; %60 to (120 or Nyquist, whichever is lower). -1 Hz to stay away from actual Nyquist (filter fails)
    cfg_detector.ft.hilbert      = 'abs';
    cfg_detector.ft.boxcar       = 0.2; %default taken from ft

    %st
    cfg_detector.st.method='threshold';
    cfg_detector.st.zscore='yes';
    cfg_detector.st.thresholddirection='above';
    cfg_detector.st.thresholdvalue =3; %artifact threshold: > maxZ.
    cfg_detector.st.paddingduration = 0.1; %extend/pad (in seconds). minor extension before/after event
    cfg_detector.st.mergeduration=1; %merge artifacts if within this range (in seconds)

    all_detectors_cfg{end+1}=cfg_detector;
else
    fprintf('creation of high-frequency detector failed: sample rate should be at least %i Hz\n',minHigherBand*2)
end
%%
%--jump detector
cfg_detector=[];
cfg_detector.label='jump';

%ft
cfg_detector.ft.medianfilter  =  'yes'; %will attempt to use compiled C function (faster)
cfg_detector.ft.medianfiltord = 9; %default from ft (to do: consider adjusting to sample rate)
cfg_detector.ft.absdiff= 'yes'; % computes abs(diff(data))

%st
cfg_detector.st.method='threshold';
cfg_detector.st.zscore='yes';
cfg_detector.st.thresholddirection='above';
cfg_detector.st.thresholdvalue =25; %artifact threshold: > maxZ.
cfg_detector.st.paddingduration = 0.1; %extend/pad (in seconds). minor extension before/after event
cfg_detector.st.mergeduration=1; %merge artifacts if within this range (in seconds)

all_detectors_cfg{end+1}=cfg_detector;
%%
%--flatline detector (also smooth)
cfg_detector=[];
cfg_detector.label='flatline';

%st
cfg_detector.st.method='threshold';
cfg_detector.st.diff='yes';
cfg_detector.st.abs='yes';
cfg_detector.st.thresholddirection='below';
cfg_detector.st.thresholdvalue = 1*250/fsample; %artifact threshold: < thresholdvalue. For a sample rate of 250 Hz, max allowed sample diff is 1.Assuming data in microV (was 0.5. changed to 1 for zmax)
cfg_detector.st.minduration = 1; %min artifact duration (in seconds)
cfg_detector.st.paddingduration = 0.1; %extend/pad (in seconds). minor extension before/after event
cfg_detector.st.mergeduration=1; %merge artifacts if within this range (in seconds)

all_detectors_cfg{end+1}=cfg_detector;
%%
%--low amplitude detector
cfg_detector=[];
cfg_detector.label='lowamp';

%st
cfg_detector.st.method='threshold';
cfg_detector.st.abs='yes';
cfg_detector.st.thresholddirection='below';
cfg_detector.st.thresholdvalue = 10; %artifact threshold: < thresholdvalue. Assuming data in microV
cfg_detector.st.minduration = 30; %min artifact duration (in seconds)
cfg_detector.st.paddingduration = 0.1; %extend/pad (in seconds). minor extension before/after event
cfg_detector.st.mergeduration=1; %merge artifacts if within this range (in seconds)

all_detectors_cfg{end+1}=cfg_detector;
%%
%--deviation from neighbors detector
cfg_detector=[];
cfg_detector.label='deviant';

%st
cfg_detector.st.method='neighbours';
cfg_detector.st.metric='correlation';
cfg_detector.st.elec=cfg.elec;
neighbours=st_get_default_neighbours(cfg);
cfg_detector.st.neighbours=neighbours;


cfg_detector.st.thresholddirection='below';
cfg_detector.st.thresholdvalue  = 0.3; % artifact threshold: < corrthresh
cfg_detector.st.channelthreshold =0.5; % artifact threshold: > chanpropthresh

cfg_detector.st.mergeduration=60; %merge artifacts if within this range (in seconds). Note: long merge window to avoid channels switching between good/bad every window.

if ~isempty(neighbours)
    all_detectors_cfg{end+1}=cfg_detector;
else
    ft_warning('creation of detector deviant failed: no neighbourhood structure\n')

end


%%
%--similar to neighbors detector
cfg_detector=[];
cfg_detector.label='similar';

%st
cfg_detector.st.method='neighbours';
cfg_detector.st.metric='maxabsdiff';
cfg_detector.st.elec=cfg.elec;
neighbours=st_get_default_neighbours(cfg);
cfg_detector.st.neighbours=neighbours;

cfg_detector.st.thresholddirection='below';
cfg_detector.st.thresholdvalue  = 0.5; %artifact if diff wave < 0.5 microV
cfg_detector.st.channelthreshold =0; % artifact if fraction of neighboring channels meeting thresholdvalue > 0 (= single channel pair suffices)

cfg_detector.st.mergeduration=60; %merge artifacts if within this range (in seconds). Note: long merge window to avoid channels switching between good/bad every window.

if ~isempty(neighbours)
    all_detectors_cfg{end+1}=cfg_detector;
else
    ft_warning('creation of detector similar failed: no neighbourhood structure\n')

end

%%

%select requested (and possible) detectors
if ~strcmp(cfg.include,'all') 
   
    %determine which ones we have
    [detectorAvailable,idx]=ismember(cfg.include,cellfun(@(X) X.label,all_detectors_cfg,'UniformOutput',false));
    all_detectors_cfg=all_detectors_cfg(idx(detectorAvailable));
end

%function for combining individual detector cfgs
cfg_detector_set=st_combine_detectors(all_detectors_cfg);

%add elec struct
cfg_detector_set.elec=cfg.elec;

fprintf('created detector set of %i detectors:\n%s\n',cfg_detector_set.number,strjoin(cfg_detector_set.label))
%--core function end---

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