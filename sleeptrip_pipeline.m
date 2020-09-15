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


%% minimal example pipeline for a sleeptrip analysis

%%% PLEASE MAKE SURE:
%%% 1. YOU HAVE Matlab2013b or later
%%% 2. YOU HAVE ALL REQUIRED TOOLBOXES (typically the case)
%%% 3. YOUR SLEEPTRIP FOLDER DIRECTLY HAS THE FILES IN THEM, 
%%%    e.g. while extracting a zip from github maybe you need to refer to 
%%%    D:/sleeptrip-master/sleeptrip-master instead of D:/sleeptrip-master
%%% 4. NO FieldTrip NOR related toolboxes like EEGLAB are loaded or in the 
%%%    PATH of your MATLAB
%%% 5. You have downloaded the example datasets here:
%%%    https://drive.google.com/open?id=160i_FG3zSZrBf9mr2zxZQX-u30OnLDMf
%%% 6. ...and the files from the example datasets are in the 
%%%    same folder as this script or somewhere in the matlab path
   


%% initialization

%%% add sleeptrip to the path
%addpath('D:/sleeptrip')
pathToSleepTrip = uigetdir('','choose spisop path, e.g. D:\sleeptrip');
addpath(pathToSleepTrip)

%%% disable some toolboxes if necessary and check if 
%%% 'signal_toolbox', 'signal_blocks' are available 
%%% because they are helpful to have.
% toggleToolbox('names')
% toggleToolbox('dsp','off')
% toggleToolbox('all','query')
% license('inuse')

%%% load all the defaults of SleepTrip and FieldTrip
st_defaults


%% prepare the subject data, e.g. zmax data in a compressed format

% here you can prepare everyting you would need later on again
subject = [];
subject.name               = 'example_zmax';
subject.dataset            = 'example_zmax.zip';
%subject.dataset            = 'example_zmax/BATT.edf';
subject.scoringfile        = 'example_zmax.csv';
subject.scoringformat      = 'zmax'; % e.g. 'zmax' or 'spisop'
subject.standard           = 'aasm'; % 'aasm' or 'rk'
% does the scoring start 
%   at the beginning of data (=0) or 
%   before (<0) or
%   after (>0)
subject.scoring_dataoffset = 0; % in seconds
% at which time in seconds was the lights off
% with respect to the beginning of scoring,
% only relevant for
subject.lightsoff          = 0; 
subject.eegchannels        = {'EEG L', 'EEG R'};
subject.ppgchannels        = {'OXY_IR_AC'};
subject.oxychannels        = {'OXY_IR_DC'};



save('subject-1','subject');

%% alternative subject brainvision

% subject = [];
% subject.name               = 'example_brainvision';
% subject.dataset            = 'test_lindev_sleep_score_preprocout_datanum_1.eeg';
% subject.scoringfile        = 'test_lindev_sleep_score_preprocout_datanum_1.txt';
% subject.scoringformat      = 'spisop'; % e.g. 'zmax' or 'spisop'
% subject.standard           = 'rk'; % 'aasm' or 'rk'
% % does the scoring start 
% %   at the beginning of data (=0) or 
% %   before (<0) or
% %   after (>0)
% subject.scoring_dataoffset = 0; % in seconds
% % at which time in seconds was the lights off
% % with respect to the beginning of scoring,
% % only relevant for
% subject.lightsoff          = 108.922;
% subject.eegchannels        = {'C3:A2', 'C4:A1'};
% % 
% save('subject-2','subject');
%
%% alternative subject somnomedics edf

% subject = [];
% subject.name               = 'example_somnomedics_edf';
% subject.dataset            = 'example_somnomedics.edf';
% subject.scoringfile        = 'example_zmax.txt';
% subject.scoringformat      = 'spisop'; % e.g. 'zmax' or 'spisop'
% subject.standard           = 'rk'; % 'aasm' or 'rk'
% % does the scoring start 
% %   at the beginning of data (=0) or 
% %   before (<0) or
% %   after (>0)
% subject.scoring_dataoffset = 0; % in seconds
% % at which time in seconds was the lights off
% % with respect to the beginning of scoring,
% % only relevant for
% subject.lightsoff          = 108.922;
% subject.eegchannels        = {'C3:A2', 'C4:A1'};
%
% % save('subject-3','subject');

%% read in the scoring, files 

cfg = [];
cfg.scoringfile   = subject.scoringfile;
cfg.scoringformat = subject.scoringformat;
cfg.standard      = subject.standard; % 'aasm' or 'rk'

[scoring] = st_read_scoring(cfg);

% do we know the lights-off moment, if yes then add it
% this is good for SleepTrip to know, e.g. for the st_sleepdescriptive
% function
scoring.lightsoff = subject.lightsoff;

% do we know the lights-on moment, in this case not, but if this would be
% good to determine more acurately where things at the end of the recoring
% are
%scoring.lightson = subject.lightson;

%practice: read in another format, maybe a custom format.

%% check when is sleep onset/offset

cfg = [];
cfg.sleeponsetdef  = 'AASM'; % also try 'AASM' and many more
[onsetnumber, lastsleepstagenumber, onsetepoch, lastsleepstage] = st_sleeponset(cfg,scoring);
onsetnumber
lastsleepstagenumber
onsetepoch
lastsleepstage

% practice: how does the sleep onset change according to the different
% definitions?

%% plot the scoring

cfg = [];
cfg.title           = subject.name;
cfg.plottype        = 'colorbar'; %'classic' 'colorblocks' 'colorbar'
% cfg.yaxdisteqi      = 'yes';
% cfg.timeticksdiff   = 60;
cfg.plotunknown     = 'no'; % 'yes' or 'no' 
% cfg.plotsleeponset  = 'no'; % 'yes' or 'no' 
% cfg.plotsleepoffset = 'no'; % 'yes' or 'no' 
% cfg.timemin         = 600 % in minutes, e.g. plot on a 10-hour time axis.   
% cfg.timerange          = [50 100];
% cfg.considerdataoffset = 'no';
% cfg.plotexcluded       = 'no';     

%%% if you want to export the figure immediately add these parameters
cfg.figureoutputfile   = [subject.name '.pdf'];
cfg.figureoutputformat = 'pdf';

figure_handle = st_hypnoplot(cfg, scoring);
% close(figure_handle) close the figure automatically e.g. after exporting

%practice: play around with the plotting paramters. Export the file as a
%pdf or a png, change the figure dimensions and resolution, manipulate the
%figure by using the figure_handle, what happens when you use an RK or AASM
%based scoring?

%% create a sleep table with the scoring parameters

cfg = [];
res_sleepdescriptive = st_scoringdescriptives(cfg, scoring);

cfg = [];
cfg. cycle = 'all';
res_sleepdescriptive_cycle = st_scoringdescriptives(cfg, scoring);


% lets take a look what we got in there
res_sleepdescriptive.table

% practice: how is sleep onset and lights-off moment related? 
% When is the end of the Total sleep time?
% Why are there several versions of the columns?

% export the results
cfg = [];
cfg.prefix = 'example';
cfg.infix  = subject.name;
cfg.posfix = '';
filelist_res_sleepdescriptives = st_write_res(cfg, res_sleepdescriptive); 

% note: you can write mutliple results with st_write_res(cfg, res_sleep1, res_sleep2, ...)
%       or if they are in a cell array with st_write_res(cfg, res_all{:})

%practice: were are the files stored, can this be changed, and how to use
%the prefix, infix and postfix best to differentiate your analysis from other runs/studies/subjects/recordings?


%% read in the data

cfg = [];
cfg.dataset    = subject.dataset;
cfg.continuous = 'yes';
cfg.channel    = subject.eegchannels;
%  possible channels for zmax are: 
% {'BATT' 'BODY TEMP' ...
%  'dX' 'dY' 'dZ' ...
%  'EEG L' 'EEG R' 
%  'LIGHT' 'NASAL L' 'NASAL R' 'NOISE' ...
%  'OXY_DARK_AC' 'OXY_DARK_DC' 'OXY_IR_AC' 'OXY_IR_DC' ...
%  'OXY_R_AC' 'OXY_R_DC' 'RSSI'}
cfg.lpfilter = 'yes';
cfg.lpfreq   = 35;
data = st_preprocessing(cfg);
% plot(data.trial{1}(1,1:(data.fsample*10)))

%practice, do some more preprocessing, how would you read in your own data
% search on the web how fieldtrip can help to read in data http://www.fieldtriptoolbox.org/faq/dataformat/
% Note you can read in .zip data directly, this saves typically one third
% of space your hard drive at least but needs longer to open

%% take a look at the data, MIGHT NOT WORK ON ALL MATLAB VERSIONS!

% get the indices of the respective epochs for later marking
epochs_W    = find(strcmp(scoring.epochs, 'W'));
epochs_N1   = find(strcmp(scoring.epochs, 'N1'));
epochs_N2   = find(strcmp(scoring.epochs, 'N2'));
epochs_N3   = find(strcmp(scoring.epochs, 'N3'));
epochs_R    = find(strcmp(scoring.epochs, 'R'));

% in the original data there are 30 seconds times the sample_rate samples per epoch
epochlenght_samples = scoring.epochlength*data.fsample;

% the first epoch is from sample 1 to sample 7680 in 256 Hz sample rate, etc.
artfctdef                  = [];
artfctdef.W.artifact  = [(epochs_W(:)  -1)*epochlenght_samples+1 (epochs_W(:)  +0)*epochlenght_samples];
artfctdef.N1.artifact = [(epochs_N1(:) -1)*epochlenght_samples+1 (epochs_N1(:) +0)*epochlenght_samples];
artfctdef.N2.artifact = [(epochs_N2(:) -1)*epochlenght_samples+1 (epochs_N2(:) +0)*epochlenght_samples];
artfctdef.N3.artifact = [(epochs_N3(:) -1)*epochlenght_samples+1 (epochs_N3(:) +0)*epochlenght_samples];
artfctdef.R.artifact  = [(epochs_R(:)  -1)*epochlenght_samples+1 (epochs_R(:)  +0)*epochlenght_samples];

cfg               = [];
cfg.continuous    = 'yes';
cfg.artfctdef     = artfctdef;
cfg.blocksize     = scoring.epochlength;
cfg.viewmode      = 'vertical'; % 'vertical' or 'butterfly'
cfg.artifactalpha = 0.7;
cfg.renderer      = 'opengl'; % 'painters' or 'opengl' or 'zbuffer'
ft_databrowser(cfg, data);


%% calculate power density

% power values in some sleep stages, e.g. non-REM
cfg = [];
cfg.scoring     = scoring;
cfg.stages      = {'N2', 'N3'}; % {'R'};
%cfg.stages      = {'S2', 'S3', 'S4'}; % {'R'};
cfg.channel     = subject.eegchannels;
% cfg.foilim      = [0.01 30];
% cfg.bands       = ...
%    {{'SO' ,0.5, 1},...
%     {'SWA', 0.5, 4},...
%     {'delta', 1, 4},...
%     {'theta', 4, 8},...
%     {'alpha' ,8, 11},...
%     {'slow_spindle', 9, 12},...
%     {'fast_spindle', 12, 15},...
%     {'spindle', 9, 15},...
%     {'beta', 15, 30}};

% this saves RAM, but is slower, especially with zmax data that is
% compressed (e.g. in a .zip file)
% cfg.dataset     = subject.dataset;
% [res_power_bins, res_power_bands] = st_power(cfg); 

[res_power_bin, res_power_band] = st_power(cfg, data);

%let us take a look
res_power_band.table

%export the data
cfg = [];
cfg.prefix = 'example';
cfg.infix  = subject.name;
cfg.posfix = '';
filelist_res_power = st_write_res(cfg, res_power_bin, res_power_band); % write mutliple results with  st_write_res(cfg, res_sleep1, res_sleep2)

%practice: use different frequency bands and ranges, what do you notice,
%what happens if the sleep stage is not present in the data? use the
%function without pre-loaded data structure and cfg.dataset

%% plot power density, example, only one channel

%which channel would you like to select, we use the first eeg channel
indChannel = strcmp(res_power_bin.table.channel,{subject.eegchannels{1}});
%indFreq = (res_power_bins.table.freq > 8) & (res_power_bins.table.freq < 20);
power_plot = figure;
freq_x     = res_power_bin.table.freq(indChannel);
%freq_x     = res_power_bins.table.freq(indChannel & indFreq);
power_y    = res_power_bin.table.mean_powerDensity_over_segments(indChannel);
%power_y    = res_power_bins.table.mean_powerDensity_over_segments(indChannel & indFreq);
% scale the power values on a logarithmic scale, in dB
% to avoid negative values add 1 to all values to have them >= 1;
power_y_logscale = 10*log10(power_y+1);
%plot(freq_x,power_y)
plot(freq_x,power_y_logscale)
plot(freq_x,power_y_logscale.*freq_x)


% practice: Try different scaling and normalization, how would this change
% the results? Can you plot muliple channels at the same time.


%% determine the spindle(s) frequency peak(s) from the power spectrum

cfg = [];
cfg.peaknum = 1; % either 1 or 2 (default)
[res_freqpeaks] = st_freqpeak(cfg,res_power_bin);

res_freqpeaks

%practice: find two peaks by default by changing the configuration of the function,
%would limiting the frequency band help, are there any artifacts (e.g. sharp peaks in the power spectra?)

%% detect sleep spindles

% define the frequency peak to use from the previous peak determination
% we will take the first peak we defined and the first value of that
freqpeak = res_freqpeaks.table.freqpeak1(1);

cfg = [];
cfg.scoring          =  scoring;
cfg.stages           = {'N2','N3'};
% cfg.stages           = {'S2', 'S3', 'S4'}; % {'R'};
cfg.channel          = subject.eegchannels;
% freqpeak1 for slow spindles
% freqpeak2 for fast spindles
cfg.centerfrequency  = freqpeak; % e.g. 12 for a frontal channel 13.3 for central channel

% this saves RAM, but is slower, especially with zmax data that is
% compressed (e.g. in a .zip file)
% cfg.dataset     = subject.dataset;
% [res_spindles_channel res_spindles_event res_spindles_filter] = st_spindles(cfg); 

[res_spindles_channel, res_spindles_event, res_spindles_filter] = st_spindles(cfg, data);
res_spindles_channel.table

%write out the results
cfg = [];
cfg.prefix = 'example';
cfg.infix  = subject.name;
cfg.posfix = '';
filelist_res_spindles = st_write_res(cfg, res_spindles_channel, res_spindles_event, res_spindles_filter); % write mutliple results with  st_write_res(cfg, res_sleep1, res_sleep2)

%practice: what if you change the threshold parameters, do you think
%filtering of the data matters? Search in different sleep stages, does the
%threshold change depending on the sleep stages used? if so, what to do
%then?

%% define events for further analysis

% begin and end
spindle_begins_subject.eegchannels{1}  = res_spindles_event.table.seconds_begin(strcmp(res_spindles_event.table.channel,{subject.eegchannels{1}}));
spindle_ends_subject.eegchannels{1}    = res_spindles_event.table.seconds_end(strcmp(res_spindles_event.table.channel,{subject.eegchannels{1}}));
spindle_begins_subject.eegchannels{2}  = res_spindles_event.table.seconds_begin(strcmp(res_spindles_event.table.channel,{subject.eegchannels{2}}));
spindle_ends_subject.eegchannels{2}    = res_spindles_event.table.seconds_end(strcmp(res_spindles_event.table.channel,{subject.eegchannels{2}}));

% the trough with the largest amplitude of the spindle shall give its time
% point, this is important for time-locked event related potentials
spindle_troughs_subject.eegchannels{1} = res_spindles_event.table.seconds_trough_max(strcmp(res_spindles_event.table.channel,{subject.eegchannels{1}}));
spindle_troughs_subject.eegchannels{2} = res_spindles_event.table.seconds_trough_max(strcmp(res_spindles_event.table.channel,{subject.eegchannels{2}}));

% we can also get the amplitudes
spindle_amplitude_subject.eegchannels{1} = res_spindles_event.table.amplitude_peak2trough_max(strcmp(res_spindles_event.table.channel,{subject.eegchannels{1}}));
spindle_amplitude_subject.eegchannels{2} = res_spindles_event.table.amplitude_peak2trough_max(strcmp(res_spindles_event.table.channel,{subject.eegchannels{2}}));

% ... or the frequency of each sleep spindle
spindle_frequency_subject.eegchannels{1} = res_spindles_event.table.frequency_by_mean_pk_trgh_cnt_per_dur(strcmp(res_spindles_event.table.channel,{subject.eegchannels{1}}));
spindle_frequency_subject.eegchannels{2} = res_spindles_event.table.frequency_by_mean_pk_trgh_cnt_per_dur(strcmp(res_spindles_event.table.channel,{subject.eegchannels{2}}));


%% plot spindle events on a hypnogram
cfg = [];
cfg.plotunknown  = 'no'; 
cfg.eventtimes = {spindle_troughs_subject.eegchannels{1}';...
                  spindle_troughs_subject.eegchannels{2}'};
cfg.eventlabels = {['spd ' subject.eegchannels{1}], ['spd ' subject.eegchannels{2}]};
figure_handle = st_hypnoplot(cfg, scoring);

% also plot event properties like amplitude and frequency for each event
cfg = [];
cfg.plotunknown  = 'no'; 
cfg.eventtimes  = {spindle_troughs_subject.eegchannels{1}';...
                   spindle_troughs_subject.eegchannels{2}';...
                   spindle_troughs_subject.eegchannels{1}';...
                   spindle_troughs_subject.eegchannels{2}'};
cfg.eventvalues = {spindle_amplitude_subject.eegchannels{1}';...
                   spindle_amplitude_subject.eegchannels{2}';...
                   spindle_frequency_subject.eegchannels{1}';...
                   spindle_frequency_subject.eegchannels{2}'};
cfg.eventranges = {[min(spindle_amplitude_subject.eegchannels{1}), max(spindle_amplitude_subject.eegchannels{1})];...
                   [min(spindle_amplitude_subject.eegchannels{2}), max(spindle_amplitude_subject.eegchannels{2})];...
                   [min(spindle_frequency_subject.eegchannels{1}), max(spindle_frequency_subject.eegchannels{1})];...
                   [min(spindle_frequency_subject.eegchannels{2}), max(spindle_frequency_subject.eegchannels{2})]};
cfg.eventlabels = {['spd ' 'ampl ' subject.eegchannels{1}], ['spd ' 'ampl ' subject.eegchannels{2}], ...
                   ['spd ' 'freq ' subject.eegchannels{1}], ['spd ' 'freq ' subject.eegchannels{2}]};
figure_handle = st_hypnoplot(cfg, scoring);

%% check events in the data browser, MIGHT NOT WORK ON YOUR MATLAB VERSION

cfg                                 = [];
% we will only look at one channel
cfg.channel                         = subject.eegchannels;
cfg.continuous                      = 'yes';
cfg.viewmode                        = 'vertical'; % 'vertical' or 'butterfly'
cfg.blocksize                       = scoring.epochlength; %view the data in 30-s blocks
%cfg.event                           = struct('type', {}, 'sample', {});


artfctdef                  = [];
artfctdef.spindles_ch1.artifact   = fix(data.fsample*[spindle_begins_subject.eegchannels{1}, spindle_ends_subject.eegchannels{1}]);
%artfctdef.spindles_ch1_trough.artfctpeak = fix(data.fsample*spindle_troughs_subject.eegchannels{1});
artfctdef.spindles_ch2.artifact   = fix(data.fsample*[spindle_begins_subject.eegchannels{2}, spindle_ends_subject.eegchannels{2}]);
%artfctdef.spindles_ch2_trough.artfctpeak = fix(data.fsample*spindle_troughs_subject.eegchannels{2});
% cfg.selectmode                        =  'markpeakevent'; %'markartifact', 'markpeakevent', 'marktroughevent' (default = 'markartifact')
cfg.artfctdef = artfctdef;
cfg.plotevents    = 'yes';
cfg.renderer      = 'opengl'; % 'painters' or 'opengl' or 'zbuffer'
ft_databrowser(cfg, data);


%% plot a event related potentials (ERP) and frequency (ERF) of e.g. spindles
% a buffer we need to have padding left and right to make nice
% time-frequency graph later on
% we look only in the first EEG channel
event_minimum_samples_ch1 = spindle_troughs_subject.eegchannels{1}*data.fsample;
event_minimum_samples_ch2 = spindle_troughs_subject.eegchannels{2}*data.fsample;

padding_buffer = 4*data.fsample; % 8 seconds

cfg     = [];

cfg.trl = [event_minimum_samples_ch1-data.fsample-padding_buffer,...
           event_minimum_samples_ch1+data.fsample+padding_buffer,...
           repmat(-(data.fsample+padding_buffer),numel(event_minimum_samples_ch1),1)];

% round trials to integers
cfg.trl = round(cfg.trl);
       
% remove trials that overlap with the beginning of the file
sel = cfg.trl(:,1)>1;
cfg.trl = cfg.trl(sel,:);
       
% remove trials that overlap with the end of the file
sel = cfg.trl(:,2)<data.sampleinfo(2);
cfg.trl = cfg.trl(sel,:);


data_events_ch1 = ft_redefinetrial(cfg, data);

cfg.trl = [event_minimum_samples_ch2-data.fsample-padding_buffer,...
           event_minimum_samples_ch2+data.fsample+padding_buffer,...
           repmat(-(data.fsample+padding_buffer),numel(event_minimum_samples_ch2),1)];
       
% round trials to integers
cfg.trl = round(cfg.trl);
       
% remove trials that overlap with the beginning of the file
sel = cfg.trl(:,1)>1;
cfg.trl = cfg.trl(sel,:);
       
% remove trials that overlap with the end of the file
sel = cfg.trl(:,2)<data.sampleinfo(2);
cfg.trl = cfg.trl(sel,:);

data_events_ch2 = ft_redefinetrial(cfg, data);


%View the event average signal timelocked to the trough.
cfg        = [];
[timelock_ch1] = ft_timelockanalysis(cfg, data_events_ch1);
[timelock_ch2] = ft_timelockanalysis(cfg, data_events_ch2);

%filter the timelocked data in the slow band to make underlying slow wave
%components more visible.
cfg = [];
% cfg.hpfilter = 'yes';
% cfg.hpfreq = 0.5;
cfg.lpfilter = 'yes';
cfg.lpfreq = 3;
[timelock_ch1_SWA_fitlered] = st_preprocessing(cfg, timelock_ch1);
[timelock_ch2_SWA_fitlered] = st_preprocessing(cfg, timelock_ch2);

% view the timelocked signals for both channels
figure
cfg        = [];
cfg.xlim   = [-1.5 1.5];
cfg.title  = 'Non-REM event ERP time-locked to trough';
cfg.graphcolor = 'brbr';
cfg.linestyle = {'-','-','-.','-.'};
cfg.linewidth = 1;
ft_singleplotER(cfg,timelock_ch1,timelock_ch2,timelock_ch1_SWA_fitlered,timelock_ch2_SWA_fitlered)

% ERF
cfg           = [];
cfg.channel   = subject.eegchannels{1};
cfg.method    = 'wavelet';
cfg.length    = 4;
cfg.foi       = 1:0.5:16; % 0.5 Hz steps
cfg.toi       = (-padding_buffer-1.5):0.1:(1.5+padding_buffer); % 0.1 s steps
event_freq_ch1 = ft_freqanalysis(cfg, data_events_ch1);
event_freq_ch2 = ft_freqanalysis(cfg, data_events_ch2);


% view the time-frequency of a slow wave or spindle event, for each channel
cfg                = [];
cfg.baseline       = [-1.5 1.5]; % a 3 s baseline around the event as it has no clear start or end.
cfg.baselinetype   = 'normchange';
cfg.zlim           = [-0.2 0.2];
cfg.xlim           = [-1.5 1.5];
cfg.title          = 'Event, time-frequency';
figure
ft_singleplotTFR(cfg,event_freq_ch1);
figure
ft_singleplotTFR(cfg,event_freq_ch2);

%practice, why did we extract more data then we needed? Does this change if
% we detect sleep spindles in another way? are spindles riding on slow
% waves?

%% cut the scoring in fixed length parts

cfg = [];
cfg.start = 1;
cfg.end = st_sleeponset(cfg,scoring)-1;
scorings = st_cutscoring(cfg,scoring);

cfg = [];
cfg.timemin = numel(scoring.epochs)*scoring.epochlength/60 + scoring.dataoffset/60;
cfg.plotsleeponset = 'no' ;  
cfg.plotsleepoffset = 'no' ;      
figure_handle = st_hypnoplot(cfg, scorings{:});

cfg = [];
cfg.length = 120;
[scorings] = st_cutscoring(cfg,scoring);

cfg = [];
cfg.timemin = numel(scoring.epochs)*scoring.epochlength/60 + scoring.dataoffset/60;
cfg.plotsleeponset = 'no' ;  
cfg.plotsleepoffset = 'no' ;      
for s = scorings
    figure_handle = st_hypnoplot(cfg, s{:});
end

%practice, make overlapping parts, can this be useful for futher analysis?

%% cut by sleep cycles

%full sleep cycles
cfg = [];
res_cycle = st_sleepcycles([],scoring);
cfg.start = res_cycle.table.startepoch;
cfg.end = res_cycle.table.endepoch;
scoring_cycles = st_cutscoring(cfg,scoring);

cfg = [];
cfg.timemin = numel(scoring.epochs)*scoring.epochlength/60 + scoring.dataoffset/60;
cfg.plotsleeponset = 'no' ;  
cfg.plotsleepoffset = 'no' ;      
for s = scoring_cycles
    figure_handle = st_hypnoplot(cfg, s{:});
end

%what remains? non-REM?
cfg = [];
res_cycle = st_sleepcycles([],scoring);
cfg.start = res_cycle.table.NRstartepoch;
cfg.end = res_cycle.table.NRendepoch;
[scoring_remains] = st_cutscoring(cfg,scoring);

cfg = [];
cfg.timemin = numel(scoring.epochs)*scoring.epochlength/60 + scoring.dataoffset/60;
cfg.plotsleeponset = 'no' ;  
cfg.plotsleepoffset = 'no' ;      
for s = scoring_remains
    fh = st_hypnoplot(cfg, s{:});
end

%practice: use the first sleep cycle to run a st_power on it, why would
%such an analysis be useful?

%% append results of the same type

[res_power_bins_appended] = st_append_res(res_power_bin,res_power_bin);
%note there is a new column in the table!
res_power_bins_appended.table

%practice: append some more results, how would we append results from a
%structure, could 
% res = {res1, res2, ...} be used with 
% res_appended = st_append_res(res{:}) ???


%% with many subjects, do an analysis, first prepare the subjects

datasets = {...
    'test_lindev_sleep_score_preprocout_datanum_1.eeg',...
    'test_lindev_sleep_score_preprocout_datanum_2.eeg'};
scoringfiles = {...
    'test_lindev_sleep_score_preprocout_datanum_1.txt',...
    'test_lindev_sleep_score_preprocout_datanum_2.txt'};
lightsoffs = {...
    108.922,...
    73.102};

% create some subjects

nDatasets = numel(datasets);
subjects = cell(1,nDatasets);

for iDataset = 1:nDatasets
    
    subject = [];
    % take the numer as subject name
    subject.name               = num2str(iDataset);
    subject.dataset            = datasets{iDataset};
    subject.scoringfile        = scoringfiles{iDataset};
    subject.lightsoff          = lightsoffs{iDataset};
    
    %these things are the same in all subjects
    subject.scoringformat      = 'spisop'; % e.g. 'zmax' or 'spisop'
    subject.standard           = 'rk'; % 'aasm' or 'rk'
    subject.scoring_dataoffset = 0;
    subject.eegchannels        = {'C3:A2', 'C4:A1'};
    
    %save subject individually
    save(['subject-' num2str(iDataset)],'subject');
    subjects{iDataset} = subject;
    
end
%save all the subjects in a cell array
save('subjects','subjects');


%% get the first sleep cycle

scoring_cycles_firsts = cell(1,numel(subjects));
scorings = cell(1,numel(subjects));
for iSubject = 1:numel(subjects)
    subject = subjects{iSubject};
    
    % read scoring
    cfg = [];
    cfg.scoringfile   = subject.scoringfile;
    cfg.scoringformat = subject.scoringformat;
    cfg.standard      = subject.standard; % 'aasm' or 'rk'
    [scoring] = st_read_scoring(cfg);
    scoring.lightsoff = subject.lightsoff;
    
    scorings{iSubject} = scoring;

    % cut into sleep cycles
    cfg = [];
    res_cycle = st_sleepcycles([],scoring);
    cfg.start = res_cycle.table.startepoch;
    cfg.end = res_cycle.table.endepoch;

    scoring_cycles = st_cutscoring(cfg,scoring);
    
    scoring_cycles_firsts{iSubject} = scoring_cycles{1};
    
end

save('scorings', 'scorings');
save('scoring_cycles_firsts', 'scoring_cycles_firsts');

%% caluculate the power

res_power_bins = cell(1,numel(subjects));
res_power_bands = cell(1,numel(subjects));

for iSubject = 1:numel(subjects)
    subject = subjects{iSubject};
    % non-rem power in the first sleep cycle
    cfg = [];
    cfg.scoring     = scoring_cycles_firsts{iSubject};
    cfg.stages      = {'S2', 'S3', 'S4'}; % {'R'};
    cfg.channel     = subject.eegchannels;
    cfg.dataset     = subject.dataset;
    [res_power_bin, res_power_band] = st_power(cfg);
    
    res_power_bins{iSubject} = res_power_bin;
    res_power_bands{iSubject} = res_power_band;
end

% put the results together and write them out
[res_power_bins_appended] = st_append_res(res_power_bins{:});
[res_power_bands_appended] = st_append_res(res_power_bands{:});
cfg = [];
cfg.prefix = 'example_multi';
cfg.infix  = subject.name;
cfg.posfix = '';
filelist_res_power_appended = st_write_res(cfg, res_power_bins_appended, res_power_bands_appended);

%% find the frequency power peaks, only one, e.g. for the 'fast spindles'
cfg = [];
cfg.peaknum = 1; %only one peak, to keep things simple
[res_freqpeaks] = st_freqpeak(cfg,res_power_bins_appended);

%results are stored in the table take a look
res_freqpeaks.table

%now lets exlude falsely marked peaks in a 
%non-spindle band range of below 9 and above 16 Hz
res_freqpeaks.table.freqpeak1((res_freqpeaks.table.freqpeak1 < 9) | (res_freqpeaks.table.freqpeak1 > 16)) = NaN;

% write our results
cfg = [];
cfg.prefix = 'example';
cfg.infix  = subject.name;
cfg.posfix = '';
filelist_res_freqpeaks = st_write_res(cfg, res_freqpeaks);

%keep the results for later
save('spindle_power_peak_freqs', 'res_freqpeaks');

%% detect spindle on the frequency power peaks
res_spindles_channels = cell(1,numel(subjects));
res_spindles_events = cell(1,numel(subjects));

for iSubject = 1:numel(subjects)
    subject = subjects{iSubject};
    
    freqpeak = res_freqpeaks.table.freqpeak1(iSubject);
    
    cfg = [];
    cfg.scoring          = scoring_cycles_firsts{iSubject};
    cfg.stages           = {'S2', 'S3', 'S4'}; % {'R'};
    cfg.channel          = subject.eegchannels;
    cfg.centerfrequency  = freqpeak;
    cfg.dataset          = subject.dataset;
    [res_spindles_channel, res_spindles_event, res_spindles_filter] = st_spindles(cfg);
    
    res_spindles_channels{iSubject} = res_spindles_channel;
    res_spindles_events{iSubject} = res_spindles_event;
end

% put the results together and write them out
[res_spindles_channels_appended] = st_append_res(res_spindles_channels{:});
[res_spindles_events_appended] = st_append_res(res_spindles_events{:});
cfg = [];
cfg.prefix = 'example_multi';
cfg.infix  = subject.name;
cfg.posfix = '';
filelist_res_spindles_appended = st_write_res(cfg, res_spindles_channels_appended, res_spindles_events_appended);

%% detect slow waves
res_slowwaves_channels = cell(1,numel(subjects));
res_slowwaves_events = cell(1,numel(subjects));

for iSubject = 1:numel(subjects)
    subject = subjects{iSubject};
    cfg = [];
    cfg.scoring          = scoring_cycles_firsts{iSubject};
    cfg.stages           = {'S2', 'S3', 'S4'}; % {'R'};
    %cfg.thresholdstages  = {'S2'}; % if you want threshold on another sleep stage
    cfg.channel          = subject.eegchannels;
    cfg.dataset          = subject.dataset;
    [res_slowwaves_channel, res_slowwaves_event, res_slowwaves_filter] = st_slowwaves(cfg);
    
    res_slowwaves_channels{iSubject} = res_slowwaves_channel;
    res_slowwaves_events{iSubject} = res_slowwaves_event;
end

% put the results together and write them out
[res_slowwaves_channels_appended] = st_append_res(res_slowwaves_channels{:});
[res_slowwaves_events_appended] = st_append_res(res_slowwaves_events{:});
cfg = [];
cfg.prefix = 'example_multi';
cfg.infix  = subject.name;
cfg.posfix = '';
filelist_res_slowwaves_appended = st_write_res(cfg, res_slowwaves_channels_appended, res_slowwaves_events_appended);

%% plot spindles and slow waves on a hypnogram per subject and save as pdf

for iSubject = 1:numel(subjects)
    subject = subjects{iSubject};

    res_spindles_event = res_spindles_events{iSubject};
    res_slowwaves_event = res_slowwaves_events{iSubject};

    
    scoring = scorings{iSubject};
    
    % the trough with the larges amplitude of the spindle shall give its time
    % point, this is important for time-locked event related potentials
    spindle_troughs_subject.eegchannels{1} = res_spindles_event.table.seconds_trough_max(strcmp(res_spindles_event.table.channel,{subject.eegchannels{1}}));
    spindle_troughs_subject.eegchannels{2} = res_spindles_event.table.seconds_trough_max(strcmp(res_spindles_event.table.channel,{subject.eegchannels{2}}));
    slowwave_troughs_subject.eegchannels{1} = res_slowwaves_event.table.seconds_trough_max(strcmp(res_slowwaves_event.table.channel,{subject.eegchannels{1}}));
    slowwave_troughs_subject.eegchannels{2} = res_slowwaves_event.table.seconds_trough_max(strcmp(res_slowwaves_event.table.channel,{subject.eegchannels{2}}));



    % we can also get the amplitudes
    spindle_amplitude_subject.eegchannels{1} = res_spindles_event.table.amplitude_peak2trough_max(strcmp(res_spindles_event.table.channel,{subject.eegchannels{1}}));
    spindle_amplitude_subject.eegchannels{2} = res_spindles_event.table.amplitude_peak2trough_max(strcmp(res_spindles_event.table.channel,{subject.eegchannels{2}}));
    slowwave_amplitude_subject.eegchannels{1} = res_slowwaves_event.table.amplitude_peak2trough_max(strcmp(res_slowwaves_event.table.channel,{subject.eegchannels{1}}));
    slowwave_amplitude_subject.eegchannels{2} = res_slowwaves_event.table.amplitude_peak2trough_max(strcmp(res_slowwaves_event.table.channel,{subject.eegchannels{2}}));

    % ... or the frequency of each sleep spindle
    spindle_frequency_subject.eegchannels{1} = res_spindles_event.table.frequency_by_mean_pk_trgh_cnt_per_dur(strcmp(res_spindles_event.table.channel,{subject.eegchannels{1}}));
    spindle_frequency_subject.eegchannels{2} = res_spindles_event.table.frequency_by_mean_pk_trgh_cnt_per_dur(strcmp(res_spindles_event.table.channel,{subject.eegchannels{2}}));

    % also plot event properties like amplitude and frequency for each event
    cfg = [];
    cfg.plotunknown        = 'no'; 
    cfg.figureoutputfile   = [subject.name '.pdf'];
    cfg.figureoutputformat = 'pdf';
    cfg.eventtimes  = {spindle_troughs_subject.eegchannels{1}';...
                       spindle_troughs_subject.eegchannels{2}';...
                       spindle_troughs_subject.eegchannels{1}';...
                       spindle_troughs_subject.eegchannels{2}';...
                       slowwave_troughs_subject.eegchannels{1}';...
                       slowwave_troughs_subject.eegchannels{2}'};
    cfg.eventvalues = {spindle_amplitude_subject.eegchannels{1}';...
                       spindle_amplitude_subject.eegchannels{2}';...
                       spindle_frequency_subject.eegchannels{1}';...
                       spindle_frequency_subject.eegchannels{2}';...
                       slowwave_amplitude_subject.eegchannels{1}';...
                       slowwave_amplitude_subject.eegchannels{2}'};
    cfg.eventranges = {[min(spindle_amplitude_subject.eegchannels{1}), max(spindle_amplitude_subject.eegchannels{1})];...
                       [min(spindle_amplitude_subject.eegchannels{2}), max(spindle_amplitude_subject.eegchannels{2})];...
                       [min(spindle_frequency_subject.eegchannels{1}), max(spindle_frequency_subject.eegchannels{1})];...
                       [min(spindle_frequency_subject.eegchannels{2}), max(spindle_frequency_subject.eegchannels{2})];...
                       [min(slowwave_amplitude_subject.eegchannels{1}), max(slowwave_amplitude_subject.eegchannels{1})];...
                       [min(slowwave_amplitude_subject.eegchannels{2}), max(slowwave_amplitude_subject.eegchannels{2})]};
    cfg.eventranges = cellfun(@(x) round(x,2), cfg.eventranges,'UniformOutput',false);               
    cfg.eventlabels = {['spd ' 'ampl ' subject.eegchannels{1}], ['spd ' 'ampl ' subject.eegchannels{2}], ...
                       ['spd ' 'freq ' subject.eegchannels{1}], ['spd ' 'freq ' subject.eegchannels{2}],...
                       ['sw ' 'ampl ' subject.eegchannels{1}], ['sw ' 'ampl ' subject.eegchannels{2}]};
    figure_handle = st_hypnoplot(cfg, scoring);
end


