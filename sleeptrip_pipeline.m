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
pathToSleepTrip = uigetdir('','choose SleepTrip path, e.g. D:\sleeptrip-main');
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


%% create a montage

%minimal montage example for zmax data
montage_zmax          = [];
montage_zmax.labelold = {'EEG L','EEG R'};
montage_zmax.labelnew = {'EEG L-R'};
%transition matrix 
montage_zmax.tra = [
%labelold: EEG L   EEG R  % labelnew:
             1       -1   % EEG L-R 
];

% other examples of an average mastoid montage
%         montage          = [];
%         montage.labelold = {
%             [mastoidChannelPrefixes{iDataset} '1'],[mastoidChannelPrefixes{iDataset} '2'],...
%             'F3','F4',...
%             'C3','C4',...
%             'O1','O2'};
%         montage.labelnew = {
%             'F3' 'F4'...
%             'C3' 'C4'...
%             'O1' 'O2'};
%         montage.tra      = [
%             %  A1   A2   F3  F4  C3  C4  O1  O2 
%               -0.5 -0.5  factor   0   0   0   0   0 % F3-A1A2
%               -0.5 -0.5  0   1   0   0   0   0 % F4-A1A2
%               -0.5 -0.5  0   0   1   0   0   0 % C3-A1A2
%               -0.5 -0.5  0   0   0   1   0   0 % C4-A1A2
%               -0.5 -0.5  0   0   0   0   1   0 % O1-A1A2
%               -0.5 -0.5  0   0   0   0   0   1 % O2-A1A2
%               ];

%% write/read/handle montages with standard settings
cfg = [];
cfg.filename = 'montage_zmax.txt';
montage_filepath = st_write_montage(cfg,montage_zmax);

cfg = [];
montage_zmax = st_read_montage(cfg,montage_filepath);

%%% combine multiple montages together to get one resulting
montage_zmax_level2_example = ft_apply_montage(montage_zmax, montage_zmax);

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
% only relevant for assessing sleep (stage) onset delays
subject.lightsoff          = 0; 
subject.eegchannels        = {'EEG L', 'EEG R'};
subject.ppgchannels        = {'OXY_IR_AC'};
subject.oxychannels        = {'OXY_IR_DC'};
%  possible other channels for zmax are: 
% {'BATT' 'BODY TEMP' ...
%  'dX' 'dY' 'dZ' ...
%  'EEG L' 'EEG R' 
%  'LIGHT' 'NASAL L' 'NASAL R' 'NOISE' ...
%  'OXY_DARK_AC' 'OXY_DARK_DC' 'OXY_IR_AC' 'OXY_IR_DC' ...
%  'OXY_R_AC' 'OXY_R_DC' 'RSSI'}

subject.datagrammerfile = 'datagrammer.txt';
subject.montage = montage_zmax;

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
% 
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
% cfg.to = 'aasm' % immediately convert it to a standard format if possible
scoring = st_read_scoring(cfg);

% do we know the lights-off moment, if yes then add it
% this is good for SleepTrip to know, e.g. for the st_sleepdescriptive
% function
scoring.lightsoff = subject.lightsoff;

% do we know the lights-on moment, in this case not, but if this would be
% good to determine more acurately where things at the end of the recoring
% are
%scoring.lightson = subject.lightson;

%practice: read in another format, maybe a custom format.

%%% exclude arousals (with respect to data start time = 0, if present
%%% exclude those epochs for further analysis
if isfield(scoring,'arousals')
    [scoring] = st_exclude_events_scoring(cfg, scoring, scoring.arousals.start, scoring.arousals.stop);
end


%% write out the scoring again and re-read it in.
cfg = [];
cfg.filename = 'scoring_test';
cfg.datatype = 'spisop';% 'sleeptrip' 'numbersincolumns'
%cfg.to = 'rk'
scoring_filepath = st_write_scoring(cfg, scoring);

%%% note that writing in cfg.scoringformat = 'numbersincolumns' is the same
%%% as to cfg.scoringformat = 'spisop', just the file endings and delimiters might differ
cfg = [];
cfg.scoringfile = scoring_filepath;
cfg.scoringformat = 'spisop'; %'sleeptrip' 'spisop' 'numbersincolumns'
%cfg.to = 'rk'
scoring_reread = st_read_scoring(cfg);

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

%% calculate a consensus scoring and stats (e.g. Fleiss's kappa)

scoring_ref = scoring;
scoring_alt = scoring;

%alter something, here just one epoch and one excluded
scoring_alt.epochs{100} = 'R';
scoring_alt.excluded(101) = 1;
scoring_alt.excluded(102) = 1;
scoring_ref.excluded(102) = 1;

cfg = [];
%cfg.align = 'sleeponset';
%cfg.reference = 'all'; % 'first' 'consensus' 'all'
%cfg.stat_alpha = 0.05;
%cfg.agreementthres  = 0.5;
%cfg.agreementthres_excl = 0.5;
[res_comp_fleiss_stat res_comp_cohen_stat res_contingency res_contingency_excluded scoring_consensus contingency_tables contingency_excluded_tables] = st_scoringcomp(cfg, scoring_ref, scoring_alt, scoring_alt);

res_comp_fleiss_stat.table
res_comp_cohen_stat.table
res_contingency.table
res_contingency_excluded.table
scoring_consensus
contingency_tables
contingency_excluded_tables

%practice: play around with cfg.agreementthres and see how the
%scoring_consensus introduces unknown scoring (i.e. '?') into the
%scoring_consensus.epochs
%Try to align on the sleep onset, when would this matter?
%What if not all periods are scored?



%% plot the scoring

cfg = [];
cfg.title           = subject.name;
cfg.plottype        = 'deepcolorblocks'; %'classic' 'colorblocks' 'colorbar' 'deepcolorblocks'
% cfg.yaxdisteqi      = 'yes';
% cfg.timeticksdiff   = 60;
cfg.plotunknown     = 'no'; % 'yes' or 'no' 
% cfg.plotsleeponset  = 'no'; % 'yes' or 'no' 
% cfg.plotsleepoffset = 'no'; % 'yes' or 'no' 
% cfg.timemin         = 600 % in minutes, e.g. plot on a 10-hour time axis.   
% cfg.timerange          = [50 100];
% cfg.considerdataoffset = 'no';
% cfg.plotexcluded       = 'no'; 
% cfg.colorblocksconnect = 'yes';
% cfg.colorscheme        = 'restless'; % 'bright' or 'dark' or 'restless'

%%% if you want to export the figure immediately add these parameters
cfg.figureoutputfile   = [subject.name '.pdf'];
cfg.figureoutputformat = 'pdf';

[figure_handle figure_axis_handle] = st_hypnoplot(cfg, scoring);
% close(figure_handle) close the figure automatically e.g. after exporting

%practice: play around with the plotting paramters. Export the file as a
%pdf or a png, change the figure dimensions and resolution, manipulate the
%figure by using the figure_handle, what happens when you use an RK or AASM
%based scoring?

%% create a sleep table with the scoring parameters

cfg = [];
res_scoringdescriptive = st_scoringdescriptives(cfg, scoring);

% lets take a look what we got in there
res_scoringdescriptive.table

% or do so sleep-cycle wise
cfg = [];
cfg. cycle = 'all';
res_sleepdescriptive_cycle = st_scoringdescriptives(cfg, scoring);

% lets take a look what we got in there
res_sleepdescriptive_cycle.table

% practice: how is sleep onset and lights-off moment related? 
% When is the end of the Total sleep time?
% Why are there several versions of the columns?

% export the results
cfg = [];
cfg.prefix = 'example';
cfg.infix  = subject.name;
cfg.posfix = '';
filelist_res_sleepdescriptives = st_write_res(cfg, res_scoringdescriptive); 

% note: you can write mutliple results with st_write_res(cfg, res_sleep1, res_sleep2, ...)
%       or if they are in a cell array with st_write_res(cfg, res_all{:})

%practice: were are the files stored, can this be changed, and how to use
%the prefix, infix and postfix best to differentiate your analysis from other runs/studies/subjects/recordings?


%% read in the data

cfg = [];
cfg.dataset    = subject.dataset;
cfg.continuous = 'yes';
cfg.channel    = subject.eegchannels;

%low pass filter to filter out signals higher than 35 Hz
cfg.lpfilter = 'yes';
cfg.lpfreq   = 35;

data = st_preprocessing(cfg);
% plot(data.trial{1}(1,1:(data.fsample*10)))

%% using and applying a montage
%if there is a montage defined you might whant to use it or not
if isfield(subject,'montage')
    cfg.montage = subject.montage;
end

data_with_montage = st_preprocessing(cfg);

%alternatively if the data is already read in
data_with_montage = ft_apply_montage(data,subject.montage);


%% extract the date from zmax/edf data header. 
% this will not work for other data formats than edf
hdr = ft_read_header(subject.dataset);
hdr = data.hdr;

try
    year = hdr.orig.T0(1);
    month = hdr.orig.T0(2);
    day = hdr.orig.T0(3);
    hour = hdr.orig.T0(4);
    minute = hdr.orig.T0(5);
    second = hdr.orig.T0(6);
    
    dt = datetime(year,month,day,hour,minute,second);
    dt
catch e
end



%% resample the data, e.g. before exporting/convert it to your hard drive
% not necessary as for many analysis the data 
% is resampled to values that are sufficient
cfg = [];
cfg.resamplefs = 100;
data_100Hz = ft_resampledata(cfg, data);

%% write the data out in a BIDS compatible dataformat (but zipped) by default

% chose a file name (without extension or with, does not matter)
% data is stored in 16bit brainvision format and compressd in a zip to save
% a lot of space on the hard drive.
cfg = [];
cfg.filename = [subject.name '_rs100'];
%cfg.format = 'edf';
%cfg.compress = 'no';
%cfg.posmarker = 'no';
[cfg filepaths] = st_write_continuous_data(cfg, data_100Hz);

% Note: All in allwe have downsampled from 256 to 100 Hz, stored the data in 16 bit and
% zipped it, if the original data would have been a typical brainvision file with the float32
% data format the we only use (100/256) * (16/32) * 0.7 < 14% of space
% that is we saved more than 85% of space!

%practice, do some more preprocessing, how would you read in your own data
% search on the web how fieldtrip can help to read in data http://www.fieldtriptoolbox.org/faq/dataformat/
% Note you can read in .zip data directly with ST_PREPROCESSING, this saves typically one third
% of space your hard drive at least, ... but it needs longer to open


%% convert the scoring to a homogenous or well known format 

% to AASM format
cfg = [];
cfg.to = 'aasm';
scoring_aasm = st_scoringconvert(cfg, scoring);

% to Rechtschaffen&Kales format (i.e. having sleep Stage 4 preserved)
cfg = [];
cfg.to = 'rk';
scoring_rk = st_scoringconvert(cfg, scoring);

% you can also convert to your custom scoring format, but keep in mind, 
% that using custom scorings will be limited, however you can use this to
% convert your custom format back into 'aasm' or 'rk' which work here.

%% use the aasm scoring as a basis
scoring = scoring_aasm;

%% reorder channels
cfg = [];
cfg.order = 'alphabetical';
data_reord = st_reorderdata(cfg, data);

%% applying grammers to data

%%% one grammer to multiple or signgle channels
cfg = [];
%cfg.grammer = 'hp 0.3 lp 35';
cfg.grammer = '( ( bp 12 14 mult 4 ) + ( hp 0.5 lp 2 ) ) + ( bp 6 8 mult 2 )';
cfg.channel = 1;
data_new_single_dg = st_datagrammer(cfg, data);

%%% read a data grammer
cfg = [];
datagrammer = st_read_datagrammer(cfg, subject.datagrammerfile);

datagrammer.channel
datagrammer.grammer

%%% multiple grammers maped by channels
data_new_multiple_dg = st_apply_datagrammer(data, datagrammer);


%% whole night (multitapper) spectrogram

cfg = [];
cfg.approach = 'spectrogram'; % 'spectrogram' 'mtmfft_segments' 'mtmconvol_memeff'
%cfg.approach = 'mtmfft_segments';
cfg.taper  = 'dpss'; % 'hanning' 'hanning_proportion' 'dpss'
cfg.transform  = 'log10'; % 'none' 'db' 'db(p+1)' 'log10' 'log10(p+1)'
cfg.channel = subject.eegchannels;
cfg.powvalue = 'power';
freq_continous = st_tfr_continuous(cfg, data);

cfg = [];
cfg.zlim           = [-3 0];
cfg.colormap       = jet(256); % still a good colormap, however not for some forms of colorblindness
% cfg.layout         = ...
% ft_multiplotTFR(cfg, freq_continous);

cfg.channel        = subject.eegchannels{1};
ft_singleplotTFR(cfg, freq_continous);

%%% to make this a bit prettier and in hours as units
set(gca,'TickDir','out');
time = freq_continous.time;
timeticksdiff = 3600;
xTick = [0:timeticksdiff:(max([max(time)]))];
set(gca, 'xTick', xTick);
set(gca, 'xTickLabel', arrayfun(@num2str,round(xTick/3600,2),'UniformOutput',false)); 
timeunit = 'h';
set(gca, 'box', 'off')
xlabel(['Time [' timeunit ']']);
ylabel('Frequency [Hz]');

%practice: compare the results of cfg.approach = 'mtmfft_segments' with
%cfg.approach = 'spectrogram', the latter computation is about 10 times
%faster as it uses matlab internal functions. Also play around with the
%cfg.taper options and see how that influences the speed of calculation vs.
% 'prettyness'. Finally check the the cfg.transform options to tansform 
% the power values and see how to change z axis scaling cfg.zlim (after commenting it out
% first to check the scale) to adapt to the different transformations.

%% looking at the data in the "score" browser STILL ALPHA version
cfg = [];
cfg.renderer = 'opengl'; % to get things plotted a little faster
cfg.scoring = scoring;
cfg.bgcolor = 'dark'; % 'white' or 'dark'
% Instead of specifiying a dataset, data header, a montage or a scoring you
% can also dot this interactively by setting
% cfg.datainteractive = 'yes';
cfg_postscore = st_scorebrowser(cfg ,data);

% st_scorebrowser(cfg, data); % without an output argument the browser will give the handling back to Matlab immediately. Not recommended!
% There are many features, it is recommended to never press Ctrl+C while
% using it (as to not crash it) and to see the button "shortcuts" to get
% some help on what can be done.
   

%% take a look at the data, MIGHT NOT WORK ON ALL MATLAB VERSIONS!

% get the indices of the respective epochs for later marking
epochs_W    = find(strcmp(scoring.epochs, 'W'));
epochs_N1   = find(strcmp(scoring.epochs, 'N1'));
epochs_N2   = find(strcmp(scoring.epochs, 'N2'));
epochs_N3   = find(strcmp(scoring.epochs, 'N3'));
epochs_R    = find(strcmp(scoring.epochs, 'R'));

% take the epochs and transform them into sample values, st_times2samples
% helps to match this to the samples to the data structure
artfctdef                  = [];
artfctdef.W.artifact  = st_times2samples(data,[(epochs_W(:)  -1)*scoring.epochlength epochs_W(:)*scoring.epochlength],'rmnanrows');
artfctdef.N1.artifact = st_times2samples(data,[(epochs_N1(:) -1)*scoring.epochlength epochs_N1(:)*scoring.epochlength],'rmnanrows');
artfctdef.N2.artifact = st_times2samples(data,[(epochs_N2(:) -1)*scoring.epochlength epochs_N2(:)*scoring.epochlength],'rmnanrows');
artfctdef.N3.artifact = st_times2samples(data,[(epochs_N3(:) -1)*scoring.epochlength epochs_N3(:)*scoring.epochlength],'rmnanrows');
artfctdef.R.artifact  = st_times2samples(data,[(epochs_R(:)  -1)*scoring.epochlength epochs_R(:)*scoring.epochlength],'rmnanrows');

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

% this saves RAM, but is slower, especially with zmax data because it is
% compressed (e.g. in a .zip file).
% thus better use this option only when the EEG files are not compressed,
% otherwise might be slow.
% cfg.dataset     = subject.dataset;
% [res_power_bins, res_power_bands] = st_power(cfg); 

cfg.quick = 'yes'; % if you want to do quick and dirty determination
% cfg.windowproportion = 1/scoring.epochlength;
% cfg.segmentlength = scoring.epochlength;
% cfg.segmentoverlap = 0;
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
figure;
%power_y = sqrt(power_y)
plot(freq_x,power_y) % raw
figure;
power_y_logscale = 10*log10(power_y+1);
plot(freq_x,power_y_logscale) % db scaled
figure;
plot(freq_x,power_y_logscale.*freq_x) % db scaled, 1/freq corrected

% practice: Try different scaling and normalization, how would this change
% the results? Can you plot muliple channels at the same time.

%% detect slow waves

cfg = [];
cfg.scoring          = scoring;
cfg.stages           = {'N2', 'N3'};
%cfg.thresholdstages  = {'N3'}; % if you want threshold on another sleep stage
cfg.channel          = subject.eegchannels;

% this saves RAM, but is slower, especially with zmax data because it is
% compressed (e.g. in a .zip file).
% thus better use this option only when the EEG files are not compressed,
% otherwise might be slow.
% cfg.dataset     = subject.dataset;
% [res_spindles_channel res_spindles_event res_spindles_filter] = st_spindles(cfg);

[res_slowwaves_channel, res_slowwaves_event, res_slowwaves_filter] = st_slowwaves(cfg, data);

%export the data
cfg = [];
cfg.prefix = 'example';
cfg.infix  = subject.name;
cfg.posfix = '';
filelist_res_slowwave = st_write_res(cfg, res_slowwaves_channel, res_slowwaves_event); % write mutliple results with  st_write_res(cfg, res_sleep1, res_sleep2)

%% plot a event related potentials (ERP) of slow waves

%%% set the time limits for plotting
time_min = -1.5;
time_max = 1.5;

cfg = [];
%cfg.maxevents = 50; % only use random 50 events per channel
[timelock reschannelcolumnname] = st_channel_event_erp(cfg, res_slowwaves_event, data);
%[timelock reschannelcolumnname] = st_channel_event_erp(cfg, res_slowwaves_event, data_100Hz);


% split by the ERP data between the channels to plot in single figure for
% one channel
cfg = [];
cfg.channel = subject.eegchannels{1};
timelock_ch1 = ft_selectdata(cfg,timelock);
cfg.channel = subject.eegchannels{2};
timelock_ch2 = ft_selectdata(cfg,timelock);

%overwrite the labels as to making plotting with ft_singleplotER feasible
timelock_ch1.label = {'chan'};
timelock_ch2.label = {'chan'};


% view the timelocked signals for both channels
figure
cfg        = [];
cfg.xlim   = [time_min time_max];
cfg.title  = 'Non-REM event ERP time-locked to trough';
cfg.graphcolor = 'br';
cfg.linestyle = {'-','-'};
cfg.linewidth = 1;
fh = ft_singleplotER(cfg,timelock_ch1,timelock_ch2);


%%% to make this a bit prettier and in hours as units
time_min = -1.5;
time_max = 1.5;
time_ticks = time_min:0.25:time_max;
time_tickLabels = arrayfun(@(t) num2str(t), time_ticks(:), 'UniformOutput', false)';
time_tickLabels{time_ticks == 0} = 'Trough';
time_tickLabels(2:2:(end-1)) = {' '};

ylabel(gca,'Signal [µV]');
xlabel(gca,'Time [s]');
set(gca, 'xTick', time_ticks);
set(gca, 'xTickLabel', time_tickLabels);
%set(gca, 'xMinorTick', 'on');
set(gca, 'TickDir','out');
set(gca, 'box', 'off')
set(gca, 'LineWidth',1)
set(gca, 'TickLength',[0.01 0.01]);

%%% practice: if the data is in sampling rate of 256 Hz but the events have
%%% been detected in 100 Hz note that the time-point 0 is off due to
%%% rounding issues. You can resolve by using a downsampled data e.g.
%%% data_100Hz from further above to create the ERPs


%% plot a time-frequency power (TFR power) of slow waves

%%% set the time limits for plotting
time_min = -1.5;
time_max = 1.5;

cfg = [];
cfg.maxevents = Inf; % only use random 50 events per channel
cfg.baselinecorrect = 'no';

% cfg.baselinecorrect = 'yes'; 
% cfg.baseline = [time_min time_max];    
% cfg.baselinetype = 'zscore'; 


%%% the following lines are taken from the end of ft_freqbaseline and show
%%% how the baseline correction calculations are based on and the options
% if (strcmp(baselinetype, 'absolute'))
%   data = data - meanVals;
% elseif (strcmp(baselinetype, 'relative'))
%   data = data ./ meanVals;
% elseif (strcmp(baselinetype, 'relchange'))
%   data = (data - meanVals) ./ meanVals;
% elseif (strcmp(baselinetype, 'normchange')) || (strcmp(baselinetype, 'vssum'))
%   data = (data - meanVals) ./ (data + meanVals);
% elseif (strcmp(baselinetype, 'db'))
%   data = 10*log10(data ./ meanVals);
% elseif (strcmp(baselinetype, 'db_absolute'))
%   data = 10*log10(data - meanVals);
% elseif (strcmp(baselinetype, 'db(data+1)'))
%   data = 10*log10((data ./ meanVals)+1);
% elseif (strcmp(baselinetype, 'db(data+1)_absolute'))
%   data = 10*log10((data - meanVals)+1);
% elseif (strcmp(baselinetype, 'log10'))
%   data = log10(data ./ meanVals);
% elseif (strcmp(baselinetype, 'log10_absolute'))
%   data = log10(data - meanVals);
% elseif (strcmp(baselinetype, 'log10(data+1)'))
%   data = log10((data ./ meanVals)+1);
% elseif (strcmp(baselinetype, 'log10(data+1)_absolute'))
%   data = log10((data - meanVals)+1);
% elseif (strcmp(baselinetype,'zscore'))
%     stdVals = repmat(nanstd(data(:,:,baselineTimes),1, 3), [1 1 size(data, 3)]);
%     data=(data-meanVals)./stdVals;

% cfg.method = 'mtmconvol';
% cfg.taper = 'dpss';
% cfg.foi = 4:0.5:30;
% cfg.tapsmofrq  = 0.4*cfg.foi;
% cfg.t_ftimwin  = 5./cfg.foi;

cfg.foi = 4:0.1:24;
cfg.timestep = 0.05;
cfg.bounds = [time_min time_max];
[event_freq reschannelcolumnname] = st_channel_event_tfr(cfg, res_slowwaves_event, data);

cfg = [];
cfg.baseline = [time_min time_max];    
cfg.baselinetype = 'zscore';    
[event_freq_baselined] = ft_freqbaseline(cfg, event_freq);

%%% localize the fast spindle delay (of the power maximum) 
%  and primary frequency of the slow-wave-locked spindle activity. 
cfg = [];
cfg.timewin = [0 0.75];
cfg.freqwin = [11 16];
[time_freq_pow_centroids] = st_tfr_localize(cfg, event_freq_baselined);
time_freq_pow_centroids

% split by the TFR data between the channels to plot in single figure for
% one channel
cfg = [];
cfg.channel = subject.eegchannels{1};
event_freq_ch1 = ft_selectdata(cfg,event_freq_baselined);
cfg.channel = subject.eegchannels{2};
event_freq_ch2 = ft_selectdata(cfg,event_freq_baselined);

% view the time-frequency of a slow wave or spindle event, for each channel
cfg                = [];
%cfg.baseline       = [time_min time_max]; % a 3 s baseline around the event as it has no clear start or end.
%cfg.baselinetype   = 'normchange';
%cfg.zlim           = [-0.2 0.2];
cfg.xlim           = [time_min time_max];
cfg.title          = 'Event, time-frequency';
fh = ft_singleplotTFR(cfg,event_freq_ch1);
fh = ft_singleplotTFR(cfg,event_freq_ch2);



%% determine the spindle(s) frequency peak(s) from the power spectrum

cfg = [];
cfg.peaknum = 1; % either 1 or 2 (default)
%cfg.powlawnorm = 'yes';
[res_freqpeaks] = st_freqpeak(cfg,res_power_bin);

res_freqpeaks.table

%practice: find two peaks by default by changing the configuration of the function,
%would limiting the frequency band help, are there any artifacts (e.g. sharp peaks in the power spectra?)

%% detect sleep spindles

% define the frequency peak to use from the previous peak determination
% we will take the first peak we defined and the first value of that
freqpeak = res_freqpeaks.table.freqpeak1(1);

cfg = [];
cfg.scoring          =  scoring;
cfg.stages           = {'N2','N3'};
%cfg.stages           = {'S2', 'S3', 'S4'}; % {'R'};
cfg.channel          = subject.eegchannels;
% freqpeak1 for slow spindles
% freqpeak2 for fast spindles
cfg.centerfrequency  = freqpeak; % e.g. 12 for a frontal channel 13.3 for central channel

% this saves RAM, but is slower, especially with zmax data because it is
% compressed (e.g. in a .zip file).
% thus better use this option only when the EEG files are not compressed,
% otherwise might be slow.
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

%% plot a event related potentials (ERP) of spindles

%%% set the time limits for plotting
time_min = -1.5;
time_max = 1.5;

cfg = [];
%cfg.maxevents = 50; % only use random 50 events per channel
[timelock reschannelcolumnname] = st_channel_event_erp(cfg, res_spindles_event, data);
%[timelock reschannelcolumnname] = st_channel_event_erp(cfg, res_spindles_event, data_100Hz);


% split by the ERP data between the channels to plot in single figure for
% one channel
cfg = [];
cfg.channel = subject.eegchannels{1};
timelock_ch1 = ft_selectdata(cfg,timelock);
cfg.channel = subject.eegchannels{2};
timelock_ch2 = ft_selectdata(cfg,timelock);

%overwrite the labels as to making plotting with ft_singleplotER feasible
timelock_ch1.label = {'chan'};
timelock_ch2.label = {'chan'};


% view the timelocked signals for both channels
figure
cfg        = [];
cfg.xlim   = [time_min time_max];
cfg.title  = 'Non-REM event ERP time-locked to trough';
cfg.graphcolor = 'br';
cfg.linestyle = {'-','-'};
cfg.linewidth = 1;
fh = ft_singleplotER(cfg,timelock_ch1,timelock_ch2);


%%% to make this a bit prettier and in hours as units
time_min = -1.5;
time_max = 1.5;
time_ticks = time_min:0.25:time_max;
time_tickLabels = arrayfun(@(t) num2str(t), time_ticks(:), 'UniformOutput', false)';
time_tickLabels{time_ticks == 0} = 'Trough';
time_tickLabels(2:2:(end-1)) = {' '};

ylabel(gca,'Signal [µV]');
xlabel(gca,'Time [s]');
set(gca, 'xTick', time_ticks);
set(gca, 'xTickLabel', time_tickLabels);
%set(gca, 'xMinorTick', 'on');
set(gca, 'TickDir','out');
set(gca, 'box', 'off')
set(gca, 'LineWidth',1)
set(gca, 'TickLength',[0.01 0.01]);

%%% practice: if the data is in sampling rate of 256 Hz but the events have
%%% been detected in 100 Hz note that the time-point 0 is off due to
%%% rounding issues. You can resolve by using a downsampled data e.g.
%%% data_100Hz from further above to create the ERPs


%% plot a time-frequency power (TFR power) of spindles

%%% set the time limits for plotting
time_min = -1.5;
time_max = 1.5;

cfg = [];
cfg.maxevents = Inf; % only use random 50 events per channel
cfg.baselinecorrect = 'no';

% cfg.baselinecorrect = 'yes'; 
% cfg.baseline = [time_min time_max];    
% cfg.baselinetype = 'zscore'; 


%%% the following lines are taken from the end of ft_freqbaseline and show
%%% how the baseline correction calculations are based on and the options
% if (strcmp(baselinetype, 'absolute'))
%   data = data - meanVals;
% elseif (strcmp(baselinetype, 'relative'))
%   data = data ./ meanVals;
% elseif (strcmp(baselinetype, 'relchange'))
%   data = (data - meanVals) ./ meanVals;
% elseif (strcmp(baselinetype, 'normchange')) || (strcmp(baselinetype, 'vssum'))
%   data = (data - meanVals) ./ (data + meanVals);
% elseif (strcmp(baselinetype, 'db'))
%   data = 10*log10(data ./ meanVals);
% elseif (strcmp(baselinetype, 'db_absolute'))
%   data = 10*log10(data - meanVals);
% elseif (strcmp(baselinetype, 'db(data+1)'))
%   data = 10*log10((data ./ meanVals)+1);
% elseif (strcmp(baselinetype, 'db(data+1)_absolute'))
%   data = 10*log10((data - meanVals)+1);
% elseif (strcmp(baselinetype, 'log10'))
%   data = log10(data ./ meanVals);
% elseif (strcmp(baselinetype, 'log10_absolute'))
%   data = log10(data - meanVals);
% elseif (strcmp(baselinetype, 'log10(data+1)'))
%   data = log10((data ./ meanVals)+1);
% elseif (strcmp(baselinetype, 'log10(data+1)_absolute'))
%   data = log10((data - meanVals)+1);
% elseif (strcmp(baselinetype,'zscore'))
%     stdVals = repmat(nanstd(data(:,:,baselineTimes),1, 3), [1 1 size(data, 3)]);
%     data=(data-meanVals)./stdVals;

% cfg.method = 'mtmconvol';
% cfg.taper = 'dpss';
% cfg.foi = 4:0.5:30;
% cfg.tapsmofrq  = 0.4*cfg.foi;
% cfg.t_ftimwin  = 5./cfg.foi;

cfg.foi = 4:0.1:24;
cfg.timestep = 0.05;
cfg.bounds = [time_min time_max];
[event_freq reschannelcolumnname] = st_channel_event_tfr(cfg, res_spindles_event, data);

cfg = [];
cfg.baseline = [time_min time_max];    
cfg.baselinetype = 'zscore';    
[event_freq_baselined] = ft_freqbaseline(cfg, event_freq);

%%% localize the fast spindle delay (of the power maximum) 
%  and primary frequency of the slow-wave-locked spindle activity. 
cfg = [];
cfg.timewin = [-0.5 0.5];
cfg.freqwin = [11 16];
[time_freq_pow_centroids] = st_tfr_localize(cfg, event_freq_baselined);
time_freq_pow_centroids

% split by the TFR data between the channels to plot in single figure for
% one channel
cfg = [];
cfg.channel = subject.eegchannels{1};
event_freq_ch1 = ft_selectdata(cfg,event_freq_baselined);
cfg.channel = subject.eegchannels{2};
event_freq_ch2 = ft_selectdata(cfg,event_freq_baselined);

% view the time-frequency of a slow wave or spindle event, for each channel
cfg                = [];
%cfg.baseline       = [time_min time_max]; % a 3 s baseline around the event as it has no clear start or end.
%cfg.baselinetype   = 'normchange';
%cfg.zlim           = [-0.2 0.2];
cfg.xlim           = [time_min time_max];
cfg.title          = 'Event, time-frequency';
fh = ft_singleplotTFR(cfg,event_freq_ch1);
fh = ft_singleplotTFR(cfg,event_freq_ch2);

%% SW-(fast)spindles/SO-spindles detection from cooccurrence of both events

%%% by default the time window in a slow wave that a spindle 
%%% is considered a SW-spindle is if its maximal trough of that spindles 
%%% falls withing the time window of the trough and the end of the slow wave.
%%% This is mostly true for fast spindles and mastoid/neutral referenced
%%% channels. st_swsp is a wrapper function for st_cooc

cfg = [];
[res_swsp_summary, res_swsp_channel_stat, res_nonswsp_channel_stat, res_nonspsw_channel_stat, res_swsp_event, res_nonswsp_event, res_nonspsw_event, res_excluded_sp_event, res_excluded_sw_event] ...
    = st_swsp(cfg, res_spindles_event, res_slowwaves_event);

res_swsp_summary.table
res_swsp_channel_stat.table

%write out the results
cfg = [];
cfg.prefix = 'example';
cfg.infix  = subject.name;
cfg.posfix = '';
filelist_res_swsp = st_write_res(cfg, res_swsp_summary, res_swsp_channel_stat, res_nonswsp_channel_stat, res_nonspsw_channel_stat, res_swsp_event, res_nonswsp_event, res_nonspsw_event, res_excluded_sp_event, res_excluded_sw_event); % write mutliple results with  st_write_res(cfg, res_sleep1, res_sleep2)

%% plot a event related potentials (ERP) of  SW-spindles

%%% set the time limits for plotting
time_min = -1.5;
time_max = 1.5;

cfg = [];
%cfg.maxevents = 50; % only use random 50 events per channel
cfg.eventtimecolumn = 'sp_seconds_trough_max';
%cfg.eventtimecolumn = 'sw_seconds_trough_max';
[timelock reschannelcolumnname] = st_channel_event_erp(cfg, res_swsp_event, data);
%[timelock reschannelcolumnname] = st_channel_event_erp(cfg, res_swsp_event, data_100Hz);


% split by the ERP data between the channels to plot in single figure for
% one channel
cfg = [];
cfg.channel = subject.eegchannels{1};
timelock_ch1 = ft_selectdata(cfg,timelock);
cfg.channel = subject.eegchannels{2};
timelock_ch2 = ft_selectdata(cfg,timelock);

%overwrite the labels as to making plotting with ft_singleplotER feasible
timelock_ch1.label = {'chan'};
timelock_ch2.label = {'chan'};


% view the timelocked signals for both channels
figure
cfg        = [];
cfg.xlim   = [time_min time_max];
cfg.title  = 'Non-REM event ERP time-locked to trough';
cfg.graphcolor = 'br';
cfg.linestyle = {'-','-'};
cfg.linewidth = 1;
fh = ft_singleplotER(cfg,timelock_ch1,timelock_ch2);


%%% to make this a bit prettier and in hours as units
time_min = -1.5;
time_max = 1.5;
time_ticks = time_min:0.25:time_max;
time_tickLabels = arrayfun(@(t) num2str(t), time_ticks(:), 'UniformOutput', false)';
time_tickLabels{time_ticks == 0} = 'Trough';
time_tickLabels(2:2:(end-1)) = {' '};

ylabel(gca,'Signal [µV]');
xlabel(gca,'Time [s]');
set(gca, 'xTick', time_ticks);
set(gca, 'xTickLabel', time_tickLabels);
%set(gca, 'xMinorTick', 'on');
set(gca, 'TickDir','out');
set(gca, 'box', 'off')
set(gca, 'LineWidth',1)
set(gca, 'TickLength',[0.01 0.01]);

%%% practice: if the data is in sampling rate of 256 Hz but the events have
%%% been detected in 100 Hz note that the time-point 0 is off due to
%%% rounding issues. You can resolve by using a downsampled data e.g.
%%% data_100Hz from further above to create the ERPs


%% plot a time-frequency power (TFR power) of SW-spindles

%%% set the time limits for plotting
time_min = -1.5;
time_max = 1.5;

cfg = [];
cfg.maxevents = Inf; % only use random 50 events per channel
cfg.baselinecorrect = 'no';

% cfg.baselinecorrect = 'yes'; 
% cfg.baseline = [time_min time_max];    
% cfg.baselinetype = 'zscore'; 


%%% the following lines are taken from the end of ft_freqbaseline and show
%%% how the baseline correction calculations are based on and the options
% if (strcmp(baselinetype, 'absolute'))
%   data = data - meanVals;
% elseif (strcmp(baselinetype, 'relative'))
%   data = data ./ meanVals;
% elseif (strcmp(baselinetype, 'relchange'))
%   data = (data - meanVals) ./ meanVals;
% elseif (strcmp(baselinetype, 'normchange')) || (strcmp(baselinetype, 'vssum'))
%   data = (data - meanVals) ./ (data + meanVals);
% elseif (strcmp(baselinetype, 'db'))
%   data = 10*log10(data ./ meanVals);
% elseif (strcmp(baselinetype, 'db_absolute'))
%   data = 10*log10(data - meanVals);
% elseif (strcmp(baselinetype, 'db(data+1)'))
%   data = 10*log10((data ./ meanVals)+1);
% elseif (strcmp(baselinetype, 'db(data+1)_absolute'))
%   data = 10*log10((data - meanVals)+1);
% elseif (strcmp(baselinetype, 'log10'))
%   data = log10(data ./ meanVals);
% elseif (strcmp(baselinetype, 'log10_absolute'))
%   data = log10(data - meanVals);
% elseif (strcmp(baselinetype, 'log10(data+1)'))
%   data = log10((data ./ meanVals)+1);
% elseif (strcmp(baselinetype, 'log10(data+1)_absolute'))
%   data = log10((data - meanVals)+1);
% elseif (strcmp(baselinetype,'zscore'))
%     stdVals = repmat(nanstd(data(:,:,baselineTimes),1, 3), [1 1 size(data, 3)]);
%     data=(data-meanVals)./stdVals;

% cfg.method = 'mtmconvol';
% cfg.taper = 'dpss';
% cfg.foi = 4:0.5:30;
% cfg.tapsmofrq  = 0.4*cfg.foi;
% cfg.t_ftimwin  = 5./cfg.foi;

cfg.foi = 4:0.1:24;
cfg.timestep = 0.05;
cfg.bounds = [time_min time_max];
cfg.eventtimecolumn = 'sp_seconds_trough_max';
%cfg.eventtimecolumn = 'sw_seconds_trough_max';
[event_freq reschannelcolumnname] = st_channel_event_tfr(cfg, res_swsp_event, data);

cfg = [];
cfg.baseline = [time_min time_max];    
cfg.baselinetype = 'zscore';    
[event_freq_baselined] = ft_freqbaseline(cfg, event_freq);

%%% localize the fast spindle delay (of the power maximum) 
%  and primary frequency of the slow-wave-locked spindle activity. 
cfg = [];
cfg.timewin = [-0.5 0.5];
cfg.freqwin = [11 16];
[time_freq_pow_centroids] = st_tfr_localize(cfg, event_freq_baselined);
time_freq_pow_centroids

% split by the TFR data between the channels to plot in single figure for
% one channel
cfg = [];
cfg.channel = subject.eegchannels{1};
event_freq_ch1 = ft_selectdata(cfg,event_freq_baselined);
cfg.channel = subject.eegchannels{2};
event_freq_ch2 = ft_selectdata(cfg,event_freq_baselined);

% view the time-frequency of a slow wave or spindle event, for each channel
cfg                = [];
%cfg.baseline       = [time_min time_max]; % a 3 s baseline around the event as it has no clear start or end.
%cfg.baselinetype   = 'normchange';
%cfg.zlim           = [-0.2 0.2];
cfg.xlim           = [time_min time_max];
cfg.title          = 'Event, time-frequency';
fh = ft_singleplotTFR(cfg,event_freq_ch1);
fh = ft_singleplotTFR(cfg,event_freq_ch2);


%% stack all the results together 

res_all_stacked = st_stack_res(...
res_sleepdescriptive_cycle,...
res_scoringdescriptive,...
res_slowwaves_channel, res_slowwaves_event,...
res_spindles_channel, res_spindles_event, res_spindles_filter,...
res_swsp_summary, res_swsp_channel_stat, res_nonswsp_channel_stat, res_nonspsw_channel_stat, res_swsp_event, res_nonswsp_event, res_nonspsw_event, res_excluded_sp_event, res_excluded_sw_event...
);


%take a look
res_all_stacked.table

%write out the results
cfg = [];
cfg.prefix = 'example';
cfg.infix  = subject.name;
cfg.posfix = '';
filelist_res_all_stacked = st_write_res(cfg, res_all_stacked); % write mutliple results with  st_write_res(cfg, res_sleep1, res_sleep2)




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
cfg.eventvalueranges = {[min(spindle_amplitude_subject.eegchannels{1}), max(spindle_amplitude_subject.eegchannels{1})];...
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
artfctdef.spindles_ch1.artifact   = st_times2samples(data,[spindle_begins_subject.eegchannels{1}, spindle_ends_subject.eegchannels{1}]);
%artfctdef.spindles_ch1.artifact   = fix(data.fsample*[spindle_begins_subject.eegchannels{1}, spindle_ends_subject.eegchannels{1}]);
%artfctdef.spindles_ch1_trough.artfctpeak = fix(data.fsample*spindle_troughs_subject.eegchannels{1});
artfctdef.spindles_ch2.artifact   = st_times2samples(data,[spindle_begins_subject.eegchannels{2}, spindle_ends_subject.eegchannels{2}]);
%artfctdef.spindles_ch2.artifact   = fix(data.fsample*[spindle_begins_subject.eegchannels{2}, spindle_ends_subject.eegchannels{2}]);
%artfctdef.spindles_ch2_trough.artfctpeak = fix(data.fsample*spindle_troughs_subject.eegchannels{2});
% cfg.selectmode                        =  'markpeakevent'; %'markartifact', 'markpeakevent', 'marktroughevent' (default = 'markartifact')
cfg.artfctdef     = artfctdef;
cfg.plotevents    = 'yes';
cfg.renderer      = 'opengl'; % 'painters' or 'opengl' or 'zbuffer'
ft_databrowser(cfg, data);


%% plot a event related potentials (ERP) and frequency (ERF) of e.g. spindles
% a buffer we need to have padding left and right to make nice
% time-frequency graph later on, so +-5s should be enough buffer
% for the +- 1s time window proper

cfg = [];
cfg.seconds = spindle_troughs_subject.eegchannels{1};
cfg.bounds = [-5 5]; % 10 seconds around the events
data_events_ch1 = ft_redefinetrial(cfg, data);
cfg.seconds = spindle_troughs_subject.eegchannels{2};
data_events_ch2 = ft_redefinetrial(cfg, data);

%View the event average signal timelocked to the trough.
cfg        = [];
timelock_ch1 = ft_timelockanalysis(cfg, data_events_ch1);
timelock_ch2 = ft_timelockanalysis(cfg, data_events_ch2);

%filter the timelocked data in the slow band to make underlying slow wave
%components more visible.
cfg = [];
% cfg.hpfilter = 'yes';
% cfg.hpfreq = 0.5;
cfg.lpfilter = 'yes';
cfg.lpfreq = 3;
timelock_ch1_SWA_fitlered = st_preprocessing(cfg, timelock_ch1);
timelock_ch2_SWA_fitlered = st_preprocessing(cfg, timelock_ch2);

% view the timelocked signals for both channels
figure
cfg        = [];
cfg.xlim   = [-1.5 1.5];
cfg.title  = 'Non-REM event ERP time-locked to trough';
cfg.graphcolor = 'brbr';
cfg.linestyle = {'-','-','-.','-.'};
cfg.linewidth = 1;
fh = ft_singleplotER(cfg,timelock_ch1,timelock_ch2,timelock_ch1_SWA_fitlered,timelock_ch2_SWA_fitlered);

% ERF
padding_buffer = 1;
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
fh = ft_singleplotTFR(cfg,event_freq_ch1);
figure
fh = ft_singleplotTFR(cfg,event_freq_ch2);

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
%cfg.plottype = 'colorbar';
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
res_freqpeaks = st_freqpeak(cfg,res_power_bins_appended);

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
    cfg.plottype = 'colorblocks';
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
    cfg.eventvalueranges = {[min(spindle_amplitude_subject.eegchannels{1}), max(spindle_amplitude_subject.eegchannels{1})];...
                       [min(spindle_amplitude_subject.eegchannels{2}), max(spindle_amplitude_subject.eegchannels{2})];...
                       [min(spindle_frequency_subject.eegchannels{1}), max(spindle_frequency_subject.eegchannels{1})];...
                       [min(spindle_frequency_subject.eegchannels{2}), max(spindle_frequency_subject.eegchannels{2})];...
                       [min(slowwave_amplitude_subject.eegchannels{1}), max(slowwave_amplitude_subject.eegchannels{1})];...
                       [min(slowwave_amplitude_subject.eegchannels{2}), max(slowwave_amplitude_subject.eegchannels{2})]};
    cfg.eventvalueranges = cellfun(@(x) round(x,2), cfg.eventvalueranges,'UniformOutput',false);               
    cfg.eventlabels = {['spd ' 'ampl ' subject.eegchannels{1}], ['spd ' 'ampl ' subject.eegchannels{2}], ...
                       ['spd ' 'freq ' subject.eegchannels{1}], ['spd ' 'freq ' subject.eegchannels{2}],...
                       ['sw ' 'ampl ' subject.eegchannels{1}], ['sw ' 'ampl ' subject.eegchannels{2}]};
    figure_handle = st_hypnoplot(cfg, scoring);
end

%% Topoplots example

% lets create an artificial result structure that is similar to what a
% st_power function would return but has different column names. basically
% any table with values associated to a channel will work.
res = [];
res.ori = 'st_power';
res.type = 'power_band';
res.cfg = [];
% here an example to read from a file
% res.table = readtable('power_values-N2_0.5-4Hz.tsv','FileType','text','Delimiter','\t');
% 'rel_difference' with the values and one that is called 'channel'%
% but we define it here explicitly 
res.table = table({'C3'; 'C4'; 'F3'; 'F4'; 'P3'; 'P4'; 'O1'; 'O2'}, ...
[-0.12 -0.14 -0.11 -0.14 -0.30 -0.19 -0.1 -0.08]',...
[1 0 1 0 0.5 0.5 1 1]',...
[1 1 1 1 1 1 0 0]'...
,'VariableNames',{'channel','rel_difference','mask','use'});
%{'yes'; 'yes'; 'yes'; 'yes'; 'yes'; 'yes'; 'no'; 'no'}...


%{
channel rel_difference  mask  use
C3      -0.12           1     1
C4      -0.14           0     1
F3      -0.11           1     1
F4      -0.14           0     1 
P3      -0.30           0.5   1
P4      -0.19           0.5   1
O1      -0.10           0.5   0
O2      -0.08           0.5   0
%}

% here another colorscheme than the standard one from Matlab
load('colormap_bluewhitered.mat')

figure
cfg = [];
cfg.property = 'rel_difference';

%see http://www.fieldtriptoolbox.org/tutorial/layout/
cfg.layout = 'easycapM25';

%cfg.channel = {'all','-C4'};
cfg.renderer = 'painters';
cfg.maskproperty   = 'mask';
cfg.gridscale = 256;
%cfg.interplimits = 'electrodes';
cfg.highlightchannel = {'C3', 'C4'};
cfg.highlightcolor   = [1 0 0];
cfg.highlight = 'on';
cfg.colorbar = 'yes';
cfg.colormap = colormap_bluewhitered;
cfg.zlim = [-0.4 0];
% FIXME: the next two lines wont work.
%cfg.filtercolumns =  {'use'};
%cfg.filtervalues  = {[1]};
cfg = st_topoplotres(cfg, res);
