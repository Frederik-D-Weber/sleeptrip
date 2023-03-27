%%
%Load the included mat file, containing:
% -data: ~2 h of 11-channel data starting right before sleep onset, referenced to linked mastoids, 0.3 Hz
% high-pass and 50 Hz notch-filtered.
% -elec: corresponding channel locations
% -scoring: corresponding sleep scores

load('tutorial_cleaning.mat')

data=ft_struct2double(data);
%%
%Let's first inspect the channel locations

%build a cfg (see st_plot_elec)
cfg=[];
cfg.elec=elec;

%call st_plot_elec, and check that channels are oriented the correct way.
st_plot_elec(cfg);

%(Note: channels here are named according to their EGI labels, not 10-20 standard.)
%%
%Let's also check the scoring by creating a hypnogram

%build a cfg (see st_hypnoplot)
cfg=[];
cfg.plottype        = 'colorblocks';
cfg.colorblocksconnect = 'yes';
cfg.colorscheme = 'restless';
cfg.timeticksdiff   = 60;
cfg.plotunknown     = 'no';
cfg.plotexcluded='no';
cfg.plotlegend='no';

%call st_hypnoplot
st_hypnoplot(cfg, scoring);

%Note: the hypnogram starts later because the first 50 minutes of the scoring (and data) have been removed for this
%tutorial
%%
%-----------------------DETECTOR SETS--------------------
%We can now set up things to clean our data.
%
% First, we need to specify a "detector set" (= a collection of individual
%detectors). SleepTrip offers different ways of doing so, ranging from taking the default
%detector set, to building a set of fully customized detectors.
%%
%-------------default detector sets---------
%1a) The simplest approach is to call st_get_default_detector_set. This will
%return the default detector set:

%build a cfg
cfg=[];
cfg.elec=elec; %elec influences which individual detectors to include (e.g., with only 2 channels, meaningful comparisons to neigboring channels cannot be performed)
cfg.fsample=data.fsample; %fsample influences which individual detectors to include (considering Nyquist etc.)

%get the default detector set
detector_set_default_all=st_get_default_detector_set(cfg);

%%
%what does a detector set look like?

%inspect the contents of the default detector set...
detector_set_default_all %default contains 7 individual detectors

%their names are contained in the field 'label':
detector_set_default_all.label

%inspect a particular detector (e.g., lowfreq)...
detector_set_default_all.detectors{2}

%...and the individual fields inside.
detector_set_default_all.detectors{2}.ft %these are details that will be provided to FieldTrip's ft_preprocessing
detector_set_default_all.detectors{2}.st %these are details for SleepTrip


%%
%1b) You can also specify which of the default detectors you want to
%include in your detector set:

cfg=[];
cfg.elec=elec;
cfg.fsample=data.fsample;
cfg.include={'absamp','flatline','deviant'}; %use the names of desired detectors

detector_set_default_custom=st_get_default_detector_set(cfg);

%Note: contents of cfg.elec and cfg.fsample may override which individual
%detectors are returned!
%%
% ----------custom detector sets------

%2a) A simple approach to customize limited aspects of your detectors is by starting from the
%default detector set.

%first make a copy of the default set
my_detector_set=detector_set_default_all;

%suppose we want to change the bandpass filter frequencies of the low-pass
%filter (detector 2). As it turns out, they're set here:
my_detector_set.detectors{2}.ft.bpfreq %vector of length 2

%simply change to new values and you're done
my_detector_set.detectors{2}.ft.bpfreq = [1 10]

%%
%2b) We can also specify individual detectors, and then combine them
%into a detector set.

%To simplify, we'll start by taking the 'absamp' detector from the default set:
cfg=[];
cfg.elec=elec;
cfg.fsample=data.fsample;
cfg.include={'absamp'}; %even with one detector, name should be inside a cell
detector_set_absamp=st_get_default_detector_set(cfg);

%note this is still a detector SET.
%we need to extract the (first) detector from the field 'detectors':
dtct_absamp=detector_set_absamp.detectors{1}

%the field 'st' contains the details used by SleepTrip for artifact
%detection. We won't go over every field right now, but 'maxAmp' indicates
%the maximum absolute amplitude the signal is allowed to have.
dtct_absamp.st

%%
%suppose we want evaluate the effect of different 'maxAmp' settings on our
%artifact detection: for example 200/300/400 microV. We make copies of our detector and adjust setting(s)
%as needed (also remember to provide a unique name).

dtct_absamp_200=dtct_absamp;
dtct_absamp_200.label='absamp_200';
dtct_absamp_200.st.thresholdvalue=200;

dtct_absamp_300=dtct_absamp;
dtct_absamp_300.label='absamp_300';
dtct_absamp_300.st.thresholdvalue=300;

dtct_absamp_400=dtct_absamp;
dtct_absamp_400.label='absamp_400';
dtct_absamp_400.st.thresholdvalue=400;

%now we can use st_combine_detectors to turn these into a custom detector
%set. The detector cfgs should be provided inside a cell:
detector_set_different_amplitudes=st_combine_detectors({dtct_absamp_200,dtct_absamp_300,dtct_absamp_400});

%lastly, we need to manually add elec:
detector_set_different_amplitudes.elec=elec

%%
%2c) It's not required to use default detectors as starting points:
%detectors can be built from scratch.
%
%internally, basic detection proceeds in two steps:
%1) [optionally] process the data using ft_preprocessing (FieldTrip
%function). Typically used for filtering, Hilbert transform, etc.
%2) apply thresholds, durations, merge intervals etc. to the data

% We'll postpone listing all the options till later, but
%here's a brief example. Let's assume we're specifically interested in 60-80
%Hz artifacts

my_custom_detector=[]; %intialize (will become a structure)
my_custom_detector.label='mydetector'; %provide a name

%we'll perform filtering using ft_preprocessing. Create
%a substructure "ft" with valid options for ft_preprocessing
my_custom_detector.ft.bpfilter='yes'; %specify (default) bandpass filter
my_custom_detector.ft.bpfreq=[60 80]; %frequency range
my_custom_detector.ft.hilbert='yes'; %extract the amplitude envelope of the filtered signal
my_custom_detector.ft.boxcar = 0.2; %smooth the amplitude envelope using specified window size

%now create a substructure "st"
my_custom_detector.st.method='threshold'; %we want to threshold the provided signal (here, amplitude envelope)
my_custom_detector.st.zscore='yes'; %zscore each channel
my_custom_detector.st.thresholddirection='above';
my_custom_detector.st.thresholdvalue = 2.5; %we consider an artifact when Z > 2.5
my_custom_detector.st.minduration=1; %only consider artifact if duration > 1 s
my_custom_detector.st.artpadding = 1; %extend detected artifact by 1 s (on both sides)
my_custom_detector.st.mergeduration = 0.5; %merge artifacts if inter-artifact gap is < 0.5 s

%finally, remember that we need a detector SET (even when we only have one detector in there)
my_custom_detector_set=st_combine_detectors({my_custom_detector});

%lastly, we need to manually add elec:
my_custom_detector_set.elec=elec

%%
%-----------------------RUN ARTIFACT DETECTION--------------------
%Performing artifact detection is relatively straightforward:
%call st_run_detector_set with 2 arguments:
%-1) a detector set cfg
%-2) data
cfg_artifacts=st_run_detector_set(detector_set_default_all,data)

%cfg_artifacts is an output configuration, with the field "continuous" containing detailed artifact
%information (start and end time, duration, artifact type, channel)

%Uncomment below to run one of the other detector sets we built:
%cfg_artifacts=st_run_detector_set(detector_set_different_amplitudes,data)
%cfg_artifacts=st_run_detector_set(my_custom_detector_set,data)
%%
%---------------------VISUALIZE (CONTINUOUS) ARTIFACTS--------------------------

%-----CODE BELOW DOES NOT WORK!!----
% first 3000 s of data and scoring have been cut, and events have correct timing i.e., start times >3000),
% but events are not shown at correct time. no problems if data and scoring start at t=0.
% to do with timing within st_scorebrowser?
data_start_time=data.time{1}(1); %3000
scoring_offset= scoring.dataoffset; %3000


eventTable=vertcat(cfg_artifacts.continuous.artifacts{:});%concatenate all artifact tables
eventTable = sortrows(eventTable,'start','ascend');
event_one_start=eventTable{1,'start'}; %at or after 3000 (depending on detector set)

%%

eventTable=vertcat(cfg_artifacts.continuous.artifacts{:});


numChan=length(data.label);
yScale=50;

cfg=[];
cfg.signallinewidth=1;
cfg.channeldisplayed=1:numChan;
cfg.chanyrange =repmat([-yScale yScale],[numChan 1]);
cfg.colorgroups='allblack';
cfg.channelcolormap=[0.5 0.5 0.5];

cfg.highlightscoringchannels='no';

cfg.bgcolor = 'white'; % 'white' or 'dark'
cfg.renderer = 'opengl'; % to get things plotted a little faster

numEventTypes=length(unique(eventTable.event));
cfg.events = eventTable;
cfg.eventchannelmatching = 'exact';
cfg.eventhighlighting = 'box';
cfg.eventcolors=colorcube(numEventTypes);
cfg.eventcolorsalpha= 0.3;

cfg.scoring = scoring;

st_scorebrowser(cfg, data);