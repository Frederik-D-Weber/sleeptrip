function scoring = st_autoscore(cfg, data)

% ST_AUTOSCORE does autoscoring on the data through several APIs or
% methods.
%
% Use as
%   scoring = st_autoscore(cfg, data)
%
% Required configuration parameters are:
%   cfg.method      = currently only 'z3scoreV1' using EEG and EOG is supported
%   cfg.interactive = either 'yes' or 'no' (default = 'yes')
%   cfg.confirm     = either 'yes' or 'no' or 'yesifpaid' (default = 'yesifpaid')
%   cfg.chanorder   = maunal overwrite when not in interactive mode e.g. 
%   cfg.dataoffset  = offset in seconds of the scoring to the data (default = 0)
%   cfg.duration    = the duration to maximally use in epochs since data offset
%                     (default = Inf)
%
% Optional configuration parameters are:
%   cfg.channel       = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION
%   cfg.downsamplefs  = downsample the data to this frequency in Hz before doing the anlysis (default = 100)
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
% METHODS
%    'z3scoreV1' and 'z3scoreV2', see https://z3score.com for details on
%                                 the algorithm
%        optional parameters are 
%            cfg.chanorder = a Nx1 cell string ordered in the with N = 4/5 channels in the data for 
%                    C3:A2
%                    C4:A1
%                    EOGl:A1/EOGl
%                    EOGr:A2/EOGr
%                    EMG (for V2 only)
%            cfg.licensefile = string, file path with a license file, haveing two
%                    lines. First line is username and second the api key,
%                    e.g,:
%                    example@email.com
%                    1234567890abcdefghijklmnopqrstuvwxyz12345678
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
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% set defaults
downsamplefsNotSet = false;
if ~isfield(cfg, 'downsamplefs')
    downsamplefsNotSet = true;
end
cfg.downsamplefs     = ft_getopt(cfg, 'downsamplefs', 100);
cfg.channel  = ft_getopt(cfg, 'channel', 'all', 1);
cfg.interactive  = ft_getopt(cfg, 'interactive', 'yes');
cfg.confirm  = ft_getopt(cfg, 'interactive', 'yesifpaid');

cfg.dataoffset  = ft_getopt(cfg, 'dataoffset', 0);
cfg.duration  = ft_getopt(cfg, 'duration', Inf);

ispaidservice = false;
switch cfg.method
    case {'z3scoreV1' 'z3scoreV2'}
        ispaidservice = true;
end

doconfirm = false;               
if (strcmp(cfg.confirm, 'yesifpaid') && ispaidservice) || strcmp(cfg.confirm,'yes')
    doconfirm = true;
end

fprintf([functionname ' function initialized\n']);

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


cfg_int = cfg;
if isfield(cfg, 'reref'),       cfg_int.reref = cfg.reref;             end
if isfield(cfg, 'refchannel'),  cfg_int.refchannel = cfg.refchannel;   end
if isfield(cfg, 'refmethod'),   cfg_int.refmethod = cfg.refmethod;     end
if isfield(cfg, 'implicitref'), cfg_int.implicitref = cfg.implicitref; end
if isfield(cfg, 'montage'),     cfg_int.montage = cfg.montage;         end

if hasdata
    cfg_int.channel = ft_channelselection(cfg.channel, data.label);
else
    cfg_int.channel = ft_channelselection(cfg.channel, hdr.label);
end

if hasdata
    data_t = st_preprocessing(cfg_int, data);
    if isempty(data_t.trial)
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
    cfg_int.continuous   = 'yes';
    if isempty(cfg_int.trl)
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

cfg.chanorder   = ft_getopt(cfg, 'chanorder', 1:numel(data.label));

if (data.fsample == 128) && downsamplefsNotSet
    ft_warning('leaving 128 Hz sampling rate as default sampling rate')
    cfg.downsamplefs = 128;
end

if (cfg.downsamplefs < data.fsample) && ~downsamplefsNotSet
    fprintf('resample data from %i to %i Hz\n',data.fsample,cfg.downsamplefs);
    cfg_resample = [];
    cfg_resample.resamplefs = cfg.downsamplefs;
    cfg_resample.detrend = 'no';
    %cfg_resample.feedback = core_cfg.feedback;
    data = ft_resampledata(cfg_resample,data);
elseif (cfg.downsamplefs > data.fsample)
    ft_error('upsampling prohibited, choose cfg.downsamplefs (now as %d Hz)  smaller or equal to the sampling frequency of provided data %d Hz', cfg.downsamplefs, data.fsample);
end

fsample = data.fsample;

cfg_sd = [];
cfg_sd.latency = [min(data.time{1})+cfg.dataoffset min(data.time{1})+min(cfg.dataoffset+cfg.duration,max(cfg.dataoffset+data.time{1}(1:end-1))-min(data.time{1}))];
data = ft_selectdata(cfg_sd, data);


channelorder = 1:numel(data.label);
if isfield(cfg,'chanorder')
    channelorder = [];
    for iCh = 1:numel(cfg.chanorder)
        [dummy chidx] = ismember(cfg.chanorder{iCh},data.label);
        channelorder = [channelorder chidx];
    end
end
channelorder = channelorder(:);

fprintf([functionname ' data prepared\n']);

scoring = [];
switch cfg.method
    case {'z3scoreV1' 'z3scoreV2'}
        if cfg.duration < 12
            ft_error('Z3score needs at least 6 minutes or 12 30-second epochs of EEG/EOG data for scoring. \n chose a higher cfg.duration')
        end
        if (max(data.time{1}) - min(data.time{1})) < 60*6
            ft_error('Z3score needs at least 6 minutes of EEG/EOG data for scoring. \n chose more data or change cfg.dataoffset')
        end
        switch cfg.method
            case 'z3scoreV1'
                if data.fsample < 100
                    ft_error('Z3scoreV1 needs at least a sample rate of 100 Hz for EEG/EOG data.')
                end
            case 'z3scoreV2'
                if data.fsample < 200
                    ft_error('Z3scoreV2 needs at least a sample rate of 200 Hz for EMG (and thus also for EEG and EOG) data.')
                end
        end
        ask_again = true;
        while ask_again
            cfg_tmp = [];
            %FIXME: in version 2 also EMG requires 200 Hz sampling rate if
            %used, EOG and EEG can be lower at 100 Hz, so resampling would
            %need to be considered or data provided in 200 Hz fully.
            cfg_tmp.version = cfg.method;
            if isfield(cfg,'channel')
                cfg_tmp.channel = channelorder(1:5);
            end
            if isfield(cfg,'licensefile')
                cfg_tmp.licensefile = cfg.licensefile;
            end
            sca = AutoSleepScoringZ3score(cfg_tmp,data);
            answer_temp = questdlg(sca.status(),'Automatic Scoring Status', ...
                'OK','OK');
            
            %%sca = sca.update_channel();
            answer_sca = questdlg('Perform Scoring?', ...
                'Sure to try Automatic Sleep Scoring now?', ...
                'Yes','Change','Cancel','Cancel');
            switch answer_sca
                case 'Yes'
                    sca = sca.score();
                    scoring = sca.scoring;
                    scoring.cfg = cfg;
                    import_success = true;
                    ask_again = false;
                case 'Change'
                    ask_again = true;
                case 'Cancel'
                    ask_again = false;
                    return
            end
        end
    otherwise
        ft_error('cfg.method = ''%s'' is unknown', cfg.method);
end



scoring.dataoffset = cfg.dataoffset;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous data
% data = datanew
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data

fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end
