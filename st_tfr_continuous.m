function freq = st_tfr_continuous(cfg, data)

% ST_CONTINUOUS calculates a time-frequency plot of continous data
%
% Use as
%   freq = st_tfr_continuous(cfg, data)
%
% Optional configuration parameters are:
%
% cfg.channel  = select channels;
% cfg.approach = the approach to calculate either 'mtmconvol_memeff', 'mtmfft_segments', 'spectrogram' (default = 'spectrogram');
% cfg.length   = timewindow/segment length in seconds (default = 30);
% cfg.overlap  = the amount of overlap for each segment (default = 0.85);
% cfg.taper    = choose the tapper either 'hanning', 'hanning_proportion'
%                or 'dpss' (see ft_freqanalysis for details)
% cfg.windowproportion = if cfg.taper = 'hanning_proportion' then the window proportion can be set to here (default = 0.2 i.e. 10% left and right of each segment is 'hanninged' leaving the ceter 80% not tapered or attenuated);
% cfg.foi      = the frequencies of interrest with their granularity (default = 0.5:0.5:30);
% cfg.transform = power value transformations in a different scaleing 'none' 'db' 'db(p+1)' 'log10' (like in Helfrich et al. 2018) or 'log10(p+1)' (default = 'none');
% cfg.powvalue  = non yet relevant modifier for cfg.approach = 'spectromgram'  with either 'power' or 'psd' (default = 'psd');
% cfg.feedback = feedback flag (default = 'text');
%
% See also FT_PREPROCESSING, FT_FREQANALYSIS

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
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
    return
end

data = ft_checkdata(data, 'datatype', {'raw', 'comp'}, 'feedback', 'yes', 'hassampleinfo', 'yes');

% set defaults
cfg.feedback = ft_getopt(cfg, 'feedback', 'text');
cfg.channel  = ft_getopt(cfg, 'channel', 'all', 1);
cfg.approach = ft_getopt(cfg, 'approach', 'spectrogram');
cfg.length  = ft_getopt(cfg, 'length', 30);
cfg.overlap  = ft_getopt(cfg, 'overlap', 0.85);
cfg.taper  = ft_getopt(cfg, 'taper', 'dpss');
cfg.foi   = ft_getopt(cfg, 'foi', 0.5:0.5:30);
switch cfg.approach
    case {'mtmconvol_memeff', 'mtmconvol'}
        cfg.tapsmofrq = ft_getopt(cfg, 'tapsmofrq', repmat(0.5,1,numel(cfg.foi)));
    otherwise
        cfg.tapsmofrq = ft_getopt(cfg, 'tapsmofrq', 0.5);
end
cfg.transform = ft_getopt(cfg, 'transform', 'none');
cfg.windowproportion = ft_getopt(cfg, 'windowproportion', 0.2);
cfg.powvalue = ft_getopt(cfg, 'powvalue', 'psd');


fprintf([functionname ' function initialized\n']);


cfg_sd = [];
cfg_sd.channel = cfg.channel;
data = ft_selectdata(cfg_sd, data);

switch cfg.approach
    case 'mtmconvol_memeff'
        
        cfg_fa           = [];
        %cfg.channel   = ch1;
        cfg_fa.method    = 'mtmconvol_memeff';
        cfg_fa.foi       =  cfg.foi; % 0.5 Hz steps
        cfg_fa.tapsmofrq =  cfg.tapsmofrq;
        cfg_fa.taper     = 'dpss'; % 'hanning' or 'dpss'
        cfg_fa.t_ftimwin = repmat(cfg.length,1,numel(cfg_fa.foi));
        %cfg_fa.toi       = time(1):(cfg_fa.t_ftimwin*(1-0.85)):time(end);
        cfg_fa.toi       = [num2str(cfg.overlap*100) '%']; % 85% overlap of cfg.t_ftimwin
        freq        = ft_freqanalysis(cfg_fa, data);
        
        %
        % % to save memory for extensive con
        % timeStep = (30*0.85);
        % padding = timeStep*1.5;
        % time_junck = timeStep*50;
        % data_freqs = {};
        % for iJunk = 1:(ceil(max(data.time{1}/time_junck)))
        %     cfg = [];
        %     cfg.channel = [ch '*'];
        %     cfg.latency = [max((iJunk-1)*time_junck - padding,min(data.time{1})),
        %         min((iJunk)*time_junck+padding,max(data.time{1}))];
        %     data_temp = ft_selectdata(cfg, data);
        %
        %
        %     cfg           = [];
        %     %cfg.channel   = ch1;
        %     cfg.method    = 'mtmconvol';
        %     cfg.foi       = 0.5:0.5:1; % 0.5 Hz steps
        %     cfg.tapsmofrq = repmat(0.5,1,numel(cfg.foi));
        %     cfg.taper     = 'dpss'; % 'hanning' or 'dpss'
        %     cfg.t_ftimwin = repmat(30,1,numel(cfg.foi));
        %     %cfg.toi       = time(1):(cfg.t_ftimwin*(1-0.85)):time(end);
        %     cfg.toi       = '85%'; % 85% overlap of cfg.t_ftimwin
        %     data_freq_temp = ft_freqanalysis(cfg, data_temp);
        %
        %     data_temp = [];
        %     data_freqs{iJunk} = data_freq_temp;
        %     data_freq_temp = [];
        % end
        %
        % save(['data_freqs-' num2str(iDataset)], 'data_freqs')
        % load(['data_freqs-' num2str(iDataset)])
        %
        % for iJunk = 1:(ceil(max(data.time{1}/time_junck)))
        %     cfg = [];
        %     cfg.latency = [max((iJunk-1)*time_junck-0.01 ,min(data.time{1})),
        %         min((iJunk)*time_junck+0.01,max(data.time{1}))];
        %     data_freqs{iJunk} = ft_selectdata(cfg, data_freqs{iJunk});
        % end
        %
        % cfg = [];
        % cfg.appenddim = 'time';
        % cfg.tolerance = 0.01;
        % cfg.parameter  = 'powspctrm';
        % data_freq = ft_appendfreq(cfg,data_freqs{:});
        % unique(diff(data_freq.time))
        %
        % cfg = [];
        % cfg.appenddim = 'time';
        % cfg.tolerance = 1;
        % cfg.parameter  = 'powspctrm';
        % data_freq_tt = ft_appendfreq(cfg,data_freq_temp3,data_freq_temp);
        %
        
        %% alternative
    case 'mtmfft_segments'
        
        cfg_rd = [];
        cfg_rd.length = cfg.length;
        cfg_rd.overlap = cfg.overlap;
        data_segmented = ft_redefinetrial(cfg_rd, data);
        
        %data = [];

        switch cfg.taper
            case 'dpss'
                cfg_fa           = [];
                %cfg.channel   = ch1;
                cfg_fa.method     = 'mtmfft';
                cfg_fa.taper      = 'dpss';
                cfg_fa.foi       =  cfg.foi; % 0.5 Hz steps
                cfg_fa.tapsmofrq = cfg.tapsmofrq;
                %cfg.toi       = time(1):(cfg.t_ftimwin*(1-0.85)):time(end);
                %cfg.toi       = '85%'; % 85% overlap of cfg.t_ftimwin
                cfg_fa.keeptrials = 'yes';
                data_freq_segmented = ft_freqanalysis(cfg_fa, data_segmented);
                
            case {'hanning' 'hanning_proportion'}
                cfg_fa = [];
                cfg_fa.method     = 'mtmfft';
                cfg_fa.taper      = cfg.taper;
                cfg_fa.foi        = cfg.foi;
                cfg_fa.keeptrials = 'yes';
                data_freq_segmented = ft_freqanalysis(cfg_fa, data_segmented);
                
            otherwise
                ft_error('cfg.taper = ''%s'' not supported.',cfg.taper)
        end
        
        begsample = data_segmented.sampleinfo(:,1);
        endsample = data_segmented.sampleinfo(:,2);
        
        freq           = data_freq_segmented;
        freq.powspctrm = permute(data_freq_segmented.powspctrm, [2, 3, 1]);
        freq.dimord    = 'chan_freq_time'; % it used to be 'rpt_chan_freq'
        freq.time      = ((begsample+endsample)/2) / data_segmented.fsample; % add the description of the time dimension
        
        data_freq_segmented = [];
    case 'spectrogram'
        ft_progress('init', cfg.feedback, 'calculating spectrograms channel and tapper wise');

        windowsegment = data.fsample*cfg.length;
        switch cfg.taper
            case 'hanning'
                windowsegment = hanning(windowsegment);
                windowsegment = windowsegment./norm(windowsegment, 'fro');

            case 'hanning_proportion'
                
                temp_windowFunctionValues = window('hanning', floor(cfg.windowproportion*windowsegment));                
                windowFunctionValuesLeft = temp_windowFunctionValues(1:floor(end/2));
                windowFunctionValuesRight = temp_windowFunctionValues(floor(end/2)+1:end);
                windowFunctionValues = ones(windowsegment,1);
                windowFunctionValues(1:length(windowFunctionValuesLeft)) = windowFunctionValuesLeft;
                windowFunctionValues(end-length(windowFunctionValuesRight)+1:end) = windowFunctionValuesRight;
                
                windowsegment = windowFunctionValues;
                windowsegment = windowsegment./norm(windowsegment, 'fro');
            case 'dpss'
                
                %num_seq = numel(cfg.foi);
                time_halfbandwidth = double(windowsegment*(cfg.tapsmofrq./data.fsample));
                [windowsegment, lambda] = dpss(double(windowsegment),time_halfbandwidth(1));
                windowsegment = windowsegment';
                % remove the last unideal taper 
                windowsegment = windowsegment(1:(end-1), :);
        
            otherwise
        end
        nf = data.fsample/2;
        if max(cfg.foi > nf)
            ft_error('maximal frequency can only be %d with sampling rate of %d. you can limit the frequency in cfg.foi = [xxx %d]',floor(nf),data.fsample,floor(nf))
        end
        f = round(nf/median(diff(cfg.foi)))*2;
        
        nChannel = numel(data.label);
        freq = [];
        freq.label = data.label;
        freq.dimord    = 'chan_freq_time'; % it used to be 'rpt_chan_freq'
        freq.freq = [];
        freq.powspctrm = [];
        freq.time = [];
        for iCh = 1:nChannel
            %nfft_points = fix(diff(freq_range)/0.05);
            switch cfg.taper
                case 'dpss'
                    nTap = size(windowsegment,1);
                    lenTap = size(windowsegment,2);
                    for iTap = 1:nTap
                        ft_progress((iCh*iTap)/(nChannel*nTap), 'spectrogram: channelnum %d tap %d',iCh,iTap);
                        [s,freqs,t_seconds] = spectrogram(data.trial{1}(iCh,:),windowsegment(iTap,:),round(data.fsample*cfg.length*cfg.overlap),f,data.fsample,'yaxis',cfg.powvalue);
                        if iTap == 1
                            timefreq = zeros(size(s));
                        end
                        s = s .* sqrt(2 ./ lenTap);
                        timefreq = timefreq + abs(s).^2;
                        %timefreq = timefreq + ps;
                    end
                    timefreq = timefreq./nTap;
                otherwise

                    ft_progress(iCh/nChannel,'spectrogram: channelnum %d tap %d',iCh,1);
                    
                    [s,freqs,t_seconds] = spectrogram(data.trial{1}(iCh,:),windowsegment,round(data.fsample*cfg.length*cfg.overlap),f,data.fsample,'yaxis',cfg.powvalue);
                    s = s .* sqrt(2 ./ (data.fsample*cfg.length));
                    timefreq = abs(s) .^2;
                    %timefreq =  ps;

            end
            nfreqpoints = numel(freqs);
            %ffreqs = flip(freqs);
            %timfreq = flip(10*log10(abs(s)+1));
            if iCh == 1
                freq.freq = freqs';
                freq.powspctrm = nan(nChannel,nfreqpoints,numel(t_seconds));
                freq.time      = t_seconds; % add the description of the time dimension
            end
            freq.powspctrm(iCh,:,:) = timefreq;
        end
        cfg_sd = [];
        cfg_sd.frequency = [min(cfg.foi) max(cfg.foi)];
        if ~isempty(freq.powspctrm)
            freq = ft_selectdata(cfg_sd,freq);
        end
        ft_progress('close');
    otherwise
        ft_error('cfg.approach = ''%s'' is unknown',cfg.approach)
end




switch cfg.transform
    case 'none'
        
    case 'db'
        freq.powspctrm(:) = 10*log10(freq.powspctrm(:));
    case 'db(p+1)'
        freq.powspctrm(:) = freq.powspctrm(:)+1;
        freq.powspctrm(:) = 10*log10(freq.powspctrm(:));
    case  'log10' %helfrich et al.
        freq.powspctrm(:) = log10(freq.powspctrm(:));
    case  'log10(p+1)'
        freq.powspctrm(:) = freq.powspctrm(:)+1;
        freq.powspctrm(:) = log10(freq.powspctrm(:));
    otherwise
        ft_error('cfg.transform = ''%s'' unknown.',cfg.transform)
        % for iCh = 1:size(data_freq.powspctrm,1)
        %     for iFreq = 1:size(data_freq.powspctrm,2)
        %         data_freq.powspctrm(iCh,iFreq,:) = data_freq.powspctrm(iCh,iFreq,:) * data_freq.freq(iFreq);
        %     end
        % end
end



% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous data
ft_postamble provenance data_freq
ft_postamble history    data_freq
ft_postamble savevar    data_freq


fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end
